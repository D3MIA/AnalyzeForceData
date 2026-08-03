#!/usr/bin/env python3
"""Generate a small, synthetic ``analysis_data.json`` sample.

The real analysis export (produced by ``surgical_force_analysis.ipynb``) is a
single JSON blob that bundles, for every trial of the study:

* per-frame time series (position, velocity, force, inter-instrument distance …)
* per-trial and per-participant summary scalars
* two flat summary tables ready to drop into a ``DataFrame``

A full export is large (tens of MB) because every trial keeps thousands of
frames per signal. This script rebuilds the *exact same schema* with a handful
of short, physically plausible trials so the file stays a few hundred KB and can
be committed to the repo as a fixture for developing / testing plotting and
reporting code.

The numbers here are **fake** – they are drawn from smooth random processes with
roughly the same magnitudes as the real recordings. They are not real surgical
measurements and should only be used as a structural example.

Run:  ``python Examples/generate_sample_analysis_data.py``
Output: ``Examples/sample_analysis_data.json``
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

# numpy 1.x / 2.x compatibility (trapz -> trapezoid, ndarray.ptp -> np.ptp)
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
_ptp = getattr(np, "ptp", None)

# ---------------------------------------------------------------------------
# Study configuration (mirrors the `meta` block of a real export)
# ---------------------------------------------------------------------------
INSTRUMENTS = ["Bipolar", "Cavitron", "Scissors"]
PAIRS = [["Bipolar", "Cavitron"], ["Bipolar", "Scissors"]]
PAIR_KEYS = ["Bipolar-Cavitron", "Bipolar-Scissors"]
PARTICIPANTS = ["Bianca", "Matheus", "Mohammed"]

SMOOTH_WINDOW = 11        # frames, matches the notebook's smoothing kernel
IDLE_SPEED = 5.0          # mm/s below which a tip is considered idle
VOXEL_SIZE = 2.0          # mm, occupancy-grid voxel edge
PLANE_DIST_MAX = 50.0     # mm, in-use gate for the paired instruments
FS = 20.0                 # Hz sampling rate used to fabricate timestamps

# Keep it small: 3 participants x 2 trials, ~60 frames each.
TRIALS_PER_PARTICIPANT = 2
FRAMES_RANGE = (55, 68)

# Rough tip "home" positions (mm) in the common registered frame, so the fake
# trajectories look like the real ones (README §Data layout).
CENTERS = {
    "Bipolar": np.array([185.0, 246.0, -522.0]),
    "Cavitron": np.array([172.0, 250.0, -515.0]),
    "Scissors": np.array([205.0, 233.0, -530.0]),
}
# Per-instrument motion "energy" – Cavitron sweeps most, Scissors is used least.
MOTION_SCALE = {"Bipolar": 9.0, "Cavitron": 16.0, "Scissors": 6.0}
# Fraction of the trial each instrument is actually engaged with tissue.
INUSE_TARGET = {"Bipolar": 0.78, "Cavitron": 0.75, "Scissors": 0.11}

# The 4 registration fiducials (BipolarCollectedPoint0..3) in the common frame — the
# quadrilateral working area the tips move within (README §4). Ordered around the
# rectangle; roughly encloses the tip extent above.
FIDUCIALS = [
    [155.0, 226.0, -520.0],
    [213.0, 226.0, -520.0],
    [213.0, 268.0, -520.0],
    [155.0, 268.0, -520.0],
]

rng = np.random.default_rng(20260728)  # deterministic output


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def smooth(a: np.ndarray, w: int = SMOOTH_WINDOW) -> np.ndarray:
    """Centred moving average, edge-padded (cheap stand-in for the notebook)."""
    if w <= 1:
        return a
    k = np.ones(w) / w
    pad = w // 2
    ap = np.pad(a, (pad, pad), mode="edge")
    return np.convolve(ap, k, mode="valid")[: len(a)]


def smooth_walk(n: int, scale: float) -> np.ndarray:
    """A smooth 1-D random walk with zero mean, scaled to ~`scale` mm amplitude."""
    steps = rng.standard_normal(n)
    walk = np.cumsum(steps)
    walk = smooth(walk, 9)
    walk -= walk.mean()
    m = np.abs(walk).max() or 1.0
    return walk / m * scale


def nan_list(a: np.ndarray) -> list:
    """numpy array -> plain python list, keeping NaN (JSON emits `NaN`)."""
    return [float(x) for x in a]


def trajectory(center: np.ndarray, n: int, scale: float) -> np.ndarray:
    """(n, 3) smooth tip path around `center`."""
    return np.column_stack([center[i] + smooth_walk(n, scale) for i in range(3)])


def tracking_status(n: int) -> np.ndarray:
    """Mostly tracked (1) with a couple of short dropout runs (0)."""
    st = np.ones(n, dtype=int)
    for _ in range(rng.integers(0, 3)):
        start = rng.integers(0, n - 4)
        st[start : start + rng.integers(1, 4)] = 0
    return st


def inuse_mask(n: int, frac: float) -> np.ndarray:
    """Boolean per-frame engagement mask covering ~`frac` of the trial."""
    base = smooth(rng.standard_normal(n), 15)
    thr = np.quantile(base, 1.0 - frac)
    return base >= thr


# ---------------------------------------------------------------------------
# Per-instrument kinematics for one trial
# ---------------------------------------------------------------------------
def instrument_signals(name: str, n: int, dt: float):
    pos = trajectory(CENTERS[name], n, MOTION_SCALE[name])
    track = tracking_status(n)
    use = inuse_mask(n, INUSE_TARGET[name]) & (track == 1)

    # speed / accel / jerk from finite differences, in mm/s, mm/s^2, mm/s^3
    vel_vec = np.gradient(pos, dt, axis=0)
    speed = smooth(np.linalg.norm(vel_vec, axis=1))
    accel = np.abs(np.gradient(speed, dt))
    jerk = np.abs(np.gradient(accel, dt))

    # angular signals (deg/s, deg/s^2) – small, smooth
    ang_speed = np.abs(smooth(rng.standard_normal(n))) * (5 + 8 * rng.random())
    ang_accel = np.abs(np.gradient(ang_speed, dt))

    # mask everything not in use (the notebook sets these to NaN)
    def mask(a):
        out = a.astype(float).copy()
        out[~use] = np.nan
        return out

    pos_masked = pos.astype(float).copy()
    pos_masked[~use] = np.nan

    return {
        "pos": pos,                     # raw, for distances
        "pos_plot": pos_masked,         # NaN where not in use
        "track": track,
        "use": use,
        "velocity": mask(speed),
        "acceleration": mask(accel),
        "jerk": mask(jerk),
        "ang_speed": mask(ang_speed),
        "ang_accel": mask(ang_accel),
    }


def occupancy_effective(pos: np.ndarray, use: np.ndarray) -> tuple[float, int]:
    """exp(H) of the dwell distribution over VOXEL_SIZE voxels, + voxel count."""
    p = pos[use]
    if len(p) == 0:
        return float("nan"), 0
    keys = np.floor(p / VOXEL_SIZE).astype(int)
    _, counts = np.unique(keys, axis=0, return_counts=True)
    prob = counts / counts.sum()
    H = -np.sum(prob * np.log(prob))
    return float(np.exp(H)), int(len(counts))


def ellipsoid(pos: np.ndarray, use: np.ndarray) -> tuple[float, float]:
    """Covariance-ellipsoid volume (mm^3) and anisotropy lambda1/lambda3."""
    p = pos[use]
    if len(p) < 4:
        return float("nan"), float("nan")
    cov = np.cov(p.T)
    lam = np.sort(np.abs(np.linalg.eigvalsh(cov)))[::-1]
    vol = (4.0 / 3.0) * math.pi * math.sqrt(max(lam.prod(), 1e-9))
    aniso = float(lam[0] / max(lam[2], 1e-9))
    return float(vol), aniso


# ---------------------------------------------------------------------------
# One trial
# ---------------------------------------------------------------------------
def make_trial(participant: str, trial_idx: int) -> dict:
    n = int(rng.integers(*FRAMES_RANGE))
    dt = 1.0 / FS
    duration = n * dt
    tnorm = np.linspace(0.0, 1.0, n)

    sig = {name: instrument_signals(name, n, dt) for name in INSTRUMENTS}

    # --- force (separate sensor, never masked) -----------------------------
    fmag = np.abs(smooth(rng.standard_normal(n))) * 0.02 + 0.004
    dFdt = np.abs(np.gradient(fmag, dt))
    torque = np.abs(smooth(rng.standard_normal(n))) * 0.4 + 0.05

    # --- inter-instrument distance / angle for each pair -------------------
    dist, angle = {}, {}
    for (a, b), key in zip(PAIRS, PAIR_KEYS):
        d = np.linalg.norm(sig[a]["pos"] - sig[b]["pos"], axis=1)
        # gate: not-in-use frames -> NaN (mirrors PLANE_DIST_MAX behaviour)
        both_use = sig[a]["use"] & sig[b]["use"]
        d = np.where(both_use, d, np.nan)
        ang = smooth(rng.standard_normal(n)) * 20 + (60 if key.endswith("Cavitron") else 100)
        ang = np.where(both_use, ang, np.nan)
        dist[key] = d
        angle[key] = ang

    # --- per-instrument summary scalars ------------------------------------
    def scal(fn):
        return {name: fn(name) for name in INSTRUMENTS}

    inuse_frac = scal(lambda nm: float(sig[nm]["use"].mean()))
    pathlen = scal(lambda nm: float(np.nansum(np.abs(np.diff(
        np.where(sig[nm]["use"][:, None], sig[nm]["pos"], np.nan), axis=0))).sum()))
    netdisp = scal(lambda nm: float(np.linalg.norm(sig[nm]["pos"][-1] - sig[nm]["pos"][0])))
    straightness = scal(lambda nm: float(netdisp[nm] / pathlen[nm]) if pathlen[nm] else float("nan"))
    bbox_vol = scal(lambda nm: float(np.prod(
        _ptp(sig[nm]["pos"][sig[nm]["use"]], axis=0))) if sig[nm]["use"].any() else float("nan"))
    idle_frac = scal(lambda nm: float(np.mean(
        np.nan_to_num(sig[nm]["velocity"], nan=0.0) < IDLE_SPEED)))
    ell = {nm: ellipsoid(sig[nm]["pos"], sig[nm]["use"]) for nm in INSTRUMENTS}
    occ = {nm: occupancy_effective(sig[nm]["pos"], sig[nm]["use"]) for nm in INSTRUMENTS}

    name = f"SurgicalTrainingRecording_2026072{trial_idx}_1353{50 + trial_idx:02d}"

    return {
        "participant": participant,
        "trial": trial_idx,
        "name": name,
        "label": f"{participant}·T{trial_idx}",
        "n_frames": n,
        "duration": duration,
        "rmse": float(rng.random() * 1e-13),   # registration is ~perfect on fakes
        "present_instruments": list(INSTRUMENTS),
        "tnorm": nan_list(tnorm),
        "fmag": nan_list(fmag),
        "dFdt": nan_list(dFdt),
        "kin": {nm: {
            "velocity": nan_list(sig[nm]["velocity"]),
            "acceleration": nan_list(sig[nm]["acceleration"]),
            "jerk": nan_list(sig[nm]["jerk"]),
        } for nm in INSTRUMENTS},
        "pos_plot": {nm: [[float(c) for c in row] for row in sig[nm]["pos_plot"]]
                     for nm in INSTRUMENTS},
        "ang_speed": {nm: nan_list(sig[nm]["ang_speed"]) for nm in INSTRUMENTS},
        "ang_accel": {nm: nan_list(sig[nm]["ang_accel"]) for nm in INSTRUMENTS},
        "tracking_status": {nm: [int(x) for x in sig[nm]["track"]] for nm in INSTRUMENTS},
        "dist": {k: nan_list(v) for k, v in dist.items()},
        "angle": {k: nan_list(v) for k, v in angle.items()},
        "force_mean": float(np.mean(fmag)),
        "force_peak": float(np.max(fmag)),
        "force_cov": float(np.std(fmag) / np.mean(fmag) * 100),
        "impulse": float(_trapz(fmag, dx=dt)),
        "dFdt_abs_mean": float(np.mean(dFdt)),
        "torque_mean": float(np.mean(torque)),
        "inuse_frac": inuse_frac,
        "pathlen": pathlen,
        "netdisp": netdisp,
        "straightness": straightness,
        "bbox_vol": bbox_vol,
        "idle_frac": idle_frac,
        "ell_vol": {nm: ell[nm][0] for nm in INSTRUMENTS},
        "aniso": {nm: ell[nm][1] for nm in INSTRUMENTS},
        "occ_eff": {nm: occ[nm][0] for nm in INSTRUMENTS},
        "occ_voxels": {nm: occ[nm][1] for nm in INSTRUMENTS},
        # keep the per-frame masks around so the table builder can reuse them
        "_sig": sig,
        "_dist": dist,
        "_angle": angle,
    }


# ---------------------------------------------------------------------------
# Summary tables (flat rows, one per trial / participant)
# ---------------------------------------------------------------------------
def nanmean(a):
    a = np.asarray(a, float)
    return float(np.nanmean(a)) if np.any(~np.isnan(a)) else float("nan")


def per_trial_row(t: dict) -> dict:
    sig = t["_sig"]
    row = {
        "participant": t["participant"],
        "trial": t["trial"],
        "file": t["name"],
        "n_frames": t["n_frames"],
        "duration_s": round(t["duration"], 3),
        "reg_rmse_mm": round(t["rmse"], 3),
        "force_mean_N": round(t["force_mean"], 5),
        "force_peak_N": round(t["force_peak"], 5),
        "force_cov_pct": round(t["force_cov"], 3),
        "force_impulse": round(t["impulse"], 4),
        "force_dFdt_Nps": round(t["dFdt_abs_mean"], 4),
        "torque_mean": round(t["torque_mean"], 5),
    }
    for nm in INSTRUMENTS:
        row[f"{nm}_tracked_pct"] = round(float(np.mean(sig[nm]["track"])) * 100, 2)
        row[f"{nm}_inuse_pct"] = round(t["inuse_frac"][nm] * 100, 2)
        row[f"{nm}_path_mm"] = round(t["pathlen"][nm], 1)
        row[f"{nm}_speed_mean"] = round(nanmean(sig[nm]["velocity"]), 2)
        row[f"{nm}_speed_peak"] = round(float(np.nanmax(sig[nm]["velocity"])), 2)
        row[f"{nm}_accel_mean"] = round(nanmean(sig[nm]["acceleration"]), 2)
        row[f"{nm}_jerk_mean"] = round(nanmean(sig[nm]["jerk"]), 1)
        row[f"{nm}_straightness"] = round(t["straightness"][nm], 3)
        row[f"{nm}_workvol_mm3"] = round(t["bbox_vol"][nm], 1)
        row[f"{nm}_ellvol_mm3"] = round(t["ell_vol"][nm], 1)
        row[f"{nm}_aniso"] = round(t["aniso"][nm], 2)
        row[f"{nm}_occ_eff"] = round(t["occ_eff"][nm], 1)
        row[f"{nm}_idle_frac"] = round(t["idle_frac"][nm], 3)
        row[f"{nm}_angspeed_dps"] = round(nanmean(sig[nm]["ang_speed"]), 2)
        row[f"{nm}_angaccel_dps2"] = round(nanmean(sig[nm]["ang_accel"]), 2)
    for key in PAIR_KEYS:
        d, a = t["_dist"][key], t["_angle"][key]
        row[f"dist_{key}_mean"] = round(nanmean(d), 1)
        row[f"dist_{key}_min"] = round(float(np.nanmin(d)), 1) if np.any(~np.isnan(d)) else float("nan")
        row[f"dist_{key}_max"] = round(float(np.nanmax(d)), 1) if np.any(~np.isnan(d)) else float("nan")
        row[f"angle_{key}_mean"] = round(nanmean(a), 1)
    return row


def per_participant_row(participant: str, rows: list) -> dict:
    cols = [k for k in rows[0] if k not in ("participant", "trial", "file", "n_frames")]
    agg = {"participant": participant, "n_trials": len(rows)}
    for c in cols:
        vals = [r[c] for r in rows if isinstance(r[c], (int, float)) and not (
            isinstance(r[c], float) and math.isnan(r[c]))]
        agg[c] = round(float(np.mean(vals)), 3) if vals else float("nan")
    return agg


# ---------------------------------------------------------------------------
# Build & write
# ---------------------------------------------------------------------------
def main() -> None:
    trials = []
    for p in PARTICIPANTS:
        for k in range(1, TRIALS_PER_PARTICIPANT + 1):
            trials.append(make_trial(p, k))

    per_trial = [per_trial_row(t) for t in trials]
    per_participant = []
    for p in PARTICIPANTS:
        rows = [r for r in per_trial if r["participant"] == p]
        per_participant.append(per_participant_row(p, rows))

    # strip the private helper keys before serialising
    for t in trials:
        for k in ("_sig", "_dist", "_angle"):
            t.pop(k, None)

    out = {
        "meta": {
            "instruments": INSTRUMENTS,
            "pairs": PAIRS,
            "pair_keys": PAIR_KEYS,
            "participants": PARTICIPANTS,
            "smooth_window": SMOOTH_WINDOW,
            "idle_speed": IDLE_SPEED,
            "voxel_size": VOXEL_SIZE,
            "plane_dist_max": PLANE_DIST_MAX,
            "time_unit": "s",
            "fiducials": FIDUCIALS,
        },
        "trials": trials,
        "tables": {"per_trial": per_trial, "per_participant": per_participant},
    }

    dest = os.path.join(os.path.dirname(__file__), "sample_analysis_data.json")
    with open(dest, "w") as f:
        json.dump(out, f)  # default allow_nan=True -> emits `NaN` like the real export
    size_kb = os.path.getsize(dest) / 1024
    print(f"Wrote {dest} ({size_kb:.1f} KB, {len(trials)} trials)")


if __name__ == "__main__":
    main()
