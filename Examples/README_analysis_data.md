# `analysis_data.json` — schema & sample

`surgical_force_analysis.ipynb` bundles everything it computes for a study into a
single JSON export, `analysis_data.json`. A real export is large (tens of MB)
because every trial keeps a full per-frame time series for each signal.

This folder ships a **small synthetic copy** so plotting / reporting / dashboard
code can be developed and tested without the real recordings:

| File | What it is |
|------|------------|
| [`sample_analysis_data.json`](sample_analysis_data.json) | ~180 KB fake export: 3 participants × 2 short trials, identical schema |
| [`generate_sample_analysis_data.py`](generate_sample_analysis_data.py) | Reproducible, seeded generator that produced it |

> ⚠️ **The sample is fake.** Values come from smooth random processes scaled to
> roughly the same magnitudes as real recordings. They are *not* surgical
> measurements — use the file only as a structural example. To keep it small the
> trials are ~60 frames (~3 s @ 20 Hz) instead of the real ~2–4 min; the schema
> is byte-for-byte the same, only the arrays are shorter.

Regenerate at any time:

```bash
python Examples/generate_sample_analysis_data.py
```

## Loading

The export is written with `json.dump(..., allow_nan=True)`, so missing samples
appear as the JavaScript-style token `NaN` (not `null`). Python's `json` and
`pandas` read this out of the box; strict/standard JSON parsers do not.

```python
import json
data = json.load(open("Examples/sample_analysis_data.json"))
data["meta"]          # study-level config
data["trials"]        # list of per-trial records (time series + scalars)
data["tables"]        # flat summary tables, ready for pandas

import pandas as pd
per_trial = pd.DataFrame(data["tables"]["per_trial"])
per_participant = pd.DataFrame(data["tables"]["per_participant"])
```

## Top-level structure

```
{
  "meta":   { … study configuration … },
  "trials": [ { …one record per trial… }, … ],
  "tables": { "per_trial": [ … ], "per_participant": [ … ] }
}
```

### `meta`

| Key | Type | Meaning |
|-----|------|---------|
| `instruments` | `[str]` | Tracked tips, in order: `Bipolar`, `Cavitron`, `Scissors`. |
| `pairs` | `[[str,str]]` | Instrument pairs analysed together. |
| `pair_keys` | `[str]` | `"A-B"` string key for each pair (used in `dist`/`angle`). |
| `participants` | `[str]` | Participant ids in the study. |
| `smooth_window` | `int` | Frame width of the smoothing kernel. |
| `idle_speed` | `float` | mm/s below which a tip counts as idle. |
| `voxel_size` | `float` | mm edge of the occupancy-grid voxel. |
| `plane_dist_max` | `float` | mm gate: paired instruments farther apart are "not in use". |
| `time_unit` | `str` | Unit of `duration` (`"s"`). |

### `trials[i]` — one trial

Every list below has length `n_frames` and is aligned to the same frame index;
`tnorm[k]` is the normalised time of frame `k`. Per-instrument entries are keyed
by instrument name; per-pair entries by `pair_key`.

**Identity / sizing**

| Key | Type | Meaning |
|-----|------|---------|
| `participant`, `trial` | `str`, `int` | Who / which trial (1-based). |
| `name` | `str` | Source recording filename (no extension). |
| `label` | `str` | Short plot label, e.g. `"Bianca·T1"`. |
| `n_frames` | `int` | Number of frames. |
| `duration` | `float` | Trial length in `meta.time_unit`. |
| `rmse` | `float` | Fiducial-registration RMSE (mm); ~0 = good fit. |
| `present_instruments` | `[str]` | Instruments actually tracked in this trial. |

**Per-frame time series** (length `n_frames`)

| Key | Shape | Meaning |
|-----|-------|---------|
| `tnorm` | `[float]` | Time normalised to `[0, 1]`. |
| `fmag` | `[float]` | Force magnitude `√(fx²+fy²+fz²)` (N). Force is never masked. |
| `dFdt` | `[float]` | Absolute rate of force change (N/s). |
| `kin[inst].velocity` | `[float]` | Tip speed (mm/s). |
| `kin[inst].acceleration` | `[float]` | Tip acceleration (mm/s²). |
| `kin[inst].jerk` | `[float]` | Tip jerk (mm/s³). |
| `pos_plot[inst]` | `[[x,y,z]]` | Registered tip position (mm) for 3-D trajectory plots. |
| `ang_speed[inst]` | `[float]` | Long-axis angular speed (deg/s). |
| `ang_accel[inst]` | `[float]` | Long-axis angular acceleration (deg/s²). |
| `tracking_status[inst]` | `[int]` | `1` tracked, `0` missing (frozen pose). |
| `dist[pair]` | `[float]` | Inter-instrument tip distance (mm). |
| `angle[pair]` | `[float]` | Angle between long axes (deg). |

> **NaN convention.** Tracking-derived signals (`kin.*`, `pos_plot`,
> `ang_speed`, `ang_accel`, `dist`, `angle`) are `NaN` on frames where the
> instrument is untracked or **not in use** (a paired tip more than
> `plane_dist_max` from the Bipolar). Smoothing also leaves a few `NaN` at the
> edges. `fmag`/`dFdt` come from a separate sensor and are never masked.

**Per-trial scalars** (force)

| Key | Meaning |
|-----|---------|
| `force_mean`, `force_peak` | Mean / peak force magnitude (N). |
| `force_cov` | Force coefficient of variation `SD/mean × 100%`. |
| `impulse` | Force impulse `∫|F| dt`. |
| `dFdt_abs_mean` | Mean `|dF/dt|` (N/s). |
| `torque_mean` | Mean torque magnitude. |

**Per-instrument scalars** (each is `{inst: value}`)

| Key | Meaning |
|-----|---------|
| `inuse_frac` | Fraction of frames the tip is in use. |
| `pathlen` | Total path length travelled (mm). |
| `netdisp` | Net displacement start→end (mm). |
| `straightness` | `netdisp / pathlen` (economy of motion). |
| `bbox_vol` | Axis-aligned bounding-box working volume (mm³). |
| `idle_frac` | Fraction of in-use frames below `idle_speed`. |
| `ell_vol` | Covariance-ellipsoid volume `(4/3)π√(λ₁λ₂λ₃)` (mm³). |
| `aniso` | Anisotropy `λ₁/λ₃` (line/plane vs. isotropic spread). |
| `occ_eff` | Effective occupied voxels `exp(H)` (spatial concentration). |
| `occ_voxels` | Count of distinct occupied voxels. |

### `tables`

Flat, denormalised rows ready for `pandas.DataFrame` and CSV export — the same
data as `trials[*]` scalars, plus per-instrument aggregates of the time series
(`*_speed_mean`, `*_speed_peak`, `*_accel_mean`, `*_jerk_mean`,
`*_angspeed_dps`, `*_angaccel_dps2`) and per-pair distance/angle stats
(`dist_<pair>_{mean,min,max}`, `angle_<pair>_mean`). Column suffixes carry the
unit (`_mm`, `_N`, `_pct`, `_dps`, `_mm3`, …).

* `per_trial` — one row per trial (matches `data/metrics_per_trial.csv`).
* `per_participant` — one row per participant, averaged across their trials
  (matches `data/metrics_per_participant.csv`).
