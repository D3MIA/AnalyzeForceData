# AnalyzeForceData

Analysis of surgical instrument tracking + force recordings stored as IGSIO/PLUS
sequence metafiles (`*.igs.mha`).

Each file is one trial and contains, per frame:

- `*TipToWorldTransform` — 4×4 pose of each instrument tip (**Bipolar**, **Cavitron**, **Scissors**)
- `Force` — `fx fy fz tx ty tz` (3D force + 3D torque)
- `BipolarCollectedPoint0..3` — 4 fiducial points used to register trials into a common frame
- `Timestamp` — frame time in seconds

## Data layout

Package the whole study as a **single zip** whose contents are one subfolder per
participant; each participant runs several trials:

```
study.zip
 └─ P01/  trial1.igs.mha  trial2.igs.mha  trial3.igs.mha
 └─ P02/  trial1.igs.mha  trial2.igs.mha  trial3.igs.mha
 └─ ...
```

The subfolder name is used as the participant id, and trials are numbered by
filename order within each folder. An extra wrapper folder inside the zip (e.g.
everything under a top-level `data/`) is fine — participants are detected by the
folder that directly contains the `.igs.mha` files, at any depth.

## Usage

1. Zip your participant folders into a single `.zip` (default path `data/study.zip`).
2. Open [`surgical_force_analysis.ipynb`](surgical_force_analysis.ipynb) in
   **Google Colab** (or Jupyter), set `DATA_ZIP`, and run all cells. In Colab you can
   upload the zip from the notebook (see §2) or mount Google Drive. The zip is
   extracted into `data_unzipped/` on each run.

## What the notebook produces

- **Point-wise rigid registration** of every trial (all participants) onto one common
  reference frame using the 4 fiducial points (SVD / Kabsch).
- **Time normalized to [0, 1]** per trial.
- **Force magnitude** `√(fx²+fy²+fz²)` — one panel per participant.
- **Velocity, acceleration, jerk** per instrument per trial — one figure per participant.
- **3D instrument trajectories** colored by normalized time — one figure per participant.
- **Path length** traveled per instrument.
- **Inter-instrument distance** for the pairs used together (Bipolar–Cavitron,
  Bipolar–Scissors), both along-trial and summative. Distances above 100 mm are set to
  NaN (instruments not interacting) and excluded.
- **Instrument orientation**: each instrument's long axis is the pivot (shaft) direction
  `normalize(−Rᵀd)` from its `*TipTo*Transform`; the angle between Bipolar–Cavitron and
  Bipolar–Scissors long axes is reported along-trial and summatively, plus per-instrument
  angular speed.
- **Summative figures**: cross-participant comparison of average force / velocity /
  acceleration / jerk / path length / straightness.
- **Missing-tracking detection**: a per-instrument `tracking_status` array (0 = missing,
  1 = tracked) flagging frames whose 4×4 pose is frozen (identical to the previous
  frame), with **% tracked per instrument** reported in the summary and table.
- **In-use masking**: every tracking-derived signal (velocity, acceleration, jerk,
  angular speed, 3D trajectory, inter-instrument distance/angle, path length, etc.) is
  set to NaN — and hidden from plots and averages — on frames where the instrument is
  untracked, or where Cavitron/Scissors is more than 100 mm from Bipolar (treated as not
  in use). Force is from a separate sensor and is not masked. `% in use` is reported per
  instrument.
- **Additional metrics** (suggested surgical-dexterity indicators): net displacement &
  straightness (economy of motion), working volume, idle fraction, force impulse,
  **force coefficient of variation** (`SD/mean × 100%` — normalised force variability;
  lower = steadier control) and **mean |dF/dt|** (N/s — mean absolute rate of force
  change; lower = more controlled, gradual force application, distinct from impulse).
- **Instrument-use localization** (how tightly clustered an instrument's use is): alongside
  the axis-aligned bounding-box working volume (kept as-is), two orientation-invariant
  measures per instrument — the **covariance-ellipsoid volume** `(4/3)π·√(λ₁λ₂λ₃)` from the
  3×3 position covariance (overall spread, frame-independent), its **anisotropy** `λ₁/λ₃`
  (motion confined to a line/plane vs. an isotropic blob), and a dwell-time-weighted
  **occupancy entropy** reported as the effective number of occupied 2 mm voxels `exp(H)`
  (spatial concentration — low means the tip revisits a small core). Computed over in-use
  frames; smaller ellipsoid volume / occupancy ⇒ more localized use.
- **Per-trial statistics table** with every metric above — displayed and exported to
  `data/metrics_per_trial.csv` and `data/metrics_per_participant.csv`.

`.igs.mha` data files are git-ignored (see `.gitignore`); only the notebook and
docs are tracked.
