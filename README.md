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
- **Summative figures**: per-participant averages and a cross-participant comparison
  of average force / velocity / acceleration / jerk.
- **Per-trial statistics table**: duration, registration RMSE, mean/peak force, and
  per-instrument path length, mean/peak speed, mean acceleration/jerk — displayed and
  exported to `data/metrics_per_trial.csv` and `data/metrics_per_participant.csv`.

`.igs.mha` data files are git-ignored (see `.gitignore`); only the notebook and
docs are tracked.
