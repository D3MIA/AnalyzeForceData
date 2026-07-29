# AnalyzeForceData

Analysis of surgical instrument tracking + force recordings stored as IGSIO/PLUS
sequence metafiles (`*.igs.mha`).

Each file is one trial and contains, per frame:

- `*TipToWorldTransform` — 4×4 pose of each instrument tip (**Bipolar**, **Cavitron**, **Scissors**)
- `Force` — `fx fy fz tx ty tz` (3D force + 3D torque)
- `BipolarCollectedPoint0..3` — 4 fiducial points used to register trials into a common frame
- `Timestamp` — frame time in seconds

## Data layout

One subfolder per participant; each participant runs several trials:

```
data/
  P01/  trial1.igs.mha  trial2.igs.mha  trial3.igs.mha
  P02/  trial1.igs.mha  trial2.igs.mha  trial3.igs.mha
  ...
```

The subfolder name is used as the participant id, and trials are numbered by
filename order within each folder. (Files placed directly in `data/` still work —
they are grouped under a single participant named after the folder.)

## Usage

1. Put each participant's `.igs.mha` files in `data/<participant>/`.
2. Open [`surgical_force_analysis.ipynb`](surgical_force_analysis.ipynb) in
   **Google Colab** (or Jupyter) and run all cells. In Colab you can upload a zip of
   the `data/` tree from the notebook (see §2) or mount Google Drive.

## What the notebook produces

- **Point-wise rigid registration** of every trial (all participants) onto one common
  reference frame using the 4 fiducial points (SVD / Kabsch).
- **Time normalized to [0, 1]** per trial.
- **Force magnitude** `√(fx²+fy²+fz²)` — one panel per participant.
- **Velocity, acceleration, jerk** per instrument per trial — one figure per participant.
- **3D instrument trajectories** colored by normalized time — one figure per participant.
- **Summative figures**: per-participant averages and a cross-participant comparison
  of average force / velocity / acceleration / jerk.

`.igs.mha` data files are git-ignored (see `.gitignore`); only the notebook and
docs are tracked.
