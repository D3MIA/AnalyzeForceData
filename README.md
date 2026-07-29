# AnalyzeForceData

Analysis of surgical instrument tracking + force recordings stored as IGSIO/PLUS
sequence metafiles (`*.igs.mha`).

Each trial is one `.igs.mha` file and contains, per frame:

- `*TipToWorldTransform` — 4×4 pose of each instrument tip (**Bipolar**, **Cavitron**, **Scissors**)
- `Force` — `fx fy fz tx ty tz` (3D force + 3D torque)
- `BipolarCollectedPoint0..3` — 4 fiducial points used to register trials into a common frame
- `Timestamp` — frame time in seconds

## Usage

1. Put every trial's `.igs.mha` file in the [`data/`](data/) folder.
2. Open [`surgical_force_analysis.ipynb`](surgical_force_analysis.ipynb) in
   **Google Colab** (or Jupyter) and run all cells. In Colab you can also upload
   the files directly from the notebook (see §2) or mount Google Drive.

## What the notebook produces

- **Point-wise rigid registration** of every trial onto a common reference frame
  using the 4 fiducial points (SVD / Kabsch), so all reported tracking shares one frame.
- **Time normalized to [0, 1]** per trial.
- **Force magnitude** `√(fx²+fy²+fz²)` per trial.
- **Velocity, acceleration, jerk** per instrument per trial.
- **3D instrument trajectories** as lines colored by normalized time.
- **Summative figures**: average force per trial and average velocity /
  acceleration / jerk per instrument.

`.igs.mha` data files are git-ignored (see `.gitignore`); only the notebook and
docs are tracked.
