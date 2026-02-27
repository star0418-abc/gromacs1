# 08_RESUME_AND_CHECKPOINTS.md — Resume/Skip & Checkpoint Resolution

## A) Sentinels
- Each stage/substage writes done.ok on success.
- Skip when done.ok exists and --resume true and --force false.

## B) Checkpoint logic
1) Same-stage resume (production continuing):
- Use: mdrun -cpi latest.cpt -append
- Do not regenerate .tpr unless explicitly requested.

2) Cross-stage transition (NPT → PROD):
- Use grompp with:
  - -c npt.gro
  - -t npt.cpt
  then mdrun.

## C) Auto-detect latest checkpoint
- Prefer checkpoints within the same RUN_ID:
  - For gmx_prod: check OUT_GMX/<RUN_ID>/03_gromacs/md/*.cpt first
  - else fall back to OUT_GMX/<RUN_ID>/03_gromacs/npt/*.cpt
- Record chosen checkpoint path in manifest.
