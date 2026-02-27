# 00_CONTRACT.md — Pipeline Contract & Acceptance

## A) Hard Contract (MUST)
1) Inputs-only-from-IN:
- PACKMOL/HTPOLYNET/GROMACS/analysis must not read inputs outside `IN/` (except tool binaries).

2) Outputs-only-to-OUT:
- All `.tpr/.xtc/.trr/.edr/.log/.cpt` and analysis outputs must be written under `OUT_GMX/<RUN_ID>/`.

3) Publish is explicit + atomic:
- Only reusable assets may be published back into `IN/` (e.g., `.gro/.pdb/.itp/.top/.mdp/.ndx/.yaml`).
- Publish only on success.
- Publish must be atomic: write temp → rename.

4) HTPOLYNET sandbox:
- HTPOLYNET must run with `workdir = OUT_GMX/<RUN_ID>/02_htpolynet/workdir/`.
- No intermediate outputs may land in repo root or `IN/`.

5) Fail fast:
- Ambiguous asset selection (multiple candidates) must stop with a clear error.

## B) Acceptance Checklist (MUST PASS)
1) A full run creates:
- `OUT_GMX/<RUN_ID>/manifest.json`
- stage directories exist and contain expected outputs

2) No runtime outputs under `IN/`:
- `IN/` must not contain `.tpr/.xtc/.edr/.cpt/.log` etc.

3) HTPOLYNET intermediates are confined:
- `OUT_GMX/<RUN_ID>/02_htpolynet/workdir/` contains htpolynet intermediate files
- `IN/` contains only published final `system.{gro,itp,top}` outputs

4) Sanitizer correctness:
- `IN/systems/<SYSTEM_ID>/gromacs/combined_atomtypes_current.itp` exists
- `IN/systems/<SYSTEM_ID>/gromacs/itp_sanitized_current/` exists
- `system.top` includes combined atomtypes first, then sanitized itps

5) MDP patching correctness:
- Final mdp used in OUT reflects CLI T/P
- tc-grps mode matches config, ndx exists if split mode

6) Resume/skip:
- `--resume` skips stages with `done.ok`
- `--stage gmx_prod` can resume from latest `.cpt` if present

7) Failure behavior:
- Failure prints last 50 log lines, absolute workdir, and command
- `crash_report.txt` exists in failing stage directory
- No cleanup is performed on failure
