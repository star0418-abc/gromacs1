# 02_CLI_AND_STAGES.md — CLI, Stages, Sentinels, Dependencies

## A) Required CLI flags
- --ff {GAFF2, OPLS-AA}
- --charge {RESP, CM5}
- --system <SYSTEM_ID>
- --run-id <RUN_ID> (auto-generate only with explicit allow flag)
- --stage {packmol, htpolynet, sanitizer, gmx_em, gmx_eq, gmx_prod, analysis}
- --resume (default true)
- --force (default false)

## B) Resource flags
- --gmx-nt <int>           (mdrun -nt)
- --gmx-gpu-id <string>    (mdrun -gpu_id or equivalent)
- --omp-num-threads <int>  (set OMP_NUM_THREADS)

## C) Stage definitions & expected outputs
1) packmol:
- OUT: OUT_GMX/<RUN_ID>/01_packmol/...
- PUBLISH: IN/systems/<SYSTEM_ID>/packmol/initial.{gro,pdb}

2) htpolynet:
- OUT: OUT_GMX/<RUN_ID>/02_htpolynet/{system_config.yaml,workdir/}
- PUBLISH: IN/systems/<SYSTEM_ID>/htpolynet/system_<RUN_ID>.{gro,itp,top} + current pointers
- STAGING: copy/symlink to IN/systems/<SYSTEM_ID>/gromacs/{system.gro,system.itp,system.top}

3) sanitizer:
- OUT: OUT_GMX/<RUN_ID>/02_5_sanitizer/
- PUBLISH: IN/systems/<SYSTEM_ID>/gromacs/combined_atomtypes_<RUN_ID>.itp + itp_sanitized_<RUN_ID>/
- MUST run before any grompp stage

4) gmx_em:
- OUT: OUT_GMX/<RUN_ID>/03_gromacs/em/
- Sentinel: done.ok

5) gmx_eq:
- OUT: OUT_GMX/<RUN_ID>/03_gromacs/nvt/ and npt/
- Sentinel: done.ok for each substage and a combined done.ok

6) gmx_prod:
- OUT: OUT_GMX/<RUN_ID>/03_gromacs/md/
- Must resume from latest checkpoint if exists

7) analysis:
- OUT: OUT_GMX/<RUN_ID>/analysis/{rdf,cn}/
- Reads analysis_pairs.yaml from IN only

## D) Skip policy
- If --resume true and stage has done.ok and --force false: skip.
- If --force true: rerun stage (do not destroy prior outputs unless explicitly requested).

## E) Dependency policy
- If user requests a later stage, verify prerequisites exist:
  - For sanitizer: staged gromacs inputs from htpolynet must exist.
  - For gmx_* stages: sanitizer outputs must exist (unless explicit unsafe override is enabled).
  - For analysis: trajectory exists (or error).
