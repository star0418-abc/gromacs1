# 12_MANIFEST_SCHEMA.md — manifest.json Required Fields

## A) Identity
- run_id, system_id, timestamp
- selected ff, charge
- selected stage and substage status

## B) Tool versions
- packmol version
- gmx version
- htpolynet version

## C) Inputs provenance
- For every IN asset used:
  - absolute path
  - hash (sha256 recommended)
- Record chosen rules file path

## D) Composition & box
- N_i per component
- MW_i used (especially OEGMA/PEGDMA representative choice)
- rho0 and box length L
- actual wt% after rounding and deviation

## E) Sanitizer outputs
- combined_atomtypes path (versioned + current)
- sanitized itps directory (versioned + current)
- conflict checks results
- allow-override flag if used

## F) Charge neutrality
- computed total Q before fix
- correction applied (yes/no)
- patched itp path if applied

## G) MDP patch
- template mdp path
- patched mdp path
- final T/P values in mdp
- tc-grps mode + ndx used

## H) Resources
- gmx-nt, gmx-gpu-id, omp-num-threads
- how applied to htpolynet internal relax (config vs env/wrapper)

## I) Resume/checkpoints
- selected checkpoint paths and method (-cpi/-t)
- stages skipped due to resume

## J) Analysis
- analysis_pairs.yaml path
- pair list, rdf params, cn cutoffs used

## K) Commands executed
- full commands for each stage (with args)
- working directories

## L) Errors
- if failure: pointer to crash_report.txt
