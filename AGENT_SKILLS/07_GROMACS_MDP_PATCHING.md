# 07_GROMACS_MDP_PATCHING.md — MDP Templates & Dynamic Patching

## A) Templates live in IN, execution uses OUT copies
- IN/systems/<SYSTEM_ID>/gromacs/mdp/*.mdp are templates.
- Copy to OUT_GMX/<RUN_ID>/03_gromacs/<stage>/<stage>.mdp before grompp.
- Never patch IN templates.

## B) CLI must override mdp defaults
- CLI temperature T must override ref_t (and gen_temp if used).
- CLI pressure P must override ref_p.
- Record final values in manifest.

## C) Thermostat grouping guardrails (GPE-safe defaults)
- Default mode is auto:
  - if polymer exists and both Polymer/NonPolymer groups are valid and not tiny -> split mode,
  - otherwise fallback to tc-grps=System with warning.
- Safety criteria for split:
  - exact ndx groups Polymer and NonPolymer exist,
  - each group passes absolute atom threshold (tc_grps_min_atoms),
  - each group passes atom-fraction threshold (tc_grps_min_fraction).

## D) tc-grps vectors and alignment
- If split is active, tau_t and ref_t must be vectors aligned with tc-grps.
- Different tau_t values per group are supported and recommended for heterogeneous polymer/solvent systems.
- If safe vector alignment cannot be guaranteed in auto mode, pipeline must fallback to System and warn.

## E) ndx group naming contract
- Auto split requires exact group names: Polymer and NonPolymer.
- Pipeline must validate ndx presence before grompp whenever custom groups are requested.
- Missing/ambiguous ndx must fail fast with actionable error.

## F) Manifest audit requirements
- Record thermostat_type (tcoupl), tc-grps, tau_t, ref_t, and per-group sizes.
- Record the mode decision (auto/split/system), including fallback reason/warnings.

## G) Prevent config drift
- Fail if CLI T/P disagrees with final mdp (post-patch).
