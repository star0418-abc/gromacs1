# 10_ANALYSIS_PAIRS.md — RDF/CN Pair Definitions & Output Contract

## A) Analysis intent must be config-driven
- Read pair definitions from:
  IN/systems/<SYSTEM_ID>/gromacs/analysis_pairs.yaml

## B) Minimal schema (recommended)
pairs:
  - name: "Li-TFSI"
    A: "group Li"                  # ndx group or gmx selection
    B: "group TFSI"
    rdf:
      rmax: 1.2
      dr: 0.002
    cn:
      cutoff: "first_minimum_smoothed"      # or numeric
      search_rmax: 1.2

## C) CN cutoff rules (robust)
- Default robust mode: first_minimum_smoothed (first_minimum accepted as alias).
- Procedure requirements:
  - smooth RDF with explicit window settings,
  - detect first peak with minimum prominence/height threshold,
  - search minimum after peak with separation/depth constraints (plus optional local basin check),
  - do **not** use raw noisy argmin refinement by default (`refine_points` is explicit opt-in only).
- Recovery/rise tests are advisory by default:
  - they add reliability flags/warnings (`no_recovery_after_minimum`) but do not block candidate selection.
- Shoulder handling:
  - no global-minimum fallback by default,
  - if no first local minimum exists after the first peak, return `no_first_minimum` and mark `shoulder_no_minimum`.
- If no reliable minimum:
  - strict mode -> fail and require explicit numeric cutoff, or
  - non-strict mode -> use configured conservative fallback cutoff and mark fallback_used=true.

## D) rmax capping rule (triclinic-safe)
- `safe_rmax` must be based on minimum box altitude (face-to-face repeat distance), not vector norms:
  - `safe_rmax_tpr = 0.5 * min_box_height_tpr - margin`
  - `safe_rmax_traj = 0.5 * min_box_height_traj - margin`
  - `safe_rmax_used = min(...)` when both exist
- Record basis explicitly: `safe_rmax_basis = "min_altitude"`.

## E) Required cutoff audit fields
- cutoff_method
- cutoff_nm
- smoothing_params
- fallback_used
- quality_flags

## F) Reproducibility guidance
- For production studies, prefer explicit numeric CN cutoffs in config.
- If heuristic cutoff is used, keep all cutoff diagnostics in manifest/logs for audit.

## G) Strict-mode semantics
- If `analysis_settings.strict=true`, the stage fails when **any** pair status is `failed` or `skipped`.
- If `strict=false`, stage may finish with `overall_status=completed_with_errors` and per-pair error counts in manifest.

## H) Outputs
- RDF to: OUT_GMX/<RUN_ID>/analysis/rdf/
- CN to:  OUT_GMX/<RUN_ID>/analysis/cn/
- Do not write analysis outputs into IN/.
