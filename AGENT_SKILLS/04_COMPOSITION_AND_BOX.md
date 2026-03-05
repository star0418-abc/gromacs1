# 04_COMPOSITION_AND_BOX.md — wt% Composition, Rounding, Box Sizing

## A) wt% is allowed
- Composition may be defined in wt%
- For polydisperse components (for example OEGMA, PEGDMA), simulation must use a chosen representative MW / representative structure
- The chosen representation must be recorded in the manifest

## B) wt% -> integer `N_i` (recommended algorithm)
- Compute relative mole ratios from wt% and `MW_i`
- Anchor one component with a target count (for example solvent anchor)
- Scale all other components proportionally
- Enforce minimum counts for key species if configured
- Round to integer counts deterministically
- Recompute actual wt% after rounding
- Record per-component deviation in the manifest
- If deviation exceeds threshold, scale up and retry

## C) Distinguish final target density from PACKMOL fill density
The pipeline must distinguish:

1. **final physical target density**
   - the density expected after relaxation / crosslinking / NPT equilibration

2. **initial PACKMOL fill density**
   - the density used to generate the pre-reaction / pre-relaxation initial box

These two densities are not always equal.

## D) Polymerization / crosslinking shrinkage contract
For polymerizing or crosslinking systems, the initial PACKMOL box must optionally account for expected volume shrinkage.

Introduce:
- `polymerization_shrinkage_vol_frac`

Definition:
- `polymerization_shrinkage_vol_frac = (V_initial - V_final) / V_initial`

If mass is conserved, then:
- `rho_packmol = rho_target * (1 - polymerization_shrinkage_vol_frac)`

Where:
- `rho_target` = final intended physical density
- `rho_packmol` = initial density used for PACKMOL packing

### Rules
- For clearly non-reactive / non-crosslinking systems, default `polymerization_shrinkage_vol_frac = 0.0`
- For HTPOLYNET / crosslinking workflows, the value should come from:
  - system config, or
  - documented repo default, or
  - explicit CLI / context override
- It must never be guessed silently

## E) Box sizing (PACKMOL)
- Compute cubic box length `L` from total mass and `rho_packmol`
- If `polymerization_shrinkage_vol_frac = 0`, then `rho_packmol = rho_target`
- If shrinkage is nonzero, use the shrinkage-adjusted `rho_packmol`
- PACKMOL should target the **initial** box, not the already-shrunk final box

This avoids:
- unrealistic initial overpacking
- excessive steric clashes during packing
- downstream LINCS / NPT instability caused by starting too dense

## F) Density reduction is NOT the same as shrinkage allowance
Do not confuse:
- **planned shrinkage-aware initial packing**, with
- **unsafe density reduction fallback after PACKMOL failure**

These are different mechanisms.

### Shrinkage-aware packing
- normal, physics-motivated workflow
- not an unsafe override by itself

### Density reduction fallback
- only enabled with explicit allow flag: `--allow-density-reduction`
- must never go below configured `rho_min`
- must enforce hard caps for:
  - total density reduction
  - implied box expansion
- all retry steps and caps must be logged

## G) Template-fit diagnostic contract
On PACKMOL failure, report:
- intended box length
- per-template span
- minimum margin to box faces
- whether any template is intrinsically too extended for the requested box

Rules:
- if a polymer-like template cannot fit, fail fast and instruct user to provide a more compact / pre-relaxed conformer
- if conversion or fitting would require box expansion beyond configured small-fraction cap, fail fast
- do not silently accept catastrophic density collapse

## H) Post-crosslink / post-build relaxation expectation
For shrinkage-aware workflows:
- PACKMOL generates the larger initial box
- HTPOLYNET / build stages create the connected network
- downstream NPT relaxation is responsible for compressing toward the final physical density

The barostat, not PACKMOL, is the default mechanism for the final approach to `rho_target`.

## I) Mandatory manifest fields
Record at minimum:
- `target_density_g_cm3`
- `packmol_density_g_cm3`
- `polymerization_shrinkage_vol_frac`
- `achieved_density_g_cm3`
- `achieved_over_target_ratio`
- `box_length_target_nm`
- `box_length_packmol_nm`
- `box_length_final_nm`
- `box_expand_fraction`
- `density_reduction_used`
- `density_reduction_reason`
- `count_rounding_summary`
- `wtpct_deviation_summary`