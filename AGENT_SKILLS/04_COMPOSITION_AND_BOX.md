# 04_COMPOSITION_AND_BOX.md — wt% Composition, Rounding, Box Sizing

## A) wt% is allowed
- Composition may be defined as wt%.
- For polydisperse components (OEGMA, PEGDMA), simulation must use a chosen representative MW/structure.

## B) wt% -> integer N_i (recommended algorithm)
- Compute relative mole ratios from wt% and MW_i.
- Anchor one component (e.g., PC) with N_PC target and scale others proportionally.
- Enforce minimum counts for key species (configurable). If below threshold, scale up.
- Recompute actual wt% after rounding; record deviation in manifest. If deviation > threshold, scale up and retry.

## C) Box sizing (PACKMOL)
- Compute cubic box length L from total mass and initial density rho0.
- Default rho0 in [0.95, 1.05] g/cm^3.
- Default retry policy is fail-fast on PACKMOL failure.
- Density reduction is NOT a default remedy.

## D) Optional density reduction (unsafe/explicit)
- Only enabled with explicit allow flag: --allow-density-reduction.
- Must never go below configured rho_min.
- Must enforce hard caps for total density reduction and implied total box expansion.
- Cap/floor decisions must be written to manifest retry metadata.

## E) Template-fit diagnostic contract
- On PACKMOL failure, report per-template span against intended box length and margin.
- If polymer-like template span cannot fit in the intended box, fail fast and instruct user to provide compact/pre-relaxed conformer.
- If PDB->GRO conversion needs box expansion above configured small-fraction cap, fail fast (do not silently accept density collapse).

## F) Mandatory manifest fields
- target_density_g_cm3
- achieved_density_g_cm3
- achieved_over_target_ratio
- box_length_target_nm
- box_length_final_nm
- box_expand_fraction
- density_reduction_used (bool)
- density_reduction_reason
