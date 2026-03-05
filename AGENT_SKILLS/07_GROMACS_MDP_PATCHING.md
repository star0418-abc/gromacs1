# 07_GROMACS_MDP_PATCHING.md — MDP Templates & Dynamic Patching

## A) Templates live in IN, execution uses OUT copies
- `IN/systems/<SYSTEM_ID>/gromacs/mdp/*.mdp` are templates
- Copy to `OUT_GMX/<RUN_ID>/03_gromacs/<stage>/<stage>.mdp` before `grompp`
- Never patch IN templates

## B) CLI must override mdp defaults
- CLI temperature `T` must override `ref_t` (and `gen_temp` if used)
- CLI pressure `P` must override `ref_p`
- Record final values in the manifest

## C) Thermostat grouping guardrails (`tc-grps`)
Default mode is `auto`.

### Auto decision rule
- if polymer exists and both `Polymer` / `NonPolymer` groups are valid and not tiny -> split thermostat mode
- otherwise fallback to `tc-grps = System` with warning

### Safety criteria for split
- exact ndx groups `Polymer` and `NonPolymer` exist
- each group passes absolute atom threshold: `tc_grps_min_atoms`
- each group passes atom-fraction threshold: `tc_grps_min_fraction`

If safe vector alignment cannot be guaranteed, fallback to `System`.

## D) `tc-grps` vectors and alignment
If split thermostat mode is active:
- `tau_t` and `ref_t` must be vectors aligned with `tc-grps`
- different `tau_t` values per group are allowed
- post-patch validation must confirm exact length alignment

If alignment is unsafe or ambiguous:
- fallback to `tc-grps = System`
- log the reason loudly

## E) COM removal guardrails (`comm-mode` / `comm-grps`)
`tc-grps` and `comm-grps` are different controls and must **not** be coupled mechanically.

### Default rule
Use:
- `comm-mode = Linear`
- `comm-grps = System`

This is the default for all pipeline runs unless an explicit policy says otherwise.

### Network / gel / percolated-system safety rule
If the system is any of the following:
- crosslinked polymer network
- gel-like percolated system
- solid-like connected framework
- a system where the polymer phase spans PBC or is expected to behave as one connected body

then the pipeline must **force**:
- `comm-mode = Linear`
- `comm-grps = System`

Do **not** set:
- `comm-grps = Polymer`
- `comm-grps = Polymer NonPolymer`
- or any split COM-removal grouping for such systems

Reason:
- split COM removal in connected / percolated systems can accumulate momentum imbalance and trigger flying-ice-cube-like artifacts

## F) Detection rule for network-sensitive systems
The pipeline should treat the system as network-sensitive if any of the following is true:
- HTPOLYNET / crosslinking stage is active or completed
- system metadata marks the system as gel / network / solid-like
- topology or build metadata indicates polymer connectivity across the box
- the pipeline cannot confidently disprove a network-like state

If uncertain, choose the safer default:
- `comm-grps = System`

## G) Unsafe override policy for COM removal
If a user explicitly requests non-`System` `comm-grps` in a network-sensitive system:
- hard-fail by default

Only allow with explicit unsafe override, for example:
- `--allow-unsafe-comm-grps`
- `--unsafe-comm-grps-reason "<why>"`

Override behavior is mandatory:
- loud warning
- manifest entry with:
  - requested `comm-grps`
  - detected network-sensitive status
  - override reason
  - explicit statement that the run is outside the pipeline’s safe default regime

## H) ndx group naming contract
Auto split thermostat mode requires exact group names:
- `Polymer`
- `NonPolymer`

The pipeline must validate ndx presence before `grompp` whenever custom groups are requested.

Missing or ambiguous ndx groups must fail fast with actionable diagnostics.

## I) Manifest audit requirements
Record at minimum:
- `tcoupl`
- `tc-grps`
- `tau_t`
- `ref_t`
- per-group sizes
- thermostat mode decision (`auto` / `split` / `system`)
- thermostat fallback reason
- `comm-mode`
- `comm-grps`
- whether network-sensitive protection forced `System`
- whether unsafe override was used

## J) Prevent config drift
Post-patch validation must fail if:
- CLI `T` disagrees with final mdp
- CLI `P` disagrees with final mdp
- final `tc-grps` vectors are misaligned
- final `comm-grps` violates network-sensitive safety policy without explicit unsafe override