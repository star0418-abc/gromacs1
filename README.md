# MD Simulation Pipeline

**PACKMOL 鈫?HTPOLYNET 鈫?SANITIZER 鈫?GROMACS 鈫?RDF/CN**

Automated pipeline for molecular dynamics simulations of gel polymer electrolytes.

## Features

- 馃敀 **Strict path isolation**: Inputs from `IN/` only, outputs to `OUT_GMX/` only
- 馃摝 **Manifest tracking**: Full reproducibility with SHA256 hashes and command logs
- 鈴革笍 **Resume support**: Sentinel-based stage completion tracking with checkpoint resume
- 馃敡 **Modular stages**: Easy to extend and customize
- 馃И **Composition model**: wt% 鈫?integer molecule count conversion with auto-scaling
- 馃搻 **Box sizing**: Density-based cubic box calculation with retry on overlap
- 馃敆 **HTPolyNet sandboxing**: All intermediates confined to workdir
- 馃摛 **Atomic publishing**: Versioned outputs with current symlinks
- 馃К **ITP Sanitizer**: Atomtype extraction, conflict detection, defaults/comb-rule validation, charge neutrality
- 馃尅锔?**MDP Patching**: Dynamic temperature/pressure/nsteps overrides from CLI
- 鈿欙笍 **Resource Pinning**: Control CPU threads, GPU, and OMP_NUM_THREADS
- 馃挜 **Crash Context**: Detailed failure reports with log excerpts
- 馃摑 **mdrun Log Streaming**: stdout/stderr streamed to stage logs with signal forwarding for graceful checkpoints
- 馃Ь **Stage Fingerprints**: `stage_state.json` captures inputs, hashes, commands, and exact git/GROMACS executable provenance
- 鉁?**Strict Resume/Force**: Resume only on exact fingerprint match; `--force` archives old artifacts
- 馃И **Quality Gates (Dispatcher)**: pipeline success can include thermodynamic sanity checks (temperature/pressure/density) with `off|warn|error` policy
- 馃К **Extended Provenance**: raw tool probes with parsed GROMACS build features, runtime/platform details, git code state, and resolved executable/token metadata
- 馃И **Audit Artifacts**: per-stage `repro.json`, dispatcher `run_report.json`, standalone `provenance.txt`, and SI-ready `computational_details.md`

## Quick Start

```bash
# Run all stages with GAFF2 forcefield (explicit run-id required)
python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage all

# Run with custom temperature and pressure
python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage all -T 350 -P 1.0

# Run GROMACS stages with total threads (mdrun -nt)
python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage gmx_prod \
    --gmx-nt 8

# Run GROMACS stages with OpenMP threads per rank (-ntomp + OMP_NUM_THREADS)
python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage gmx_prod \
    --omp-num-threads 8

# Dry-run (inspect generated inputs without execution)
python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage gmx_em --dry-run

# Force re-run a completed stage (archives prior outputs)
python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage packmol --force

# Skip completed stages (default behavior with resume=True)
python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage gmx_prod

# Allow auto-generated run ID (for quick tests)
python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --allow-auto-run-id --stage all

# Allow velocity reset if checkpoint is missing (explicit opt-in)
python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage gmx_prod --allow-velocity-reset

# Allow grompp warnings explicitly (default maxwarn=0)
python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage gmx_em --grompp-maxwarn 1

# Export SI-ready computational details (partial pipelines supported by default)
python tools/export_computational_details.py --run-root OUT_GMX/<RUN_ID>

# Strict publication mode: require canonical EM/NVT/NPT/MD and custom output path
python tools/export_computational_details.py --run-root OUT_GMX/<RUN_ID> \
    --require-stages em,nvt,npt,md --out computational_details_publication.md
```

## Directory Layout

```
project/
鈹溾攢鈹€ IN/                              # ALL INPUTS (read-only for pipeline)
鈹?  鈹溾攢鈹€ molecules/<NAME>/            # Molecular structure files
鈹?  鈹?  鈹溾攢鈹€ pdb/                     # PDB templates (PACKMOL)
鈹?  鈹?  鈹溾攢鈹€ gro/, itp/, top/, mol2/  # Other formats
鈹?  鈹溾攢鈹€ systems/<SYSTEM_ID>/
鈹?  鈹?  鈹溾攢鈹€ packmol/{pdb,gro}/       # PACKMOL published outputs
鈹?  鈹?  鈹溾攢鈹€ htpolynet/               # HTPolyNet inputs/outputs
鈹?  鈹?  鈹?  鈹溾攢鈹€ gro/, itp/, top/     # Published crosslinked system
鈹?  鈹?  鈹?  鈹斺攢鈹€ rules/{GAFF2,OPLSAA}/ # Reaction rules
鈹?  鈹?  鈹斺攢鈹€ gromacs/                 # GROMACS staged inputs
鈹?  鈹?      鈹溾攢鈹€ gro/, itp/, top/     # Staged from HTPolyNet
鈹?  鈹?      鈹溾攢鈹€ mdp/                 # MDP templates (em, nvt, npt, md)
鈹?  鈹?      鈹溾攢鈹€ ndx/                 # Index files
鈹?  鈹?      鈹溾攢鈹€ combined_atomtypes_*.itp   # Sanitizer output
鈹?  鈹?      鈹溾攢鈹€ itp_sanitized_*/           # Sanitized ITPs
鈹?  鈹?      鈹斺攢鈹€ analysis_pairs.yaml
鈹?  鈹斺攢鈹€ forcefield/{oplsaa,gaff2}/
鈹?  鈹斺攢鈹€ htpolynet/rules/{GAFF2,OPLSAA}/  # Shared HTPolyNet rules (optional)
鈹?
鈹溾攢鈹€ OUT_GMX/<RUN_ID>/                # ALL OUTPUTS (write-only)
鈹?  鈹溾攢鈹€ 01_packmol/                  # PACKMOL output
鈹?  鈹?  鈹溾攢鈹€ packmol.inp, initial.pdb
鈹?  鈹?  鈹斺攢鈹€ done.ok
鈹?  鈹溾攢鈹€ 02_htpolynet/                # HTPolyNet output
鈹?  鈹?  鈹溾攢鈹€ system_config.yaml       # Generated config
鈹?  鈹?  鈹溾攢鈹€ workdir/                 # Sandbox (all intermediates here)
鈹?  鈹?  鈹溾攢鈹€ system.{gro,itp,top}     # Final outputs
鈹?  鈹?  鈹斺攢鈹€ done.ok
鈹?  鈹溾攢鈹€ 02_5_sanitizer/              # ITP Sanitizer logs
鈹?  鈹溾攢鈹€ 03_gromacs/
鈹?  鈹?  鈹溾攢鈹€ em/                      # Energy minimization
鈹?  鈹?  鈹?  鈹溾攢鈹€ em.mdp, em.tpr, em.gro, em.log
鈹?  鈹?  鈹?  鈹溾攢鈹€ stage_state.json
鈹?  鈹?  鈹?  鈹溾攢鈹€ repro.json
鈹?  鈹?  鈹?  鈹斺攢鈹€ done.ok
鈹?  鈹?  鈹溾攢鈹€ nvt/                     # NVT equilibration
鈹?  鈹?  鈹?  鈹溾攢鈹€ nvt.mdp, nvt.tpr, nvt.gro, nvt.cpt
鈹?  鈹?  鈹?  鈹溾攢鈹€ stage_state.json
鈹?  鈹?  鈹?  鈹溾攢鈹€ repro.json
鈹?  鈹?  鈹?  鈹斺攢鈹€ done.ok
鈹?  鈹?  鈹溾攢鈹€ npt/                     # NPT equilibration
鈹?  鈹?  鈹?  鈹溾攢鈹€ npt.mdp, npt.tpr, npt.gro, npt.cpt
鈹?  鈹?  鈹?  鈹溾攢鈹€ stage_state.json
鈹?  鈹?  鈹?  鈹溾攢鈹€ repro.json
鈹?  鈹?  鈹?  鈹斺攢鈹€ done.ok
鈹?  鈹?  鈹斺攢鈹€ md/                      # Production MD
鈹?  鈹?      鈹溾攢鈹€ md.mdp, md.tpr, md.xtc, md.cpt
鈹?  鈹?      鈹溾攢鈹€ logs/mdrun.stdout.log, logs/mdrun.stderr.log
鈹?  鈹?      鈹溾攢鈹€ stage_state.json
鈹?  鈹?      鈹溾攢鈹€ repro.json
鈹?  鈹?      鈹斺攢鈹€ done.ok
鈹?  鈹溾攢鈹€ analysis/{rdf,cn}/           # RDF and coordination numbers
鈹?  鈹溾攢鈹€ logs/
鈹?  鈹斺攢鈹€ manifest.json
鈹?  鈹斺攢鈹€ run_report.json
鈹?  鈹斺攢鈹€ computational_details.md
鈹?
鈹溾攢鈹€ pipeline/                        # Python package
鈹?  鈹溾攢鈹€ mdp_patcher.py               # MDP template patching
鈹?  鈹溾攢鈹€ resource_scheduler.py        # CPU/GPU pinning
鈹?  鈹溾攢鈹€ crash_context.py             # Failure diagnostics
鈹?  鈹斺攢鈹€ stages/
鈹?      鈹溾攢鈹€ packmol.py               # 鉁?Implemented
鈹?      鈹溾攢鈹€ htpolynet.py             # 鉁?Implemented
鈹?      鈹溾攢鈹€ sanitizer_stage.py       # Sanitizer stage orchestrator (run/wiring)
鈹?      鈹溾攢鈹€ topology_sanitizer.py    # Topology/ITP/include sanitization helpers
鈹?      鈹溾攢鈹€ spatial_checker.py       # GRO/PBC/dipole/grompp consistency helpers
鈹?      鈹溾攢鈹€ sanitizer.py             # Backward-compatible facade (re-exports)
鈹?      鈹溾攢鈹€ gromacs.py               # 鉁?Implemented
鈹?      鈹斺攢鈹€ analysis.py              # 鉁?Implemented
鈹?
鈹溾攢鈹€ run_pipeline.py
鈹斺攢鈹€ README.md
```

## Stage Execution Order

```
1. packmol     鈫?Generate initial packed system      鉁?IMPLEMENTED
2. htpolynet   鈫?Crosslink polymer network           鉁?IMPLEMENTED
3. sanitizer   鈫?ITP atomtype sanitization           鉁?IMPLEMENTED
4. gmx_em      鈫?Energy minimization                 鉁?IMPLEMENTED
5. gmx_eq      鈫?NVT + NPT equilibration             鉁?IMPLEMENTED
6. gmx_prod    鈫?Production MD                       鉁?IMPLEMENTED
7. analysis    鈫?RDF/CN post-processing              鉁?IMPLEMENTED
```

## Forcefield/Charge Compatibility

Default-supported parameter-set pairs are enforced at CLI validation time:
- `GAFF2 + CM5`
- `OPLS-AA + RESP`

Unsupported pairs hard-fail by default.  
Unsafe override requires both:
- `--allow-charge-model-mismatch`
- `--charge-model-mismatch-reason \"...\"`

When unsafe override is used, the pipeline emits a loud warning and records:
- `unsafe_override=true`
- `reason`
- `ff`
- `charge_model`
- `validated_parameter_set=false` (explicitly marks unvalidated FF/charge mixing)

## HTPolyNet Stage - Rules, Box, and Outputs

### Rules Resolution (Fail-Fast)
Search order for reaction rules (first match wins):
1) `IN/systems/<SYSTEM_ID>/htpolynet/rules/<FF>/`
2) `IN/htpolynet/rules/<FF>/`
3) `IN/forcefield/<ff>/htpolynet/rules/<FF>/`
4) `IN/forcefield/<ff>/rules/<FF>/`

Each rules directory must include:
- A rules YAML (e.g., `reactions.yaml` or `rules.yaml`)
- A forcefield marker (`forcefield.itp` or `forcefield.top`, optionally under `ff/`)

Missing or ambiguous rules trigger a hard error with actionable guidance.

### Box Specification (Real, Not Comments)
- `.gro` inputs: box vectors are read directly and written into `system_config.yaml`.
- `.pdb` inputs: you must provide a box length via composition (`box_size_nm`) or an explicit context override.
- Triclinic boxes are rejected by default (require re-boxing to orthorhombic/cubic).

### Optional Multi-Step Curing
If `ctx.htpolynet_curing_steps` is provided (list of dicts with `temperature_K`, `duration_ps`, optional `ramp_ps`),
the steps are injected into `system_config.yaml` and recorded in the manifest. Absent 鈫?no behavior change.

### Output Selection and Staging
- Final output directory priority: `systems/final-results` 鈫?`systems/final` 鈫?`output/final` 鈫?other (only if completion marker exists).
- Outputs are accepted only if `system.gro + system.itp + system.top` exist together.
- `index.ndx` is staged to `IN/systems/<SYSTEM_ID>/gromacs/ndx/` if present alongside outputs.

### Placeholder Policy
- If `--allow-placeholder` is used and HTPolyNet is unavailable, placeholder outputs are written to `OUT_GMX` only.
- The stage returns **failure** and does **not** publish/stage placeholders to `IN/`.

## PACKMOL Stage - Reproducibility Policies

The PACKMOL stage enforces strict reproducibility with explicit unit handling, density floor limits, and versioned publishing.

### Unit Policy

| Context | Unit | Rationale |
|---------|------|-----------|
| Internal (pipeline) | **nm** | GROMACS native unit |
| PACKMOL input | **脜** | PACKMOL native unit (auto-converted: nm 脳 10) |
| PDB output (from PACKMOL) | **脜** | PACKMOL writes 脜 |
| GRO output | **nm** | PDB coordinates are scaled to nm in Python; `editconf` runs without `-scale` |

**Fail-fast**: After PACKMOL generates `initial.pdb`, the pipeline validates and infers units using multiple signals (winsorized extent + local nearest-neighbor distance sampling). This reduces false positives from long polymer chains that appear to span PBC in unwrapped templates.

**Template units**:
- `.pdb` templates: unit inference is applied (脜 vs nm) and may require explicit `--packmol-pdb-scale` if signals are ambiguous or conflicting.
- `.gro` templates: treated as **nm** inputs; unit inference is bypassed and `--packmol-pdb-scale` must be `1.0`.

### Retry and Density Floor (Hardened)

Default behavior is fail-fast on PACKMOL failure.
No density reduction is applied unless `--allow-density-reduction` is explicitly set.

**Single source of truth for box**: the intended PACKMOL box is authoritative for density. PDB鈫扜RO conversion preserves the intended box by default; if a minimal expansion is required to fit the structure, the delta is recorded and a warning is emitted.

**Density floor policy (always enforced):**
- `rho0` = initial target density from recipe.
- `rho_effective = rho0 / (1 - polymerization_shrinkage_vol_frac)` is used for initial box sizing when shrinkage is configured.
- `effective_rho_min = max(rho0_min_effective, rho_effective * 0.85, rho_effective * (1 - max_density_reduction_fraction))`
- Default attempt policy: one PACKMOL attempt at `rho_i = rho_effective` (fail-fast)
- If `--allow-density-reduction`: `rho_i = max(rho_effective * (0.95 ** attempt_index), effective_rho_min)`
- The absolute floor (`rho0_min_g_cm3`, mapped to effective density) is never bypassed.
- Total density reduction and implied total box expansion are capped and audited in manifest retry metadata.
- Achieved density is checked after conversion; if below floor, stage fails by default.

**Safeguard B2: Versioned attempt outputs**
- Each attempt writes: `packmol_try{N}.inp` AND `packmol_try{N}.pdb`
- The generated PACKMOL input always uses the explicit `output_file` path for `output ...` (no ambiguous fallback name)
- On success, copies to canonical names (`packmol.inp`, `initial.pdb`)
- Manifest records: `successful_attempt_index`, `successful_attempt_files`
- Failed attempts retain versioned logs: `packmol_try{N}.log`
- PACKMOL seeds are deterministic by default (`seed>0` uses deterministic attempt offsets; `seed=-1` is deterministically derived from context+attempt unless explicit random mode is enabled).
- Explicit random mode keeps `seed=-1`, emits a reproducibility warning, and records the non-deterministic seed mode in retry metadata.
- Retry policy metadata is recorded in `composition.retry_policy`.

**Template-fit diagnostics before density reduction**
- On failure, the stage computes and prints a per-template fit table:
  - `name`, `polymer_like`, `span_nm`, `box_nm`, `margin_nm`, `fit_ok`
- The same table is recorded in manifest under `composition.retry_policy.template_fit`.
- If polymer-like template span is the limiting factor, the stage fails fast with remediation:
  - provide a compact/coiled pre-relaxed polymer conformer, or intentionally increase box via density settings.
- If PDB鈫扜RO conversion requires box expansion above the configured cap, the stage hard-fails (no silent density collapse acceptance).

Manifest records:
- `target_density_g_cm3`, `achieved_density_g_cm3`, `achieved_over_target_ratio`
- `box_length_target_nm`, `box_length_final_nm`, `box_expand_fraction`
- `density_reduction_used`, `density_reduction_reason`

```bash
# Default behavior: single fail-fast attempt at target density:
python run_pipeline.py ... --stage packmol

# Explicit opt-in to allow density-backoff retries (still floor-capped):
python run_pipeline.py ... --stage packmol --allow-density-reduction
```

**Retry diagnostics recorded in manifest** (`retry_attempts[]`):
```json
{
  "attempt_index": 0,
  "rho_target": 1.0,
  "box_length_nm": 5.5,
  "tolerance": 2.0,
  "seed": 74293115,
  "composition_changed": false,
  "success": true,
  "inp_file": "OUT_GMX/.../packmol_try1.inp",
  "pdb_file": "OUT_GMX/.../packmol_try1.pdb"
}
```

**Density retry summary** (`density_retry_summary`):
```json
{
  "total_attempts": 3,
  "successful_attempt": 2,
  "final_density_ratio": 0.9025,
  "final_attempt_target_density_ratio": 0.95,
  "floor_enforced": false,
  "all_failed": false
}
```

**Type hints (RetryAttemptInfo)**

To avoid runtime import cycles while keeping IDE/mypy support, `RetryAttemptInfo` is imported under `TYPE_CHECKING` and referenced via forward references in annotations (e.g., `List["RetryAttemptInfo"]`).

### Strict Publish Mode (Default)

By default, publish failures (writing to `IN/systems/<SYS_ID>/packmol/`) fail the stage.

```bash
# Default strict mode - publish failure = stage failure
python run_pipeline.py ... --stage packmol

# Non-strict mode (legacy) - publish failure is a warning
python run_pipeline.py ... --stage packmol --allow-partial-publish
```

Manifest records `publish_status` and explicit `publish_details` (`pdb_published`, `gro_published`, `both_published`).
`packmol_publish_ok` is only true when **both** PDB and GRO were updated.
If GRO conversion is not publishable, PACKMOL publish is blocked to prevent `new PDB + stale GRO` drift and the stage is marked non-complete.

### Strict GRO Conversion Mode

When `--strict-gro-conversion` is enabled (recommended for production), the pipeline fails immediately if:
- the resolved GROMACS command is missing, or
- `editconf` execution fails / returns nonzero.

```bash
# Strict mode (recommended) - fail if gmx missing
python run_pipeline.py ... --stage packmol --strict-gro-conversion

# Non-strict mode (default) - warn and continue with PDB only
python run_pipeline.py ... --stage packmol --no-strict-gro-conversion
```

**Unit diagnostics recorded in manifest** (`unit_diagnostics`):
```json
{
  "box_length_nm": 5.5,
  "scale_factor_used": 0.1,
  "decision_source": "extent+distance",
  "raw_extent_A": [45.2, 44.8, 45.1],
  "winsor_extent_A": [44.9, 44.4, 44.8],
  "sampled_nn_distance_stats": {
    "count": 1200,
    "min": 0.98,
    "median": 1.42,
    "max": 1.87,
    "units": "raw"
  },
  "extent_nm": 4.52,
  "required_box_nm": 5.5,
  "box_adjusted": false,
  "box_change_reason": "none",
  "translation_nm": [0.02, -0.01, 0.00],
  "margin_nm": 0.2,
  "margin_nm_effective": 0.11,
  "max_margin_fraction": 0.02,
  "epsilon_nm": 0.001,
  "intended_box_nm": 5.5,
  "box_delta_nm": 0.0,
  "box_delta_fraction": 0.0,
  "density_change_expected": false,
  "gro_conversion_status": "success",
  "strict_packmol_units": true
}
```

| `gro_conversion_status` | Meaning |
|-------------------------|---------|
| `success` | GRO file created successfully |
| `skipped_no_gmx` | gmx not found, non-strict mode, continued with PDB |
| `skipped_strict_fail` | strict mode blocked conversion due missing/invalid gmx command |
| `failed` | command executed but conversion failed |

### Versioned Publishing

Before overwriting `IN/systems/<SYS_ID>/packmol/{pdb,gro}/initial.*`:
1. Previous files archived to `_archive/YYYYMMDD-HHMMSS/` (with collision-safe `-NN` suffix if needed)
2. SHA256 hashes written to `current.sha256` via atomic temp-file + `os.replace`
3. `os.replace()` used for atomic overwrites (safe on Windows/Linux)
4. Temp files use collision-safe naming: `.{name}.tmp.{pid}_{uuid8}` (safe for parallel runs)

### Per-System Composition Config

If `IN/systems/<SYS_ID>/composition.yaml` exists, it is used instead of the default recipe:

```yaml
molecules:
  - name: LiTFSI
    mw_g_mol: 287.09
    wt_pct: 15.0
    min_count: 10
    role: salt
  - name: PC
    mw_g_mol: 102.09
    wt_pct: 25.0
    role: solvent
  - name: OEGMA
    mw_g_mol: 475.0
    wt_pct: 20.0
    role: monomer
    mw_note: representative MW for polydisperse material
    pdi: 1.10        # optional, recorded only
  - name: PEGDMA
    mw_g_mol: 750.0
    wt_pct: 8.0
    role: crosslinker
    min_count: 3
    mw_note: representative MW for polydisperse material
target_total_molecules: 500
gel_min_crosslinker_count: 12
gel_policy: warn                 # off | warn | error
auto_scale_to_gel_min_crosslinker: false
rho0_g_cm3: 1.0
rho0_min_g_cm3: 0.95
rho0_max_g_cm3: 1.05
polymerization_shrinkage_vol_frac: 0.06
```

The composition source (path + SHA256) is recorded in manifest for full provenance.

Role-aware fields are optional and backward-compatible:
- `role`: freeform string, with recommended values `monomer`, `crosslinker`, `solvent`, `salt`, `additive`.
- `mw_note` / `pdi`: metadata only (counts are still computed from fixed `mw_g_mol`).

### Polydispersity Notes (Explicit Approximation)

For materials like OEGMA/PEGDMA, fixed MW is an approximation unless you explicitly model multiple bins.
- If `pdi` is provided or `mw_note` indicates representative MW, this metadata is written to manifest.
- If a molecule is marked `role: monomer` or `role: crosslinker`, name contains `OEGMA`/`PEGDMA`, and `mw_note` is missing, a warning is emitted and recorded (`composition.warnings.polydispersity_approximation`).
- No automatic count redistribution by `pdi` is performed in this model.

### Dry-Run Semantics

`--dry-run` generates PACKMOL input but **does not execute PACKMOL**:
- Missing templates 鈫?returns **False** with actionable message (fail-fast)
- Manifest marked `preview_only: true`
- No `done.ok` written

### Graceful Failure Handling

The PACKMOL stage catches and handles errors gracefully instead of crashing with tracebacks:

| Error Condition | Behavior |
|-----------------|----------|
| Missing `composition.yaml` (no `--allow-default-recipe`) | `mark_failed()` + actionable error message + return False |
| Malformed/invalid `composition.yaml` schema | `mark_failed()` + schema/type/value error message + return False |
| Malformed/Empty PDB (no atoms parsed) | `mark_failed()` + error message + return False |
| Publish failure (strict mode) | `mark_failed()` + error message + return False |

All failures record context in the manifest for debugging.

### Polymer Percolation Guardrails

For polymer/gel electrolyte (GPE) systems, insufficient molecule counts can prevent HTPolyNet from achieving network percolation and the resulting structure will be non-physical.

**Default behavior (strict mode):**
- Missing `composition.yaml` 鈫?**HARD FAILURE** with message to create one or use `--allow-default-recipe`
- `target_total_molecules=500` (default) with polymer-like species 鈫?**WARNING** (too low for percolation)
- `total_molecules < min_total_molecules_polymer` (5000) 鈫?**HARD FAILURE**
- `polymer_chains < min_polymer_chains` (100) 鈫?**HARD FAILURE**

**Polymer-likeness detection:**
- MW >= 1000 g/mol 鈫?polymer-like
- Name contains: PEO, PEG, PVDF, PAM, PEGDMA, OEGMA, PMMA, PVA, POLYMER, etc.

**Override flags:**
- `--allow-default-recipe`: Permit demo mode (not for production)
- `--min-polymer-chains N`: Override minimum (default: 100)
- `--min-total-molecules-polymer N`: Override minimum (default: 5000)
- `--no-strict-polymer-check`: Convert failures to warnings
- `--allow-composition-changing-retries`: Allow retries to reduce molecule counts (DANGEROUS)

### Crosslinker Gel Guardrails

In addition to the polymer-chain checks above, composition conversion supports role-based gel checks:
- Set `role: crosslinker` on relevant molecules.
- Set `gel_min_crosslinker_count` to enforce a minimum expected crosslinker count from mole-fraction scaling.
- `gel_policy` controls behavior when expected count is below threshold:
  - `off`: disabled
  - `warn`: emit warning and record `composition.warnings.gel_point_risk` in manifest
  - `error`: fail fast with `ValueError`
- `auto_scale_to_gel_min_crosslinker: true` can increase `scale_factor` deterministically (same buffer style as `min_count`) until the threshold is satisfied.

Warning/error messages include crosslinker name and expected count and explicitly note the non-percolating microgel risk under PBC.

### Unit Handling (脜 vs nm) 鈥?Hardened

PACKMOL outputs coordinates in 脜ngstr枚ms (脜) but GROMACS expects nanometers (nm). The pipeline now combines:
- **winsorized extent signal** (0.5鈥?9.5% per-axis span) to reduce false positives from chain/PBC spanning outliers
- **residue-local nearest-neighbor signal**: atoms are grouped by `(chain_id, resname, resid)` and nearest-neighbor distances are sampled within residues (bounded, deterministic)

Decision policy:
- If extent and distance agree: accept inferred units.
- If only one signal is informative: use that signal.
- If signals disagree or both are inconclusive: **ambiguous**.
  - Ambiguity always fails fast and requires explicit `--packmol-pdb-scale`.
  - `strict_packmol_units=False` does **not** allow ambiguous scale fallback; it only relaxes non-scale checks.

Override for ambiguous cases:
- `--packmol-pdb-scale 0.1`: Force 脜鈫抧m conversion
- `--packmol-pdb-scale 1.0`: Force nm鈫抧m (no scaling)

Extent guard thresholds are configurable via `tolerance_factor`:
- Upper bound multiplier: `1.0 + 0.2 * tolerance_factor`
- Lower bound multiplier: `0.1 * tolerance_factor` (or explicit `min_fraction` override)
- Default (`tolerance_factor=0.5`) remains compatible with prior behavior (`1.10x` upper / `0.05x` lower).

Manifest keys (recorded for reproducibility):
```json
"unit_handling": {
  "raw_extent_A": [120.5, 118.3, 119.8],
  "winsor_extent_A": [118.1, 117.4, 118.0],
  "box_length_nm": 12.0,
  "ratio_extent": 9.84,
  "sampled_nn_distance_stats": { "count": 4000, "min": 0.95, "median": 1.41, "max": 1.90, "units": "raw" },
  "decision_source": "extent+distance",
  "unit_inferred": "A",
  "ambiguity_reason": null,
  "scale_factor_used": 0.1,
  "extent_nm": 12.05,
  "strict_packmol_units": true,
  "required_box_nm": 12.0,
  "box_adjusted": false,
  "box_change_reason": "none",
  "translation_nm": [0.003, -0.002, 0.001],
  "clearance_low_nm": 0.004,
  "clearance_high_nm": 0.004,
  "margin_nm": 0.2,
  "margin_nm_effective": 0.2,
  "max_margin_fraction": 0.02,
  "epsilon_nm": 0.001,
  "intended_box_nm": 12.0,
  "box_delta_nm": 0.0,
  "box_delta_fraction": 0.0,
  "density_change_expected": false,
  "centering_strategy": "center_in_box_via_translate"
}
```

### Boundary Margin and Centering 鈥?Hardened

Atoms placed exactly at box edges can cause PBC wrap clashes during GROMACS preprocessing.

**Safeguard B4: Density-preserving centering**

Fixed margins can dilute density for small systems if applied by expanding the box during conversion.
The pipeline now preserves the intended PACKMOL box and centers the structure within it.

The conversion uses **explicit `-translate`** (equivalent to `gmx editconf -c`) with the intended `-box`:
- Centers the structure in the existing box (no default expansion)
- Applies a tiny epsilon shift (<< 0.01 nm) to avoid exact boundary hits
- Expands only if the structure cannot fit in the intended box (minimal expansion)

**Adaptive margin target (no box inflation):**
- `--packmol-edge-margin NM`: Requested margin in nm (default: 0.2)
- Effective margin is capped: `min(requested, max_margin_fraction * min(box_dim))`
- `max_margin_fraction` defaults to 0.02 (2% of the smallest box dimension)
- Used only as a *clearance target*; it does **not** enlarge the box by default
- Fail-fast: if `box_length_ang - 2*margin` leaves no valid inside-box span, PACKMOL input generation aborts early with computed values (`box_length_ang`, `margin`, `usable_span`).

**Post-run boundary handling:**
- Conversion enforces epsilon-clearance on both sides; if centering cannot satisfy it, box is minimally expanded.
- Auto-expansion is hard-capped by `packmol_max_box_expand_fraction` (default `0.03`).
- If required expansion exceeds the cap, stage fails fast with actionable guidance.
- Manifest records `box_adjusted`, `box_delta_fraction`, `expected_density_drop_fraction`, and final achieved density.

### Density Floor Logic

The retry strategy respects per-system configuration:
- `polymerization_shrinkage_vol_frac` adjusts packing density before PACKMOL:
  - `rho_effective = rho0_g_cm3 / (1 - shrinkage_vol_frac)`
  - valid range: `0 <= shrinkage_vol_frac < 0.3`
- `rho0_min_g_cm3`/`rho0_max_g_cm3` are converted to effective-density candidates for retry logic and manifest audit.
- Effective floor uses the effective-density domain:
  - `effective_min = max(rho0_min_effective, initial_rho_effective * density_floor_fraction)`
- Console logs: `rho0`, `effective_rho_min`, `rho_i` for each attempt
- Final density check uses `achieved_density >= max(rho0_min_effective, density_floor_fraction * target_density)`
- `achieved_density` is computed from the **final GRO box**
- `--allow-density-reduction` only enables density-backoff attempts; it does **not** bypass floor enforcement.

### Optional Preassembly Hints (GPE-focused)

Default behavior remains unchanged: `packmol_preassembly_mode = none`.

Optional modes:
- `li_near_polymer`: bias Li placement toward the polymer-rich region.
- `li_solvent_shell`: bias Li toward solvent-rich region and reduce initial Li-TFSI close contacts.

How it works:
- Adds lightweight PACKMOL region hints (sub-box/slice placement) for Li/polymer/solvent/TFSI species.
- Applies an inner-contact safety floor via PACKMOL tolerance for biased modes (`li_inner_exclusion_nm`, default `0.19 nm`).
- Runs a cheap post-pack quality check on `initial.pdb` (no MD):
  - fraction of Li within `r_li_polymer` of polymer oxygen
  - fraction of Li within `r_li_solvent` of solvent oxygen
  - fraction of Li within `r_li_tfsi` of any TFSI atom (close-pair proxy)
- If polymer/solvent oxygen separation is unreliable (weak/ambiguous resname metadata), metrics are marked degraded instead of silently reusing one oxygen pool for both.
- Computes Li-to-nearest-heavy-atom distances for diagnostics in all modes; hard-fail enforcement remains preassembly-bias-only (`li_near_polymer` / `li_solvent_shell`) by design.
- If thresholds are not met, runs deterministic preassembly repair rounds controlled by `packmol_preassembly_retry_count` (independent of `--allow-density-reduction`).

Default metric cutoffs and thresholds:
- `r_li_polymer = 0.35 nm`
- `r_li_solvent = 0.35 nm`
- `r_li_tfsi = 0.40 nm`
- `target_li_polymer_fraction = 0.40` (for `li_near_polymer`)
- `target_li_solvent_fraction = 0.50` (for `li_solvent_shell`)
- `max_li_tfsi_close_fraction = 0.30` (both modes)
- `li_inner_exclusion_nm = 0.19`
- `packmol_preassembly_retry_count = 2`

Manifest (`composition.preassembly`) records:
- selected mode and settings
- per-round metrics, thresholds, violations, seed, applied hints
- contact safety diagnostics (Li-heavy min/median)
- whether repair was attempted and how many retries were used
- final preassembly report for downstream correlation with equilibration effort

`composition.yaml` example:
```yaml
packmol_preassembly:
  mode: li_solvent_shell
```



## ITP Sanitizer Stage - Hardening Policies

The sanitization stage enforces correctness, safety, and auditability for ITP files.

Module layout:
- `pipeline/stages/sanitizer_stage.py`: `SanitizerStage` orchestration entrypoint (`run()` and stage wiring).
- `pipeline/stages/topology_sanitizer.py`: topology/ITP/include/defaults handling, sanitizer-managed `system.top` updates, and final output writing helpers.
- `pipeline/stages/spatial_checker.py`: GRO parsing, PBC unwrapping, dipole checks, and grompp-based GRO/TOP consistency helpers.
- `pipeline/stages/sanitizer.py`: compatibility facade that re-exports the stage and legacy symbols for existing imports.

### Stage 1 Internal Parsing Refactor

`topology_sanitizer.py` now includes a small internal IR layer used only for the most fragile parsing paths:
- `AtomRecord`: one supported `[ atoms ]` row
- `MoleculeTypeIR`: one target moleculetype with parsed atom rows
- `TopMoleculesEntry`: one `[ molecules ]` entry with both raw count token and optional resolved integer count

Stage 1 is intentionally narrow:
- existing sanitizer entrypoints and `TopologySanitizerMixin` call shapes stay the same
- same-name moleculetype signature/charge comparison now consumes typed `[ atoms ]` rows internally
- ordered `[ molecules ]` extraction now delegates low-level row parsing to a dedicated helper and preserves unresolved count tokens explicitly
- only unambiguous supported 7/8-column `[ atoms ]` rows are accepted in this path
- unsupported/ambiguous row formats are surfaced as explicit parser states instead of being silently treated as valid data

### Stage 2 IR Migration Scope

Stage 2 keeps the same public sanitizer facade, but moves two safety-critical internals onto the Stage 1 IR:
- same-name moleculetype conflict detection now compares explicit IR parse results instead of falling through raw string heuristics
- deterministic moleculetype signatures now include typed atom fields and any locally resolvable `[ atomtypes ]` LJ definitions used by the target moleculetype
- uncertain same-name comparisons are never treated as equivalent:
  - strict mode fails closed
  - non-strict mode records a degraded comparison and treats the definitions as conflicting
- charge-fix classification/protection now uses the same typed moleculetype parse result
- protected/degraded molecules are removed from the checker's default auto-target pool before any charge-fix write occurs; explicit `--charge-fix-target-allowlist` is now the intended opt-in path for those targets


### Include Resolution
- **Shadowing Detection**: If an included file (e.g., `posre.itp`) exists in multiple search paths, the pipeline detects it.
  - Default: Warnings printed to console and audit log.
  - Strict mode (`--strict-include-resolution`): Hard error unless `--allow-include-shadowing` is present.
- **Same-basename policy (winner-only)**:
  - Raw ITP collection no longer hard-fails just because two files share a basename (for example, custom `PC.itp` shadowing forcefield `PC.itp`).
  - Ambiguity is handled at actual include resolution time using configured search order.
  - Sanitized output materializes only the resolved include winner for a given include target; shadowed losers are not written.
- **Coverage**: The same resolver/shadowing policy is applied to forcefield and non-forcefield include-closure expansion (not just top-level forcefield files).
- **Single enforcement path**: include policy checks (roots/escape/shadowing/strictness) are enforced through one include-closure expansion path; forcefield and non-forcefield files use identical policy semantics.
- **Priority**: Search order is configurable via `--itp-include-dirs-priority`.
  - `forcefield_first` (default; legacy alias: `last`): Local > System > FF > User.
  - `sanitized_first` (legacy alias: `first`): Local > User > System > FF.
- **Recursive include context**: nested `#include` parsing re-derives the search stack from the currently included file, so deeper includes do not accidentally inherit the original parent file's local precedence.
- **Escape protection**: include targets are blocked if they resolve outside allowed roots (including file parent, configured include search dirs, project root). Use `--allow-unsafe-include-escape` only when intentional.

### Charge Neutrality Protection
Solvents, ions, and protected polymers are guarded against accidental charge correction modification.
- **Active-molecule gate**: charge-neutrality checks are skipped only when `[molecules]` has no active species (`count <= 0` for all rows). If active molecules exist but sanitized molecule-to-ITP mapping is missing/empty, sanitizer fails fast with actionable diagnostics.
- **Rounding-only default policy**:
  - Auto-correction is allowed only when all active moleculetype charges satisfy `d = |Q_mol - round(Q_mol)| <= tol`.
  - `tol` is controlled by `--charge-fix-moleculetype-rounding-tol` (default `1e-4`).
  - If any active moleculetype exceeds `tol`, sanitizer fails fast by default (treated as topology/FF stitching error), unless `--charge-fix-allow-non-rounding` is explicitly set.
- **Diagnostics**: Charge audit prints a per-moleculetype table (`Q_mol`, `round(Q_mol)`, `d`, count, `Q_system` contribution) and dominant contributors; the same data is written to manifest.
- **Auto-protection**: 
  - Molecules matching `PROTECTED_SOLVENT_PATTERNS` (e.g., SOL, HOH, TIP3P).
  - Molecules matching `PROTECTED_ION_PATTERNS` (e.g., LI/LI+, TFSI, CL, NA).
  - Simple ions (monatomic) and known polymers.
  - Molecules matching `--charge-fix-protect-resnames <CSV>` are classified as user-protected before any charge-based correction targeting.
- **Classification reliability**:
  - If an active molecule cannot be classified reliably because `[ atoms ]` parsing is incomplete, preprocessed, or partially unparseable, strict mode fails closed.
  - Non-strict mode records the molecule as degraded/unknown and refuses correction on that target unless it is explicitly allowlisted via `--charge-fix-target-allowlist`.
  - The checker's default target allowlist is prefiltered through this audit, so degraded/protected molecules do not reach auto-patching just because they happen to share a legacy solvent token.
- **Overrides**:
  - `--charge-fix-allow-solvents`: Allow modification of protected solvents.
  - `--charge-fix-allow-ions`: Allow modification of protected ions.
  - `--charge-fix-target-allowlist <MOL>`: Explicitly allow specific molecules (case-insensitive).
  - `--charge-fix-protect-resnames <CSV>`: Additional protected residue/molecule names. Matching stays boundary-aware after separator canonicalization, so names like `PEG-400`, `PEG_400`, `TFSI-1`, `TFSI_1`, `TFSI-alpha`, or `G4-ALT` can be protected without broad short-token false positives.
- **Protected polymer net-charge policy**:
  - `--polymer-net-charge-tol` (default `1e-3`): if a protected polymer has `|q| <= tol`, sanitizer accepts drift and skips correction with manifest audit.
  - If `|q| > tol`: correction requires explicit target allowlist; non-strict mode warns+skips, strict mode fails.
  - `--charge-fix-polymer-method spread_safe`: when correction is allowlisted, distribute across many C/H atoms while excluding hetero atoms and their bond-neighborhood.
  - Polymer correction still requires rounding-only gate pass (or explicit `--charge-fix-allow-non-rounding` override).
- **Hetero safeguard**: any correction touching hetero atoms (`O/N/S/B/F/P`) is refused unless those atom names are explicitly allowlisted via `--charge-fix-target-atomnames`.
- **Strict Mode**: If strict charge neutrality is enabled, modifying any protected species requires explicit opt-in.
  - In strict mode, if any protected molecules exist, `--charge-fix-target-allowlist` is required.
  - Protected polymer correction is allowed only with explicit allowlist and explicit `--charge-fix-polymer-method spread_safe`.

### Parameter Validation
- **LJ Validations**:
  - Use `--lj-outlier-policy {off,warn,error}` with `--lj-bounds-profile {aa,coarse_grained,polarizable,custom}`.
  - For `--lj-bounds-profile custom`, provide `--lj-sigma-max-nm` and `--lj-epsilon-max-kj-mol`.
  - Dummy/virtual (`ptype D/V`) LJ outlier checks are skipped.
  - Negative sigma is warning-only by default because some legacy/special GROMACS encodings use it intentionally.
  - Negative sigma paired with `epsilon <= 0` is still warning-only by default for compatibility, but it is diagnosed as more suspicious than the legacy/special-encoding case.
  - If you want negative-sigma or other non-fatal LJ findings to hard-fail, use `--lj-outlier-policy error` (or legacy `--strict-lj-validation`).
  - For comb-rule `1`, near-zero positive `p1/p2` values that are more than 12 orders of magnitude below the median positive magnitude are excluded from log-domain MAD statistics and reported explicitly instead of silently distorting the outlier score.
  - Unit-screaming values (for example sigma > 10 nm or epsilon > 1e6 kJ/mol) always hard-fail.
- **Triclinic Support**:
  - Charges are validated on rectangular boxes.
  - Dipole correction for triclinic boxes is NOT supported (safe failure) to prevent artifacts.

### Robustness
- **Molecule-count source-of-truth**:
  - If staged `system.top` `[molecules]` parsing is uncertain/empty, sanitizer requires `grompp -pp` preprocessed `[molecules]` for trusted counts.
  - For staged/post-reaction topologies, sanitizer no longer silently falls back to manifest pre-reaction counts when trusted topology counts are unavailable.
  - `system.top` preserves existing `[molecules]` only when trusted raw staged `[molecules]` is used; otherwise it is regenerated from the trusted source used for validation.
- **Grompp Timeout**: Preprocessing (finding atoms) has an adaptive timeout logic.
  - Retries once with extended duration.
  - Configurable via `GMX_GROMPP_PP_TIMEOUT_S` env var or `--grompp-preprocess-timeout-s`.
  - Complexity metrics (includes, lines) are logged to manifest.
- **Context arg auditing**: `PipelineContext.from_args()` now coercively normalizes wrapper/orchestrator string inputs and tracks consumed keys.
  - Enable `--strict-context-args` to hard-fail on unconsumed pipeline-related keys.
  - Without strict mode, unconsumed pipeline-related keys are warned to `stderr` when `--verbose` is enabled.

## Strict Resume/Force Semantics (GROMACS)

- **Resume (default)**: `resume=True` by default. With no flags, completed stages are skipped only when `done.ok` and key stage artifacts/state are coherent.
- **No-resume**: Use `--no-resume` to run all stages even if they have `done.ok` (e.g., for validation).
- **Force**: `--force` archives prior outputs to `<stage_dir>/_archive/YYYYmmdd-HHMMSS/` before rerun.
- **Archive naming is collision-safe**: repeated force reruns in the same second use `-NN` suffixes to avoid collisions/overwrites.
- **Conflict**: `--force` and `--no-resume` cannot both be set (they are redundant/ambiguous).
- **MD resume is gated**: production resume requires `md.cpt`, `md.mdp`, and `stage_state.json` (plus `md.tpr` unless repair is explicitly enabled) with a compatibility-gate fingerprint match (MDP, recorded inputs, overrides, commands, GROMACS/git provenance). Mutable runtime artifacts (`md.cpt`/`md.tpr` hashes) remain in audit metadata but are excluded from strict resume compatibility gating.
- **Resume fingerprint uses recorded inputs first**: MD resume comparison resolves `fingerprint_payload` recorded paths from prior `stage_state.json`; if unresolved, it falls back to current inputs with warnings.
- **Repairable partial state**: if `md.cpt` exists but `md.tpr` is missing, pipeline fails fast by default and tells you to use `--allow-resume-repair` (or restore `md.tpr`). With `--allow-resume-repair`, compatibility is validated first, then `md.tpr` is rebuilt via `grompp`, then continuation runs with `mdrun -cpi md.cpt -append`.
- **Non-repairable partial state**: if `md.tpr` exists but `md.cpt` is missing, continuation is blocked with actionable guidance (restore `md.cpt` or start fresh with `--force`).
- **No destructive 鈥渞epair鈥?*: resume-repair never archives/deletes prior outputs.
- **Resume state refresh**: after successful production resume, `stage_state.json` is rewritten so subsequent resumes stay consistent with latest checkpoint/audit state.
- **`done.ok` is a final marker**: GROMACS stage sentinels are only written after required outputs and metadata are persisted.
- **No unsafe implicit restart**: missing `nvt.cpt` or `npt.cpt` is a hard error unless an explicit recovery path is selected, or non-strict auto-fallback can validate `.gro` velocities.
- **Velocity reset is explicit**: when `--allow-velocity-reset` is used, the MDP is forced to `continuation=no`, `gen_vel=yes`, and `gen_temp` is set (from `--temperature` or template/ref_t).
- **Velocity reset MDP validation is strict**: for both NPT and production MD velocity-reset paths, final patched MDP must still contain `continuation=no`, `gen_vel=yes`, and matching `gen_temp`.
- **POSRES detection is robust**: POSRES is detected from any `define` token containing `POSRES` (e.g., `-DPOSRES_HEAVY`, `-DPOSRES_FC=1000`, `-D POSRES_LIGAND`, quoted values).
- **POSRES reference is explicit or sticky-pinned**: if an MDP enables POSRES, `grompp` always receives `-r`. Use `--posres-reference` in strict mode; with `--no-strict-posres-reference`, the first encountered structure is pinned to `OUT_GMX/<RUN_ID>/03_gromacs/_shared/posres_reference.gro` and reused across stages/resume.
- **Restart strictness is decoupled**: checkpoint/restart strictness can be controlled with `--strict-restart-policy`. If unset, restart strictness inherits `--strict-mdp-validation` for backward compatibility.
- **Runtime-only fingerprint ignore (opt-in)**: `--resume-ignore-runtime` ignores runtime-only resume fingerprint keys (currently `mdrun_command_base`, e.g., thread/GPU/pinning changes) while preserving full command recording in stage audit files.
- **Sanitizer is required**: `grompp` fails if `combined_atomtypes_current.itp` or `itp_sanitized_current/` is missing unless `--allow-unsanitized-grompp` is set (not recommended).
- **grompp warnings are explicit**: default `--grompp-maxwarn=0`; override only when needed.
- **Dry-run** never writes `done.ok` or `stage_state.json`, and never archives outputs.
- **Index selection**: use `IN/.../gromacs/ndx/index.ndx` if present; otherwise allow exactly one `*.ndx`; multiple `.ndx` files are a hard error.
- **Index requirement**: if the final MDP requests custom groups (tc-grps/comm-grps/energygrps/etc.) and no `.ndx` is found, the pipeline fails before grompp with a clear error.

## Dispatcher Quality Gates (Thermodynamic Sanity)

Program/tool success (`exit_code == 0`) is necessary but not sufficient for physically meaningful sampling.  
Dispatcher QC gates add lightweight sanity checks on machine-readable stage metrics (currently emitted by `gmx_eq` after successful NPT):

- **Density plausibility** (when `expected_density_g_cm3` is available):
  - `rel_err = abs(density - expected) / expected`
  - warn/error if `rel_err > --qc-density-rel-tol`
- **Temperature sanity**: warn/error if average is NaN or outside `[50 K, 1000 K]`
- **Pressure sanity**: warn/error if average is NaN or `|P| > 5000 bar`

Policy controls:
- `--qc-policy {off,warn,error}`  
  - default: `warn` for `--stage all`, otherwise `off`
- `--qc-enable-after-stages` (comma list, default `gmx_eq`)
- `--qc-density-rel-tol` (default `0.05`)
- `--all-qc-stop-on-warn` (default `False`, applies to `--stage all`)

All-mode review points:
- In `--stage all` with `--qc-policy warn`, dispatcher emits highlighted review points after `htpolynet` and `gmx_eq`.
- If `--all-qc-stop-on-warn` is enabled (or policy is `error` and a QC check fails), dispatcher stops early after writing `run_report.json`.

## Audit Artifacts & SI Exporter

- **Per-stage**: `OUT_GMX/<RUN_ID>/03_gromacs/<stage>/repro.json` with hashes, commands, host/git/GROMACS provenance (launcher vs binary: `gmx_launcher_token`, `gmx_binary_token`, resolved paths, and version).
- **Run-level (dispatcher-owned)**: `OUT_GMX/<RUN_ID>/run_report.json` is always written, including on stage failure/interruption. It includes:
  - run identity + selected stage/resume/force
  - deterministic `stages_to_run`
  - per-stage `started_at_utc` / `ended_at_utc` / `duration_s` / `status` / error summary / traceback log path
  - top-level `warnings`, `qc.events` decisions, and `provenance.{tools,env,runtime,code}`
  - backward-compatible `stages` aggregate of GROMACS `repro.json` plus parse errors
- **Provenance text artifact**: `OUT_GMX/<RUN_ID>/provenance.txt` stores raw tool probe outputs (`gmx --version`, `gmx mdrun -version`, etc.), runtime/platform fields (`python_version`, `platform`, `uname`, `hostname`, `cwd`), git code state (`git_commit`, `git_branch`, `git_dirty`, or structured git error), and resolved executable/token metadata used for each probe.
- **Atomic JSON writes**: `stage_state.json`, `repro.json`, and dispatcher `run_report.json` are written via temp-file+fsync+replace to avoid truncated JSON on interrupts.
- **Atomic provenance text writes**: `provenance.txt` is written via `provenance.txt.tmp.<pid>` in the same directory, then flush+fsync+replace, with best-effort parent-directory fsync for durability.
- **Bounded probe capture**: stored probe `stdout`/`stderr` are capped to a fixed preview length and annotated with truncation flags/original lengths so artifacts stay bounded.
- **SI-ready snippet**: `OUT_GMX/<RUN_ID>/computational_details.md` generated by:
  `python tools/export_computational_details.py --run-root OUT_GMX/<RUN_ID>`
- **MDP source-of-truth**: exporter parses each stage's final `mdp.path` (absolute state), verifies file SHA256 against `run_report.json`, and fails fast on mismatch.
- **Partial pipelines**: exporter reports available stages by default; expected EM/NVT/NPT/MD sequence is still listed with missing stages marked `not run`.
- **Strict publication gate (optional)**: enforce canonical stages with `--require-stages em,nvt,npt,md`.
- **Output override (optional)**: use `--out <path>`; default remains `<run_root>/computational_details.md`.

## Manifest Schema

The manifest keeps backward-compatible keys while enforcing canonical composition semantics and crash-safe persistence.

```json
{
  "schema_version": 2,
  "version": "1.1",
  "composition": {
    "initial_counts": {"PEO": 50, "PEGDMA": 20, "LI": 10},   // immutable PACKMOL recipe
    "counts": {"PEO": 50, "PEGDMA": 20, "LI": 10},            // backward-compat alias of initial_counts
    "molecule_counts": {"PEO": 50, "PEGDMA": 20, "LI": 10},   // legacy alias of initial_counts
    "post_counts": {"CrosslinkedNetwork": 1, "LI": 10},       // post-crosslink state
    "post_counts_source": "htpolynet_system_top",
    "post_counts_timestamp_utc": "2026-01-26T11:00:00+00:00",
    "box_size_nm": 5.5,                                        // legacy scalar alias (when cubic)
    "box": {
      "box_vectors_nm": [5.5, 5.5, 5.5],                      // accepts scalar/[3]/[9]/3x3 on input
      "box_type": "orthorhombic"
    },
    "history": [
      {
        "stage": "packmol",
        "counts": {"PEO": 50, "PEGDMA": 20, "LI": 10},
        "reason": "set_initial_composition",
        "at_utc": "2026-01-26T10:00:00+00:00"
      }
    ],
    "conflicts": []
  },
  "box": {
    "packmol_box_vectors_nm": [5.5, 5.5, 5.5],
    "packmol_box_is_triclinic": false,
    "post_htpolynet_box_vectors_nm": [[5.4, 0.0, 0.0], [0.1, 5.3, 0.0], [0.0, 0.0, 5.2]],
    "post_htpolynet_box_is_triclinic": true
  },
  "commands": [
    {
      "argv": ["gmx", "mdrun", "-nt", "8"],
      "cmd": "gmx mdrun -nt 8",
      "returncode": 0,
      "stdout_preview": "... capped ...",
      "stderr_preview": "... capped ..."
    }
  ],
  "charge_corrections": [
    {
      "molecule": "PEO",
      "method": "uniform_over_atoms",
      "target_file": "IN/systems/SYS_001/gromacs/system.itp",
      "target_section": "atoms",
      "total_delta": 0.001,
      "affected_atoms": [{"atom_index": 1, "before": 0.12, "after": 0.121, "delta": 0.001}],
      "affected_atoms_count": 1,
      "validation": {"sum_matches_total_delta": true}
    }
  ]
}
```

> **Manifest immutability and compatibility rules**
> - `composition.initial_counts` is first-writer-wins and never overwritten.
> - `composition.counts` and `composition.molecule_counts` are compatibility aliases of `initial_counts`.
> - `composition.post_counts` is explicitly post-crosslink data; consumers must opt in to it.
> - `set_composition()` remains supported as a legacy wrapper: first call sets initial composition, later calls are interpreted as post composition.
> - Legacy manifests with only `composition.box_size_nm` are auto-migrated to `composition.box.box_vectors_nm`.

## HTPolyNet Stage Safety Policies

### Fail-Fast Behavior

**Without htpolynet installed**, the pipeline **fails immediately**. This prevents publishing non-physical placeholder outputs.

```bash
# This will FAIL if htpolynet is not installed:
python run_pipeline.py --ff GAFF2 --charge CM5 --stage htpolynet

# For demos/testing only, use --allow-placeholder:
python run_pipeline.py --ff GAFF2 --charge CM5 --stage htpolynet --allow-placeholder
# 鈿狅笍 Placeholders stay in OUT_GMX only (never published to IN/)
# 鈿狅笍 Outputs are marked placeholder=true and cannot proceed to GROMACS by default
# 鈿狅笍 Demo publish requires --allow-placeholder-propagate
# 鈿狅笍 Demo staging to IN/.../gromacs requires --allow-placeholder-stage-to-gromacs
# 鈿狅笍 grompp remains blocked by poison pill unless --allow-placeholder-gromacs-compile is set
```

### Input Provenance (Strict Mode)

The initial structure MUST come from `IN/systems/<SYS_ID>/packmol/` (published inputs). OUT fallback is **disabled by default** to prevent "stale artifact" contamination.

```bash
# If PACKMOL outputs are not published to IN/, this fails with actionable instructions
python run_pipeline.py --ff GAFF2 --charge CM5 --stage htpolynet

# Override with explicit flag (NOT RECOMMENDED for reproducibility):
python run_pipeline.py ... --unsafe-allow-out-fallback
```

### Cure Temperature Boiling-Risk Guard Rail (GPE-safe defaults)

`HTPolyNetConfig` now defaults to a safer curing profile for solvent-containing GPE systems:
- `cure_temperature_K = 350.0` (was 400.0)
- `cure_temp_policy = "warn"`
- `cure_temp_safety_margin_K = 10.0`

The config generator detects common volatile co-solvents from `molecule_counts` names (case-insensitive substring match):
- `DMC` / `DIMETHYL CARBONATE` (`bp 鈮?363.15 K`)
- `DME` / `DIMETHOXYETHANE` / `1,2-DIMETHOXYETHANE` (`bp 鈮?358.15 K`)

If requested `cure_temperature_K > (min_bp - margin)`:
- `cure_temp_policy="warn"`: keep requested value, emit prominent YAML warnings
- `cure_temp_policy="cap"`: cap to `(min_bp - margin)` in generated YAML
- `cure_temp_policy="error"`: fail fast with actionable error
- `cure_temp_policy="off"`: no automatic cap/error action (risk warning still documented in YAML comments when detected)

Unsafe override is explicit:
- `allow_cure_above_solvent_bp=true`
- `allow_cure_above_solvent_bp_reason="<why this is physically justified>"`

If override is enabled without a reason, config generation fails.

### Max Conversion Guard Rail (stress/LINCS risk)

`HTPolyNetConfig` now defaults to:
- `max_conversion = 0.85` (was 0.95)
- `recommended_max_conversion = 0.85`
- `conversion_policy = "warn"`

Rationale: forcing very high late-stage conversion can inject large internal stress
in vitrifying networks and later destabilize GROMACS (e.g., LINCS failures).

If `max_conversion > recommended_max_conversion`:
- `conversion_policy="warn"`: keep requested value, emit YAML warning comments
- `conversion_policy="cap"`: cap to `recommended_max_conversion`
- `conversion_policy="error"`: fail fast
- `conversion_policy="off"`: no automatic cap/error action

Unsafe high-conversion override is explicit:
- `allow_high_conversion=true`
- `allow_high_conversion_reason="<why high conversion is required>"`

If override is enabled without a reason, config generation fails.

### Search-Radius Strategy (deterministic + staged recipe)

Default single-radius behavior remains backward-compatible, with a safer default:
- `search_radius_nm = 0.35` (was 0.50)

Optional `search_radius_schedule` can be provided in config, but this pipeline does
**not** emit schedule keys into `system_config.yaml` unless HTPolyNet schema support
is confirmed. This avoids introducing unknown keys that may break parsing.

Instead, generated YAML includes a staged-run recipe comment for dense GPE systems:
- Stage 1: run to conversion `0.30` with `search-radius = 0.35 nm`
- Stage 2: run to conversion `0.60` with `search-radius = 0.45 nm`
- Stage 3: run to conversion `0.80` with `search-radius = 0.60 nm`

### Output Selection (Deterministic)

If multiple candidate outputs exist (e.g., `system.gro` and `final.gro`), the stage **fails fast** with a candidate list. No silent mtime-based guessing.

Resolution: clean the workdir, use unique run_ids, or specify outputs explicitly.

### Box Vector Handling

Both pre- and post-HTPolyNet boxes are recorded:

| Field | Description |
|-------|-------------|
| `box.packmol_box_vectors_nm` | Box from initial structure (before crosslinking) |
| `box.post_htpolynet_box_vectors_nm` | Box from staged GRO (after HTPolyNet/NPT relaxation) |

**Triclinic boxes (default fail-fast)**: If a 9-float GRO box line is detected, HTPolyNet stage fails by default because orthorhombic/cubic assumptions would otherwise require unsafe collapsing that can distort volume/density and PBC geometry.

Use `gmx editconf` to re-box to orthorhombic/cubic, then rerun.

**Anisotropy warning**: If box dimensions differ by >2%, a warning is emitted.

**Unsafe escape hatch**: `--allow-triclinic-unsafe-diagonal-approx`
- Enables the old diagonal approximation explicitly (not recommended)
- Emits loud warnings
- Records:
  - `htpolynet.box_scalar_policy = "UNSAFE_DIAGONAL_APPROX"`
  - `box.original_box_matrix_nm` (full 3脳3 for audit)
  - approximated diagonal vector used for config

### Configurable Timeout

HTPolyNet subprocess timeout is configurable (default: 43200s = 12 hours):

```bash
# For large/high-crosslink GPE systems:
python run_pipeline.py ... --htpolynet-timeout 86400  # 24 hours

# Disable timeout (emits warning, use with caution):
python run_pipeline.py ... --htpolynet-timeout 0  # None = no limit
```

### Gelation Feasibility Precheck

Before launching HTPolyNet, the stage computes conservative degree moments over reactive species (`reactive_site_count > 0`), weighted by molecule counts:
- `mean_f = 危(N_i f_i)/危(N_i)`
- `mean_f2 = 危(N_i f_i^2)/危(N_i)`
- `p_c_est = mean_f / (mean_f2 - mean_f)`

Default failure conditions:
- No species with `f_i >= 3`
- `mean_f2 - mean_f <= 0`
- `p_c_est > 1` (with small epsilon)

Overrides:
- `--skip-gelation-precheck`
- `--allow-gelation-precheck-fail` (downgrade fail to warning)

Manifest always records full precheck summary (`htpolynet.gelation_precheck`), including inputs, computed moments, `p_c_est`, pass/fail status, and override use.

### Crosslinking QC

HTPolyNet stage parses stdout/stderr for conversion and gel fraction using anchored, label-specific patterns only:
- `^(total )?conversion[:=] ... (%|fraction)`
- `^gel fraction[:=] ... (%|fraction)`

Generic `%` scraping is intentionally rejected to avoid false positives from unrelated progress lines.

To enforce a minimum conversion threshold:
```bash
python run_pipeline.py ... --min-conversion 0.5
# Fails if conversion < 50%
```

**QC strictness**:
- Default: if neither conversion nor gel fraction is parsed, stage fails (`allow-missing` is off)
- Override: `--allow-missing-htpolynet-qc` (unsafe)
- If multiple conflicting anchored matches are found, stage warns, picks the last anchored match deterministically, and stores matched lines in QC warnings for audit
- `--min-conversion` behavior is unchanged: missing conversion with threshold set still fails

### Placeholder File Sentinels

When `--allow-placeholder-propagate` is used, placeholder files include explicit sentinels for downstream detection:

| File Type | Sentinel Format |
|-----------|-----------------|
| .gro | Title line: `PLACEHOLDER=1 run_id=<RUN_ID> ts=<ISO8601>` |
| .itp/.top | Header comment: `;;; PLACEHOLDER=1` + run_id + timestamp |

Downstream stages can detect placeholders by checking:
- Manifest: `htpolynet.placeholder_propagated_downstream_unsafe == true`
- File content: `PLACEHOLDER=1` sentinel in first lines

Default propagation safeguards:
- Placeholder publish to `IN/.../htpolynet` is explicit via `--allow-placeholder-propagate`
- Placeholder staging to `IN/.../gromacs` is OFF by default; explicit opt-in: `--allow-placeholder-stage-to-gromacs`
- Poison pill is ON by default for propagated placeholders (intentional missing include in `system.top`) to block accidental grompp
- To disable poison pill (unsafe): `--allow-placeholder-gromacs-compile`


## ITP Sanitizer

The ITP Sanitizer prepares topology files for GROMACS by extracting atomtypes, detecting conflicts, and ensuring consistent `[ defaults ]` settings.

### What It Does

1. **Atomtype Extraction**: Scans all ITP files and extracts `[ atomtypes ]` definitions
2. **Conflict Detection**: Detects atomtypes with same name but different physical parameters
3. **Defaults Validation**: Ensures all `[ defaults ]` tuples are consistent across all reachable sources (`nbfunc/comb-rule/gen-pairs/fudgeLJ/fudgeQQ`)
4. **Combined Atomtypes**: Generates a single `combined_atomtypes_current.itp` with all unique atomtypes
5. **Sanitized ITPs**: Removes `[ atomtypes ]` and `[ defaults ]` from individual ITPs
6. **system.top Generation**: Creates topology with correct include order (defaults 鈫?atomtypes 鈫?extra includes 鈫?molecules), writes one deterministic sanitizer-managed block in `system.top`, deduplicates managed includes, and avoids duplicate `[ defaults ]` injection when existing `system.top` already has defaults/forcefield include; existing `[molecules]` is preserved only when raw staged `[molecules]` was the trusted source
7. **Molecule Scope**: `system.top` includes only molecule ITPs for species present in `[molecules]` (include-aware moleculetype detection)
8. **Include Integrity**: Resolves `#include` targets like GROMACS `-I`, detects shadowing in strict mode, sanitizes include-closure files with `[ defaults ]`/`[ atomtypes ]`, validates paths on disk, and writes winner-only sanitized files for include-target collisions

### Include Resolution (GROMACS -I)

For `#include` lines, the sanitizer rebuilds the search stack for each currently included file. Priority semantics are defined in one place (`_build_include_search_paths()`), and the per-file order is:
1. Including file's directory
2. User-provided `--itp-include-dirs` immediately after local dir only when `--itp-include-dirs-priority sanitized_first`
3. `IN/systems/<SYSTEM_ID>/gromacs/top` and `.../itp`
4. `IN/systems/<SYSTEM_ID>/htpolynet/itp`
5. Forcefield directory (`IN/forcefield/<FF>/gromacs`)
6. Sanitized output dir (`.../itp_sanitized_current`, if present)
7. Configured GROMACS share dirs (if provided via context)
8. User-provided `--itp-include-dirs` at the end when priority is `forcefield_first`

Relative paths in `--itp-include-dirs` are resolved against `--project-root` for reproducible behavior across working directories.
Compatibility-only extra include roots are appended after the configured order; they do not create a second priority mechanism.

If an include cannot be resolved, the sanitizer reports attempted paths (or fails fast if `--strict-include-resolution` is set).
`system.top` validation uses the same resolver and parses both `"file.itp"` and `<file.itp>` includes.
`system.top` include validation and managed-block updates now use the same conservative preprocessor policy: strict mode fails on unresolved conditional or macro-like include content, while non-strict mode logs an explicit degraded warning and skips conditionally ambiguous includes instead of treating them as active.

**Shadowing**: If multiple include candidates exist, shadowing is reported. In strict mode this is a hard error unless `--allow-include-shadowing` is enabled.
For same-basename shadows, sanitizer output is winner-only per include target (no deterministic renaming of losers).
Collection no longer enforces global basename uniqueness across all scanned ITP inputs; ambiguity is handled at include-resolution time using the configured search order.
**Escape guard**: includes that resolve outside configured include roots are blocked by default; use `--allow-unsafe-include-escape` only as an explicit unsafe override.
**Conditional/macro safety**: Python include/name discovery does not try to evaluate GROMACS preprocessor conditionals. If it sees unresolved macro/conditional directives in touched include-discovery paths, strict mode fails closed; non-strict mode logs explicit degraded fallback and skips conditionally ambiguous content instead of pretending it is active.

**Include-closure sanitization**: Any reachable file containing `[ defaults ]` or `[ atomtypes ]` is sanitized and written into `itp_sanitized_current/` at its include path (relative includes only; `..` is not allowed for sanitization).

### Section Header Parsing

Section headers with trailing comments are now correctly parsed:
```
[ atomtypes ] ; Define my atoms here  鈫? parsed as "atomtypes"
```

### Conflict Signature and Mass Tolerance

Atomtype conflicts are detected based on **physical parameters only**:
- `mass` - uses **robust abs+rel tolerance** (configurable):
  - `ITP_MASS_ABS_TOL` (default: 1e-3) - absolute tolerance
  - `ITP_MASS_REL_TOL` (default: 1e-4) - relative tolerance
  - Formula: `max(abs_tol, rel_tol * max(|a|, |b|))`
  - Example: Carbon mass 12.0110 vs 12.0107 **does not** conflict by default
- `sigma` and `epsilon` (relative tolerance 1e-5 with absolute floor)
- `ptype` (validated against A/D/S/V, missing ptype allowed)
- `comb-rule` (must match for comparison to be valid)

**Comb-Rule Awareness**: The sanitizer resolves the full include-closure first, then parses `[ defaults ]` to determine the `comb-rule` (1=C6/C12, 2/3=sigma/epsilon) before parsing atomtypes. This ensures correct LJ parameter interpretation and comparison.

**Mixed Defaults Policy (hard error by default)**:
- Distinct `[ defaults ]` tuples across sources fail fast by default because GROMACS compiles a single global `[ defaults ]`.
- The report classifies mismatch type and records all source files.
- Deterministic primary-defaults selection policy:
  1. Prefer `forcefield.itp` / `*.ff` defaults when present
  2. Else first encountered in include order
  3. Else previous selected defaults
- Unsafe mixed-default continuation requires:
  - `--allow-mixed-defaults`
  - `--allow-mixed-defaults-reason "TEXT"` (required audit string)
- Representation mismatch is never overridable:
  - If comb-rule `1` appears together with `2/3`, sanitizer hard-fails (`representation-mismatch`) even when `--allow-mixed-defaults` is set.
- For comb-rule `2` vs `3` mixes, use cross-group policy:
  - `--mixed-defaults-cross-group-policy {off,warn,generate}`
  - default is `warn` when mixed defaults override is enabled, else `off`
  - `--mixed-defaults-cross-group-rule {lorentz-berthelot,geometric}` (auto default from canonical comb-rule: `2->LB`, `3->geometric`)
  - `--mixed-defaults-cross-group-max-pairs N` (default: `20000`, hard cap)
  - `--mixed-defaults-cross-group-reason "TEXT"` (required when policy=`generate`)
- `generate` is implemented for mixed-defaults classification `comb-rule-only` with canonical `nbfunc=1` and canonical comb-rule `2` or `3`.
- `generate` writes versioned/current sanitizer-managed `nonbond_params_cross_group_current.itp`, and `topology_sanitizer.py` includes it deterministically when generating a new `system.top` or updating an existing managed block.
- When policy=`generate`, sanitizer verifies the expected cross-group artifact path exists before continuing; missing remediation artifacts are treated as hard failures.
- Existing explicit `pairtypes` / `nonbond_params` pairs are never overwritten; generated pairs skip preexisting entries and are reported.

### Limitations: Mixed Forcefields & Combination Rules

GROMACS uses one global `[ defaults ]` tuple for the compiled topology. If GAFF2-like (`comb-rule=2`) and OPLS-like (`comb-rule=3`) defaults are mixed, silently forcing a single comb-rule can distort LJ mixing and cause artifacts (for example phase separation or unstable diffusion behavior).  
This pipeline therefore fails by default on mixed defaults, and only allows continuation with explicit unsafe override + reason. For `2/3` mixed defaults, `--mixed-defaults-cross-group-policy generate` provides an explicit auditable override path via `[ nonbond_params ]`.

**Comb-Rule Safety**: If atomtypes are found but no `[ defaults ]` section exists, the sanitizer fails fast with an actionable error. Use `--allow-default-defaults` to permit a fallback defaults block (recorded in the manifest and used consistently for `system.top`).

**Relative Tolerances**: LJ parameters are compared with relative tolerances to avoid both false mismatches (very small C12 values) and false positives (large sigma differences):
- sigma/epsilon: relative tolerance 1e-5, absolute floor 1e-9/1e-12
- C6/C12: relative tolerance 1e-4, absolute floor 1e-18/1e-24

**Charge is excluded** from conflict detection by default because in many forcefields the `[ atomtypes ]` charge column is a placeholder (often 0.0), with real charges defined in `[ atoms ]` sections. Charge mismatches generate warnings but don't fail unless `--strict-charge` is enabled. For conflicting atomtype charges, the sanitizer now keeps deterministic source values (no silent zeroing) and reports the mismatch loudly.

### Virtual Sites / Dummy Atoms

TIP4P-like water models and other systems with virtual sites (ptype D or V) are handled correctly:
- **No LJ warnings** for dummy/virtual atoms with 蟽=0 and 蔚=0 (this is expected)
- Warnings still issued if ptype=D/V but LJ params are nonzero

### Robust Atomtype Parsing

The LJ-validation parser is intentionally conservative:
- Supported tails are `... mass charge ptype sigma epsilon` and `... mass charge sigma epsilon` (missing `ptype`)
- `ptype` is recognized only in the canonical third token from the tail and accepts `A/D/S/V`
- Bare five-/six-field rows whose first numeric token is integer-like are accepted only when the following token is lexically an explicit charge literal (signed or zero-like); otherwise they are rejected as ambiguous between `mass charge` and `at_num mass`
- Middle tokens ahead of `mass/charge` may exist, but rows with additional ptype-like middle tokens are rejected as ambiguous instead of guessed
- Files with unusual or shifted layouts produce actionable errors rather than silently reinterpreting columns

### Deterministic Ordering

For reproducibility:
- Files are processed in sorted path order (case-insensitive for cross-platform stability)
- Case-insensitive output path collisions (including basenames) are treated as errors (to avoid OS-dependent behavior)
- Override winners are selected by: **active forcefield match** (highest), then priority, then first file
- Same-priority ties go to the first file (deterministic)

### Override Logic and Active Forcefield

When conflicts are allowed via `--allow-override`:
1. Prefer entries matching the active forcefield (`--ff` value)
2. Fall back to file priority if no FF match
3. Emit clear warning showing winner/loser with source files and forcefields

### Prefix Safety

If `--atomtype-prefix` is used and implicit/library parameter lookup is detected, use:
- `--prefix-implicit-topology-policy {error,warn,bake,inject}` (default: `error`)
- `--prefix-implicit-topology-reason "TEXT"` (required for `bake` and `inject`)

Policy behavior:
- `error`: fail with actionable `PrefixUnsafeError` (safe default)
- `warn`: continue without prefixing (loud warning about namespace-collision risk)
- `bake` (preferred): materialize missing bonded parameters (`bonds`/`angles` and exact-match `dihedrals`) into explicit topology lines before prefixing
- `inject`: generate sanitizer-managed prefixed `[ bondtypes ]/[ angletypes ]/[ dihedraltypes ]` include entries for used system types (exact entries only by default)

All bake/inject actions are reported in sanitizer output/manifest with counts and unresolved items.

### Extra Includes for system.top

The `generate_system_top()` function supports an `extra_includes` parameter for injecting additional includes (e.g., POSRES files, water model ITPs) after atomtypes but before molecule definitions:

```python
from pipeline.itp_sanitizer import generate_system_top

content = generate_system_top(
    combined_atomtypes_path=...,
    sanitized_itp_paths=...,
    extra_includes=["posre_polymer.itp", "tip3p.itp"],
    ...
)
```

Generated output order:
1. `[ defaults ]`
2. `#include "<relative path to combined_atomtypes_current.itp>"`
3. Extra includes (POSRES, water models, sanitizer-generated sidecars), deduplicated in first-seen order
4. Molecule ITPs
5. `[ system ]` and `[ molecules ]`

When the sanitizer itself builds extra includes, the managed order is explicit:
1. `prefix_injected_types_current.itp` (if present)
2. `nonbond_params_cross_group_current.itp` (if present)
3. `nonbond_params_secondary_current.itp` (if present)

### Instance Reuse Safety

The sanitizer resets all internal state at the start of `run()`. This prevents state leakage if an instance is reused across multiple runs.

### Deduplicated Includes

`resolve_includes()` uses a visited set to ensure each file appears at most once in the resolved include list, preventing duplicates when multiple parents include the same file.

### Output Locations

**Note**: The sanitizer writes to `IN/systems/<SYSTEM_ID>/gromacs/` because these are "published inputs" for GROMACS鈥攔eusable assets that future runs can reference. Runtime outputs (logs, stage metadata) go to `OUT_GMX/<RUN_ID>/02_5_sanitizer/`.

GROMACS stages require the sanitizer outputs:
- `IN/systems/<SYSTEM_ID>/gromacs/combined_atomtypes_current.itp`
- `IN/systems/<SYSTEM_ID>/gromacs/itp_sanitized_current/`

Missing outputs cause an early failure unless `--allow-unsanitized-grompp` is explicitly set (not recommended).

### Defaults Handling

If no consistent `[ defaults ]` section can be determined, the sanitizer fails by default (to avoid guessing physics). Use `--allow-default-defaults` to permit a fallback defaults block.

**Float Tolerance**: Defaults are compared using `math.isclose(rel_tol=1e-4)` so minor formatting differences like 0.8333 vs 0.833333 do **not** trigger false conflicts.

### LJ Parameter Outlier Warnings (Configurable)

Configurable thresholds for atomtype LJ parameter warnings via environment variables:

```bash
# Mass tolerance
export ITP_MASS_ABS_TOL=1e-3
export ITP_MASS_REL_TOL=1e-4

# Sigma bounds (nm)
export ITP_SIGMA_MIN_NM=0.01
export ITP_SIGMA_MAX_NM=2.0

# Epsilon bounds (kJ/mol)
export ITP_EPSILON_MIN_KJ_MOL=0.0
export ITP_EPSILON_MAX_KJ_MOL=1000.0

# C6/C12 bounds (for comb-rule 1)
export ITP_C6_MIN=0.0
export ITP_C6_MAX=1e-2
export ITP_C12_MIN=0.0
export ITP_C12_MAX=1e-4
```

Also supported:
- `PIPELINE_LJ_SIGMA_MIN_NM`, `PIPELINE_LJ_SIGMA_MAX_NM`
- `PIPELINE_LJ_EPSILON_MAX_KJ_MOL`
- aggregate override `LJ_OUTLIER_THRESHOLDS` (for example `sigma_max_nm=1.5,epsilon_max_kj_mol=80`)

Warnings include explicit units and the file location for easy debugging.

### HTPolyNet [molecules] Preservation

After HTPolyNet crosslinks, the `[molecules]` section may differ from the pre-reaction manifest (e.g., multiple monomers become a single crosslinked network molecule). The Sanitizer detects this by:

1. Checking if `IN/systems/<SYSTEM_ID>/gromacs/top/system.top` exists (staged by HTPolyNet)
2. Parsing the `[molecules]` section from that file
3. Using those molecule counts instead of manifest pre-reaction counts

This ensures `grompp` receives a consistent `.gro/.itp/.top` set.

### GRO/TOP Consistency Validation

The Sanitizer validates that the coordinate file (`.gro`) matches the topology's `[molecules]` section by comparing atom counts.

- If an already available preprocessed topology exists (for example the sanitizer's cached `grompp -pp` artifact, or an explicitly provided preprocessed topology path in context), sanitizer prefers that as the topology truth source for semantic reads before falling back to raw Python parsing.
- If raw `system.top` `[molecules]` parsing encounters conditional-control or macro directives in that region, the raw parse is marked **uncertain** and conditional blocks are not treated as active by the fallback parser.
- When no trusted preprocessed topology is already available and the raw parse is uncertain/empty, sanitizer prefers `grompp -pp` preprocessed topology as the source of truth for molecule counts.
- If preprocessing is unavailable, sanitizer falls back to weaker sources (manifest/ITP fallback) and keeps uncertainty metadata in sanitizer report/manifest fields.
- Manifest QC/report fields include the uncertainty status and source choice (`topology_truth_source`, `raw_top_uncertain`, `raw_top_uncertainty_reasons`, `source`, `fallback_used`).

For atom-count validation, if GROMACS is available, sanitizer preprocesses topology via `gmx grompp -pp` and parses the fully expanded output. Otherwise it falls back to a strict ITP parser that ignores preprocessor lines and skips the check if ambiguity is detected (unless `--strict-gro-top-check` is set).

If a mismatch is detected, the pipeline fails fast with a clear error message:

```
[ERROR] GRO/TOP atom count mismatch!
        GRO file: 15000 atoms
        Expected from [molecules]: 12000 atoms
```

This catches problems early (before `grompp`) when HTPolyNet crosslinking changed the molecule composition but the coordinate file wasn't updated.

### Configurable grompp Preprocessing Timeout

The Sanitizer uses `grompp -pp` to preprocess topology for validation. For large crosslinked systems on slow filesystems (NFS), you can increase the timeout:

```bash
# Via environment variable (preferred for scripts)
export GMX_GROMPP_PP_TIMEOUT_S=600
python run_pipeline.py ... --stage sanitizer

# Default: 300 seconds, with adaptive scaling for large systems (>10k atoms)
```

If preprocessing times out under `--strict-gro-top-check`, an actionable error is provided with instructions to increase the timeout.

Executable resolution follows project-standard precedence: context override 鈫?`GMX_CMD` 鈫?default `gmx`.

`grompp -pp` stdout/stderr are streamed to debug log files under `OUT_GMX/<RUN_ID>/02_5_sanitizer/grompp_preprocess/` (no large in-memory capture). These debug logs are troubleshooting artifacts only and are excluded from reproducibility-critical fingerprinting.

`_read_gro_coords()` is safe by default: when called without `retain_resnames`, it does not materialize atom coordinates. Full-coordinate parsing requires `allow_full_parse=True` and is hard-capped by `GRO_ABSURD_ATOMS` with an actionable error on oversized GRO files.

### v4 Hardening: Dipole Check and Gating

The dipole shift check now runs in **heuristic mode by default**鈥攊t warns but never hard-fails unless `--strict-dipole-check` is enabled AND all preconditions pass.

`--strict-charge-physics` is not enforced inside `charge_neutrality.py`; dipole enforcement lives in sanitizer-stage dipole checks (`--strict-dipole-check` with `--charge-fix-max-dipole-shift`).

**Gating conditions** (enforcement skipped if any apply):
| Condition | Reason |
|-----------|--------|
| Giant molecule | Atom count > `dipole_unwrap_max_atoms` (default 512) |
| Giant residue group | Residue group atom count > `dipole_group_max_atoms` (default 2048) |
| Charged molecule | Net charge 鈮?0 (dipole is origin-dependent) |
| Percolated topology (bond-graph) | Alternate bond paths imply non-zero image/lattice shift (`percolated_topology`) |
| Missing bond graph | Dipole check is skipped as unreliable (`bond_graph_missing`) |
| Triclinic box | Dipole check skipped or failed per triclinic policy |

Wrapped span ratio is now a **suspicion metric**, not the decision-maker, when bond topology is available.

Dipole output now carries explicit reliability metadata in existing sanitizer charge-neutrality report fields:
- `dipole.reliable` (bool)
- `dipole.reliability_reasons` (list)
- `dipole.skip_reasons` (list)

Reliability rules:
- Missing bond graph skips dipole computation in non-strict mode and records explicit reliability/skip reasons.
- In strict dipole mode, missing bond graph fails with a clear reliability error (no silent downgrade).
- Span ratio alone never marks percolation.

Loop-consistency tolerance in bond-graph unwrapping is precision-aware:
- Default: `max(base_resolution_nm * 10, 1e-6)` where `base_resolution_nm` is inferred from GRO coordinate decimal precision.
- Override env knob: `GPE_UNWRAP_LOOP_TOL_NM`.

**CLI Options:**
| Flag | Default | Description |
|------|---------|-------------|
| `--strict-dipole-check` | false | Hard-fail on dipole violations if preconditions pass |
| `--dipole-unwrap-max-atoms` | 512 | Skip dipole enforcement for molecules with more atoms |
| `--dipole-group-max-atoms` | 2048 | Skip per-group dipole work for oversized residue groups |
| `--dipole-percolation-threshold` | 0.45 | Fallback span-ratio threshold for non-bonded unwrap path (and suspicion logging) |

**Example output (heuristic mode):**
```
[WARN] Dipole shift |螖渭|=0.0234 D exceeds limit 0.0100 D (heuristic mode)
[INFO] Dipole enforcement gated: giant_molecule (1024 > 512 atoms)
```

### v4/v7 Hardening: Include Resolution

Include resolution now uses a restricted search path list:
1. Current file's parent directory
2. `system_gromacs_dir/top` and `/itp`
3. `system_htpolynet_dir/itp`
4. Forcefield directory
5. Sanitized output directory
6. User-provided `--itp-include-dirs`

**v7 Hardening Additions:**

- **Full include-closure expansion**: All ITPs (not just forcefield) have their `#include` chains fully expanded via BFS, discovering moleculetypes in nested includes
- **Unified resolution helper**: `_resolve_include_in_dirs()` implements GROMACS `-I` style search order consistently across all include resolution
- **Context propagation**: Include directories are built per-file using `_build_include_search_paths()`, respecting the file's location
- **Class/FF map inheritance**: Included files inherit classification and forcefield from their parent ITP

**v7 Hardening: Protected Pattern Matching:**

Pattern matching for protected solvents/ions now uses boundary-aware matching:
- Normalization removes **all non-alphanumerics** (e.g., `+`, `-`, whitespace, parentheses) so `Li+` matches `LI`
- Short patterns (鈮? chars like "PC", "EC") require **exact match** after normalization
- Longer patterns can match as substrings but require word boundaries on the **raw name**
- Prevents false positives like "PC" matching "POLYMER_CHAIN"

### v4 Hardening: gen-pairs Consistency

After defaults are validated, the sanitizer scans ITPs for `[ pairs ]` sections and warns on potential conflicts:

- **gen-pairs=yes + [pairs] present**: May duplicate 1-4 interactions
- **gen-pairs=no + no [pairs]**: 1-4 interactions may be missing

This is warning-only since `#ifdef` conditional compilation may change actual behavior at `grompp` time.

### v4 Hardening: Comb-Rule Aware LJ Outlier Check

LJ outlier checks now use a **policy + profile** model:
- `--lj-outlier-policy {off,warn,error}` (default: `warn`)
- `--lj-bounds-profile {aa,coarse_grained,polarizable,custom}` (default: `aa`)
- `--lj-sigma-max-nm`, `--lj-epsilon-max-kj-mol` (required only for `profile=custom`)

Behavior:
- Comb-rule `2/3`: interpret LJ columns as sigma/epsilon and apply profile bounds
- Comb-rule `1`: treat LJ columns as raw `p1/p2` and use robust statistical outlier checks (MAD on log-scale) after excluding explicitly reported near-zero values that would otherwise pollute the log-domain statistics
- Dummy/virtual atomtypes (`ptype D/V`) are skipped to reduce false positives
- Negative sigma with `epsilon > 0` is treated as a compatibility warning/review case, not auto-classified as unit corruption
- Negative sigma with `epsilon <= 0` stays warning-only by default but is reported as more suspicious than the legacy/special-encoding case
- A second-level unit-screaming guard always hard-fails for absurd values (e.g., sigma > 10 nm or epsilon > 1e6 kJ/mol)

A structured LJ report is always emitted (manifest `sanitizer.lj_validation`) with active profile/policy, thresholds, parse issues, and top offenders including provenance and source locations.

### v4 Hardening: Sentinel Block system.top Updates

The Sanitizer now uses **sentinel-block based updates** to preserve custom content in `system.top`:

```
; >>> sanitizer managed block begin
[ defaults ]
...
#include "../combined_atomtypes_current.itp"
#include "molecule.itp"
; <<< sanitizer managed block end
```

**Preservation guarantees:**
- Content **outside** sentinel markers is never modified
- `#define`, `#ifdef`, extra includes, and custom sections are preserved
- If markers already exist, sanitizer always uses that managed block location (preferred placeholder mechanism)
- Managed-block body generation is deterministic and idempotent: repeated runs rewrite the same block body, dedupe repeated managed includes, and keep molecule includes in the sanitizer's trusted `[molecules]` order
- If markers don't exist, insertion is deterministic: immediately after the last already-owned defaults prefix before the first non-defaults topology section (direct `[ defaults ]`, direct pre-section include with direct defaults, or unconditional `forcefield.itp` / `*.ff` include); if no such ownership prefix exists, after leading comment/blank/`#define`/`#undef` preamble and before other pre-section includes
- Sanitizer does not attempt broad macro-block parsing to find a "smart" insertion point
- If **multiple blocks** are found, the sanitizer replaces all managed blocks with one clean block (strict mode via `--strict-top-update` or `--strict-include-resolution` fails fast)
- If markers are **malformed** (unmatched or misordered), strict mode (`--strict-top-update` or `--strict-include-resolution`) fails fast; non-strict mode only rebuilds when the damaged payload is still recognizable as managed-block content, otherwise it fails closed instead of deleting user topology sections
- Existing-top defaults ownership is checked conservatively: direct `[ defaults ]`, unconditional `forcefield.itp` / `*.ff` includes, and direct pre-section includes whose target file directly defines `[ defaults ]` are treated as owned defaults prefixes; nested pre-section include chains remain ambiguous, fail closed in strict mode, and are flagged as degraded warnings otherwise
- Atomic writes prevent partial files on crash

### v4 Hardening: Moleculetype Exact-Name Policy

Moleculetype names are now tracked with **exact case**:

- **Case collisions** (e.g., `UEMA` vs `uema`) trigger a **fail-fast error**
  - These will break grompp on case-insensitive filesystems
- **Same-name conflicts** from different sources are resolved by priority:
  - `htpolynet` > `staged` > `molecules-lib` > `forcefield`
- If same-name signatures differ:
  - strict mode (`--strict-forcefield-consistency`, `--strict-include-resolution`, or strict sanitizer flow) fails fast
  - non-strict mode logs a clear warning and picks by the priority rule above
- **Signature-based detection**: streaming `sha256` over full canonical `[ atoms ]` rows (`atomtype|atomname|charge_quantized`) plus section flags
  - Quantization step is configurable by environment variable `GPE_MOLECULE_SIG_QUANT_STEP` (default `1e-6`)
  - Secondary per-atom charge check runs before accepting "same signature", with tolerance `quant_step / 2`
- Debug stats attached to signature: `qmin`, `qmax`, `sum_abs_q`, `sum_q2`
- Identical definitions from different sources are silently de-duplicated

### v4 Hardening: Dipole Vector Delta

The dipole shift calculation now computes the correct **per-group vector delta**:

- For each molecule group: `螖渭_i = ||渭_after_vec - 渭_before_vec||`
- Reports and persists both `avg(螖渭_i)` and `max(螖渭_i)` in sanitizer manifest output
- Avoids direction-cancellation artifacts from averaging vector directions
- **Wrapped-span pre-gate**: skips obvious spanning groups before any unwrap allocation
- **Streaming unwrap span check**: can early-exit when unwrapped span ratio crosses threshold

### Duplicate [molecules] Entries

GROMACS allows the same molecule name to appear on multiple lines in `[molecules]`. The Sanitizer handles this correctly:

- Duplicate entries are **collapsed** by summing counts into the first occurrence
- Each molecule ITP is included at most once in `system.top`
- No assertion errors for repeated molecule lines


### Charge Neutrality

The Sanitizer performs charge neutrality preflight checking:

**Fail-Fast Behavior**:
- Fails if any molecule in `[molecules]` with count > 0 has no corresponding ITP file
- Fails if an active molecule ITP cannot be parsed or has empty `[ atoms ]` section
- Fails on case-insensitive alias/moleculetype collisions that make active molecule resolution ambiguous
- In strict mode, fails if any active molecule has preprocessor directives in/around `[ atoms ]`, or unparsed/ambiguous `[ atoms ]` data lines
- Fails if system charge |Q| exceeds auto-correction threshold (default: 1e-4)
- Optional unsafe escape hatch exists only as an explicit checker flag (`allow_unparsed_atoms_lines=True`), skipped lines are recorded in charge-neutrality audit output, and outcomes are marked degraded in non-strict mode
- The charge checker now uses the same fail-closed row semantics as the Stage 1 `[ atoms ]` parser, so ambiguous raw 7-column, glued charge/mass, or negative-mass rows are treated as unparseable rather than silently counted

**Alias vs Real ITP Moleculetype**:
- Molecule counts remain keyed by pipeline alias for backward compatibility
- Charge patching, post-write verification, and bond-graph parsing always use the real `[ moleculetype ]` name discovered in the ITP
- Target `[ moleculetype ]` matching for patching is case-insensitive to avoid alias/casing drift bugs

**Strict vs Non-Strict Mode** (`--strict-charge-neutrality`):
- **Strict mode**: Enforces reliability at system-validation scope (all active molecules), not only correction target selection. If any active source is unreliable, the stage fails even when total `Q 鈮?0`.
- **Non-strict mode** (default): Returns explicit degraded outcomes (`degraded_unreliable_charge_sources`, `degraded_rounding_drift_non_strict`) instead of reporting a clean neutral pass.

**Safe-Subset Correction** (default, `--charge-fix-method safe_subset`):
- **Rounding-only gate (default)**:
  - Auto-correction requires `|Q_system| <= --charge-neutrality-correct`
  - and all active moleculetypes must satisfy `|Q_mol - round(Q_mol)| <= --charge-fix-moleculetype-rounding-tol` (default `1e-4`).
  - Non-rounding drift is treated as topology/FF stitching error and fails fast unless `--charge-fix-allow-non-rounding` is set.
- Adjusts ONLY atoms with low |q| and conservative name matching:
  - Hydrogens: `H`, `H1/H2/...`, or common H-class tokens (HC/HA/HP/HGA...)
  - Carbons: `C1/C2/...` or common carbon class tokens (CT/CB/CG/...)
  - Explicitly excludes element/ion-like names (CL/BR/F/I/NA/LI/K/MG/CA/Zn/Fe/etc.)
  - Prefer false negatives over false positives; override with `--charge-fix-target-atomnames` if needed
- Target must be in the configured allowlist (explicit target selection does not bypass allowlist safety)
- Auto target selection is deterministic and ranks by estimated per-atom `|螖q|` safety, with exclusion diagnostics when no candidate is acceptable
- Uniform smearing available via `--charge-fix-method uniform_all`
- Enforces max per-atom 螖q and max total 螖Q guardrails
- Per-atom limit diagnostics include `|Q| / total_adjustable_atoms`; small systems may fail this guard while larger systems pass for the same `|Q|`
- If no safe atoms found, correction is refused with actionable message
- Correction is refused if selected atoms include hetero atoms unless explicitly allowlisted atom names are provided
- Disallowed atomtypes are enforced case-insensitively (`Cl`/`cl`/`CL` are equivalent)

**Protected Polymer Policy**:
- `--polymer-net-charge-tol` (default `1e-3`) defines acceptable per-molecule polymer drift.
- If protected polymer `|q| <= tol`: correction is skipped as acceptable drift, with audit fields (`target`, `net_q`, `tol`, `decision`) recorded.
- If protected polymer `|q| > tol` and not explicitly allowlisted:
  - non-strict mode: warn + skip correction
  - strict mode: hard fail
- If allowlisted, polymer correction uses `spread_safe`:
  - distributes `螖Q` uniformly across many atoms (minimizes max `|螖q_i|`)
  - excludes hetero atoms (`O/N/S/B/F/P`)
  - excludes atoms within `--charge-fix-polymer-exclusion-bonds` (default `2`) of hetero atoms
  - fails with actionable guidance if too few eligible atoms remain

**Ion Detection** (Advisory):
- Molecules detected as ions (integer/half-integer charge, 鈮? atoms) are skipped by default
- Use `--charge-fix-allow-ions` to override this advisory and allow correction on ions

### Limitations: Fixed-Charge FF in Concentrated Electrolytes

For concentrated electrolytes, non-polarizable fixed-charge force fields miss electronic polarizability and may over-bind ions.  
When charge-fix is refused on protected ions/solvents or pressure shifts to polymer/backbone correction, the sanitizer emits targeted warnings and mitigation guidance (ECC-style scaling or polarizable models).

Recommended safer automation settings for GPE workflows:
- `--charge-fix-polymer-exclusion-bonds 3`
- `--polymer-net-charge-tol 0.05`
- keep `--strict-gro-top-check` enabled

**Robust ITP Parsing** (Format Variants):
- Handles both 8-token format (`nr type resnr residue atom cgnr charge mass`) and 7-token format (missing `cgnr`)
- Tolerates trailing comments on section headers (`[ atoms ] ; comment`)
- Recovers glued numeric fields in `[ atoms ]` charge extraction and `[ atomtypes ]` LJ tail parsing (e.g., `1.000-0.001`, `1.234-5.678`)
- Supports scientific notation in charge/mass tokens (e.g., `1.2e-5`)
- Uses token-based patching for reliable charge replacement
- Prevents silent charge under-counting by refusing unparseable target `[ atoms ]` data lines unless explicit unsafe override is enabled

**Post-Write Verification**:
- Writes patched ITP to a staging location first, then re-parses to verify molecule charge matches expected value
- Recomputes full system charge using prospective in-memory state before publishing patched output
- Publishes patched ITP and updates in-memory molecule charge cache only after full verification succeeds
- Only claims "Correction applied" if both verifications pass
- Uses high precision (12 decimal places) in patched ITP to avoid rounding no-ops

**System Charge Arithmetic**:
- `compute_system_charge()` uses float contributions plus `math.fsum`
- This avoids `Decimal(str(float))` conversion artifacts while remaining robust at the configured correction scales (`1e-4` and below)

**Run-Scoped Outputs**:
- Patched ITPs are named `{molecule}_corrected_{run_id}.itp` for reproducibility
- Patched ITP writes are transactional (staging write + verification + atomic publish), with UTF-8 I/O and ensured output directory creation
- Each patched file includes a comprehensive audit header with thresholds, target metadata, and explicit correction semantics:
  - `effective_method` and `method_detail` (actual algorithm used, not requested default)
  - `delta_per_atom`
  - `atoms_patched_per_molecule`
  - `target_count`
  - `per_molecule_delta`
  - `Total system charge correction` (system-level, ~`-Q`)
- Atom-level correction audit entries include locator fields (`nr`, `atomname`, `atomtype`, `line_no`, `charge_col_idx`, `line_excerpt`, `old_q`, `new_q`) with payload caps + truncation counters
- Runtime warning logs honor `--charge-neutrality-warn` (configured `warn_threshold`) consistently

**CLI Options**:
| Flag | Default | Description |
|------|---------|-------------|
| `--strict-charge-neutrality` | false | Fail-fast on unreliable active charge sources (preprocessor, ambiguous/unparsed `[ atoms ]`, or name collisions) |
| `--charge-neutrality-tol` | 1e-6 | Non-negative neutrality pass tolerance; must satisfy `tol <= warn` |
| `--charge-neutrality-warn` | 1e-4 | Non-negative warning threshold; must satisfy `tol <= warn <= correct` |
| `--charge-neutrality-correct` | 1e-4 | Non-negative auto-correction limit; must satisfy `warn <= correct` |
| `--charge-fix-moleculetype-rounding-tol` | 1e-4 | Max `|Q_mol-round(Q_mol)|` treated as rounding drift |
| `--charge-fix-allow-non-rounding` | false | Unsafe override to permit correction for non-rounding drift |
| `--polymer-net-charge-tol` | 1e-3 | Acceptable per-molecule protected-polymer net charge drift |
| `--charge-fix-polymer-method` | skip_if_small | Protected-polymer policy (`skip_if_small` or `spread_safe`) |
| `--charge-fix-polymer-exclusion-bonds` | 2 | Exclude atoms within N bonds of hetero atoms for `spread_safe` |
| `--charge-fix-target-molecules` | 鈥?| Comma-separated target molecules (must still satisfy allowlist policy) |
| `--charge-fix-target-atomnames` | 鈥?| Comma-separated atom name patterns |
| `--charge-fix-disallowed-atomtypes` | 鈥?| Comma-separated atomtypes to never adjust |
| `--charge-fix-allow-ions` | false | Allow correction on detected ions |
| `--charge-fix-source-counts` | auto | Source for molecule counts: auto/manifest/system_top |




### v5 Hardening: Triclinic Box Policy

For triclinic boxes, the `--dipole-triclinic-policy` controls behavior:
- `skip` (default): Skip dipole check
- `error`: Fail with actionable message
- `mic`: Correct triclinic MIC (future)

### v5 Hardening: Include Path Priority

`--itp-include-dirs-priority {forcefield_first,sanitized_first,first,last}` controls user path precedence:
- `forcefield_first` (default, alias `last`): User paths searched last (backward compatible)
- `sanitized_first` (alias `first`): User paths searched first (explicit overrides)
- invalid values are rejected with a clear error.

Shadowing is logged with warnings showing all alternatives.

### v5 Hardening: Ion Protection for Charge-Fix

Moleculetypes are classified by total charge:
- `monatomic_ion`: 1 atom with |q| > 0
- `ionic`: |round(q)| >= 1
- `neutral`: Others

Post-correction verification ensures ionic species weren't modified unless `--charge-fix-allow-ions`.

### v5 Hardening: Enhanced Signature & LJ Warnings

- Moleculetype signatures now hash canonicalized topology-relevant sections (`[ atoms ]`, connectivity, restraints, virtual-site blocks) to catch materially different same-name definitions
- Distinguishes same-topology neutral molecules across charge models (e.g., RESP vs CM5)
- Default charge quantization in signatures is `1e-6` (`GPE_MOLECULE_SIG_QUANT_STEP`) with secondary per-atom compare tolerance `quant_step/2`
- Conflict reporting includes deterministic winner and strict/non-strict behavior
- LJ negative values get loud warnings; zero values (virtual sites) skipped silently

Minimal signature-compare example:
```python
from pathlib import Path
from pipeline.stages.sanitizer import SanitizerStage

stage = SanitizerStage()
sig_resp = stage._compute_moleculetype_signature(Path("RESP.itp"), "MOL")
sig_cm5 = stage._compute_moleculetype_signature(Path("CM5.itp"), "MOL")
print(sig_resp != sig_cm5)  # Expected: True when internal charges differ
```

Minimal wrapped-span pre-gate example:
```python
from pipeline.stages.sanitizer import SanitizerStage, GROAtom

stage = SanitizerStage()
group = [GROAtom(1, "POL", "C1", 0.01, 0.01, 0.01), GROAtom(1, "POL", "C2", 9.95, 0.02, 0.02)]
span = stage._wrapped_span(group)
box = (10.0, 10.0, 10.0)
ratio = max(span[0] / box[0], span[1] / box[1], span[2] / box[2])
print(ratio > 0.45)  # Expected: True -> suspicion log only; final decision uses bond-graph unwrap
```


## MDP Dynamic Patching

MDP templates in `IN/systems/<SYSTEM_ID>/gromacs/mdp/` are copied to `OUT_GMX/`, patched with CLI overrides, and written atomically (UTF-8 temp file + replace) before `grompp`.

### Patching Logic

| CLI Flag | MDP Parameters Affected |
|----------|------------------------|
| `-T <K>` | `ref_t`, `gen_temp` |
| `-P <bar>` | `ref_p` |
| `--nsteps-em` | `nsteps` (EM only) |
| `--nsteps-nvt` | `nsteps` (NVT only) |
| `--nsteps-npt` | `nsteps` (NPT only) |
| `--nsteps-md` | `nsteps` (production only) |
| `--tc-grps-mode auto/system/split` | `tc-grps`, `tau_t`, `ref_t`, `comm-grps`, `comm-mode` |
| `--tau-t-polymer` | `tau_t` (first value) |
| `--tau-t-nonpolymer` | `tau_t` (second value) |
| `--tau-t-values` | `tau_t` (custom list; count must match `--tc-grps-groups`) |
| `--tc-grps-min-atoms` | Split safeguard threshold; tiny groups fallback to `System` |
| `--tc-grps-min-fraction` | Split safeguard threshold by atom fraction; tiny groups fallback to `System` |
| `--allow-small-tc-grps` | Expert opt-in to keep tiny split groups |
| `--comm-grps-policy` | `comm-grps`, `comm-mode` |
| `--allow-comm-mode-angular-with-pbc` | Expert opt-in to keep `comm-mode=Angular` under PBC |
| `--force-gen-vel` | `gen_vel` |
| `--gen-seed` | `gen_seed` (when `gen_vel=yes`, including `--force-gen-vel` in NPT/MD) |
| `--npt-barostat` | `pcoupl`, `tau_p`, `compressibility` safety injection/warnings (NPT only) |

### Deterministic Override Order

Overrides are applied in dependency-aware phases to avoid order-dependent tc-grps/ref_t mismatches:

1. Mode-changing overrides (`tc-grps-mode`, `tc-grps`, `comm-grps-policy`)
2. Scalar temperature/pressure overrides (`temperature`, `gen_temp`, `pressure`)
3. Explicit vector overrides (`ref_t`, `tau_t`, `ref_p`, `tau_p`)
4. Everything else

This guarantees that tc-grps split/unified happens before ref_t/tau_t are expanded.

### Final Vector Normalization

After all overrides, the patcher normalizes tc-grps-coupled vectors:

- If `tc-grps` has **N** groups and `ref_t` has 1 value 鈫?replicate to **N**
- If `tc-grps` has **N** groups and `tau_t` has 1 value 鈫?replicate to **N** only with `--allow-tau-t-autofill`
- Any other length mismatch (`ref_t` or `tau_t`) 鈫?**ERROR** (grompp would fail)

This ensures `--temperature 300` plus `--tc-grps-mode split` always yields `ref_t = 300 300` (and matching `tau_t` length).

### tc-grps Modes

- **auto** (default, conservative): split into `Polymer NonPolymer` only when clearly safe, otherwise fallback to `System`
- **system**: `tc-grps = System` 鈥?single thermostat group
- **split**: `tc-grps = Polymer NonPolymer` 鈥?separate groups (explicit/legacy mode; requires matching ndx)
- **split + custom groups**: `--tc-grps-groups "A B C"` 鈥?explicit group list (recommended for non-Polymer/NonPolymer labels)

**Note**: `--tc-grps-groups` is whitespace-delimited; group names cannot contain spaces.

**Auto mode contract**:
- Requires exact ndx groups `Polymer` and `NonPolymer`.
- Each group must satisfy both thresholds:
  - atom count >= `--tc-grps-min-atoms`
  - atom fraction >= `--tc-grps-min-fraction`
- If auto split is selected, applied groups are exactly `Polymer NonPolymer` from that decision (no fallback to arbitrary first-two ndx groups).
- If checks fail, patcher falls back to `tc-grps = System` with warning.

**System mode rewrite policy**:
- `tc-grps-mode=system` actively rewrites final state (not a no-op): `tc-grps=System`, `comm-grps=System`, and compatible `comm-mode`.
- For split鈫抯ystem fallback, vector collapse is safety-first:
  - auto-collapse `tau_t`/`ref_t` only when all values are identical,
  - otherwise fail clearly (unless expert override is explicitly enabled).

### tc-grps Split Safety (GPE Workflows)

> [!WARNING]
> In gel polymer electrolyte (GPE) simulations, Polymer and NonPolymer groups have different dominant timescales. Blindly duplicating a single `tau_t` value for both groups can cause non-ergodic energy partitioning ("hidden temperature partitioning").

**Default behavior**: Split mode **fails fast** if `tau_t` in the template has only one value:

```bash
# This will FAIL:
python run_pipeline.py --ff GAFF2 --charge CM5 --stage gmx_eq --tc-grps-mode split

# Solutions:
# 1. Provide explicit per-group tau_t (RECOMMENDED):
python run_pipeline.py ... --tc-grps-mode split --tau-t-polymer 0.5 --tau-t-nonpolymer 0.1

# 1b. Custom groups with explicit tau_t list:
python run_pipeline.py ... --tc-grps-mode split --tc-grps-groups "Li TFSI Solvent" \
    --tau-t-values "0.1,0.5,2.0"

# 2. Opt-in to legacy autofill (NOT RECOMMENDED):
python run_pipeline.py ... --tc-grps-mode split --allow-tau-t-autofill
```

**Template compatibility**: MDP templates ship with `tc-grps = System` and a single `tau_t`/`ref_t` value. When split mode is enabled, you must provide per-group `tau_t` values; the patcher will expand `ref_t` to match the group count.

**Tiny-group DOF safeguard**:
- Split mode checks group atom counts from `.ndx` when available.
- Default threshold: `--tc-grps-min-atoms 50`.
- Default fraction threshold: `--tc-grps-min-fraction 0.01`.
- Auto-detected tiny groups: split is disabled and patcher falls back to `tc-grps = System`.
- Explicit tiny groups: strict mode errors; non-strict mode falls back to `System` unless `--allow-small-tc-grps` is set.
- Rationale: very small thermostat groups can amplify noise (low DOF), causing temperature spikes and instability.

### COM Removal Policy (Issue B)

`comm-grps` policy is now explicit and safety-oriented for percolated polymer/gel/network systems:

| Policy | Behavior |
|--------|----------|
| `auto` | If `tc-grps` contains `polymer`/`gel`/`network`/`matrix` (case-insensitive): use `system`; otherwise use `match_tc` |
| `match_tc` (default) | `comm-grps = tc-grps` (per-group COM removal). Emits a **strong warning** for percolated-risk group names |
| `system` | `comm-grps = System`; removes overall drift only (usually safer default) |
| `none` | `comm-mode = None` and `comm-grps = System` (disables COM removal; drift may accumulate) |
| `explicit` | Use `--comm-grps <groups>` value |

If `comm-grps` resolves to multiple groups, the patcher validates those names against `.ndx` groups when available:
- strict mode: error on missing group
- non-strict mode: warning with remediation

```bash
# Override the default policy:
python run_pipeline.py ... --tc-grps-mode split --comm-grps-policy auto
python run_pipeline.py ... --tc-grps-mode split --comm-grps-policy system
python run_pipeline.py ... --tc-grps-mode split --comm-grps-policy none
python run_pipeline.py ... --tc-grps-mode split --comm-grps-policy explicit --comm-grps "Protein Solvent"
```

### comm-mode Angular + PBC Safety

- If `comm-mode=Angular` and PBC is active (`pbc` missing is treated as `xyz`), default behavior is to force `comm-mode=Linear` with warning.
- Strict mode (`--strict-mdp-validation`) fails fast instead of overriding.
- Expert escape hatch: `--allow-comm-mode-angular-with-pbc` keeps Angular but emits a high-visibility warning.

### Stage Transition Velocity Handling (Issue C / P0.1-P0.B)

The MDP patcher enforces consistent stage policies for `continuation` and `gen_vel`:

| Stage | continuation | gen_vel | Rationale |
|-------|--------------|---------|-----------|
| NVT | `no` | `yes` | Fresh start, generate velocities at gen_temp |
| NPT | `yes` | `no` | Continue from NVT checkpoint |
| Production MD | `yes` | `no` | Continue from NPT checkpoint |

When `--allow-gro-velocity-restart` is set for NPT/MD, stage defaults prefer `continuation=no` (expert restart intent). If a checkpoint is actually present, checkpoint policy still explicitly overrides to continuation mode.

**Override precedence** (highest 鈫?lowest):
1. `--force-gen-vel` 鈫?sets gen_vel=yes
2. `--allow-continuation-override` 鈫?skips auto-setting continuation
3. Stage defaults

> [!CAUTION]
> `gen_vel=yes` + `continuation=yes` is ALWAYS an error (contradictory settings).

**If `gen_vel = yes`**, the patcher validates:
- `gen_temp` is present (if missing, near-0K behavior)
- `gen_temp` is numeric (catches "300K" typos)
- `gen_temp` matches `ref_t`
- `gen_seed` is set (warns if not, as velocities will be non-reproducible)

### Expert Restart from .gro with Velocities (P0.2)

For advanced workflows restarting from a `.gro` with embedded velocities (no checkpoint):

> [!WARNING]
> **`--allow-velocity-reset` and `--allow-gro-velocity-restart` are mutually exclusive.**
> If both flags are set, CLI validation fails fast via `parser.error` (exit code `2`) with guidance on the two restart modes and a "choose ONE" message.

> [!IMPORTANT]
> **GRO Velocity Validation**: When using `--allow-gro-velocity-restart`, the pipeline now validates
> multiple sampled atom lines across the `.gro` atom block and requires parseable coordinate/velocity columns.
> If velocities are missing, it fails immediately with a clear error instead of silently starting near 0K.

If `.gro` velocity status is known missing, this is a hard error even in non-strict mode:
- thermostat startup can inject large energy
- hotspots/constraint failures become likely
- remediation: provide `.gro` with velocities, or use `gen_vel=yes + gen_temp`, or `continuation=yes + .cpt`

**Recovery Flag Decision Matrix:**

| Checkpoint | allow_velocity_reset | allow_gro_velocity_restart | Result |
|------------|---------------------|---------------------------|--------|
| Exists | 鈥?| 鈥?| Use checkpoint (normal continuation) |
| Missing | True | False | Regenerate velocities (thermal shock risk) |
| Missing | False | True | Use .gro velocities (validated) |
| Missing | True | True | **CLI error** (`parser.error`, mutually exclusive) |
| Missing | False | False | ERROR (strict) / Auto-fallback only if `.gro` velocities are validated |

**Non-strict/no-flag fallback hardening**:
- If no checkpoint and no restart flags are set, pipeline now requires `prev_gro_path` validation.
- If `.gro` lacks velocities or cannot be validated, stage fails (no silent near-0K starts).
- If no `prev_gro_path` is available, stage fails with actionable guidance to provide one or use explicit flags.

```bash
# Regenerate velocities (fresh Maxwell-Boltzmann distribution):
python run_pipeline.py ... --allow-velocity-reset

# Expert restart from .gro with velocities (GRO must contain velocities):
python run_pipeline.py ... --allow-gro-velocity-restart
```

**Checkpoint-context guard**: Even with `--allow-gro-velocity-restart`, the combo `(gen_vel=no, continuation=no)` is an **ERROR** if a checkpoint is available (prevents accidentally dropping a checkpoint).

**`velocity_mode` recording**: The chosen velocity handling mode is recorded in `stage_state.json` for reproducibility:

| `velocity_mode` | Description |
|-----------------|-------------|
| `checkpoint` | Normal continuation from checkpoint |
| `velocity_reset` | Velocities regenerated from Maxwell-Boltzmann |
| `gro_velocity_restart` | Velocities read from .gro file |
| `error` | Configuration invalid, stage failed |

### Validation and Fail-Fast Behavior (Issue D)

The patcher validates after patching:

| Check | Behavior |
|-------|----------|
| Override status is true failure/missing template key | **FAIL** with explicit diagnostics |
| Override status is skipped with warning (non-strict safety downgrade) | **WARN + continue** |
| Vector length mismatch (tc-grps vs tau_t/ref_t) | **FAIL** |
| gen_vel=yes but gen_temp 鈮?ref_t | **FAIL** |
| nsteps 鈮?0 for production stage | **FAIL** |
| Post-write value mismatch | **FAIL** |

### GPE Default Time Scales

For science-quality results in gel polymer electrolyte (GPE) simulations, template defaults are:

| Stage | Default Duration | Rationale |
|-------|------------------|----------|
| NVT | 1 ns (500,000 steps) | Polymer segment relaxation |
| NPT | 1 ns (500,000 steps) | Density equilibration under stress |
| Production | 100 ns (50M steps) | Ion transport / diffusion / conductivity |

**Thermostat defaults (templates)**:
- `tcoupl = V-rescale`
- `tc-grps = System` in templates, but patcher default mode is `--tc-grps-mode auto`:
  split to `Polymer NonPolymer` when safe, else fallback to `System`
- `tau_t = 0.1` ps for NVT/NPT
- `tau_t = 1.0` ps for production MD (reduced overdamping of transport)

**Quick smoke tests** can override with shorter runs:
```bash
python run_pipeline.py ... --nsteps-nvt 5000 --nsteps-npt 5000 --nsteps-md 50000
```

**Integrator/LINCS guidance (templates)**:
- `dt = 0.002` with `constraints = h-bonds`
- `lincs_iter = 3`, `lincs_order = 8` (moderate hardening)
- If LINCS warnings persist: try `dt = 0.001` and/or increase `lincs_iter`

**Output cadence (templates)**:
- `nstxout = 0`, `nstvout = 0` (xtc-only trajectories)
- `nstxout-compressed = 1000` (2 ps at 2 fs)
- `nstenergy = 1000`, `nstlog = 1000` (2 ps)

### Barostat Strategy for Crosslinked Gels

> [!WARNING]
> Parrinello-Rahman barostat can blow up for far-from-equilibrium or high-stress crosslinked gel systems.

- **NPT equilibration** (default): Uses **Berendsen** barostat for stable density settling
- **Production MD**: Uses **Parrinello-Rahman** for correct NPT fluctuations
  - Template defaults: `tau_p = 2.0` (NPT Berendsen) and `tau_p = 5.0` (Production PR)
  - If pressure coupling is unstable or box oscillates, increase `tau_p` and/or reduce `compressibility`

**Valid barostat values**: `Berendsen`, `Parrinello-Rahman` only. C-rescale is a thermostat concept (not a barostat) and is not supported.

**tau_p Preservation**: When switching barostats, the patcher now **preserves existing tau_p** from the template if already set, and emits a warning if the default differs from what the barostat typically uses. tau_p is only auto-set if the template lacks it.

**Compressibility safety when switching to Parrinello-Rahman**:
- If template lacks `compressibility`, the patcher injects dimensionality-matched defaults based on `pcoupltype`:
- `liquid`: `4.5e-5`
- `gel`, `solid`, or unknown: `1.0e-6` (safer for crosslinked networks)
- Dimensionality is enforced: isotropic=`1`, semiisotropic=`2`, anisotropic=`6` values
- Existing `compressibility` is preserved (never silently overwritten)
- If `system_type` is `gel`/`solid` and existing `compressibility` is water-like (`~4.5e-5`), a warning is emitted about volume ringing/LINCS risk with remediation guidance

To override NPT barostat:
```bash
# Use Parrinello-Rahman for NPT (if system is already well-equilibrated):
python run_pipeline.py ... --npt-barostat Parrinello-Rahman
```

### MDP Key Alias Policy

GROMACS treats underscore and hyphen variants equivalently (e.g., `ref_t` = `ref-t`). The MDP patcher now:

1. **Detects existing key style** in the template before patching
2. **Preserves the template's style** to avoid appending duplicates
3. **Deduplicates semantic keys** before writing if both variants exist
4. **Override respects template style**: If template has `gen-temp=0` and you override with `gen_temp=300`, the patcher updates `gen-temp` (the existing key) to 300鈥攐verride is never lost

**Semantic key groups**: `gen_vel/gen-vel`, `gen_seed/gen-seed`, `gen_temp/gen-temp`, `tc_grps/tc-grps`, `tau_t/tau-t`, `ref_t/ref-t`, `ref_p/ref-p`, `tau_p/tau-p`, `comm_grps/comm-grps`, `comm_mode/comm-mode`

**Behavior**: If template uses `ref_t`, patching will edit that line in-place, not append `ref-t` at EOF.

**Duplicate alias handling**:
- Strict mode: **ERROR** with line numbers if multiple variants of the same semantic key exist
- Non-strict: **WARNING**, earlier duplicates removed and the **last** occurrence is kept (last-wins)

### Temperature Key Injection

When user explicitly sets `--temperature`:
- `ref_t` is updated (or **injected** if missing from template)
- `gen_temp` is updated to match `ref_t` if template has it
- If `gen_vel=yes` and `gen_temp` is missing, it is **injected**
- If `tc-grps` is split, `ref_t` is expanded to match the group count

This ensures GROMACS never silently uses defaults when user intends explicit temperature control.

### comm-grps Warning

When `comm-grps` length differs from `tc-grps` length and `comm-grps != System`:

> **WARNING**: This may cause COM drift / "flying ice cube" artifacts.
> **RECOMMENDED**: Set `comm-grps = System` OR match `comm-grps` to `tc-grps`.

This mismatch check is skipped when `comm-mode = None` (COM removal disabled).

### Strict vs Non-Strict MDP Validation

The patcher has two validation modes:

| Issue | Non-Strict (default) | Strict Mode |
|-------|---------------------|-------------|
| comm-grps vs tc-grps length mismatch | Warning | Error |
| comm-mode Angular with PBC | Override to Linear + warning | Error |
| gen_seed missing when gen_vel=yes | Warning | Error |
| gen_temp/ref_t thermal gradient | Warning | Error |
| tc-grps/tau_t/ref_t length mismatch | Error | Error |
| ref_p/compressibility dimensionality mismatch vs pcoupltype | Warning | Error |
| Unknown override key | Warning (adds new line) | Error (likely typo) |
| Duplicate semantic key variants in template | Warning (last wins) | Error (line numbers) |

Enable strict mode:
```bash
python run_pipeline.py ... --strict-mdp-validation
```

> [!NOTE]
> Final validation now ALWAYS runs on the output MDP, even when no patches were applied. This catches pre-existing template inconsistencies.

**gen_vel + ref_t spread**: When `gen_vel=yes` and `ref_t` has a thermal gradient (e.g., `300 400`), strict mode errors by default (non-strict warns). Rationale: `gen_temp` is a scalar; velocities are generated at **one** temperature, then the thermostat drives each group to its own `ref_t`, which can cause thermal shock and instability.

Fixes:
- Use `tc-grps-mode system` (single thermostat group)
- Use a single `ref_t` during equilibration and split later
- Set `gen_vel=no` and continue from a checkpoint
- Advanced: explicitly set `gen_temp` and set `ctx.allow_gen_vel_with_ref_t_spread = True` to acknowledge risk (strict mode downgrades to warning)

If `ref_t` values are consistent within tolerance, `gen_temp` must still match **all** values (no average-based pass).

### Hardening v2: Checkpoint-Aware Validation

The MDP patcher now validates gen_vel/continuation consistency based on **actual checkpoint availability**:

| checkpoint_available | gen_vel=no, continuation=no | Behavior |
|---------------------|----------------------------|----------|
| `False` | Expert restart from .gro | Allowed with `--allow-gro-velocity-restart` |
| `True` | Dangerous: would drop checkpoint | **ERROR** (even with `--allow-gro-velocity-restart`) |
| `None` | Unknown | Warns that checkpoint safety could not be verified (no silent `False` fallback) |

The `prepare_stage_mdp()` function now accepts `checkpoint_available` parameter, which stage runners should set based on actual `.cpt` file detection.

### Hardening v2: ref_p Dimensionality Validation

NPT pressure coupling now validates `ref_p` value count against `pcoupltype`:

| pcoupltype | Required ref_p values | Example |
|------------|----------------------|---------|
| `isotropic` | 1 | `ref_p = 1.0` |
| `semiisotropic` | 2 | `ref_p = 1.0 1.0` (xy, z) |
| `anisotropic` | 6 | `ref_p = 1.0 1.0 1.0 0 0 0` (xx yy zz xy xz yz) |

Mismatches are warnings by default, errors in strict mode.

### Hardening v2: force_gen_vel Consistency

When `--force-gen-vel` is used with NPT/MD stages, the patcher now **automatically sets continuation=no** to prevent the contradictory combination:

```bash
# OLD behavior (could create gen_vel=yes + continuation=yes = ERROR):
python run_pipeline.py ... --stage gmx_prod --force-gen-vel
# continuation would remain 'yes' from defaults

# NEW behavior (deterministic):
# --force-gen-vel 鈫?gen_vel=yes + continuation=no + gen_temp set
```

### Hardening v2: tc-grps tau_t Existence Check

When splitting tc-grps, the patcher now validates:
1. If `tcoupl != no`, `tau_t` must exist
2. `tau_t` count must match group count
3. Clear diagnostic: `"tc-grps split with tcoupl=v-rescale but tau_t is missing"`

### MDP Key Alias Normalization

The patcher now includes a reverse mapping (`VARIANT_TO_SEMANTIC`) for robust override handling:
- Overrides using `ref-t` will correctly modify template lines using `ref_t`
- No duplicate keys created when alias styles differ

### LINCS Tuning for Gel Systems

Crosslinked polymer networks can stress bond constraints. Templates use gel-safe defaults:

| Parameter | Default | Notes |
|-----------|---------|-------|
| `lincs_iter` | 2 | Increased from 1 for stressed crosslinks |
| `lincs_order` | 6 | Higher accuracy for gel systems |
| `lincs_warnangle` | 30掳 | Early warning for problematic geometries |

**If LINCS warnings occur frequently:**
1. Reduce `dt` to 1 fs (`--nsteps-xxx` values should double)
2. Increase `lincs_iter` to 4 and `lincs_order` to 8 in templates
3. Check for overly strained crosslink geometries in input structure

### Reproducibility: Velocity Generation Seed

NVT template uses a **fixed seed** (`gen_seed = 20260123`) for reproducible velocity generation.

To override:
```bash
python run_pipeline.py ... --gen-seed 42
```

The seed is recorded in the manifest under `mdp_patches.nvt`.



## Resource Scheduling / Pinning

Control CPU and GPU allocation for GROMACS:

```bash
# Use 8 threads for mdrun
python run_pipeline.py ... --gmx-nt 8

# Use specific GPU
python run_pipeline.py ... --gmx-gpu-id 0

# Set OpenMP threads
python run_pipeline.py ... --omp-num-threads 4
```

**Threading policy (deterministic):**
- `--omp-num-threads` always sets both `mdrun -ntomp` and `OMP_NUM_THREADS` (threads per rank).
- `--gmx-nt` sets `mdrun -nt` only for non-MPI (`thread-MPI`) GROMACS commands.
- If the resolved command looks like MPI GROMACS (e.g., `gmx_mpi`), `--gmx-nt` fails fast with a clear error. Use `mpirun/srun` for ranks and `--omp-num-threads` for threads per rank.

**GROMACS executable selection:**
- If a context override exists, it is used.
- Else if `GMX_CMD` is set (e.g., `GMX_CMD=gmx_mpi`), it is used.
- Otherwise defaults to `gmx`.
This resolution order is used consistently across GROMACS stages and PACKMOL PDB鈫扜RO conversion.

**GPU rule:**
- Normal command path uses `-gpu_id`.
- `GMX_GPU_ID` is only set in an env-fallback path (not when `-gpu_id` is used).

**HTPolyNet config injection:**
- Default (`--htpolynet-inject-add-missing` not set): only inject into pre-existing `gromacs.ntomp` / `gromacs.gpu_id` keys.
- With `--htpolynet-inject-add-missing`: create missing `gromacs`/`ntomp`/`gpu_id` config keys when needed.

These settings apply to all GROMACS stages and are recorded in the manifest.

## Resume & Checkpoint Logic

### Sentinel Files

Each stage writes `done.ok` as a final completion marker after required outputs + metadata are written.  
With `--resume` (default), completed stages are skipped only when `done.ok` and key artifacts/state files are coherent.

### GROMACS Checkpoint Resume

- **Same-stage resume** (e.g., continuing production):
  - Auto-detects `md.cpt` in stage directory
  - Uses `mdrun -cpi md.cpt -append`
  - Requires `md.mdp` + `stage_state.json` compatibility validation; `md.tpr` is required unless explicit repair is enabled
  - Resume compatibility gate ignores mutable runtime hashes (`md.cpt`/`md.tpr`) while preserving them in audit metadata
  - On successful resume, `stage_state.json` is refreshed for deterministic repeated resume
  - Optional: `--resume-ignore-runtime` can ignore runtime-only command-base changes

- **Cross-stage transition** (NPT 鈫?Production):
  - Uses `grompp -c npt.gro -t npt.cpt`
  - Preserves velocities from equilibration

- **Partial-state repair (production only)**:
  - `md.cpt` present + `md.tpr` missing: repairable with `--allow-resume-repair` (validate compatibility first, then rebuild `md.tpr`, then continue)
  - `md.tpr` present + `md.cpt` missing: not auto-repairable for continuation

### Force Re-run

Use `--force` to re-run a stage even if `done.ok` exists.

## Crash Context Preservation

On stage failure:

1. **Console output (stderr only)**: smart failure-centric excerpts, working directory, exact command
2. **crash_report.txt**: Written to failing stage directory with selected logs + selection reasons, excerpts, artifact summary, and env snapshot
3. **No cleanup**: All files preserved for debugging

For long-running `mdrun`, stdout/stderr are streamed to `logs/mdrun.stdout.log` and `logs/mdrun.stderr.log`.
When streamed commands fail, the raised error now includes bounded stdout/stderr tail excerpts and the same hint quality as non-streaming command paths.
Crash excerpt extraction uses a bounded progressive tail search (5MB cap): scan tail windows `200 -> 800 -> 2000` lines for failure markers (`Fatal error`, `LINCS`, `constraint`, `nan`, etc.). If found, include:
- Primary failure excerpt around the last marker (`~80` lines before + `~120` lines after)
- Trailing tail (`last 80 lines`) for exit context
If no markers are found, the report explicitly falls back to the last tail block.

Log selection now favors physical-failure signal over wrapper noise:
- Select up to 2 logs: a **primary** physics log and an optional **secondary** MPI/wrapper stderr log (clearly labeled noisy)
- Uses filename hints plus content scoring from a bounded tail slice (256KB)
- Penalizes wrapper-only tails (`mpirun has exited`, `MPI_Abort`, `srun:`, `PMI`) unless physical-failure markers are also present

Crash report writing remains atomic and MPI-safe: each write uses a unique temp name (`pid + uuid`), flush+`fsync`, then `replace`. Fallback uses `tempfile.gettempdir()` (portable, no hard-coded `/tmp`).
Log candidate discovery includes `.log`, `.out`, `.err`, `.txt`, `.mdp`, `.mdpout`, `.top`, and `.itp` (while still excluding large binary runtime outputs like `.xtc/.trr/.tpr/.cpt/.gro/.edr`).
If the tail must be truncated (e.g., >5MB read), the report will include a truncation note.

Example crash report location:
```
OUT_GMX/<RUN_ID>/03_gromacs/em/crash_report.txt
```

### Cleanup Policy (TRR)

- TRR deletion is **opt-in only** (`--cleanup-delete-trr`) and recorded in manifest cleanup records with `stage_tag` and `deleted_files`.
- If a `.trr` unlink fails, cleanup continues and manifest records `decision: "failed_to_delete"` with the file and error message.
- Cleanup still requires `--cleanup-success`.
- Stage-aware behavior:
  - `em`, `nvt`, `npt`: TRR is deleted when cleanup is enabled (equilibration TRR is disposable).
  - `md`: TRR is kept when `--analysis-requires-velocities` is set; otherwise deleted when cleanup is enabled.
    Production note: this is generally safe for coordinate-focused analysis using `.xtc`; keep TRR when velocity/force data are needed (for example Green-Kubo transport workflows).

## Analysis Stage

Computes radial distribution functions (RDF) and coordination numbers (CN) for defined atom pairs.

Reads pair definitions from `analysis_pairs.yaml`:

```yaml
pairs:
  - name: "Li-O_ether"
    group1: { mode: group, value: "Li" }                      # Exact group name (from index.ndx)
    group2: { mode: expr,  value: "resname PEO and name O*" } # GROMACS selection expression
    cn_cutoff: first_minimum_smoothed   # or numeric value in nm
    rdf:
      max_r: 1.2
      bin_width: 0.002
    cn:
      cutoff_nm: 0.35
    
rdf_settings:
  bin_width: 0.002    # nm
  max_r: 1.5          # nm (may be capped based on box size)
  
cn_settings:
  cutoff_nm: 0.35                       # Default cutoff if not per-pair
  # Optional smoothing parameters for first_minimum_smoothed detection:
  smooth_window: 5                      # Moving average window (points, odd preferred)
  smooth_window_nm: 0.01                # Physical smoothing width; overrides smooth_window when set
  peak_height_threshold: 0.08           # Peak must exceed baseline by this additive amount
  # peak_threshold: 1.2                 # Legacy factor form still accepted (interpreted as baseline*(factor-1))
  min_sep_nm: 0.05                      # Min separation between peak and min
  max_sep_nm: 0.40                      # Max separation after peak for first-shell minimum search
  min_depth: 0.05                       # Peak-to-minimum depth threshold (noise guard)
  recovery_window_nm: 0.08              # Post-minimum rise diagnostic window (advisory)
  recovery_height: 0.03                 # Rise threshold for diagnostic recovery flag
  recovery_positive_points: 4           # Positive-slope bins for diagnostic recovery flag
  local_window_nm: 0.00                 # Optional basin check window (0 disables)
  refine_points: 2                      # Optional smooth-only refinement mode (off by default)
  search_start_nm: 0.2                  # Start searching after this r
  baseline_window_nm: 0.10              # Baseline sampling window near search_start_nm
  truncation_tail_points: 3             # Descending tail points for truncation detection

analysis_settings:
  strict: true                          # Hard-fail on ambiguous/unsafe analysis
  gmx_timeout_s: 600                    # Timeout for all analysis-stage gmx calls
  traj_box_dt_ps: 10                    # Optional subsampling for trajectory box scan (gmx traj -dt)
```

### Filename Sanitization

Output filenames are sanitized and collision-safe:
- `safe_prefix = _safe_filename(pair.name)` keeps only `[A-Za-z0-9._-]` (path traversal tokens are neutralized)
- A stable hash suffix is appended: `<safe_prefix>__<hash8>`
- Outputs use `rdf_<safe_prefix>__<hash8>.xvg`, `cn_<safe_prefix>__<hash8>.xvg`, `cn_<safe_prefix>__<hash8>.dat`
- The hash is derived from pair identity (name + selections), so different pairs cannot silently overwrite each other
- If multiple pairs share the same readable prefix, analysis logs a warning and continues safely

### Boolean Parsing for `strict`

The `analysis_settings.strict` value is parsed robustly:
- Accepts: `true`, `false`, `yes`, `no`, `on`, `off`, `1`, `0` (case-insensitive strings)
- Accepts: Python `True`/`False` booleans
- Accepts: integers (`0` = false, non-zero = true)
- Default: `true` if unspecified or unrecognized

### Index.ndx Handling (Required for GROUP mode)

The analysis stage discovers index.ndx in deterministic priority:
1) manifest-provided path (if present and exists)
2) `IN/systems/<SYSTEM_ID>/gromacs/ndx/index.ndx`
3) exactly one `*.ndx` under `IN/systems/<SYSTEM_ID>/gromacs/ndx/`

If multiple `.ndx` files exist, analysis hard-fails with a list of candidates.

**Rule**: if a selection is `mode: group`, index.ndx **must** exist and the group name must be present, or analysis fails.
When index exists, `gmx rdf` always runs with `-n <index.ndx>`.

### Selection Rules (Group vs Expression)

Canonical form uses explicit mode:
- `mode: group` 鈫?exact index group name, resolved to `group "<name>"`
- `mode: expr` 鈫?raw GROMACS selection expression (passed as-is)

Backward compatibility: string `group1/group2` values are still accepted. Resolution is deterministic:
- If index exists and the string exactly matches a group name 鈫?treated as GROUP
- Else if the string clearly looks like an expression 鈫?treated as EXPR
- Else 鈫?**hard error** (ambiguous; use explicit mode or add group to index)

Examples:
```yaml
pairs:
  - name: "Li-TFSI"
    group1: { mode: group, value: "Li" }
    group2: { mode: group, value: "TFSI" }
```
```yaml
pairs:
  - name: "Sodium-O"
    group1: { mode: group, value: "Sodium ions" }   # spaces OK
    group2: { mode: expr, value: "name O*" }
```

### Strict Analysis Mode

`analysis_settings.strict` (default true) enforces fail-fast behavior to avoid misleading outputs:
- Missing trajectory or TPR 鈫?hard error
- Missing index for `mode: group` 鈫?hard error
- Ambiguous selection strings 鈫?hard error
- Unknown box size **and** no explicit `max_r` 鈫?hard error
- `cn_cutoff: first_minimum_smoothed` fails to detect a reliable minimum (including truncation) 鈫?hard error
- Any pair ending with `status: failed` or `status: skipped` 鈫?stage fails

If `strict: false`, analysis will:
- allow uncapped rmax **only** if `max_r` is explicitly provided
- skip pairs when box is too small, with `pair_skipped_reason` recorded
- fall back to `cn_settings.cutoff_nm` when first minimum is not found (with warning)
- set manifest `overall_status: completed_with_errors` when any pair fails/skips

`--force` only re-runs analysis outputs; it does **not** bypass missing mandatory inputs (trajectory, TPR, index for group mode).

### CN Computation Method

Coordination numbers are computed via **GROMACS native** `gmx rdf -cn`:
- No hardcoded density assumptions (previous code used 33 atoms/nm鲁 鈮?water)
- Produces physically meaningful CN for gel polymer electrolytes
- Both raw CN(r) curve (`cn_<safe_prefix>__<hash8>.xvg`) and scalar summary (`cn_<safe_prefix>__<hash8>.dat`) are saved
- In `--dry-run`, commands are logged but RDF/CN files are intentionally not generated

### rmax Box Capping

To prevent periodic artifacts / double counting, `rmax` is automatically capped:
- Triclinic-safe basis: **minimum box altitude** (minimum opposite-face distance), not `min(|a|,|b|,|c|)`
- `safe_rmax_tpr = 0.5 脳 min_box_height_tpr - 0.02`
- `safe_rmax_traj = 0.5 脳 min_box_height_traj - 0.02` (from `gmx traj -ob`)
- `safe_rmax_used = min(safe_rmax_tpr, safe_rmax_traj)` when both are available
- `effective_rmax = min(requested_max_r, safe_rmax_used)`
- If capped, a warning is printed and both requested/effective values are recorded with `safe_rmax_basis: min_altitude`
- **Invariant**: effective_rmax never exceeds safe_rmax_used (no minimum override)

### Outputs

- RDF: `OUT_GMX/<RUN_ID>/analysis/rdf/rdf_<safe_prefix>__<hash8>.xvg`
- CN curve: `OUT_GMX/<RUN_ID>/analysis/cn/cn_<safe_prefix>__<hash8>.xvg`
- CN summary: `OUT_GMX/<RUN_ID>/analysis/cn/cn_<safe_prefix>__<hash8>.dat`

### Cutoff Detection

When `cn_cutoff: first_minimum_smoothed` (or legacy alias `first_minimum`), the pipeline uses a robust first-shell heuristic:
1. Apply fixed-width symmetric smoothing (reflect-style edge handling; no edge window shrink)
2. Build baseline near `search_start_nm` (if empty, fallback to 1.0 and mark `baseline_unreliable`)
3. Find first significant peak (above baseline + `peak_height_threshold`)
4. Validate minimum candidates with depth (`min_depth`), minimum separation (`min_sep_nm`), and max separation (`max_sep_nm`)
5. Optionally enforce local-basin robustness (`local_window_nm`)
6. Recovery checks (`recovery_window_nm`/`recovery_height` or `recovery_positive_points`) are advisory reliability flags (not hard gating)
7. `refine_points` is explicit opt-in; default is no refinement, and selected index remains the smoothed minimum (no raw argmin noise snapping)
8. No global-minimum fallback is used by default
9. If no local minimum exists after the first peak (shoulder/plateau), return `no_first_minimum` with `shoulder_no_minimum`
10. Detect truncation: if the post-peak curve is still descending at the RDF tail, report `truncated_before_first_minimum`
11. If still unreliable:
   - strict=true 鈫?hard error
   - strict=false 鈫?fall back to `cn_settings.cutoff_nm` and record `cutoff_unreliable: true`

If `truncated_before_first_minimum` is reported:
- if `rmax` is **not** box-capped, increase `rdf.max_r`/`search_rmax` to capture the first-shell valley.
- if `rmax` **is** box-capped, use explicit `cn_cutoff`, analyze after longer/stabler NPT, or use a larger box/system.

### Manifest Provenance

Per-pair analysis results include:
- `resolved_ref`, `resolved_sel`: actual selections sent to gmx
- `output_basename`: deterministic collision-safe basename (`<safe_prefix>__<hash8>`)
- `rdf_settings`: bin_width, requested/effective rmax, rmax cap source, `safe_rmax_basis`, `min_box_height_tpr_nm`, `min_box_height_traj_nm`, `safe_rmax_tpr`, `safe_rmax_traj`, `safe_rmax_used` (legacy `min_box_tpr/min_box_traj` retained for compatibility)
- `cn_method`: `"gmx_rdf_native"` or `"dryrun"`
- `cutoff_mode`: `"fixed"` or `"first_minimum_smoothed"`
- `cutoff_method`, `peak_r`, `peak_g`, `cutoff_r`, `r_min`, `g_rmin_raw`, `g_rmin_smooth`, `cutoff_reason`
- `cutoff_nm`, `cutoff_unreliable`, `baseline_unreliable`, `reliability_flags`, `quality_flags`, `fallback_used`, `smoothing_params`, `warnings`
- `output_files`: paths to rdf_xvg, cn_xvg, cn_dat

## CLI Options

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `--ff` | `GAFF2`, `OPLS-AA` | Required | Forcefield (OPLS-AA internally normalized to OPLSAA) |
| `--charge` | `RESP`, `CM5` | Required | Charge method |
| `--allow-charge-model-mismatch` | Flag | False | Unsafe override for unsupported ff/charge pair (requires reason) |
| `--charge-model-mismatch-reason` | String | None | Required audit reason when unsafe ff/charge override is enabled |
| `--system` | `<SYSTEM_ID>` | Required | System ID for path resolution |
| `--run-id` | `<RUN_ID>` | Required* | Run ID for output directory (*unless `--allow-auto-run-id`) |
| `--allow-auto-run-id` | Flag | False | Allow auto-generated run ID if `--run-id` not provided |
| `--project-root` | Path | Auto-detected | Override project root (highest priority). Default walks up from `pipeline/cli.py` and picks the first parent containing `run_pipeline.py`, `.git`, or `IN/`; falls back to legacy detection if none are found. |
| `--stage` | `packmol`, `htpolynet`, `sanitizer`, `gmx_em`, `gmx_eq`, `gmx_prod`, `analysis`, `all` | Required | Stage to run |
| `--qc-policy` | `off`, `warn`, `error` | Auto (`warn` for `--stage all`, else `off`) | Dispatcher QC gate behavior on stage metrics |
| `--qc-density-rel-tol` | Float (0-1) | 0.05 | Relative tolerance for density-vs-expected QC check |
| `--qc-enable-after-stages` | CSV stage list | `gmx_eq` | Stages where dispatcher attempts QC if `stage_metrics.json` exists |
| `--all-qc-stop-on-warn` | Flag | False | In `--stage all`, stop early after QC warning |
| `--dry-run` | Flag | False | Generate without execution |
| `--force` | Flag | False | Force re-run even if done.ok exists (archives prior outputs) |
| `--no-resume` | Flag | False | Do not skip any stages (default: resume=True) |
| `-T`, `--temperature` | Float (K) | None | Temperature override |
| `-P`, `--pressure` | Float (bar) | None | Pressure override |
| `--tc-grps-mode` | `auto`, `system`, `split` | `auto` | Thermostat groups mode (`auto` is conservative/system-first; explicit `system` is safest) |
| `--tau-t-polymer` | Float (ps) | None | tau_t for Polymer group in legacy Polymer/NonPolymer split mode |
| `--tau-t-nonpolymer` | Float (ps) | None | tau_t for NonPolymer group in legacy Polymer/NonPolymer split mode |
| `--tau-t-values` | CSV floats | None | tau_t list for custom split groups (count must match `--tc-grps-groups`) |
| `--tc-grps-min-atoms` | Integer | 50 | Minimum atoms per split thermostat group before fallback to `System` |
| `--tc-grps-min-fraction` | Float (0-1) | 0.01 | Minimum atom fraction per split thermostat group before fallback to `System` |
| `--allow-small-tc-grps` | Flag | False | Expert opt-in to allow tiny split thermostat groups |
| `--allow-tau-t-autofill` | Flag | False | Allow duplicating scalar tau_t (legacy, not recommended for GPE) |
| `--comm-grps-policy` | `auto`, `match_tc`, `system`, `none`, `explicit` | `match_tc` | COM group policy (auto is risk-aware; none disables COM removal) |
| `--comm-mode` | `Linear`, `Angular`, `None` | `Linear` | COM removal mode |
| `--allow-comm-mode-angular-with-pbc` | Flag | False | Expert opt-in to allow Angular COM mode under PBC |
| `--comm-grps` | String | None | Explicit comm-grps (only valid with policy=explicit) |
| `--force-gen-vel` | Flag | False | Force gen_vel=yes in production (not recommended) |
| `--nsteps-em/nvt/npt/md` | Integer | None | Per-stage nsteps override (`> 0` when provided) |
| `--gmx-eq-metrics-window-ps` | Float (ps) | 20.0 | Averaging window for `gmx_eq` stage thermodynamic metrics |
| `--gmx-nt` | Integer | None | mdrun `-nt` total threads (`> 0` when provided; non-MPI/thread-MPI mode only; MPI builds fail fast) |
| `--gmx-gpu-id` | String | None | mdrun GPU ID(s), e.g., '0' or '0,1,2' (no spaces) |
| `--omp-num-threads` | Integer | None | OpenMP threads per rank (`> 0` when provided; sets `-ntomp` and `OMP_NUM_THREADS`) |
| `--htpolynet-inject-add-missing` | Flag | False | Allow HTPolyNet config injection to create missing `gromacs.ntomp/gpu_id` keys |
| `--cleanup-success` | Flag | False | Remove .trr after success |
| `--allow-override` | Flag | False | Allow atomtype conflicts (priority-based selection) |
| `--atomtype-prefix` | String | None | Prefix atomtype names |
| `--prefix-implicit-topology-policy` | `error`, `warn`, `bake`, `inject` | `error` | Prefix safety policy for implicit/library topology lookup |
| `--prefix-implicit-topology-reason` | String | None | Required audit reason when prefix policy is `bake` or `inject` |
| `--strict-charge` | Flag | False | Fail on atomtype charge mismatches |
| `--allow-unsanitized-grompp` | Flag | False | Allow grompp without sanitizer outputs (NOT recommended) |
| `--charge-fix-method` | `safe_subset`, `uniform_all` | `safe_subset` | Correction method: safe-subset or uniform |
| `--charge-fix-target-allowlist` | String | None | Comma-separated allowed targets (e.g., "PC,EC") |
| `--charge-fix-protect-resnames` | String | None | Additional protected residue/molecule names for charge-fix classification |
| `--charge-fix-allow-solvents` | Flag | False | Allow correction on protected solvent targets |
| `--charge-fix-allow-ions` | Flag | False | Allow correction on protected ion targets |
| `--polymer-net-charge-tol` | Float | 1e-3 | Acceptable protected-polymer per-molecule drift before correction is required |
| `--charge-fix-polymer-method` | `skip_if_small`, `spread_safe` | `skip_if_small` | Protected-polymer policy and correction method |
| `--charge-fix-polymer-exclusion-bonds` | Integer | 2 | Bond-distance exclusion from hetero atoms for polymer spread-safe correction |
| `--charge-fix-max-delta-per-atom` | Float | 1e-5 | Max |螖q| per atom for charge correction |
| `--charge-fix-max-total` | Float | 1e-4 | Max |螖Q| total correction for neutrality |
| `--charge-fix-max-dipole-shift` | Float | 0.01 | Max |螖渭| dipole shift (Debye) |
| `--charge-fix-moleculetype-rounding-tol` | Float | 1e-4 | Max `|Q_mol-round(Q_mol)|` treated as rounding drift |
| `--charge-fix-allow-non-rounding` | Flag | False | Unsafe override to allow correction for non-rounding drift |
| `--strict-charge-physics` | Flag | False | Deprecated in checker path; use `--strict-dipole-check` for dipole hard-fail enforcement |
| `--strict-gro-top-check` | Flag | False | Fail when GRO/TOP check is uncertain |
| `--allow-default-defaults` | Flag | False | Allow fallback [ defaults ] when missing |
| `--allow-mixed-defaults` | Flag | False | Unsafe override to allow mixed [ defaults ] tuples |
| `--allow-mixed-defaults-reason` | String | None | Required audit reason when `--allow-mixed-defaults` is enabled |
| `--mixed-defaults-cross-group-policy` | `off`, `warn`, `generate` | `warn` with mixed-default override, else `off` | Cross-group LJ handling policy for implemented comb-rule-only mixed defaults (`nbfunc=1`, canonical comb-rule 2/3) |
| `--mixed-defaults-cross-group-rule` | `lorentz-berthelot`, `geometric` | Auto from canonical comb-rule | Sigma mixing rule used when cross-group generation is enabled |
| `--mixed-defaults-cross-group-max-pairs` | Integer | 20000 | Hard cap for generated cross-group pair count |
| `--mixed-defaults-cross-group-reason` | String | None | Required audit reason when cross-group policy is `generate` |
| `--mixed-defaults-preserve-within-group` | Flag | False | Guard rail: generate within-secondary-group `nonbond_params` for comb-rule-only mixed defaults |
| `--lj-outlier-policy` | `off`, `warn`, `error` | `warn` | LJ outlier action policy (`unit-screaming` thresholds always error) |
| `--lj-bounds-profile` | `aa`, `coarse_grained`, `polarizable`, `custom` | `aa` | LJ bounds profile for sigma/epsilon outlier checks |
| `--lj-sigma-max-nm` | Float | None | Custom sigma max (required with `--lj-bounds-profile custom`) |
| `--lj-epsilon-max-kj-mol` | Float | None | Custom epsilon max (required with `--lj-bounds-profile custom`) |
| `--strict-forcefield-consistency` | Flag | False | Fail on mixed forcefield atomtype conflicts |
| `--itp-include-dirs-priority` | `forcefield_first`, `sanitized_first`, `first`, `last` | `forcefield_first` | Include search precedence (`last`/`first` are legacy aliases) |
| `--allow-unsafe-include-escape` | Flag | False | Unsafe override: allow `#include` paths to resolve outside configured roots |
| `--strict-top-update` | Flag | False | Fail if sanitizer managed block markers in `system.top` are malformed |
| `--dipole-group-max-atoms` | Integer | 2048 | Skip dipole work for oversized residue groups before unwrap |
| `--dipole-percolation-threshold` | Float (0-1) | 0.45 | Fallback span-ratio threshold for non-bonded unwrap path |
| `--allow-placeholder` | Flag | False | Allow non-physical placeholder outputs (demo only) |
| `--allow-placeholder-propagate` | Flag | False | Allow publishing placeholder trio to `IN/.../htpolynet` (demo only) |
| `--allow-placeholder-stage-to-gromacs` | Flag | False | Explicit opt-in to stage propagated placeholders to `IN/.../gromacs` |
| `--allow-placeholder-gromacs-compile` | Flag | False | Disable placeholder poison pill that intentionally breaks grompp (unsafe) |
| `--allow-triclinic-unsafe-diagonal-approx` | Flag | False | Unsafe triclinic escape hatch: diagonal approximation instead of fail-fast |
| `--skip-gelation-precheck` | Flag | False | Skip conservative gelation feasibility precheck before HTPolyNet |
| `--allow-gelation-precheck-fail` | Flag | False | Downgrade failed gelation precheck to warning and continue |
| `--allow-missing-htpolynet-qc` | Flag | False | Allow successful run even when conversion+gel fraction cannot be parsed |
| `--min-conversion` | Float (0-1) | None | Minimum crosslinking conversion threshold |
| `--htpolynet-timeout` | Integer (s) | 43200 | HTPolyNet subprocess timeout (`0` disables timeout) |
| `--unsafe-allow-out-fallback` | Flag | False | Allow reading initial structure from `OUT_GMX` when IN/published asset is missing |
| `--allow-partial-publish` | Flag | False | Downgrade immediate strict publish raise, but stage is still non-complete when both PDB+GRO are not updated |
| `--allow-density-reduction` | Flag | False | Explicitly allow density-backoff retries (still floor/cap constrained and audited) |
| `--packmol-pdb-scale` | `0.1`, `1.0` | None | Explicit PACKMOL PDB scale override (脜鈫抧m or nm鈫抧m) |
| `--no-strict-packmol-units` | Flag | False | Relax non-scale checks; ambiguous unit scale still requires explicit `--packmol-pdb-scale` |
| `--strict-gro-conversion` / `--no-strict-gro-conversion` | Flag pair | False | Fail-fast on missing/failed `editconf` conversion command (strict) or allow non-strict conversion warnings |
| `--packmol-edge-margin` | Float (nm) | 0.2 | Placement edge margin for PACKMOL inside-box constraints |
| `--packmol-preassembly-mode` | `none`, `li_near_polymer`, `li_solvent_shell` | `none` | Optional Li preassembly bias mode |
| `--packmol-preassembly-li-fraction` | Float (0-1) | 1.0 | Fraction of Li molecules targeted by preassembly hints |
| `--packmol-preassembly-retry-count` | Integer | 2 | Max preassembly repair retries (seed/hint adjustments) |
| `--packmol-preassembly-sample-cap` | Integer | 4000 | Atom sampling cap per role for preassembly metrics |
| `--packmol-li-polymer-cutoff-nm` | Float | 0.35 | Li-polymer oxygen proximity cutoff |
| `--packmol-li-solvent-cutoff-nm` | Float | 0.35 | Li-solvent oxygen proximity cutoff |
| `--packmol-li-tfsi-cutoff-nm` | Float | 0.40 | Li-TFSI close-contact proxy cutoff |
| `--packmol-target-li-polymer-fraction` | Float (0-1) | 0.40 | Target Li-polymer fraction for `li_near_polymer` mode |
| `--packmol-target-li-solvent-fraction` | Float (0-1) | 0.50 | Target Li-solvent fraction for `li_solvent_shell` mode |
| `--packmol-max-li-tfsi-close-fraction` | Float (0-1) | 0.30 | Max tolerated Li-TFSI close-contact fraction |
| `--gen-seed` | Integer | None | Deterministic seed for velocity generation (NVT) |
| `--npt-barostat` | `Berendsen`, `Parrinello-Rahman` | `Berendsen` | Barostat (pcoupl) for NPT equilibration |
| `--grompp-maxwarn` | Integer | 0 | Maximum warnings allowed by grompp (`>= 0`) |
| `--allow-velocity-reset` | Flag | False | Allow proceeding without checkpoint (mutually exclusive with `--allow-gro-velocity-restart`) |
| `--allow-resume-repair` | Flag | False | Allow repairing production resume when `md.cpt` exists but `md.tpr` is missing |
| `--strict-restart-policy` | Flag | Inherit `--strict-mdp-validation` | Enforce strict checkpoint/restart policy independently of MDP validation strictness |
| `--resume-ignore-runtime` | Flag | False | Ignore runtime-only resume fingerprint keys (e.g., `mdrun_command_base`) |
| `--strict-mdp-validation` | Flag | False | Escalate MDP validation warnings to errors |
| `--allow-gro-velocity-restart` | Flag | False | Allow expert restart (gen_vel=no + continuation=no) only when `.gro` velocities are valid (mutually exclusive with `--allow-velocity-reset`) |
| `--allow-continuation-override` | Flag | False | Skip auto-setting continuation for stage policy |
| `--posres-reference` | Path | None | Explicit POSRES reference structure passed to grompp (-r) |
| `--no-strict-posres-reference` | Flag | False | Allow sticky auto-pinned POSRES reference under `OUT_GMX/<RUN_ID>/03_gromacs/_shared/` |

## Manifest Structure

`OUT_GMX/<RUN_ID>/manifest.json` contains:
- **Schema metadata**: `schema_version` + `version` (light migration keeps old keys readable)
- **Options**: ff, charge, system_id, run_id, stage, T, P
- **Tool versions**: packmol, gromacs, htpolynet (summary strings in manifest)
- **Extended provenance**: raw version probes + parsed GROMACS build features + runtime/platform fields + git code state + capped probe previews in `run_report.json` and `provenance.txt`
- **Composition**: immutable `initial_counts`, explicit `post_counts`, history/conflicts, MW, wt%, box vectors, density, unit diagnostics, optional preassembly metrics
- **HTPolyNet**: config, rules file, published files, staged files
- **Sanitizer**: combined_atomtypes, sanitized ITPs, conflicts, defaults used, charge warnings
- **Charge neutrality**: computed Q, correction applied, patched ITP path, `effective_method` + `method_desc`, alias target + real ITP moleculetype target, 螖q/螖Q guardrails, optional atom-level correction audit (`line_no`, `charge_col_idx`, `line_excerpt`, capped records, truncation counters), and optional unparsed `[ atoms ]` line audit when explicit unsafe override is used
- **MDP patches**: original and patched values per stage
- **Resources**: gmx-nt, gmx-gpu-id, omp-num-threads
- **Checkpoints**: paths and methods used
- **Analysis**: pair definitions, cutoffs used
- **Commands**: structured command audit with `argv` + `cmd`, return code, capped stdout/stderr previews, and optional `duration_s`/`cwd`/env subset
- **Errors**: crash_report.txt path if failure
- **Recovery metadata**: corruption recovery trace when a bad manifest is quarantined (`recovered_from_corrupt_manifest`, `corrupt_backup_path`, `recovered_at_utc`)

## Requirements

- Python 3.8+ with PyYAML
- PACKMOL
- HTPolyNet
- GROMACS (gmx, gmx_mpi). To use `gmx_mpi`, set `GMX_CMD=gmx_mpi`.

## License

Internal use.
