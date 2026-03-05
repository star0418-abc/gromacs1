# 06_ITP_SANITIZER.md — ITP Sanitizer

## Purpose
The ITP Sanitizer is the topology / forcefield sanitization layer used before `grompp` and downstream MD execution.

Its job is to make topology inputs:
- deterministic
- auditable
- fail-fast when ambiguity or physical conflict would otherwise produce silent corruption

The sanitizer must protect against:
- conflicting `[ atomtypes ]` / `[ defaults ]`
- inconsistent duplicate `[ moleculetype ]`
- unsafe include shadowing
- stale or ambiguous ITP selection
- silent topology drift
- unsafe charge-neutrality repair

---

## Module Boundaries

### `sanitizer_stage.py`
Owns orchestration only.

Responsibilities:
- stage entrypoint
- sequencing
- invoking topology and spatial checks
- aggregating logs / QC / stage results
- deciding success / failure

Must **not** accumulate low-level topology parsing or geometry logic again.

### `topology_sanitizer.py`
Owns topology / forcefield sanitization and final topology writes.

Responsibilities:
- include discovery and recursive include closure
- `[ moleculetype ]` collection and deterministic winner selection
- `[ atomtypes ]` / `[ defaults ]` parsing and conflict checking
- sanitized ITP generation
- `combined_atomtypes.itp` generation
- corrected ITP syncing / atomic replacement / backup handling
- sanitizer-managed `system.top` block generation
- final topology-side output writing
- charge-neutrality correction planning and write-side application

This is the **single writer of final topology artifacts**.

### `spatial_checker.py`
Owns spatial / geometry / PBC / dipole / `grompp`-assisted consistency checks.

Responsibilities:
- GRO parsing
- filtered / streaming coordinate access
- PBC unwrap and bond-graph consistency checks
- dipole / charge-distribution validation
- `grompp -pp`-based topology truth extraction
- GRO/TOP consistency checks
- returning check results, reliability flags, and patch plans

This module must **not** become the final writer of topology artifacts.

### `sanitizer.py`
Compatibility facade only.

Responsibilities:
- preserve old import paths
- re-export public sanitizer symbols

Do not put new primary logic here.

---

## Inputs
Expected inputs include:
- molecule ITPs staged into the sanitizer work area
- forcefield ITPs and required include trees
- `system.top`
- pipeline/context include search paths
- sanitizer policy flags from context / CLI / environment
- optional charge-neutrality protection metadata

The sanitizer must treat the full include-search context as authoritative.

---

## Outputs
The sanitizer produces deterministic sanitized topology artifacts, including:
- sanitized molecule ITPs
- `combined_atomtypes.itp`
- updated `system.top`
- versioned backups for corrected ITPs when applicable
- manifest / QC / log metadata describing what changed and why
- charge-neutrality provenance output when charge correction is applied

### Winner-Only Output Rule
If multiple candidate files compete for the same include target or logical molecule definition, the sanitizer must resolve a deterministic winner and materialize **winner-only** sanitized outputs.

Provenance belongs in logs / manifests, not as multiple competing live outputs.

---

## Core Policy

### A) Fail fast on physical ambiguity
If topology or forcefield definitions conflict in a way that makes simulation behavior ambiguous, fail loudly.

Examples:
- incompatible `[ defaults ]`
- conflicting physical LJ parameters for the same atomtype
- incompatible duplicate `[ moleculetype ]` definitions
- unresolved charge-neutrality repair with no safe eligible atoms

### B) Deterministic winner selection
Any tolerated override must be:
- deterministic
- explainable
- logged

Never depend on incidental filesystem order.

### C) Loud provenance
Whenever a winner is chosen or a conflict is tolerated, record:
- competing sources
- chosen winner
- reason / policy path
- whether the decision was `exact`, `tolerated`, `uncertain`, or `unreliable`

### D) No silent downgrade
If the sanitizer cannot safely prove equivalence, it must not silently pretend equivalence exists.

---

## Include Discovery and Resolution
The sanitizer must recursively resolve `#include` directives using the full context-derived include search path.

Rules:
1. build include search paths from pipeline/context-defined roots
2. resolve includes recursively from all active topology inputs
3. preserve deterministic search order
4. if include shadowing is allowed, resolve it explicitly by policy and log the winner
5. if include shadowing is not allowed, fail on actual ambiguity

Do not reject same-basename files prematurely if policy-defined resolution is expected to choose a winner later.

---

## `[ moleculetype ]` Conflict Handling
The sanitizer must identify all `[ moleculetype ]` definitions and determine whether duplicates are:
- exactly equivalent
- safely tolerable under policy
- physically conflicting
- uncertain due to incomplete parsing

### Signature policy
Use deterministic signatures as evidence, not as infallible truth.

If signatures differ, determine whether the difference is:
- physically real
- namespace-only
- caused by uncertain parsing
- caused by unsupported formatting or preprocessing

### Uncertain charge extraction
If charge extraction is affected by preprocessing directives or ambiguous formatting:
- mark the comparison as uncertain
- try stronger secondary checks
- declare a hard difference only when evidence supports it

---

## `[ atomtypes ]` and `[ defaults ]`

### `[ defaults ]`
- must remain globally coherent
- incompatible values are a hard error unless an explicit policy says otherwise

### `[ atomtypes ]`
- physical LJ conflicts must fail loudly
- exact duplicates may be tolerated
- tolerated override paths, if allowed, must be deterministic and fully logged

### Parsing rule
- prefer explicit tokenization first
- remain robust against nonstandard numeric formatting
- if parsing remains ambiguous:
  - raise in strict mode
  - emit explicit warning and uncertainty in tolerant mode

Do not silently read the wrong sigma / epsilon columns.

---

## Charge Neutrality Repair Policy

### A) When charge repair is allowed
If a residual total charge `Δq_total` remains after topology normalization, rounding, or deterministic repair steps, the sanitizer may apply a **small deterministic redistribution** only if:
- the residual is within configured tolerance for repair
- the eligible atom pool is non-empty
- the redistribution plan is fully logged

Otherwise, fail fast.

### B) Protected atom / site policy
Charge repair must avoid chemically sensitive sites.

The sanitizer must support a **protected-site mechanism** with deterministic selectors. Preferred selector types:
- explicit protected atom indices
- protected atom names within a molecule template
- protected atom roles / metadata tags
- protected atomtypes as a coarse fallback only

Typical protected classes include:
- formally charged atoms
- redox-active centers
- reactive / crosslinking atoms
- strongly coordinating heteroatoms
- user-tagged protected sites

Do not redistribute residual charge blindly over all atoms.

### C) Eligible atom pool
The repair plan must build an eligible atom pool by excluding protected atoms / sites.

If no eligible atoms remain:
- fail fast
- report why the repair could not be applied

### D) Per-atom correction limits
The sanitizer must enforce a configurable cap on absolute per-atom repair magnitude.

If the required redistribution would exceed the cap:
- fail fast in strict mode
- optionally fail with actionable guidance in tolerant mode

Do not silently distort chemically important local charges.

### E) Provenance output
Any charge-neutrality repair must produce explicit provenance, recorded in manifest and/or sidecar report, including:
- `delta_q_total_e`
- redistribution policy name
- protected selector source
- protected atom count
- eligible atom count
- per-molecule summary
- per-atom or per-site allocation summary
- maximum absolute delta applied to any atom
- whether the result was `exact`, `tolerated`, or `uncertain`

The total applied correction must satisfy:
- `sum(delta_q_i) = -delta_q_total_e`

This must be auditable.

---

## Managed `system.top` Updates
`topology_sanitizer.py` owns final `system.top` writing.

Responsibilities:
- validate required includes
- update sanitizer-managed include blocks
- maintain deterministic include order
- preserve user content outside the managed block

Preferred behavior:
- if a managed placeholder exists, use it
- otherwise insert the managed block at the deterministic default location defined by policy

Do not scatter final topology writes across multiple modules.

---

## Corrected ITP Syncing
When corrected ITP content must replace an original file:
1. preserve the original file as a true backup first
2. atomically replace the live file
3. keep the operation crash-safe and auditable

Use temp-file + fsync + `os.replace` style semantics where supported.

The versioned backup must always represent the pre-modification original.

---

## Reliability and Uncertainty
The sanitizer must distinguish between:
- `exact`
- `tolerated`
- `uncertain`
- `unreliable`

Never collapse `uncertain` or `unreliable` into `exact`.

---

## Strict vs Tolerant Behavior
Default behavior should be **strict** unless the user explicitly enables tolerant mode.

### Strict mode
- fail on unresolved ambiguity
- fail on unsupported conflicting definitions
- fail on ambiguous physical parsing
- fail on unsafe charge-neutrality repair

### Tolerant mode
- may proceed with deterministic winner selection
- must record warnings and uncertainty
- is **not** permission for silent corruption

---

## Anti-Patterns
Do not implement or reintroduce:
- topology writing from multiple modules without a single writer
- hidden winner selection based on incidental filesystem order
- broad regex-only parsing that can silently misread physical parameters
- declaring `different` solely because parsing encountered `#ifdef` / `#include`
- materializing multiple competing winners into live sanitized outputs
- blind charge redistribution over chemically sensitive atoms
- letting `sanitizer.py` become a new primary logic home

---

## Maintenance Guidance
When modifying sanitizer logic:
1. preserve module boundaries
2. keep deterministic behavior
3. prefer explicit provenance over implicit heuristics
4. make uncertainty visible
5. keep final topology writes centralized in `topology_sanitizer.py`
6. keep `sanitizer_stage.py` thin
7. keep charge-neutrality repair policy auditable and chemistry-aware