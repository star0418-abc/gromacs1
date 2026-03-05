# 06 ITP Sanitizer

## Purpose

The ITP Sanitizer is the topology/forcefield sanitization layer used before `grompp` and downstream MD execution. Its job is to make topology inputs deterministic, auditable, and fail-fast when forcefield or molecule-definition conflicts would otherwise produce silent corruption.

The sanitizer must protect against:

- conflicting `[ atomtypes ]` / `[ defaults ]` definitions
- inconsistent duplicate `[ moleculetype ]` definitions
- unsafe include shadowing
- accidental use of stale or ambiguous ITPs
- silent topology drift between staged inputs and generated outputs

It must produce sanitized topology artifacts that are safe to consume by later pipeline stages.

---

## Module Boundaries

The sanitizer is now split into three layers. Any future edits must respect these ownership boundaries.

### `sanitizer_stage.py`
Owns orchestration only.

Responsibilities:
- stage entrypoint
- stage sequencing
- calling topology and spatial checks
- aggregating logs, QC, and stage results
- deciding stage success/failure

Must **not** become the place where low-level topology parsing, GRO parsing, or geometry logic accumulates again.

### `topology_sanitizer.py`
Owns topology/forcefield sanitization and final topology writes.

Responsibilities:
- ITP/include discovery and recursive include closure
- `[ moleculetype ]` collection and deterministic winner selection
- `[ atomtypes ]` / `[ defaults ]` parsing and conflict checking
- sanitized ITP generation
- combined atomtypes generation
- corrected ITP syncing / atomic replacement / backup handling
- `system.top` managed block generation and final include updates
- final topology-side output writing

This module is the **single writer of final topology artifacts**.

### `spatial_checker.py`
Owns spatial / geometry / PBC / dipole / `grompp`-assisted consistency checks.

Responsibilities:
- GRO parsing
- filtered/streaming coordinate access
- PBC unwrap and bond-graph consistency checks
- dipole / charge-distribution validation
- `grompp -pp`-based topology truth extraction
- GRO/TOP consistency checks
- returning check results, reliability flags, and patch plans

This module must **not** become the final writer of topology artifacts. It returns results to be applied by higher-level topology logic.

### `sanitizer.py`
Compatibility facade only.

Responsibilities:
- preserve old import paths
- re-export public sanitizer symbols required by the rest of the repo

Do not put new primary logic here.

---

## Inputs

Expected inputs include:

- molecule ITPs staged into the sanitizer work area
- forcefield ITPs and any required include trees
- `system.top`
- pipeline/context include search paths
- sanitizer policy flags from context/CLI/environment

The sanitizer must treat the full include-search context as authoritative. It is not limited to forcefield-local includes only.

---

## Outputs

The sanitizer produces deterministic sanitized topology artifacts, typically including:

- sanitized molecule ITPs
- `combined_atomtypes.itp`
- updated `system.top`
- versioned backups for corrected ITPs when applicable
- manifest/QC/log metadata describing what changed and why

### Winner-Only Output Rule

If multiple candidate files compete for the same include target or logical molecule definition, the sanitizer must resolve a deterministic winner and materialize **winner-only** sanitized outputs.

Do **not** write multiple same-target candidates into sanitized output directories merely to preserve provenance. Provenance belongs in logs/manifests, not as competing live topology files.

---

## Core Policy

### A. Fail fast on forcefield corruption
If forcefield-level definitions conflict in a way that would make simulation behavior ambiguous, fail loudly.

Examples:
- incompatible `[ defaults ]`
- conflicting physical LJ parameters for the same atomtype
- incompatible duplicate `[ moleculetype ]` definitions that cannot be safely reconciled

### B. Deterministic winner selection
When policy permits one definition to override another, the override must be:

- deterministic
- explainable
- logged in manifest/QC output

Selection must never depend on filesystem nondeterminism.

### C. Loud provenance
Whenever a winner is chosen or a conflict is tolerated by policy, record:
- competing sources
- chosen winner
- reason/policy path
- whether the decision was exact, tolerant, or uncertain

### D. No silent downgrade
When the sanitizer cannot safely prove equivalence, it must not silently pretend equivalence exists.

If evidence is incomplete or parsing is uncertain:
- raise in strict mode, or
- emit explicit uncertainty/warnings in tolerant mode

---

## Include Discovery and Resolution

The sanitizer must recursively resolve `#include` directives using the full context-derived include search path, not only the current file directory and not only forcefield-local trees.

### Rules
1. Build include search paths from the pipeline/context-defined roots.
2. Resolve includes recursively from all active ITP/topology inputs participating in the current system.
3. Preserve deterministic search order.
4. If include shadowing is allowed, resolve it explicitly by policy and record the winner.
5. If include shadowing is not allowed, fail when actual resolution becomes ambiguous.

Do not perform premature global rejection of same-basename files if actual include resolution policy is intended to decide the winner later.

---

## `[ moleculetype ]` Conflict Handling

The sanitizer must identify all `[ moleculetype ]` definitions and determine whether duplicates are:

- exactly equivalent
- safely tolerable under policy
- physically conflicting
- uncertain due to incomplete parsing

### Signature policy
Use deterministic signatures to compare definitions, but do not treat the signature as infallible truth.

If signatures indicate equality, that is evidence of equivalence.

If signatures differ, evaluate whether the difference is:
- physical and real
- namespace-only
- due to uncertain parsing
- due to unsupported/ambiguous formatting

### Uncertain charge extraction
If charge extraction is affected by preprocessing directives or ambiguous formatting, do **not** immediately classify the pair as different.

Instead:
- mark charge comparison as uncertain
- fall back to stronger/secondary checks where available
- only declare hard difference when evidence supports it

This prevents false-positive conflicts from legal topology constructs.

---

## `[ atomtypes ]` and `[ defaults ]`

### `[ defaults ]`
- Must remain globally coherent.
- Incompatible values are a hard error unless an explicit policy says otherwise.

### `[ atomtypes ]`
- Physical LJ conflicts must fail loudly.
- Exact duplicates may be tolerated.
- Tolerant override paths, if allowed, must remain deterministic and fully logged.

### Parsing rule
Prefer explicit tokenization first.

However, parsing must remain robust against nonstandard numeric formatting such as glued-number outputs or spacing anomalies.

Policy:
- use tokenization-first parsing for standard rows
- if numeric layout is nonstandard, use a conservative numeric extraction fallback
- if the result is still ambiguous:
  - raise in strict mode
  - emit explicit warning and uncertainty in tolerant mode

Do **not** silently parse the wrong sigma/epsilon columns.

---

## Managed `system.top` Updates

`topology_sanitizer.py` owns final `system.top` writing.

Responsibilities include:
- validating required includes
- updating sanitizer-managed include blocks
- maintaining deterministic include order
- preserving user content outside the managed block

Preferred behavior:
- if a managed placeholder exists, use it
- otherwise insert the managed block at the deterministic default location defined by the sanitizer policy

The sanitizer must not scatter final topology writes across multiple modules.

---

## Corrected ITP Syncing

When corrected ITP content must replace an original file:

1. preserve the original file as a true backup first
2. then atomically replace the live file
3. keep the operation crash-safe and auditable

Use temp-file + fsync + `os.replace` style semantics where supported.

The versioned backup must always represent the pre-modification original, never an already-patched copy.

---

## Spatial / Geometry Checks Interface

Spatial checks are topology-adjacent but not topology-owned.

The spatial checker may provide:

- dipole check results
- reliability flags
- skip reasons
- atom-count truth from `grompp -pp`
- patch plans / recommended corrections

But final topology artifact writes remain owned by `topology_sanitizer.py`.

---

## Reliability and Uncertainty

The sanitizer must distinguish between:

- **exact**: fully supported and confidently parsed
- **tolerated**: accepted by policy, with explicit logging
- **uncertain**: incomplete or ambiguous evidence
- **unreliable**: known fallback path not suitable for strong assertions

Never collapse `uncertain` or `unreliable` into `exact`.

---

## Strict vs Tolerant Behavior

Strict mode should:
- fail on unresolved ambiguity
- fail on unsupported conflicting definitions
- fail on ambiguous physical parsing where correctness matters

Tolerant mode may:
- proceed with deterministic winner selection
- record warnings and uncertainty
- preserve reproducibility and provenance

Tolerant mode is **not** permission for silent corruption.

---

## Advanced Namespace Prefixing

Namespace prefixing for molecule/atomtype isolation is an advanced strategy, not the default path.

Rules:
- do not enable it casually
- do not let it complicate the normal sanitizer flow
- use only when required by a clearly defined collision-isolation workflow

Default sanitizer behavior should remain as simple and deterministic as possible.

---

## Anti-Patterns

Do not implement or reintroduce the following:

- topology writing from multiple modules without a single writer
- hidden winner selection based on incidental filesystem order
- broad “just regex it” parsing that can silently misread physical parameters
- declaring `different` solely because parsing encountered `#ifdef` / `#include`
- materializing multiple competing winners into live sanitized outputs
- letting compatibility facade code become a new primary logic home

---

## Maintenance Guidance

When modifying sanitizer logic:

1. preserve module boundaries
2. keep deterministic behavior
3. prefer explicit provenance over implicit heuristics
4. make uncertainty visible
5. keep final topology writes centralized in `topology_sanitizer.py`
6. keep `sanitizer_stage.py` thin
7. keep `sanitizer.py` as a compatibility facade only

If a change crosses topology ownership and spatial-check ownership, stop and make the boundary explicit before coding.

---

## Default Expectation for Future Changes

For most future sanitizer work:
- topology-side fixes belong in `topology_sanitizer.py`
- GRO/PBC/dipole/`grompp` consistency fixes belong in `spatial_checker.py`
- orchestration-only changes belong in `sanitizer_stage.py`

Any change that blurs these boundaries must justify why.