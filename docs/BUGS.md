# Bug Knowledge Base

`docs/BUGS.md` is the human-readable summary of recurring bug knowledge for this repository. `.ai/bugs.jsonl` is the canonical append-only ledger, and `.ai/bug_schema.md` defines the JSONL record contract.

## Purpose

Use this document to capture bug knowledge that should change how future fixes are approached, not just what happened once. Organize entries by failure pattern and subsystem so the next fix starts from prior root causes, guardrails, and regression tests instead of rediscovering them.

## How to use before fixing a bug

1. Read the relevant sections of this file before editing code or workflow docs.
2. Scan `.ai/bugs.jsonl` for the same file path.
3. Scan `.ai/bugs.jsonl` for the same function, class, CLI surface, or documentation section.
4. Scan `.ai/bugs.jsonl` for matching tags.
5. Scan `.ai/bugs.jsonl` for a similar root-cause or trigger pattern even if the file differs.
6. Reuse the prior prevention rules and regression tests unless the new case proves they were insufficient.

Example searches:

```bash
grep -n '"file":"pipeline/stages/topology_sanitizer.py"' .ai/bugs.jsonl
grep -n '"function":"patch_mdp"' .ai/bugs.jsonl
grep -n '"tags":.*"atomic-write"' .ai/bugs.jsonl
```

```powershell
Select-String -Path .ai/bugs.jsonl -Pattern '"file":"pipeline/stages/topology_sanitizer.py"'
Select-String -Path .ai/bugs.jsonl -Pattern '"tags":.*"parser"'
```

## When to add a JSONL record

- Add one new record after every confirmed bug fix, duplicate classification, or `wontfix` decision that changes code, configuration, or operator-facing behavior.
- Add one new record after regression hardening that closes a real failure mode discovered by tests or review.
- In a Git repository, create the fix commit first, then append the JSONL record with the real fix commit hash, then commit the ledger and knowledge-base update.
- Do not invent historical entries.
- Do not log speculative bugs or TODOs that were not actually investigated.

## When to also update BUGS.md

- Update this file when a fix establishes a reusable rule, invariant, or search keyword.
- Update this file when the same failure pattern can recur in more than one file or stage.
- Update this file when the prevention guidance is broader than a single code location.
- Update the relevant README section at the same time when user-facing behavior or operator workflow changed.

## Recurring bug patterns

### Filesystem and Atomic Publish Rules

- Keep runtime outputs under `OUT_GMX/<RUN_ID>/`; only publish reusable assets into `IN/` after success.
- For any replacement or publish, write a temporary sibling file, flush and `fsync` it, then atomically rename it into place.
- Create backups before replacing a live topology or pointer file.
- Resume and publish logic must never leave a partially written live file behind after interruption.
- Search tags: `filesystem`, `atomic-write`, `publish`, `backup`, `current-pointer`.
- Promote a fix into this section when the bug involves truncation, partial writes, misplaced outputs, or non-atomic publish behavior.

### Parser Strictness and Semantic Normalization

- Normalize aliases, case, and spelling before semantic validation.
- Strict and non-strict modes must operate on the same canonical representation; tolerant mode may relax the outcome, not reinterpret the input.
- When tolerant mode continues, warnings and manifests must preserve the original token and the canonical key that was used.
- Detect duplicate exact keys and semantic alias collisions during the initial parse, then apply the same deterministic winner on read, patch, and write paths.
- If strict mode would reject malformed or ambiguous parser input, fail before patching rather than preserving a silently reinterpreted template.
- Treat malformed `.ndx` group-body tokens as unsafe input when thermostat auto decisions depend on group sizes; never silently undercount atoms and continue.
- Sanitize non-finite or boundary safety thresholds in one helper before auto/fallback decisions so invalid config cannot create mode-dependent behavior.
- When validation helpers stage diagnostics in a local sink, flush that sink before every early return triggered by missing dependent values; incomplete input must not suppress earlier warnings.
- Add paired tests for strict and tolerant paths using the same raw input.
- Search tags: `parser`, `strictness`, `normalization`, `semantics`, `mdp`, `itp`.
- Promote a fix into this section when one input can be interpreted differently across parsers or strictness modes.

### Stage Boundary and Resume Provenance

- Read inputs only from `IN/` and keep runtime outputs under `OUT_GMX/<RUN_ID>/`.
- Preserve deterministic sentinels, checkpoint selection, and crash context when fixing resume or failure paths.
- If a stage runner can be invoked directly or owns substage skip logic, it must write its own final `done.ok` after `stage_state.json`/`repro.json` are durable; do not rely only on an outer dispatcher wrapper to create the success sentinel.
- Same-stage resume compatibility checks must preserve recorded semantic fields such as `velocity_mode` when rebuilding `stage_state.json`, or the first normal resume can fail on self-inflicted fingerprint drift.
- Resume compatibility gates must fail closed on non-runtime schema drift; never warn and compare only shared keys unless an explicit, narrowly scoped compatibility rule reconstructs the missing field first.
- GRO velocity restart validation must require parseable and materially non-zero velocities; syntactically present all-zero velocity columns count as missing restart state.
- Under wrapper launchers such as `mpirun ... gmx_mpi`, record launcher provenance separately from the actual GROMACS binary provenance so audit output does not attribute wrapper paths to the binary.
- Do not weaken the requirement to keep the full failing workdir and the last log context on error.
- Search tags: `resume`, `checkpoint`, `sentinel`, `provenance`, `crash-context`.
- Promote a fix into this section when the bug changes skip logic, checkpoint choice, or failure preservation rules.

## Regression checklist

- Read `docs/BUGS.md` and scan `.ai/bugs.jsonl` by file, function or section, tags, and root cause before editing.
- Capture the trigger condition precisely enough to reproduce or recognize the failure.
- Record the exact primary file path and function or section in the JSONL entry.
- Fix the underlying invariant, not only the visible symptom.
- Add regression tests, or explain in the record why automated coverage was not possible.
- If the repository is under Git, commit the fix first, then append the JSONL record with the fix commit hash, then commit the ledger and documentation update.
- Update this file when the fix becomes a reusable rule.
- Update the relevant README section when behavior, workflow, or operator expectations changed.

## Example entry mapping to a JSONL record ID

The IDs below are examples only. They are templates that demonstrate how a ledger record maps into this knowledge base; they are not historical repository bugs.

| Example JSONL record ID | BUGS.md section | Why it belongs here |
| --- | --- | --- |
| `BUG-2026-03-07-901` | `Filesystem and Atomic Publish Rules` | The root cause is a non-atomic live-file replacement pattern that can recur anywhere the pipeline publishes or syncs files. |
| `BUG-2026-03-07-902` | `Parser Strictness and Semantic Normalization` | The root cause is canonicalization drift between strict and tolerant parsing, which is a reusable parser rule rather than a one-file anecdote. |
