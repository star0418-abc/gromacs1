# Bug JSONL Schema

`.ai/bugs.jsonl` is the canonical structured bug ledger for this repository. It must remain valid JSON Lines at all times: exactly one JSON object per line, no comments, no custom wrappers, and no in-place edits to older records.

## Required fields

| Field | Type | Rules |
| --- | --- | --- |
| `id` | string | Unique entry ID in the form `BUG-YYYY-MM-DD-NNN`. Never reuse an old ID. |
| `created_at` | string | RFC 3339 / ISO 8601 timestamp with timezone offset, for example `2026-03-07T14:23:11+08:00`. |
| `status` | string | One of `open`, `fixed`, `wontfix`, `duplicate`. |
| `severity` | string | One of `low`, `medium`, `high`, `critical`. |
| `project` | string | Repo-local subsystem or surface, for example `pipeline.sanitizer`, `pipeline.gromacs`, or `docs.workflow`. |
| `file` | string | Primary repo-relative path using forward slashes. Choose the file that best represents the root cause. |
| `function` | string | Function, method, class, CLI surface, or documentation section that future searches should use. |
| `lines` | string | One-based line number or inclusive line range from the fixed revision, for example `314` or `742-811`. Use `n/a` only when no stable line mapping exists. |
| `symptom` | string | Observable failure, regression, or incorrect behavior. |
| `root_cause` | string | Technical cause that made the symptom possible. |
| `trigger` | string | Concrete condition that exposes the bug. |
| `fix` | string | Concise description of the implemented fix or disposition. |
| `prevention` | array of strings | Concrete guardrails, invariants, or review checks that should prevent recurrence. |
| `tests_added` | array of strings | Added or updated regression tests. Use `[]` only when automated coverage was not possible. |
| `related_commits` | array of strings | Git commit hashes related to the fix. Leave `[]` only when the work is outside Git or no commit exists yet. |
| `related_issues` | array of strings | Issue IDs, PR IDs, URLs, or prior bug-record IDs that help cross-reference the record. |
| `tags` | array of strings | Lower-kebab-case search tags such as `atomic-write`, `parser`, or `checkpoint`. |

## Optional fields

| Field | Type | Purpose |
| --- | --- | --- |
| `run_id` | string | Link the bug to a specific `OUT_GMX/<RUN_ID>/` execution when relevant. |
| `duplicate_of` | string | Record the canonical bug or issue ID when `status` is `duplicate`. |
| `supersedes` | array of strings | Prior bug record IDs that this entry corrects or replaces in an append-only way. |
| `notes` | string | Small extra context that does not fit cleanly into the required fields. |

## Allowed enums

### `status`

- `open`: known bug recorded before a fix is available
- `fixed`: bug resolved and ready for future reuse
- `wontfix`: intentionally not fixed, with rationale preserved
- `duplicate`: no new root cause; refer to the canonical prior record

### `severity`

- `low`: localized annoyance or low-risk correctness issue
- `medium`: meaningful behavior bug with contained impact
- `high`: major correctness, data-integrity, or workflow breakage
- `critical`: severe corruption, unsafe output, or repository-wide blocker

## Writing rules

1. Write exactly one JSON object per physical line in `.ai/bugs.jsonl`.
2. Keep the ledger append-only. Never rewrite or delete older lines.
3. Do not put example or template entries into `.ai/bugs.jsonl`; examples belong only in this schema file and `docs/BUGS.md`.
4. Use concrete repo-relative paths, function or section names, and search-friendly tags.
5. Use empty arrays instead of placeholder strings like `"none"` or `"n/a"` for `tests_added`, `related_commits`, `related_issues`, and `tags`.
6. In Git repositories, preferred order is: commit the actual fix first, append the ledger record with the real fix commit hash, then commit the ledger and documentation update.
7. If a later correction is needed, append a new record with a new `id` and cross-reference the prior record through `related_issues`, `duplicate_of`, or `supersedes`.

## Example JSONL records

The following records are examples only. They are templates for future entries and do not describe actual historical repository bugs.

### Example: filesystem / atomic write bug

```json
{"id":"BUG-2026-03-07-901","created_at":"2026-03-07T15:10:00+08:00","status":"fixed","severity":"high","project":"pipeline.sanitizer","file":"pipeline/stages/topology_sanitizer.py","function":"sync_corrected_itp","lines":"1480-1544","symptom":"A corrected topology file could become truncated after an interrupted replacement step.","root_cause":"The live file was overwritten directly instead of writing a sibling temporary file and using an atomic replace after persistence barriers.","trigger":"The process exits or is interrupted after truncation begins but before the replacement write is fully durable.","fix":"Write corrected content to a temporary sibling path, fsync the file, create the backup first, replace the live file atomically, and fsync the parent directory.","prevention":["Use temp-file plus fsync plus atomic replace for every publish or in-place topology update.","Create backups before replacing a live file.","Test interruption-safe replacement paths whenever publish logic changes."],"tests_added":["tests/test_atomic_publish.py::test_atomic_replace_preserves_old_file_on_failure","tests/test_atomic_publish.py::test_backup_created_before_live_replace"],"related_commits":["example-fix-commit-1234abcd"],"related_issues":["docs/BUGS.md#filesystem-and-atomic-publish-rules"],"tags":["filesystem","atomic-write","publish","backup","sanitizer"]}
```

### Example: parser / strict-vs-nonstrict semantic bug

```json
{"id":"BUG-2026-03-07-902","created_at":"2026-03-07T15:25:00+08:00","status":"fixed","severity":"medium","project":"pipeline.mdp","file":"pipeline/mdp_patcher.py","function":"normalize_mdp_keys","lines":"210-312","symptom":"Strict and non-strict parsing produced different effective meanings for the same MDP input.","root_cause":"Key aliases and spelling variants were normalized after mode-specific validation, so each mode inspected a different logical key set.","trigger":"An input file contains alias spellings such as ref-t, ref_t, or mixed-case variants that should map to the same canonical key.","fix":"Normalize aliases before validation, preserve source-token provenance for warnings, and run strict and tolerant checks against the same canonical representation.","prevention":["Normalize parser tokens before semantic validation.","Make tolerant mode relax outcomes, not input meaning.","Add paired strict and tolerant regression tests for the same raw input."],"tests_added":["tests/test_mdp_patcher.py::test_aliases_normalized_before_strict_validation","tests/test_mdp_patcher.py::test_tolerant_mode_keeps_same_canonical_keys"],"related_commits":["example-fix-commit-5678efgh"],"related_issues":["docs/BUGS.md#parser-strictness-and-semantic-normalization"],"tags":["parser","strictness","normalization","semantics","mdp"]}
```
