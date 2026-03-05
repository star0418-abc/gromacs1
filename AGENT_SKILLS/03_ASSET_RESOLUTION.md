# 03_ASSET_RESOLUTION.md — Asset Resolution Rules

## A) Charge selection drives BOTH mol2 and charge-consistent topology assets
- If `--charge RESP`:
  - resolve molecular charge-bearing assets from `IN/molecules/<NAME>/mol2/RESP/...`
  - resolve topology assets only from a RESP-compatible topology family for the selected forcefield
- If `--charge CM5`:
  - resolve molecular charge-bearing assets from `IN/molecules/<NAME>/mol2/CM5/...`
  - resolve topology assets only from a CM5-compatible topology family for the selected forcefield

Never silently mix:
- RESP mol2 with CM5-derived topology assets
- CM5 mol2 with RESP-derived topology assets
- one molecule’s charge model with another molecule’s incompatible topology family unless the repo explicitly defines that mixed workflow

## B) Forcefield selection scope
- If `--ff GAFF2`: use only GAFF2-compatible topology assets
- If `--ff OPLS-AA`: use only OPLS-AA-compatible topology assets
- Never silently mix forcefields inside one resolved topology set

Asset existence alone is **not** proof of physical compatibility.

## C) Compatibility policy (forcefield × charge_model)
The pipeline must not hard-code a universal “validated pair set” unless the repository explicitly documents it.

Supported pairs are determined by a **repo-local validated whitelist**.

### Conservative default rule
If no repo-local whitelist is provided, the conservative fallback among the currently supported charge models is:
- `(GAFF2, RESP)`

The following must **not** be assumed validated by default merely because files exist:
- `(GAFF2, CM5)`
- `(OPLS-AA, RESP)`
- `(OPLS-AA, CM5)`

A repository may additionally allow pairs such as `(OPLS-AA, CM5)` **only if** it explicitly marks them as validated or accepted-by-project with documented rationale.

## D) Repo-local compatibility whitelist
The repository should expose one deterministic source of truth for allowed pairs, for example:
- `IN/forcefields/compatibility_matrix.yaml`, or
- an equivalent config artifact owned by the pipeline

Each allowed pair entry must contain at least:
- `ff`
- `charge_model`
- `status` (`validated` / `project-accepted` / `experimental`)
- `evidence_ref`
- `notes`

Resolution rules:
- if pair is `validated` or `project-accepted`: allow
- if pair is `experimental`: require explicit unsafe override
- if pair is absent: hard-fail before stage execution

## E) Unsafe override policy (explicit only)
Unsupported or experimental pairs are allowed only with explicit override:
- `--allow-charge-model-mismatch`
- `--charge-model-mismatch-reason "<why>"`

Override behavior is mandatory:
- loud console warning
- manifest entry with:
  - `unsafe_override=true`
  - `override_type=charge_model_mismatch`
  - `ff`
  - `charge_model`
  - `compatibility_status`
  - `reason`
- explicit statement that the run is outside the default validated parameter regime

No silent downgrade to “best effort”.

## F) Template resolution for PACKMOL
- Prefer `IN/molecules/<NAME>/gro` or `pdb` as PACKMOL template sources
- If both exist, use a deterministic documented preference (for example: `gro > pdb`)
- If neither exists, fail fast
- Do not guess, auto-infer, or fabricate geometry from unrelated assets

## G) Uniqueness guarantee
Any resolution must yield exactly one logical asset path per required role.

Error conditions:
- zero matches → missing asset
- more than one match at the same precedence level → ambiguity, require explicit disambiguation
- incompatible candidates across precedence levels → fail with actionable diagnostics

## H) Manifest audit requirements
Record at minimum:
- `ff`
- `charge_model`
- `compatibility_status`
- `compatibility_evidence_ref`
- resolved mol2 path(s)
- resolved topology asset path(s)
- whether unsafe override was used
- override reason (if any)