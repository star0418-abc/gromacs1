# 03_ASSET_RESOLUTION.md — Asset Resolution Rules

## A) Charge selection drives BOTH mol2 and itp
- If --charge RESP:
  - Use IN/molecules/<NAME>/mol2/RESP/...
  - Use an itp consistent with RESP (never mix RESP mol2 with CM5 itp)
- If --charge CM5:
  - Use IN/molecules/<NAME>/mol2/CM5/...
  - Use an itp consistent with CM5

## B) Forcefield selection scope
- If --ff GAFF2: use GAFF2-compatible itp set
- If --ff OPLS-AA: use OPLS-AA-compatible itp set
- Do not silently mix forcefields.

## C) Compatibility matrix (forcefield × charge_model)
- Default-supported (validated) pairs:
  - (GAFF2, CM5)
  - (OPLS-AA, RESP)
- Any other pair is UNSUPPORTED by default and must hard-fail before stage execution.

## D) Unsafe override policy (explicit only)
- Unsupported pairs are allowed only with explicit override:
  - --allow-charge-model-mismatch
  - --charge-model-mismatch-reason "<why>"
- Override behavior is mandatory:
  - loud console warning,
  - manifest entry with: unsafe_override=true, reason, ff, charge_model,
  - explicit statement that the run is no longer a validated parameter set.

## E) Template resolution for PACKMOL
- Prefer IN/molecules/<NAME>/gro or pdb as PACKMOL template.
- If both exist, define a deterministic preference (e.g., gro > pdb) and document it.
- If neither exists, fail fast (do not auto-OCR/guess).

## F) Uniqueness guarantee
- Any resolution must yield exactly one file path; otherwise error:
  - zero matches -> missing asset
  - >1 matches -> ambiguity, require user to disambiguate by naming/version tag
