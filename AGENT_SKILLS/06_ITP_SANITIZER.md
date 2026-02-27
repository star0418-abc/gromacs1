# 06_ITP_SANITIZER.md — Atomtype Extraction, Sanitization, Namespace Defense

## A) When to run
- Run after staging and before any grompp.

## B) Inputs
- All .itp files that will be included by the system:
  - forcefield library ITPs (from IN/forcefield/<FF>/gromacs/)
  - staged system.itp (from htpolynet)
  - molecule itps referenced by system.top
  - any additional included itps (recursive #include resolution)

## C) Outputs (versioned + current pointers)
1) Combined atomtypes:
- IN/systems/<SYSTEM_ID>/gromacs/combined_atomtypes_<RUN_ID>.itp
- IN/systems/<SYSTEM_ID>/gromacs/combined_atomtypes_current.itp (pointer)

2) Sanitized itps:
- IN/systems/<SYSTEM_ID>/gromacs/itp_sanitized_<RUN_ID>/
- IN/systems/<SYSTEM_ID>/gromacs/itp_sanitized_current/ (pointer)

3) system.top include order must be:
- [ defaults ] section (nbfunc, comb-rule, gen-pairs, fudgeLJ, fudgeQQ)
- #include "combined_atomtypes_current.itp"
- #include sanitized molecule itps (files with [ moleculetype ] only)
- [ system ] / [ molecules ]

## D) Conflict policy (fail-fast default)
- Conflict signature: name + ptype + mass + sigma + epsilon (EXCLUDES charge)
- Same atomtype name with different physical parameters -> ERROR.
- Charge mismatches -> WARNING (unless --strict-charge enabled, then ERROR).
- Charge output policy for mismatches:
  - Do NOT zero charges to 0.0.
  - Keep deterministic source values (or fail in strict mode), and report loudly in logs/manifest.

## E) [ defaults ] / comb-rule validation
- Parse [ defaults ] from all scanned files.
- If multiple defaults exist and differ in nbfunc/comb-rule/fudgeLJ/fudgeQQ -> ERROR.
- system.top MUST include [ defaults ] BEFORE atomtypes to ensure correct LJ interpretation.

## F) Namespace collision defense (cross-forcefield contamination)
- Track source_ff for each atomtype based on file origin:
  - Files from forcefield_dir/* -> ctx.ff
  - Other files -> "UNKNOWN"
- If same atomtype name appears from different source_ff with different params -> ERROR.
- Allow override only if user explicitly enables --allow-override.
- Override selection is deterministic by priority: forcefield (100) > staged (50) > htpolynet (40) > molecule (30).
- Any override must be written to manifest with a loud warning and the winner recorded.

## G) Reproducibility
- Files are processed in deterministic sorted order (Path.as_posix()).
- Override winners are selected deterministically.
- Combined atomtypes output is sorted by atomtype name.
- Logs record which file "won" for each overridden atomtype.

## H) Optional Namespace Prefixing (OFF by default)
- If enabled via --atomtype-prefix, sanitizer rewrites atomtype names with a prefix.
- SAFETY CHECK: Before prefixing, detect implicit/library topology usage:
  - Scan [ bonds ], [ angles ], [ dihedrals ] for lines lacking explicit parameters.
  - If detected -> ERROR with actionable message (prefixing not supported).
- If safe, update:
  - [ atomtypes ] definitions
  - [ atoms ] type field
- Combined atomtypes file uses prefixed names.
- Warn: htpolynet rules or any rule relying on atomtype naming may also need updates.

## I) Token-based atomtypes parser
- Uses explicit tokenization by whitespace, not fragile regex.
- Strips inline comments (;) before parsing.
- Supports formats:
  - name at_num mass charge ptype sigma epsilon (7 tokens)
  - name bond_type at_num mass charge ptype sigma epsilon (8 tokens)
  - name mass charge ptype sigma epsilon (6 tokens, simplified)
- Distinguishes format by checking if token[1] is numeric (atomic number) vs string (bond_type).
- Errors include file path and original raw line for debugging.

## J) Include dependency discovery
- Recursively resolves #include directives from forcefield files.
- Prevents infinite loops via visited set.
- Warns if include target not found.

## K) Molecule ITP classification
- Files with [ moleculetype ] section are classified as molecule ITPs.
- system.top only includes molecule ITPs, not forcefield parameter files.
