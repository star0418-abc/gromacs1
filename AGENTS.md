# AGENTS.md

Before making any changes, read **all** files under `AGENT_SKILLS/` in numeric order.

## Bug Knowledge (mandatory for bug fixes)
- Before editing a bug fix, read `docs/BUGS.md`.
- Before editing a bug fix, scan `.ai/bugs.jsonl` for the same file path, function or section, tags, and root-cause pattern.
- After a bug fix is verified, append a new record to `.ai/bugs.jsonl`.
- Update `docs/BUGS.md` when the fix exposes a recurring pattern or reusable prevention rule.
- If this repository is under Git, commit the fix and record the relevant commit hash in the bug record.

## Non-negotiable rules (summary)
1) **Read inputs only from `IN/`.**
2) **Write all GROMACS runtime outputs only to `OUT_GMX/<RUN_ID>/`.**
3) **HTPOLYNET must run in `OUT_GMX/<RUN_ID>/02_htpolynet/workdir/` (sandbox).**
4) **Publish to `IN/` is explicit, atomic, versioned (RUN_ID), and only on success.**
5) **ITP Sanitizer must run before any grompp.**
6) **CLI `T/P` must override mdp defaults (MDP patching).**
7) **Resume/skip must be deterministic with sentinels and checkpoint logic.**
8) **On failure: keep full context (logs, workdir), print last 50 lines, write crash report.**

## Definition of Done
A change is acceptable only if all acceptance checks in `AGENT_SKILLS/00_CONTRACT.md` pass.

## Housekeeping
- Keep behavior reproducible; update README when behavior changes.
- Do not leave temporary/test artifacts in the repo after finishing changes.
