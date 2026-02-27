# 11_CRASH_AND_CLEANUP.md — Crash Context Preservation & Cleanup Policy

## A) On failure (exit code != 0), MUST:
- Print last 50 lines of relevant log(s)
- Print absolute workdir path
- Print exact executed command line
- Write crash_report.txt into the failing stage directory with the above
- Preserve the entire OUT_GMX/<RUN_ID>/ subtree (never cleanup on failure)

## B) Cleanup policy
- Default: no cleanup.
- Optional flag: --cleanup-success
  - Only after successful completion may remove large files (e.g., .trr).
  - Never delete .tpr/.mdp/.log needed for provenance.

## C) Debugging pointers
- Ensure the terminal output includes where to look (workdir path).
