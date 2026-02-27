# 05_HTPOLYNET.md — Sandbox, Dynamic Config, Rules, Crosslinking

## A) Sandbox (mandatory)
- Enforce htpolynet workdir:
  OUT_GMX/<RUN_ID>/02_htpolynet/workdir/
- All intermediates must remain inside workdir.

## B) Dynamic system_config.yaml (mandatory)
- Generate OUT_GMX/<RUN_ID>/02_htpolynet/system_config.yaml at runtime using templates.
- Must include N_i, box size L, selected ff/charge, and resource settings if supported.

## C) Rules coupling to forcefield
- Select rules by FF:
  IN/systems/<SYSTEM_ID>/htpolynet/rules/GAFF2/reactions.yaml
  IN/systems/<SYSTEM_ID>/htpolynet/rules/OPLSAA/reactions.yaml
- OPLS naming must be consistent everywhere as OPLSAA in paths/rules bundles.

## D) Chemistry-driven reactive-site contract (no name hardcoding)
- Crosslinker classification is based on reactive_site_count:
  - reactive_site_count >= 2 -> crosslinker
  - reactive_site_count == 1 -> chain-stopper
  - reactive_site_count >= 3 -> multifunctional crosslinker (supported when rules define it)
- The pipeline must not infer crosslinker behavior from molecule name strings.

## E) Reactive metadata source of truth
- reactive_site_count must come from explicit metadata:
  - IN/molecules/<NAME>/meta.yaml, and/or
  - system-level recipe/config (e.g., IN/systems/<SYSTEM_ID>/composition.yaml or IN/systems/<SYSTEM_ID>/htpolynet/reactive_site_counts.yaml)
- For molecules used in HTPolyNet, missing reactive_site_count is a hard error.

## F) Publish + rename stable
- Only publish final outputs into IN/systems/<SYSTEM_ID>/htpolynet/:
  system_<RUN_ID>.gro/.itp/.top and current pointers.
- On publish, normalize names so staging does not need to guess.

## G) Staging (mandatory)
- After publish, stage into IN/systems/<SYSTEM_ID>/gromacs/:
  system.gro, system.itp, system.top (symlink or copy).
