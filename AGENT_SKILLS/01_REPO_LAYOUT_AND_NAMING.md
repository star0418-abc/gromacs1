# 01_REPO_LAYOUT_AND_NAMING.md — Authoritative Layout & Naming

## A) Layout (authoritative)
project/
  IN/
    molecules/<NAME>/{mol,mol2/{CM5,RESP},pdb,gro,itp,top}/...
    systems/<SYSTEM_ID>/
      packmol/{pdb,gro}/
      htpolynet/{gro,itp,top,rules/{GAFF2,OPLSAA}/}
      gromacs/{gro,itp,top,mdp,ndx}/
      gromacs/analysis_pairs.yaml
    forcefield/
      oplsaa/gromacs/{atomtypes.itp,nonbonded.itp}
      gaff2/{amber/,gromacs/{atomtypes.itp,nonbonded.itp}}
  OUT_GMX/<RUN_ID>/
    01_packmol/
    02_htpolynet/{system_config.yaml,workdir/}
    03_gromacs/{em,nvt,npt,md}/
    analysis/{rdf,cn}/
    logs/
    manifest.json

## B) Versioning rules
- Any published asset must be versioned by RUN_ID:
  e.g., `system_<RUN_ID>.gro` and `system_current.gro` pointer.
- Never overwrite the molecule library in-place.
- Prefer “current pointer” as symlink or CURRENT file.

## C) Naming guidance (minimal)
- Molecule folder name: stable unique `<NAME>` (e.g., PC, LiTFSI, tetraglyme, OEGMA, UEMA, MPB, PEGDMA, TPO-L).
- System folder: `<SYSTEM_ID>` (e.g., SYS_001).
- Run id: `<RUN_ID>` includes timestamp + system id recommended.

## D) Uniqueness requirement
- Asset resolution must yield exactly one file. If multiple matches exist, stop with error.
