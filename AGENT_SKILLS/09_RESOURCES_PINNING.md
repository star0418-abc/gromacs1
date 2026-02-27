# 09_RESOURCES_PINNING.md — CPU/GPU Pinning for GROMACS and HTPOLYNET

## A) Supported flags
- --gmx-nt <int>
- --gmx-gpu-id <string>
- --omp-num-threads <int> (sets OMP_NUM_THREADS)

## B) Must apply to BOTH:
1) Pipeline main gromacs calls:
- pass -nt / -gpu_id to mdrun
- set OMP_NUM_THREADS

2) HTPOLYNET internal gromacs relax calls:
- If htpolynet supports passing mdrun args in config: inject into system_config.yaml
- Else: enforce environment variables/wrapper so it cannot occupy all cores/GPUs

## C) Recordkeeping
- Record effective CPU/GPU settings in manifest.
