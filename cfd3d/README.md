# cfd3d — 3D compressible CFD solver

Unstructured finite-volume solver for the compressible Navier-Stokes equations
in modern Fortran (2008+), MPI-parallel.

**Status — Milestone 2:** serial + MPI compressible **Euler**. METIS partition
with one-layer ghost cells, persistent-request halo exchange with comm/compute
overlap (pure-interior faces computed while halos are in flight), global CFL
via `MPI_Allreduce(MIN)`, per-rank VTK output. Verified on the 3D Sod shock
tube: L1 error vs the exact Riemann solution is **2.4 %** (rho) / 2.2 % (p) at
N=200, **identical** serial and 4-rank.

M1 — serial Euler — landed previously and remains green.

## Build

From the OpenRadioss top level:

```sh
cmake -B build_cfd -Dbuild=cfd3d -DCMAKE_BUILD_TYPE=Release
cmake --build build_cfd -j
```

Dependencies: gfortran ≥ 8 (or Intel ifort/ifx), CMake ≥ 3.15,
**OpenMPI** (`libopenmpi-dev` + `openmpi-bin`) and **METIS**
(`libmetis-dev` 5.x).

## Test

```sh
ctest --test-dir build_cfd/cfd3d --output-on-failure
```

Tests:
- `test_eos` — primitive↔conservative round-trip, sound-speed sanity
- `test_riemann_hllc` — HLLC fluxes for identical / wall-mirror / rotational
- `test_mesh_topology` — face counts, ownership, volume sum on a 2-hex mesh
- `test_sod_shock_tube` — serial Sod (N=200) ≤ 6 % L1 error
- `test_sod_shock_tube_mpi4` — same Sod on 4 ranks, same tolerance

## Run

```sh
cd build_cfd/cfd3d/tests/sod3d
mpirun -np 4 ../../cfd3d_solver sod3d.nml
```

Each rank writes `sod3d_NNNNNNNN_r<rank>.vtk`; merge externally (e.g. ParaView
`Group Datasets`) for global views.

## M2 architecture notes

- `mpi_runtime` — initialize MPI, holds `t_mpi_ctx` (comm, rank, nproc).
- `metis_binding` — `iso_c_binding` interface to `METIS_PartMeshDual`.
- `partition` — every rank reads the full Gmsh mesh; rank 0 calls METIS;
  `epart` broadcast; each rank extracts local + one-layer ghost slice,
  rebuilds topology/metrics, classifies pure-interior vs partition vs
  boundary faces, builds the halo CSR via `MPI_Alltoall` then `MPI_Alltoallv`.
- `halo_exchange` — persistent `MPI_Send_init` / `MPI_Recv_init` handles;
  `halo_pack_and_start` kicks off Isend/Irecv each RK stage.
- `flux_assembly` — three loops: `residual_pure_interior` (runs while halo
  in flight), `residual_partition` (after `halo_wait`), `residual_boundary`.

## Layout

```
cfd3d/
  src/                 solver source (21 modules)
  tests/               ctest-driven unit + acceptance tests (serial + MPI4)
  cases/sod3d/         Sod shock tube input + standalone Gmsh-msh generator
  doc/
```

## Out of scope for M2 (future)

- True distributed-input ParMETIS (M2.5) — each rank reads its own chunk
  rather than the full mesh. Needed when `nv` exceeds per-rank memory.
- Viscous fluxes + MUSCL reconstruction (M3).
- Turbulence model (M4).
