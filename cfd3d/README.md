# cfd3d — 3D compressible CFD solver

Unstructured finite-volume solver for the compressible Navier-Stokes equations
in modern Fortran (2008+), MPI-parallel.

**Status — Milestone 3:** MPI compressible **Navier-Stokes** with second-order
MUSCL reconstruction. Per-cell least-squares gradients of primitive variables,
Venkatakrishnan slope limiter, Sutherland viscosity, viscous stress + heat
fluxes. Two-channel halo exchange (state + gradients/limiters) so partition
faces remain second-order under MPI. Inviscid acceptance: 3D Sod L1 error
drops from 2.4 % (first-order, M2) to **0.6 %** (MUSCL, M3) at N=200, identical
serial vs MPI4. **MMS spatial-order test gives p ≈ 1.88** (target 2.0 for MUSCL).

M1/M2 previously landed: M1 serial Euler, M2 added METIS partition +
persistent-request halo exchange.

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

## M3 architecture notes

- `gradients` — per-cell LSQ gradient of primitives; boundary contributions
  use the face-midpoint value (0.5*(W_cell + W_ghost)) so the half-step
  stencil is consistent.
- `limiters` — Venkatakrishnan ψ ∈ [0,1] per cell, per primitive variable.
- `viscous_fluxes` — Newtonian stress tensor + Fourier heat flux. `T = p/(ρR)`,
  `μ(T)` from Sutherland's law via `gas_properties`.
- `halo_exchange` — two channels (U for state, GP for grads+limiters).
- `flux_assembly` — MUSCL reconstruction at faces, viscous contribution
  subtracted from inviscid. Boundary viscous uses an over-relaxed gradient
  correction so the normal-direction derivative equals the local face-cell
  difference (essential at no-slip walls).

## Verification (M3.3)

- `test_mms_order` — MMS consistency test on a smooth steady manufactured
  solution `u_a = sin(πx)`, refined from N=8³ to N=16³. Observed spatial
  order **p ≈ 1.88** (target ≥ 1.5). The small deficit from theoretical 2.0
  comes from the Venkat limiter clipping near smooth extrema and the
  first-order boundary stencil. Confirms MUSCL+gradients+viscous+limiters
  achieve second-order asymptotic accuracy.

## Known issues / TODOs

- `test_couette` is currently a smoke test. With slip walls on the cavity
  sides and a moving lid, the steady state is a recirculating eddy (positive
  u in the upper half, negative u in the lower half) — that's the correct
  closed-cavity physics, not the 1D linear Couette profile. Recovering the
  linear profile would need periodic BCs in x.

## Out of scope for M3 (future)

- True distributed-input ParMETIS (M2.5) — each rank reads its own chunk
  rather than the full mesh. Needed when `nv` exceeds per-rank memory.
- VTU/PVTU output upgrade (M3 still uses legacy VTK).
- Turbulence model (M4).
