# cfd3d — 3D compressible CFD solver

Unstructured finite-volume solver for the compressible Navier-Stokes equations
in modern Fortran (2008+), MPI-parallel.

**Status — Milestone 5.1:** Blast-IC plumbing landed on the M4 stack —
spherical Sedov-Taylor and `blast_bubble` initial conditions, point-gauge
("probe") output to per-probe CSV files, and HLLC positivity guards for
high pressure-ratio Riemann problems. Acceptance: 3-D Sedov-Taylor blast
shock-radius matches the analytical `R(t) = ξ₀(E/ρ_∞)^{1/5} t^{2/5}` to
**5.6 %** at N=40³ (target ≤ 10 %). 17/17 ctest tests pass.

This is the first sub-phase of M5 (TNT blast in/external to a rigid 3-D
structure). M5.2 adds multi-gas / JWL EOS for TNT detonation products,
M5.3 adds external blast on a rigid cube, M5.4 adds internal blast in a
vented room.

**Earlier — Milestone 4:** MPI compressible **Navier-Stokes** with second-order
MUSCL reconstruction and an **algebraic SGS LES** model (Smagorinsky / WALE).
The SGS eddy viscosity drops into `viscous_flux_face` via `μ_eff = μ + ρ ν_t`
and `k_eff = μ cp/Pr + ρ ν_t cp/Pr_t`, with the filter width `Δ = V_cell^(1/3)`
computed per face. WALE self-damps near walls; Smagorinsky is also available
via namelist toggle. The M3 MMS framework was extended to include an
analytical SGS source (4th-order central FD of the analytical SGS stress).

Acceptance status:
- 3D Sod (MUSCL, inviscid): L1_rel(ρ) = **0.6 %**, identical serial vs MPI4.
- M3 MMS (SGS off): observed order **p ≈ 1.88**.
- M4 MMS-LES (WALE on): observed order **p ≈ 1.88** — confirms the SGS
  infrastructure preserves spatial order.

M1/M2/M3 previously landed: M1 serial Euler, M2 METIS partition +
persistent-request halo exchange, M3 MUSCL+Venkat+Sutherland+MMS.

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
