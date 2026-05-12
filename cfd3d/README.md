# cfd3d — 3D compressible CFD solver

Unstructured finite-volume solver for the compressible Navier-Stokes equations
in modern Fortran (2008+), MPI-parallel by design.

**Status — Milestone 1 (this branch):** serial compressible **Euler**.
Mesh I/O (Gmsh msh2 ASCII), HLLC Riemann flux, SSP-RK3, VTK legacy output,
verified against the exact 1D Riemann solution on a 3D Sod shock tube.

Future milestones add MPI domain decomposition (M2), viscous fluxes + MUSCL
reconstruction (M3), and a turbulence model (M4). The M1 data layout
(`nc_internal/nc_total` cell split, interior-then-boundary face ordering,
ghost-state-based BCs) is structured so that M2/M3 drop in without rework.

## Build

From the OpenRadioss top level:

```sh
cmake -B build_cfd -Dbuild=cfd3d -DCMAKE_BUILD_TYPE=Release
cmake --build build_cfd -j
```

Requires gfortran ≥ 8 (or Intel ifort/ifx) and CMake ≥ 3.15. No external
libraries in M1.

## Test

```sh
ctest --test-dir build_cfd/cfd3d --output-on-failure
```

Runs unit tests (EOS, HLLC, mesh topology) and the Sod shock tube acceptance
case (M1 ⇒ L1 error < 6 % vs. exact Riemann at N=200).

## Run the Sod case

```sh
cd build_cfd/cfd3d/tests/sod3d
./../../cfd3d_solver sod3d.nml
```

This writes `sod3d_NNNNNNNN.vtk` snapshots viewable in ParaView.

## Layout

```
cfd3d/
  src/                 solver source (15 modules)
  tests/               ctest-driven unit + acceptance tests
  cases/sod3d/         Sod shock tube input + standalone Gmsh-msh generator
  doc/
```

Source dependency DAG:

```
kinds → constants → {mesh_types, eos_ideal_gas, bc_types}
mesh_types → mesh_topology → mesh_metrics → mesh_module
mesh_io_gmsh → mesh_module
fields → {bc_apply, flux_assembly, time_integration}
eos_ideal_gas → {bc_apply, riemann_hllc, flux_assembly}
riemann_hllc → flux_assembly → time_integration → solver_driver
io_vtk_legacy, solver_control → solver_driver → main
```
