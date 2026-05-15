#
# cfd3d — 3D compressible CFD solver (unstructured FV, Fortran, MPI).
# Milestone 2: serial compressible Euler + MPI domain decomposition
#              (METIS partition, ghost-layer halos, persistent comm/compute overlap).
#
cmake_minimum_required(VERSION 3.15)

enable_language(Fortran)

set(CFD3D_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_Fortran_STANDARD 2008)

# --- External dependencies (M2) --------------------------------------------
find_package(MPI REQUIRED COMPONENTS Fortran)

# METIS (from libmetis-dev). No standard FindMETIS in CMake; locate manually.
find_path(METIS_INCLUDE_DIR metis.h)
find_library(METIS_LIBRARY metis)
if(NOT METIS_INCLUDE_DIR OR NOT METIS_LIBRARY)
   message(FATAL_ERROR "METIS not found (install libmetis-dev). "
                       "Include: ${METIS_INCLUDE_DIR}  Lib: ${METIS_LIBRARY}")
endif()
message(STATUS "METIS include: ${METIS_INCLUDE_DIR}")
message(STATUS "METIS lib:     ${METIS_LIBRARY}")

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
   set(CFD3D_FFLAGS_COMMON
       -ffree-form -ffree-line-length-none
       -fimplicit-none
       -Wall -Wno-unused-dummy-argument)
   set(CFD3D_FFLAGS_RELEASE -O2 -g)
   set(CFD3D_FFLAGS_DEBUG   -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
   set(CFD3D_FFLAGS_COMMON -free -implicitnone -warn all)
   set(CFD3D_FFLAGS_RELEASE -O2 -g)
   set(CFD3D_FFLAGS_DEBUG   -O0 -g -check all -traceback -fpe0)
else()
   set(CFD3D_FFLAGS_COMMON "")
   set(CFD3D_FFLAGS_RELEASE -O2)
   set(CFD3D_FFLAGS_DEBUG   -O0 -g)
endif()

set(CFD3D_CORE_SOURCES
   src/kinds.f90
   src/constants.f90
   src/mpi_runtime.f90
   src/metis_binding.f90
   src/mesh_types.f90
   src/mesh_io_gmsh.f90
   src/mesh_topology.f90
   src/mesh_metrics.f90
   src/mesh_module.f90
   src/gas_properties.f90
   src/eos_ideal_gas.f90
   src/fields.f90
   src/bc_types.f90
   src/bc_apply.f90
   src/gradients.f90
   src/limiters.f90
   src/riemann_hllc.f90
   src/sgs_model.f90
   src/viscous_fluxes.f90
   src/mms.f90
   src/probes.f90
   src/eos_jwl.f90
   src/partition.f90
   src/halo_exchange.f90
   src/flux_assembly.f90
   src/time_integration.f90
   src/io_vtk_legacy.f90
   src/solver_control.f90
   src/solver_driver.f90)

add_library(cfd3d_core OBJECT ${CFD3D_CORE_SOURCES})
target_compile_options(cfd3d_core PRIVATE
   ${CFD3D_FFLAGS_COMMON}
   $<$<CONFIG:Release>:${CFD3D_FFLAGS_RELEASE}>
   $<$<CONFIG:Debug>:${CFD3D_FFLAGS_DEBUG}>)
target_link_libraries(cfd3d_core PRIVATE MPI::MPI_Fortran)
set_target_properties(cfd3d_core PROPERTIES
   Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_executable(cfd3d_solver src/main.f90)
target_link_libraries(cfd3d_solver PRIVATE cfd3d_core MPI::MPI_Fortran ${METIS_LIBRARY})
target_compile_options(cfd3d_solver PRIVATE
   ${CFD3D_FFLAGS_COMMON}
   $<$<CONFIG:Release>:${CFD3D_FFLAGS_RELEASE}>
   $<$<CONFIG:Debug>:${CFD3D_FFLAGS_DEBUG}>)
set_target_properties(cfd3d_solver PROPERTIES
   Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")

# Mesh generator for the Sod 3D test case (standalone, no deps)
add_executable(gen_sod_mesh cases/sod3d/gen_sod_mesh.f90)
target_compile_options(gen_sod_mesh PRIVATE
   ${CFD3D_FFLAGS_COMMON}
   $<$<CONFIG:Release>:${CFD3D_FFLAGS_RELEASE}>
   $<$<CONFIG:Debug>:${CFD3D_FFLAGS_DEBUG}>)
set_target_properties(gen_sod_mesh PROPERTIES
   Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_gen")

# Mesh generator for the external-blast-on-cube test case (standalone, no deps)
add_executable(gen_box_with_cube cases/blast_cube/gen_box_with_cube.f90)
target_compile_options(gen_box_with_cube PRIVATE
   ${CFD3D_FFLAGS_COMMON}
   $<$<CONFIG:Release>:${CFD3D_FFLAGS_RELEASE}>
   $<$<CONFIG:Debug>:${CFD3D_FFLAGS_DEBUG}>)
set_target_properties(gen_box_with_cube PROPERTIES
   Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_gen_cube")

enable_testing()
add_subdirectory(tests)
