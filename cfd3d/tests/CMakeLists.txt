#
# cfd3d tests
#

set(TEST_FFLAGS_COMMON ${CFD3D_FFLAGS_COMMON})
set(TEST_FFLAGS_REL    ${CFD3D_FFLAGS_RELEASE})
set(TEST_FFLAGS_DBG    ${CFD3D_FFLAGS_DEBUG})

function(add_cfd3d_test name)
   add_executable(${name} ${name}.f90)
   target_link_libraries(${name} PRIVATE cfd3d_core MPI::MPI_Fortran ${METIS_LIBRARY})
   target_compile_options(${name} PRIVATE
      ${TEST_FFLAGS_COMMON}
      $<$<CONFIG:Release>:${TEST_FFLAGS_REL}>
      $<$<CONFIG:Debug>:${TEST_FFLAGS_DBG}>)
   set_target_properties(${name} PROPERTIES
      Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_${name}")
   target_include_directories(${name} PRIVATE
      "${CMAKE_BINARY_DIR}/cfd3d/mod")
endfunction()

add_cfd3d_test(test_eos)
add_cfd3d_test(test_riemann_hllc)
add_cfd3d_test(test_mesh_topology)
add_cfd3d_test(test_sod_shock_tube)
add_cfd3d_test(test_couette)
add_cfd3d_test(test_mms_order)
add_cfd3d_test(test_sedov)
add_cfd3d_test(test_eos_jwl)
add_cfd3d_test(test_cut_cell_geom)
add_test(NAME test_cut_cell_geom
         COMMAND $<TARGET_FILE:test_cut_cell_geom>)
add_cfd3d_test(test_cut_cell_build)
add_test(NAME test_cut_cell_build
         COMMAND $<TARGET_FILE:test_cut_cell_build>)
add_cfd3d_test(test_cut_cell_solver)
add_cfd3d_test(test_blast_overpressure)
add_cfd3d_test(test_blast_cube)
add_cfd3d_test(test_blast_room)
add_cfd3d_test(test_blast_l_room)

# Working directory for the integration tests — meshes + namelists are generated here.
set(CASE_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/sod3d")
file(MAKE_DIRECTORY "${CASE_WORKDIR}")

# Tiny mesh for topology test: 2 hex cells along x.
add_test(NAME generate_topo_mesh
         COMMAND $<TARGET_FILE:gen_sod_mesh> 2 1 1 2.0 1.0 1.0 "${CASE_WORKDIR}/topo2.msh"
         WORKING_DIRECTORY "${CASE_WORKDIR}")

# Sod acceptance mesh: Nx=200 along x.
add_test(NAME generate_sod_mesh
         COMMAND $<TARGET_FILE:gen_sod_mesh> 200 1 1 1.0 0.05 0.05 "${CASE_WORKDIR}/sod3d.msh"
         WORKING_DIRECTORY "${CASE_WORKDIR}")

# Copy the namelist into the case workdir so relative paths work.
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/sod3d/sod3d.nml"
   "${CASE_WORKDIR}/sod3d.nml" COPYONLY)

add_test(NAME test_eos COMMAND $<TARGET_FILE:test_eos>)
add_test(NAME test_riemann_hllc COMMAND $<TARGET_FILE:test_riemann_hllc>)

add_test(NAME test_mesh_topology
         COMMAND $<TARGET_FILE:test_mesh_topology> "${CASE_WORKDIR}/topo2.msh"
         WORKING_DIRECTORY "${CASE_WORKDIR}")
set_tests_properties(test_mesh_topology PROPERTIES DEPENDS generate_topo_mesh)

add_test(NAME test_sod_shock_tube
         COMMAND $<TARGET_FILE:test_sod_shock_tube> "${CASE_WORKDIR}/sod3d.nml"
         WORKING_DIRECTORY "${CASE_WORKDIR}")
set_tests_properties(test_sod_shock_tube PROPERTIES DEPENDS generate_sod_mesh)

# Couette flow (viscous) acceptance test.
set(COUETTE_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/couette")
file(MAKE_DIRECTORY "${COUETTE_WORKDIR}")
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/couette/couette.nml"
   "${COUETTE_WORKDIR}/couette.nml" COPYONLY)
add_test(NAME generate_couette_mesh
         COMMAND $<TARGET_FILE:gen_sod_mesh> 4 32 4 1.0 1.0 1.0 "${COUETTE_WORKDIR}/couette.msh"
         WORKING_DIRECTORY "${COUETTE_WORKDIR}")
add_test(NAME test_couette
         COMMAND $<TARGET_FILE:test_couette> "${COUETTE_WORKDIR}/couette.nml"
         WORKING_DIRECTORY "${COUETTE_WORKDIR}")
set_tests_properties(test_couette PROPERTIES DEPENDS generate_couette_mesh)

# MMS spatial-order test: refine box [0,1]^3 from N=8 to N=16; check order≈2.
set(MMS_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/mms")
file(MAKE_DIRECTORY "${MMS_WORKDIR}")
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/mms/mms.nml.coarse"
   "${MMS_WORKDIR}/mms.nml.coarse" COPYONLY)
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/mms/mms.nml.fine"
   "${MMS_WORKDIR}/mms.nml.fine" COPYONLY)
add_test(NAME generate_mms_coarse
         COMMAND $<TARGET_FILE:gen_sod_mesh> 8 8 8 1.0 1.0 1.0 "${MMS_WORKDIR}/mms_coarse.msh"
         WORKING_DIRECTORY "${MMS_WORKDIR}")
add_test(NAME generate_mms_fine
         COMMAND $<TARGET_FILE:gen_sod_mesh> 16 16 16 1.0 1.0 1.0 "${MMS_WORKDIR}/mms_fine.msh"
         WORKING_DIRECTORY "${MMS_WORKDIR}")
add_test(NAME test_mms_order
         COMMAND $<TARGET_FILE:test_mms_order> "${MMS_WORKDIR}/mms.nml"
         WORKING_DIRECTORY "${MMS_WORKDIR}")
set_tests_properties(test_mms_order PROPERTIES
   DEPENDS "generate_mms_coarse;generate_mms_fine")

# M4 MMS-LES — same MMS framework with WALE SGS turned on.
set(MMS_LES_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/mms_les")
file(MAKE_DIRECTORY "${MMS_LES_WORKDIR}")
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/mms_les/mms_les.nml.coarse"
   "${MMS_LES_WORKDIR}/mms_les.nml.coarse" COPYONLY)
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/mms_les/mms_les.nml.fine"
   "${MMS_LES_WORKDIR}/mms_les.nml.fine" COPYONLY)
add_test(NAME generate_mms_les_coarse
         COMMAND $<TARGET_FILE:gen_sod_mesh> 8 8 8 1.0 1.0 1.0 "${MMS_LES_WORKDIR}/mms_les_coarse.msh"
         WORKING_DIRECTORY "${MMS_LES_WORKDIR}")
add_test(NAME generate_mms_les_fine
         COMMAND $<TARGET_FILE:gen_sod_mesh> 16 16 16 1.0 1.0 1.0 "${MMS_LES_WORKDIR}/mms_les_fine.msh"
         WORKING_DIRECTORY "${MMS_LES_WORKDIR}")
add_test(NAME test_mms_order_les
         COMMAND $<TARGET_FILE:test_mms_order> "${MMS_LES_WORKDIR}/mms_les.nml"
         WORKING_DIRECTORY "${MMS_LES_WORKDIR}")
set_tests_properties(test_mms_order_les PROPERTIES
   DEPENDS "generate_mms_les_coarse;generate_mms_les_fine")

# M5.1 Sedov-Taylor blast (point-energy IC, strong-shock capture).
set(SEDOV_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/sedov")
file(MAKE_DIRECTORY "${SEDOV_WORKDIR}")
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/sedov/sedov.nml"
   "${SEDOV_WORKDIR}/sedov.nml" COPYONLY)
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/sedov/probes.csv"
   "${SEDOV_WORKDIR}/probes.csv" COPYONLY)
add_test(NAME generate_sedov_mesh
         COMMAND $<TARGET_FILE:gen_sod_mesh> 40 40 40 1.5 1.5 1.5 "${SEDOV_WORKDIR}/sedov.msh"
         WORKING_DIRECTORY "${SEDOV_WORKDIR}")
add_test(NAME test_sedov
         COMMAND $<TARGET_FILE:test_sedov> "${SEDOV_WORKDIR}/sedov.nml"
         WORKING_DIRECTORY "${SEDOV_WORKDIR}")
set_tests_properties(test_sedov PROPERTIES DEPENDS generate_sedov_mesh)

# M5.2 free-field TNT blast — peak overpressure vs Kingery-Bulmash.
add_test(NAME test_eos_jwl
         COMMAND $<TARGET_FILE:test_eos_jwl>)

set(BLAST_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/blast_freefield")
file(MAKE_DIRECTORY "${BLAST_WORKDIR}")
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/blast_freefield/blast_freefield.nml"
   "${BLAST_WORKDIR}/blast_freefield.nml" COPYONLY)
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/blast_freefield/probes.csv"
   "${BLAST_WORKDIR}/probes.csv" COPYONLY)
add_test(NAME generate_blast_freefield_mesh
         COMMAND $<TARGET_FILE:gen_sod_mesh> 80 80 80 0.6 0.6 0.6 "${BLAST_WORKDIR}/blast_freefield.msh"
         WORKING_DIRECTORY "${BLAST_WORKDIR}")
add_test(NAME test_blast_overpressure
         COMMAND $<TARGET_FILE:test_blast_overpressure> "${BLAST_WORKDIR}/blast_freefield.nml"
         WORKING_DIRECTORY "${BLAST_WORKDIR}")
set_tests_properties(test_blast_overpressure PROPERTIES
   DEPENDS generate_blast_freefield_mesh)

# M5.3: external blast on a rigid cube. 0.6 m domain, 30^3 hex grid with
# a 0.1 m cube cut out at [0.4,0.5]x[0.25,0.35]^2. Slip-wall on cube
# faces, farfield outside. Smoke-test scope: confirm front-face shock
# reflection and back-face shadowing — full ConWep reflected-overpressure
# comparison requires finer mesh + full-scale charge.
set(BLAST_CUBE_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/blast_cube")
file(MAKE_DIRECTORY "${BLAST_CUBE_WORKDIR}")
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/blast_cube/blast_cube.nml"
   "${BLAST_CUBE_WORKDIR}/blast_cube.nml" COPYONLY)
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/blast_cube/probes.csv"
   "${BLAST_CUBE_WORKDIR}/probes.csv" COPYONLY)
add_test(NAME generate_blast_cube_mesh
         COMMAND $<TARGET_FILE:gen_box_with_cube>
                 30 30 30 0.6 0.6 0.6
                 0.4 0.25 0.25 0.1 0.1 0.1
                 "${BLAST_CUBE_WORKDIR}/blast_cube.msh"
         WORKING_DIRECTORY "${BLAST_CUBE_WORKDIR}")
add_test(NAME test_blast_cube
         COMMAND $<TARGET_FILE:test_blast_cube> "${BLAST_CUBE_WORKDIR}/blast_cube.nml"
         WORKING_DIRECTORY "${BLAST_CUBE_WORKDIR}")
set_tests_properties(test_blast_cube PROPERTIES
   DEPENDS generate_blast_cube_mesh)

# M5.4: internal blast in a vented room. 1 m^3 hex mesh (30^3); five
# faces are slip walls, one face has a 0.25 x 0.25 m vent (farfield).
# Smoke-test scope: confirm back-wall reflection and that the vent
# does not reflect (the near-vent peak is below the back-wall peak).
set(BLAST_ROOM_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/blast_room")
file(MAKE_DIRECTORY "${BLAST_ROOM_WORKDIR}")
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/blast_room/blast_room.nml"
   "${BLAST_ROOM_WORKDIR}/blast_room.nml" COPYONLY)
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/blast_room/probes.csv"
   "${BLAST_ROOM_WORKDIR}/probes.csv" COPYONLY)
add_test(NAME generate_blast_room_mesh
         COMMAND $<TARGET_FILE:gen_box_with_vent>
                 30 30 30 1.0 1.0 1.0
                 0.375 0.375 0.25 0.25
                 "${BLAST_ROOM_WORKDIR}/blast_room.msh"
         WORKING_DIRECTORY "${BLAST_ROOM_WORKDIR}")
add_test(NAME test_blast_room
         COMMAND $<TARGET_FILE:test_blast_room> "${BLAST_ROOM_WORKDIR}/blast_room.nml"
         WORKING_DIRECTORY "${BLAST_ROOM_WORKDIR}")
set_tests_properties(test_blast_room PROPERTIES
   DEPENDS generate_blast_room_mesh)

# Arbitrary-shaped room: 1 m^3 bbox, L-shaped polygonal floor, one
# 0.15 m cube obstacle in the bottom arm, one 0.2x0.2 m vent on the
# x_max face. Driven by a spec file consumed by gen_room_mesh.
set(BLAST_L_ROOM_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/blast_l_room")
file(MAKE_DIRECTORY "${BLAST_L_ROOM_WORKDIR}")
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/blast_l_room/blast_l_room.nml"
   "${BLAST_L_ROOM_WORKDIR}/blast_l_room.nml" COPYONLY)
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/blast_l_room/probes.csv"
   "${BLAST_L_ROOM_WORKDIR}/probes.csv" COPYONLY)
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/blast_l_room/blast_l_room.spec"
   "${BLAST_L_ROOM_WORKDIR}/blast_l_room.spec" COPYONLY)
add_test(NAME generate_blast_l_room_mesh
         COMMAND $<TARGET_FILE:gen_room_mesh>
                 "${BLAST_L_ROOM_WORKDIR}/blast_l_room.spec"
                 "${BLAST_L_ROOM_WORKDIR}/blast_l_room.msh"
         WORKING_DIRECTORY "${BLAST_L_ROOM_WORKDIR}")
add_test(NAME test_blast_l_room
         COMMAND $<TARGET_FILE:test_blast_l_room> "${BLAST_L_ROOM_WORKDIR}/blast_l_room.nml"
         WORKING_DIRECTORY "${BLAST_L_ROOM_WORKDIR}")
set_tests_properties(test_blast_l_room PROPERTIES
   DEPENDS generate_blast_l_room_mesh)

# Phase 1e cut-cell conservation test: a PLAIN box mesh (all cells
# present — no center-cutout) plus a namelist cube obstacle handled
# analytically by the cut-cell engine. Closed box -> mass + energy
# conserved to machine precision.
set(CUT_CELL_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/cut_cell_box")
file(MAKE_DIRECTORY "${CUT_CELL_WORKDIR}")
configure_file(
   "${CFD3D_SOURCE_DIR}/cases/cut_cell_box/cut_cell_box.nml"
   "${CUT_CELL_WORKDIR}/cut_cell_box.nml" COPYONLY)
add_test(NAME generate_cut_cell_box_mesh
         COMMAND $<TARGET_FILE:gen_sod_mesh> 20 20 20 1.0 1.0 1.0 "${CUT_CELL_WORKDIR}/cut_cell_box.msh"
         WORKING_DIRECTORY "${CUT_CELL_WORKDIR}")
add_test(NAME test_cut_cell_solver
         COMMAND $<TARGET_FILE:test_cut_cell_solver> "${CUT_CELL_WORKDIR}/cut_cell_box.nml"
         WORKING_DIRECTORY "${CUT_CELL_WORKDIR}")
set_tests_properties(test_cut_cell_solver PROPERTIES
   DEPENDS generate_cut_cell_box_mesh)

# M2 parallel acceptance: same Sod case run on 4 ranks. The mesh is partitioned
# via METIS, halos exchanged each RK stage, and the global L1 error vs. exact
# Sod must satisfy the same 6 % threshold as serial.
find_program(MPIEXEC_EXECUTABLE NAMES mpirun mpiexec)
if(MPIEXEC_EXECUTABLE)
   add_test(NAME test_sod_shock_tube_mpi4
            COMMAND ${MPIEXEC_EXECUTABLE} --oversubscribe -np 4
                    $<TARGET_FILE:test_sod_shock_tube> "${CASE_WORKDIR}/sod3d.nml"
            WORKING_DIRECTORY "${CASE_WORKDIR}")
   set_tests_properties(test_sod_shock_tube_mpi4 PROPERTIES
      DEPENDS generate_sod_mesh
      ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
endif()
