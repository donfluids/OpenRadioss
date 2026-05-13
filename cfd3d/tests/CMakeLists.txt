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
