#
# cfd3d tests
#

set(TEST_FFLAGS_COMMON ${CFD3D_FFLAGS_COMMON})
set(TEST_FFLAGS_REL    ${CFD3D_FFLAGS_RELEASE})
set(TEST_FFLAGS_DBG    ${CFD3D_FFLAGS_DEBUG})

function(add_cfd3d_test name)
   add_executable(${name} ${name}.f90)
   target_link_libraries(${name} PRIVATE cfd3d_core)
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
