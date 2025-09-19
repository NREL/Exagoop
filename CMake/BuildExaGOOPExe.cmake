function(build_exagoop_exe exagoop_exe_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)

  add_executable(${exagoop_exe_name} "")

  target_include_directories(${exagoop_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(${exagoop_exe_name} PRIVATE ${SRC_DIR})
  target_include_directories(${exagoop_exe_name} PRIVATE ${CMAKE_BINARY_DIR})

  target_sources(${exagoop_exe_name}
     PRIVATE
       ${SRC_DIR}/constants.H
       ${SRC_DIR}/mpm_diagnostics.cpp
       ${SRC_DIR}/mpm_init.cpp
       ${SRC_DIR}/mpm_init.cpp
       ${SRC_DIR}/mpm_particle_container.H
       ${SRC_DIR}/mpm_particle_timestep.cpp
       ${SRC_DIR}/nodal_data_ops.H
       ${SRC_DIR}/constitutive_models.H
       ${SRC_DIR}/mpm_eb.cpp
       ${SRC_DIR}/mpm_kernels.H
       ${SRC_DIR}/mpm_particle_grid_ops.cpp
       ${SRC_DIR}/mpm_specs.H
       ${SRC_DIR}/interpolants.H
       ${SRC_DIR}/mpm_check_pair.H
       ${SRC_DIR}/mpm_eb.H
       ${SRC_DIR}/mpm_particle_container.cpp
       ${SRC_DIR}/mpm_particle_outputs.cpp
       ${SRC_DIR}/nodal_data_ops.cpp
       ${SRC_DIR}/main.cpp
  )


  if(EXAGOOP_ENABLE_CUDA)
    set(pctargets "${exagoop_exe_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(EXAGOOP_SOURCES ${tgt} SOURCES)
      list(FILTER EXAGOOP_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${EXAGOOP_SOURCES} PROPERTIES LANGUAGE CUDA)
    endforeach()
    set_target_properties(${exagoop_exe_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(${exagoop_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas --disable-optimizer-constants>)
  endif()

  target_link_libraries(${exagoop_exe_name} PRIVATE amrex)

  install(TARGETS ${exagoop_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
