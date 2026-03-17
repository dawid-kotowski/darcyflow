# tests
dune_add_test(SOURCES testdarcyflow.cc)
dune_add_test(SOURCES testtransport.cc)
dune_add_test(SOURCES testparametersolver.cc)
dune_add_test(SOURCES testpythonfullsolver.cc)
target_link_libraries(testpythonfullsolver PRIVATE Python3::Python)

# symlinks
dune_symlink_to_source_files(FILES parameter.ini)
