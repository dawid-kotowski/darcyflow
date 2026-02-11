# tests
dune_add_test(SOURCES testdarcyflow.cc)
dune_add_test(SOURCES testtransport.cc)
dune_add_test(SOURCES testparametersolver.cc)

# symlinks
dune_symlink_to_source_files(FILES parameter.ini)