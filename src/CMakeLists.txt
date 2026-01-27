# create custom target
add_executable(dune-darcyflow dune-darcyflow.cc)

# link target against found dune libraries (e.g. dune-common)
target_link_dune_default_libraries(dune-darcyflow)

# set include directories for the target
dune_default_include_directories(dune-darcyflow PRIVATE)

# link target against other dune registered packages (e.g. TBB, METIS, etc.)
dune_target_enable_all_packages(dune-darcyflow)
