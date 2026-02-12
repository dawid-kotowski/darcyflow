# Add the Python Module to the library
add_library(ipydarcyflow SHARED dune-darcyflow.cc)

# Link into Python Lib
target_link_libraries(ipydarcyflow ${PYTHON_LIBRARIES})

# idk why this is (removes prefix?)
set_target_properties(ipydarcyflow PROPERTIES PREFIX "")