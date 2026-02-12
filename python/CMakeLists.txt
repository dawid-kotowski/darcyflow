# Add Symlinks for convenience

# Add configuration file, that puts module in python path
configure_file(pyenv.sh.in pyenv.sh)

# Add subdir
add_subdirectory(pymor)