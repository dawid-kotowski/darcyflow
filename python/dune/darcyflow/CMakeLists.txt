# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_python_targets(darcyflow
  __init__
)
dune_add_pybind11_module(NAME _darcyflow)
target_link_libraries(_darcyflow PUBLIC Dune::DarcyFlow)
if(SKBUILD)
  install(TARGETS _darcyflow LIBRARY DESTINATION python/dune/darcyflow)
endif()
