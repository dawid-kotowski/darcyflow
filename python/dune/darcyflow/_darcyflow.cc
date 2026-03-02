// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>
#include <dune/python/pybind11/numpy.h>

#include <dune/darcyflow/pybind11/bindings.hh>

namespace py = pybind11;

PYBIND11_MODULE(_darcyflow, m)
{
  m.doc() = "pybind11 dune-darcyflow plugin";

  registerSolverIntoModule(m);
}
