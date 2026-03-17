// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <array>
#include <iostream>
#include <stdexcept>

#include <dune/pdelab.hh>
#include <dune/python/pybind11/embed.h>

#include <dune/darcyflow/pybind11/bindings.hh>

namespace py = pybind11;

template<class Vector>
double differenceNorm(const Vector& lhs, const Vector& rhs)
{
  if (lhs.size() != rhs.size())
    DUNE_THROW(Dune::Exception, "State sizes differ.");

  Vector diff(lhs);
  diff -= rhs;
  return diff.two_norm();
}

int main(int argc, char** argv)
{
  try
  {
    Dune::MPIHelper::instance(argc, argv);
    py::scoped_interpreter guard{};

    Dune::ParameterTree pTree;
    Dune::ParameterTreeParser pTreeParser;
    pTreeParser.readINITree("parameter.ini", pTree);
    pTreeParser.readOptions(argc, argv, pTree);

    // Init with py::dict
    py::dict config;
    config["reduction"] = pTree.get<double>("reduction");
    config["grid.yasp_x"] = 18;
    config["grid.yasp_y"] = 18;
    config["time.time"] = pTree.get<double>("time.time");
    config["time.solverSteps"] = 0.05;
    config["time.dt"] = 0.1;
    config["time.T"] = 0.4;
    config["problem.eta"] = pTree.get<double>("problem.eta");
    config["problem.inflowVelocity"] = pTree.get<double>("problem.inflowVelocity");
    config["problem.parametric.coatingHeight"] = pTree.get<double>("problem.parametric.coatingHeight");
    config["problem.parametric.inflowAngle"] = pTree.get<double>("problem.parametric.inflowAngle");
    config["problem.parametric.minReaction"] = pTree.get<double>("problem.parametric.minReaction");
    config["problem.parametric.coatingReaction"] = pTree.get<double>("problem.parametric.coatingReaction");
    config["problem.non-parametric.openingHeight"] = pTree.get<double>("problem.non-parametric.openingHeight");
    config["problem.non-parametric.minPermeability"] = pTree.get<double>("problem.non-parametric.minPermeability");
    config["problem.non-parametric.coatingPermeability"] = pTree.get<double>("problem.non-parametric.coatingPermeability");
    config["darcy.reduction"] = pTree.get<double>("darcy.reduction");
    config["visualization.subsampling"] = pTree.get<double>("visualization.subsampling");
    config["visualization.subsamplingVelocity"] = pTree.get<double>("visualization.subsamplingVelocity");
    config["visualization.subsamplingDG"] = pTree.get<double>("visualization.subsamplingDG");
    

    using Solver = Dune::DarcyFlow::Python::FullSolver<double, 2>;
    Solver solver(config);

    // set up two random parameters
    std::array<double, ParameterParser::parameterSize> muA{{0.10, 0.20, 0.15, 0.10}};
    std::array<double, ParameterParser::parameterSize> muB{{0.90, 0.85, 0.90, 0.95}};

    auto solA = solver.solve(muA);
    auto solB = solver.solve(muB);

    if (solA.size() < 2 || solB.size() < 2)
      DUNE_THROW(Dune::Exception, "Expected at least two saved timesteps in each trajectory.");

    // compare for solutiontrajectory
    const double paramDifference = differenceNorm(solA.back(), solB.back());
    const double timeDifferenceA = differenceNorm(solA.front(), solA.back());
    const double timeDifferenceB = differenceNorm(solB.front(), solB.back());

    std::cout << "paramDifference(final states) = " << paramDifference << std::endl;
    std::cout << "timeDifferenceA(first vs last) = " << timeDifferenceA << std::endl;
    std::cout << "timeDifferenceB(first vs last) = " << timeDifferenceB << std::endl;

    // visualize
    solver.visualize(solA, "pythontransportsolution");

    const double eps = 1e-12;
    if (paramDifference <= eps)
      DUNE_THROW(Dune::Exception, "Distinct parameters produced identical final states.");
    if (timeDifferenceA <= eps)
      DUNE_THROW(Dune::Exception, "Trajectory does not vary in time.");
    if (timeDifferenceB <= eps)
      DUNE_THROW(Dune::Exception, "Trajectory does not vary in time.");
  }
  catch (Dune::Exception& e)
  {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (std::exception& e)
  {
    std::cerr << "std::exception: " << e.what() << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }

  return 0;
}
