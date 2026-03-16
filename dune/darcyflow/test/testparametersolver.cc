// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

#include <dune/darcyflow.hh>


int main(int argc, char** argv)
{
  try
  {
    Dune::MPIHelper::instance(argc, argv);

    const int dim = 2;
    using RF = double;

    // setup parameters
    Dune::ParameterTree pTree;
    Dune::ParameterTreeParser pTreeParser;
    pTreeParser.readINITree("parameter.ini", pTree);
    pTreeParser.readOptions(argc, argv, pTree);

    // setup grid
    using Grid = Dune::YaspGrid<dim>;
    Dune::FieldVector<RF, dim> domain(1.0);
    std::array<int, dim> domainDims;
    domainDims[0] = pTree.get<int>("grid.yasp_x");
    domainDims[1] = pTree.get<int>("grid.yasp_y");
    auto grid = std::make_shared<Grid>(domain, domainDims);
    using GV = Grid::LeafGridView;
    GV gv = grid->leafGridView();

    // setup solver
    using ParameterSolverType = ParameterSolver<GV>;
    ParameterSolverType parameterSolver(gv, pTree);

    // setup some parameters
    std::array<double, 4> parameter;
    parameter[0] = 0.25;
    parameter[1] = 0.25;
    parameter[2] = 0.2;
    parameter[3] = 0.05;
    
    // solve for that choice of parameters
    auto solution = parameterSolver.solve(parameter);
    parameterSolver.visualize(solution);
  }
  catch (Dune::Exception& e)
  {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
