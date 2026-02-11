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

    // setup darcy problem and solver
    using GV = Grid::LeafGridView;
    using DarcyProblemType = DarcyProblem<GV, RF>;
    DarcyProblemType darcyProblem(pTree);
    using DarcySolverType = DarcySolver<GV, DarcyProblemType>;
    DarcySolverType darcySolver(grid->leafGridView(), darcyProblem, pTree);

    // solve darcy velocity
    darcySolver.solve();
    darcySolver.writeVTK(std::string("tmpdarcysolution"));
    auto velocityDgf = darcySolver.getDiscreteGridFunction();

    // setup transport problem
    using TransportProblemType = TransportProblem<GV, RF>;
    TransportProblemType baseTransport(pTree);
    using DGFType = DarcyTraits<GV>::DiscreteGridFunction;
    using TransportProblemAdapterType = TransportProblemAdapter<DGFType, GV, RF>;
    TransportProblemAdapterType transportProblem(baseTransport, velocityDgf);

    // setup transport solver
    using TransportSolverType = TransportSolver<GV, TransportProblemAdapterType>;
    TransportSolverType transportSolver(grid->leafGridView(), transportProblem, pTree);

    // solve transport problem with darcy velocity
    using SolutionType = TransportSolverType::SolutionType;
    SolutionType solutionTrajectory = transportSolver.solve();
    transportSolver.writeVTK(solutionTrajectory);
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