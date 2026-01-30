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

    static const int dim = 2;
    using RF = double;

    // setup default parameters
    typename Dune::ParameterTree pTree;
    typename Dune::ParameterTreeParser pTreeParser;
    pTreeParser.readINITree("parameter.ini", pTree);
    pTreeParser.readOptions(argc, argv, pTree);

    // setup grid
    using Grid = typename Dune::YaspGrid<dim>;
    Dune::FieldVector<double, dim> domain(1.0);
    std::array<int, dim> domainDims;
    domain[0] = 1.0;
    domain[1] = 1.0;
    domainDims[0] = pTree.get<int>("grid.yasp_x");
    domainDims[1] = pTree.get<int>("grid.yasp_y");
    auto grid = std::make_shared<Grid>(domain, domainDims);
    
    // setup darcy solver
    using GV = Grid::LeafGridView;
    using DarcyProblemType = DarcyProblem<GV, RF>;
    DarcyProblemType darcyProblem(pTree);
    using DarcySolver = DarcySolver<GV, DarcyProblemType>;
    DarcySolver darcySolver(grid->leafGridView(), darcyProblem, pTree);
    darcySolver.solve();
    darcySolver.writeVTK();
  }
  catch (Dune::Exception &e)
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