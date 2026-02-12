// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_UTILITY_PARAMETERSOLVER_HH
#define DUNE_DARCYFLOW_UTILITY_PARAMETERSOLVER_HH

#include <dune/pdelab.hh>

#include <dune/darcyflow.hh>

/**
 * \brief Parameterized Solver class for parametrized PDEs of Transport-DarcyFlow Type
 * 
 * \param GV GridView
 * \param ParameterTree Dune::ParameterTree object holding initial specification
 */
template<typename GV>
class ParameterSolver
{

private:
  using RF = double;

  using DarcyProblemType = DarcyProblem<GV, RF>;
  using DarcySolverType = DarcySolver<GV, DarcyProblemType>;
  using TransportProblemType = TransportProblem<GV, RF>;
  using DarcyGFS = typename DarcyTraits<GV>::GFS;
  using DarcyDGF = typename DarcyTraits<GV>::DiscreteGridFunction;
  using TransportProblemAdapterType = TransportProblemAdapter<DarcyDGF, GV, RF>;
  using TransportSolverType = TransportSolver<GV, TransportProblemAdapterType>;

public:
  using SolutionType = typename TransportSolverType::SolutionType;

public:
  ParameterSolver(const GV& gv, Dune::ParameterTree& pTree) 
    : pTree_(pTree),
      gv_(gv),
      darcyProblem_(pTree_),
      darcySolver_(gv_, darcyProblem_, pTree_),
      darcyDgf_(darcySolver_.getDiscreteGridFunction()),
      transportProblem_(pTree_),
      adaptedTransportProblem_(transportProblem_, darcyDgf_),
      transportSolver_(gv_, adaptedTransportProblem_, pTree_)
  {}

  const auto& getGfs()
  {
    return transportSolver_.getGfs();
  }

  template<typename Param>
  SolutionType solve(const Param& parameters)
  {
    // set up parameters into pTree
    ParameterParser parameterParser(pTree_);
    parameterParser.parse(parameters);
    pTree_ = parameterParser.getParameterTree();

    // update problems
    darcyProblem_.update();
    adaptedTransportProblem_.update();

    // solve problems
    darcySolver_.solve();
    adaptedTransportProblem_.updateDgf(darcyDgf_);
    return transportSolver_.solve();
  }

  void visualize(const SolutionType& solution, std::string filename = "transportsolution")
  {
    // update problem in case of gfs issues
    transportProblem_.update();

    // visualize
    transportSolver_.writeVTK(solution, filename);
  }

private:
  Dune::ParameterTree pTree_;
  const GV& gv_;
  DarcyProblemType darcyProblem_;
  DarcySolverType darcySolver_;
  DarcyDGF darcyDgf_;
  TransportProblemType transportProblem_;
  TransportProblemAdapterType adaptedTransportProblem_;
  TransportSolverType transportSolver_;
  
};

#endif // DUNE_DARCYFLOW_UTILITY_PARAMETERSOLVER_HH
