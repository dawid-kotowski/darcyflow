// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_TRANSPORT_TRANSPORTSOLVER_HH
#define DUNE_DARCYFLOW_TRANSPORT_TRANSPORTSOLVER_HH

#include <dune/pdelab.hh>

#include <dune/darcyflow/utility/traits.hh>

/**
 * \brief Transport Solver using the information of a velocity field
 * 
 * \param GV GridView
 * \param Problem TransportProblem holding specified operators A, b and bctype
 * \param ParameterTree Dune::ParameterTree object holding parametric specification
 */
template<typename GV, typename Problem, typename DGTraits = DGTraits<GV>>
class TransportSolver
{

private:
  static const int dim_ = DGTraits::dim;
  static const int degree_ = DGTraits::order_dg;
  using DF = DGTraits::DomainField;
  using RF = DGTraits::RangeField;

  using DGFEM = DGTraits::FEM;
  using DGConstraints = DGTraits::ConstraintsType;
  using DGGFS = DGTraits::GFS;
  using DGVector = DGTraits::VectorType;
  using DGDGF = DGTraits::DiscreteGridFunction;

  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;

public:
  using SolutionType = std::vector<DGVector>;

public:
  TransportSolver(const GV& gv, Problem& problem, Dune::ParameterTree& pTree)
    : gv_(gv), problem_(problem), pTree_(pTree),
      dgfem_(), dggfs_(gv_, dgfem_)
  {
    dggfs_.name("DGTransport");
  }

  void logger(std::string message, Dune::Timer& timer, const int verbose = 0)
  {
    double time = timer.elapsed();
    if (!verbose == 0 )
    {
      std::cout << "(clock:"
      << time << ") === Transport Process Info === " << message << std::endl;
    }
  }

  /**
   * \brief Generic in-place Solver for current state of pTree
   */
  SolutionType solve()
  {
    // logger setup
    const int processVerb = 1;
    Dune::Timer timer;
    timer.start();

    // stationary local operator setup
    logger(std::string("Starting Transport Problem assembly ..."), timer, processVerb);
    using RhsLocalOperator = Dune::PDELab::ConvectionDiffusionDG<Problem,
                                                              typename DGTraits::FEM>;
    RhsLocalOperator rhslop(problem_, Dune::PDELab::ConvectionDiffusionDGMethod::IIPG,
                              Dune::PDELab::ConvectionDiffusionDGWeights::weightsOff, 0.0);

    // mass local operator setup
    using LhsLocalOperator = Dune::PDELab::L2;
    LhsLocalOperator lhslop(2*degree_);

    // constraints setup
    DGConstraints cc;
    cc.clear();

    // grid function init
    dggfs_.update();

    // matrix setup
    const int upperDofBound = std::pow(2, dim_) * dggfs_.maxLocalSize();
    MBE nativeMatrix(upperDofBound);

    // stationary grid operator setup
    using RhsGO = Dune::PDELab::GridOperator<DGGFS, DGGFS,
                                          RhsLocalOperator, MBE, DF, RF, RF,
                                          DGConstraints, DGConstraints>;
    RhsGO rhsgo(dggfs_, cc, dggfs_, cc, rhslop, nativeMatrix);

    // mass grid operator setup
    using LhsGO = Dune::PDELab::GridOperator<DGGFS, DGGFS,
                                          LhsLocalOperator, MBE, DF, RF, RF,
                                          DGConstraints, DGConstraints>;
    LhsGO lhsgo(dggfs_, cc, dggfs_, cc, lhslop, nativeMatrix);

    // Full grid operator setup
    using FullGO = Dune::PDELab::OneStepGridOperator<RhsGO, LhsGO>;
    FullGO go(rhsgo, lhsgo);
    using V = typename FullGO::Traits::Domain;

    // solution vector (maybe improve generally, since this is in the wrong "space")
    DGVector vOld(dggfs_);
    vOld = 0.0;

    // employ dirichlet conditions
    using BCExtender = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>;
    BCExtender g(gv_, problem_);
    Dune::PDELab::interpolate(g, dggfs_, vOld);

    // stationary linear solver setup
    const int lsMaxIterations = 10000;
    const int lsVerb = 0;
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_CG_ILU0;
    LinearSolver ls(lsMaxIterations, lsVerb);
    
    // stationary PDE solver setup
    using PDESolver = Dune::PDELab::StationaryLinearProblemSolver<FullGO, LinearSolver, V>;
    const RF pdeReduction = 1e-10;
    const RF pdeDefect = 1e-99;
    const int pdeVerbose = 0;
    PDESolver pdesolver(go, ls, pdeReduction, pdeDefect, pdeVerbose);

    // time stepper setup
    Dune::PDELab::ExplicitEulerParameter<RF> method;
    Dune::PDELab::OneStepMethod<RF, FullGO, PDESolver, V, V> osm(method, go, pdesolver);
    osm.setVerbosityLevel(0);
    logger(std::string("Assembled Problem."), timer, processVerb);

    // time stepping
    logger(std::string("Starting Time Solver Loop ..."), timer, processVerb);
    double time = pTree_.get<double>("time.time");
    double dt = pTree_.get<double>("time.dt");
    double T = pTree_.get<double>("time.T");
    SolutionType solutionTrajectory;
    while (time < T - 1e-10)
    {
      // time step
      V vNew(dggfs_, 0.0);
      osm.apply(time, dt, vOld, vNew);

      // increment
      solutionTrajectory.push_back(vNew);
      vOld = vNew;
      time += dt;

      // assemble constraints for new time step
      problem_.setTime(time);
      Dune::PDELab::constraints(g, dggfs_, cc);
    }
    logger(std::string("Computed Time Trajectory. Saved Solution."), timer, processVerb);
    double _ = timer.stop();

    return solutionTrajectory;
  }

  /**
   * \brief writer method given a trajectory of type std::vector<DGCoefficientType>
   */
  void writeVTK(const SolutionType& solutionTrajectory,
                std::string filename = "transportsolution")
  {
    // assertion for solver
    if (solutionTrajectory.empty())
        DUNE_THROW(Dune::Exception, "Run the solver first!");

    // vtk setup
    using SVTKWriter = Dune::SubsamplingVTKWriter<GV>;
    Dune::RefinementIntervals subsampling(pTree_.get<double>("visualization.subsamplingDG"));
    auto stationaryVTKWriter = std::make_shared<SVTKWriter>(gv_, subsampling);
    using VTKWriter = Dune::VTKSequenceWriter<GV>;
    std::string vtkfile(filename);
    VTKWriter vtkwriter(stationaryVTKWriter, vtkfile);
    DGVector solutionCoefficients(dggfs_, 0.0);
    DGDGF solution(dggfs_, solutionCoefficients);
    using VTKGridAdapter = Dune::PDELab::VTKGridFunctionAdapter<DGDGF>;
    vtkwriter.addVertexData(std::make_shared<VTKGridAdapter>(solution, "uh"));

    // time stepping
    double time = pTree_.get<double>("time.time");
    double dt = pTree_.get<double>("time.dt");
    double T = pTree_.get<double>("time.T");
    for (std::size_t timeStep = 0 ; timeStep < solutionTrajectory.size() ; ++timeStep )
    {
      solutionCoefficients = solutionTrajectory[timeStep];
      vtkwriter.write(time, Dune::VTK::appendedraw);
      time += dt;
    }
  }

    
private:
  const GV& gv_;
  Problem& problem_;
  typename Dune::ParameterTree& pTree_;
  DGFEM dgfem_;
  DGGFS dggfs_;
};

#endif // DUNE_DARCYFLOW_TRANSPORT_TRANSPORTSOLVER_HH