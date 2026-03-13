// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_DARCYVELOCITY_DARCYSOLVER_HH
#define DUNE_DARCYFLOW_DARCYVELOCITY_DARCYSOLVER_HH

#include <dune/pdelab.hh>

#include <dune/darcyflow/darcyvelocity/darcyproblem.hh>
#include <dune/darcyflow/darcyvelocity/constraints.hh>
#include <dune/darcyflow/darcyvelocity/masslumping.hh>
#include <dune/darcyflow/darcyvelocity/reconstruction.hh>
#include <dune/darcyflow/darcyvelocity/schurcomplement.hh>
#include <dune/darcyflow/utility/traits.hh>

/**
 * \brief Darcy Solver using a mixed formulation
 * 
 * \param GV GridView
 * \param Problem DarcyProblem holding specified operators A, b and bctype
 * \param ParameterTree Dune::ParameterTree object holding parametric specification
 */
template<typename GV, typename Problem,
        typename DarcyTraits = DarcyTraits<GV>,
        typename DGTraits = DGTraits<GV>>
class DarcySolver
{

private:
  static const int dim = DarcyTraits::dim;
  using DF = DarcyTraits::DomainField;
  using RF = DarcyTraits::RangeField;

  using DarcyFEM = DarcyTraits::FEM;
  using DarcyConstraints = DarcyTraits::ConstraintsType;
  using DarcyGFS = DarcyTraits::GFS;
  using DarcyVector = DarcyTraits::VectorType;
  using DarcyDGF = DarcyTraits::DiscreteGridFunction;

  using DGFEM = DGTraits::FEM;
  using DGConstraints = DGTraits::ConstraintsType;
  using DGGFS = DGTraits::GFS;
  using DGVector = DGTraits::VectorType;
  using DGDGF = DGTraits::DiscreteGridFunction;

  // Tensor spaces
  using TensorVBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::bcrs>;
  using TensorGFS = Dune::PDELab::CompositeGridFunctionSpace<
                    TensorVBE, Dune::PDELab::LexicographicOrderingTag,
                    typename DarcyTraits::GFS,DGGFS>;

public:
  DarcySolver(const GV& gv, Problem& problem, Dune::ParameterTree& pTree)
    : gv_(gv), problem_(problem),
      pTree_(pTree), darcyfem_(gv_), dgfem_(),
      darcygfs_(gv_, darcyfem_),
      dggfs_(gv_, dgfem_), darcyCoefficients_(darcygfs_),
      dgCoefficients_(dggfs_)
  {
    darcygfs_.name("RT0");
    dggfs_.name("DG");
  }

  DarcyDGF getDiscreteGridFunction() const
  {
    return DarcyDGF(darcygfs_, darcyCoefficients_);
  }

  void logger(std::string message, Dune::Timer& timer, const int verbose = 0)
  {
    double time = timer.elapsed();
    if (!verbose == 0 )
    {
      std::cout << "(clock:"
      << time << ") === Darcy Process Info === " << message << std::endl;
    }
  }

  void solve()
  {
    const int processVerb = 0;
    Dune::Timer timer;
    timer.start();
    logger(std::string("Starting Darcy Problem assembly ..."), timer, processVerb);

    // make tensor gfs (requires tmp, otherwise GFS is changed after init of vectors)
    DarcyGFS darcygfsTmp(gv_, darcyfem_);
    DGGFS dggfsTmp(gv_, dgfem_);
    TensorGFS tensorgfs(darcygfsTmp, dggfsTmp);

    // make local operator
    using LocalOperator = Dune::PDELab::DiffusionMixed<Problem>;
    LocalOperator localOperator(problem_);

    // assemble constraints
    using BCAdapter = Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem>;
    BCAdapter bctype(gv_, problem_);
    using CC = typename TensorGFS::template ConstraintsContainer<RF>::Type;
    CC cc;
    cc.clear();
    Dune::PDELab::constraints(bctype, tensorgfs, cc);

    // make grid operator
    tensorgfs.update();
    const int upperDofBound = std::pow(2, dim) * tensorgfs.maxLocalSize();
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MBE matrix(upperDofBound);
    using GridOperator = Dune::PDELab::GridOperator<TensorGFS,
                                                    TensorGFS,
                                                    LocalOperator,
                                                    MBE,
                                                    DF, RF, RF, CC, CC>;
    GridOperator gridOperator(tensorgfs, cc, tensorgfs, cc, localOperator, matrix);

    // Solve matrix-based
    using XGFS = typename TensorGFS::template Child<0>::Type;
    using WrappedX = Dune::PDELab::Backend::Vector<XGFS, DF>;
    using X = Dune::PDELab::Backend::Native<WrappedX>;
    using YGFS = typename TensorGFS::template Child<1>::Type;
    using WrappedY = Dune::PDELab::Backend::Vector<YGFS, DF>;
    using Y = Dune::PDELab::Backend::Native<WrappedY>;

    // assemble rhs
    using RhsType = typename GridOperator::Range;
    using LhsType = typename GridOperator::Domain;
    LhsType zero(tensorgfs, 0.0);
    RhsType rhs(tensorgfs, 0.0);
    gridOperator.residual(zero, rhs);

    // assemble full matrix
    using JacobianType = typename GridOperator::Jacobian;
    JacobianType E(gridOperator, 0.0);
    gridOperator.jacobian(zero, E);

    // build dirichlet containers
    const auto flowFunction = [this](const auto& is, const auto& x) { 
      return problem_.j(is, x); 
    };
    const auto neutralFunction = [](const auto& x) { return 0.0; };
    auto darcyDirichlet = Dune::PDELab::makeGridFunctionFromCallable(gv_, flowFunction);
    auto dgDirichlet = Dune::PDELab::makeGridFunctionFromCallable(gv_, neutralFunction);
    const auto boundaryGridFunction = Dune::PDELab::CompositeGridFunction(darcyDirichlet, dgDirichlet);

    // assemble strong constraints
    RhsType rhsConstraints(tensorgfs, 0.0);
    Dune::PDELab::interpolate(boundaryGridFunction, tensorgfs, rhsConstraints);
    Dune::PDELab::copy_constrained_dofs(cc, rhsConstraints, rhs);
    const auto rhsInit = rhsConstraints;

    // correct rhs by strong nonzero dirichlet contributions
    using Dune::PDELab::Backend::native;
    Dune::PDELab::set_nonconstrained_dofs(cc, 0.0, rhsConstraints);
    RhsType rhsCorrection(tensorgfs, 0.0);
    native(E).mv(native(rhsConstraints), native(rhsCorrection));
    Dune::PDELab::set_constrained_dofs(cc, 0.0, rhsCorrection);
    rhs -= rhsCorrection;

    // switch to native matrices
    auto nativeE = native(E);
    auto darcyMatrix = nativeE[0][0];
    auto dgMatrix = nativeE[0][1];
    auto stabilityMatrix = nativeE[1][1];

    eliminateColumns(darcyMatrix, cc);

    logger(std::string("Assembled problem."), timer, processVerb);

    //[invert (thus solve) complete pressure problem using schurcomplement]
    logger(std::string("Starting Solver assembly ..."), timer, processVerb);
    using DarcyMatrix = decltype(darcyMatrix);
    using DGMatrix = decltype(dgMatrix);
    using StabilityMatrix = decltype(stabilityMatrix);
    Dune::MatrixAdapter<DarcyMatrix, X, X> darcyOperator(darcyMatrix);

    // invert the darcy matrix using preconditioned cg
    const RF darcyRelaxation = 1.0;
    using ILU = Dune::SeqILU<DarcyMatrix, X, X>;
    auto darcyPrec = std::make_shared<ILU>(darcyMatrix, darcyRelaxation);
    const RF darcyReduction = 1e-8;
    const int darcyMaxIterations = 10;
    const int darcyVerbose = 0;
    Dune::CGSolver<X> darcySolver(darcyOperator, *darcyPrec,
                                  darcyReduction, darcyMaxIterations, darcyVerbose);

    // perform mass lumping
    using LumpedMatrixType = typename MassLumpedMatrixType<DarcyMatrix, DGMatrix>::type;
    LumpedMatrixType approxDarcyMatrix;
    getMassLumpedMatrix(darcyMatrix, dgMatrix, approxDarcyMatrix);
    approxDarcyMatrix += stabilityMatrix;
    using LumpedOperatorType = Dune::MatrixAdapter<LumpedMatrixType, Y, Y>;
    auto approxDarcyOperator = std::make_shared<LumpedOperatorType>(approxDarcyMatrix);
    logger(std::string("Performed Mass Lumping."), timer, processVerb);

    // invert the lumped schur matrix for pressure field solution
    // set up AMG as preconditioner
    using Smoother = Dune::SeqSSOR<LumpedMatrixType, Y, Y>;
    using AMG = Dune::Amg::AMG<LumpedOperatorType, Y, Smoother>;
    using SmootherArgs = Dune::Amg::SmootherTraits<Smoother>::Arguments;
    SmootherArgs smootherArgs;
    smootherArgs.iterations = 1;
    smootherArgs.relaxationFactor = 1;
    const int maxAmgIterations = 15;
    const int coarsenTarget = 2000;
    Dune::Amg::Parameters amgParams(maxAmgIterations, coarsenTarget);
    amgParams.setDebugLevel(0);

    using Criterion = Dune::Amg::CoarsenCriterion<
                      Dune::Amg::SymmetricCriterion<LumpedMatrixType, 
                      Dune::Amg::FirstDiagonal >>;
    Criterion criterion(amgParams);

    auto dgPrec = std::make_shared<AMG>(*approxDarcyOperator, criterion, smootherArgs);

    // set up schurcomplement
    using SchurType = Schurcomplement<Y, X, DGMatrix, StabilityMatrix>;
    SchurType schurcomplement(darcySolver, dgMatrix, stabilityMatrix);
    
    // set up cg solver with schurcomplement and amg preconditioning
    const RF dgReduction = 1e-20;
    const int dgMaxIterations = 2000;
    const int dgVerbose = 0;
    Dune::CGSolver<Y> dgSolver(schurcomplement, *dgPrec,
                              dgReduction, dgMaxIterations, dgVerbose);

    // set up a second high-precision inversion of the darcy matrix
    const RF darcyPreciseReduction = 1e-12;
    const int darcyPreciseMaxIterations = 100;
    Dune::CGSolver<X> darcyPreciseSolver(darcyOperator, *darcyPrec,
                                        darcyPreciseReduction, darcyPreciseMaxIterations,
                                        darcyVerbose);
    Dune::InverseOperatorResult preciseRes;
    logger(std::string("Assembled Solvers."), timer, processVerb);

    // compute rhs term for the schurcomplement equation
    logger(std::string("Starting Inversion of Schurcomplement ..."), timer, processVerb);
    auto darcyrhs = native(rhs.block(0));
    auto dgrhs = native(rhs.block(1));
    X darcyrhsTmp(darcyMatrix.N());
    X darcyrhsCopy = darcyrhs;
    darcyPreciseSolver.apply(darcyrhsTmp, darcyrhsCopy, preciseRes);
    auto schurrhs = dgrhs;
    schurrhs *= -1.0;
    dgMatrix.umtv(darcyrhsTmp, schurrhs);

    // solve schurcomplement inversion in one go
    auto dgSolution = std::make_shared<Y>(dgMatrix.M());
    *dgSolution = 0.0;
    Dune::InverseOperatorResult res;
    dgSolver.apply(*dgSolution, schurrhs, res);
    dgCoefficients_.attach(dgSolution);
    logger(std::string("Finished Inversion."), timer, processVerb);
    //![invert (thus solve) complete pressure problem using schurcomplement]

    // reconstruct correct velocity field from pressure solution (dgSolution)
    logger(std::string("Starting Reconstruction of Velocity ..."), timer, processVerb);
    using ReconstructionType = Reconstruction<Y, X, DGMatrix>;
    ReconstructionType reconstruction(darcyPreciseSolver, dgMatrix, darcyrhs, -1);
    auto reconVector = reconstruction.getReconstructionVector();
    reconstruction.apply(*dgSolution, *reconVector);

    typename DarcyTraits::VectorType reconDarcyCoefficients(darcygfs_);
    reconDarcyCoefficients.attach(reconVector);
    darcyCoefficients_ = reconDarcyCoefficients;
    logger(std::string("Reconstructed Velocity. Saved Darcyflow."), timer, processVerb);
    double _ = timer.stop();
  }

  void writeVTK(std::string filename = "darcysolution")
  {
    using VTKWriter = Dune::SubsamplingVTKWriter<GV>;
    Dune::RefinementIntervals subsampling(pTree_.template get<double>("visualization.subsamplingVelocity"));
    VTKWriter vtkwriter(gv_, subsampling);
    std::string vtkfile(filename);

    // plot velocity
    DarcyDGF darcydgf(darcygfs_, darcyCoefficients_);
    using DarcyVTK = Dune::PDELab::VTKGridFunctionAdapter<DarcyDGF>;
    vtkwriter.addCellData(std::make_shared<DarcyVTK>(darcydgf, "Velocity"));

    // plot pressure
    DGDGF dgdgf(dggfs_, dgCoefficients_);
    using DGVTK = Dune::PDELab::VTKGridFunctionAdapter<DGDGF>;
    vtkwriter.addCellData(std::make_shared<DGVTK>(dgdgf, "Pressure"));

    // plot permeability field
    auto lambda = [this](const auto& el, const auto& x){ return this->problem_.A(el, x)[0][0]; };
    auto amag = Dune::PDELab::makeGridFunctionFromCallable(gv_,lambda);
    using VTKGridFunctionAdapter = Dune::PDELab::VTKGridFunctionAdapter<decltype(amag)>;
    vtkwriter.addCellData(std::make_shared<VTKGridFunctionAdapter>(amag, "Permeability"));

    vtkwriter.write(vtkfile, Dune::VTK::appendedraw);
  }

private:
  const GV& gv_;
  Problem& problem_;
  Dune::ParameterTree& pTree_;
  const DarcyFEM darcyfem_;
  const DGFEM dgfem_;
  DarcyGFS darcygfs_;
  DGGFS dggfs_;
  DarcyVector darcyCoefficients_;
  DGVector dgCoefficients_;
};


#endif // DUNE_DARCYFLOW_DARCYVELOCITY_DARCYSOLVER_HH
