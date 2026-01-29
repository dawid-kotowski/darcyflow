// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_TRANSPORTSOLVER_HH
#define DUNE_DARCYFLOW_TRANSPORTSOLVER_HH

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
    TransportSolver(const GV& gv, Problem& problem, Dune::ParameterTree& pTree)
      : gv_(gv), problem_(problem), pTree_(pTree),
        dgfem_(gv_), dggfs_(gv_, dgfem_), solveflag_(false)
    {
        dggfs_.name("DGTransport");
    }

    void solve()
    {
        // stationary local operator setup
        using RhsLocalOperator = Dune::PDELab::ConvectionDiffusionDG<Problem,
                                                                  typename DGTraits::FEM>;
        RhsLocalOperator rhslop(problem_, 0.0);

        // mass local operator setup
        using LhsLocalOperator = Dune::PDELab::L2;
        LhsLocalOperator lhslop(2*degree_);

        // constraints setup
        DGConstraints cc;
        cc.clear();

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
        PDESolver pdesolver(go, ls, pdeReduction);

        // time stepper setup
        Dune::PDELab::ImplicitEulerParameter<RF> method;
        Dune::PDELab::OneStepMethod<RF, FullGO, PDESolver, V, V> osm(method, go, pdesolver);
        osm.setVerbosityLevel(1);

        // time stepping
        double time = pTree_.get<double>("time.time");
        double dt = pTree_.get<double>("time.dt");
        double T = pTree_.get<double>("time.T");
        while (time < T - 1e-10)
        {
            // time step
            V vNew(dggfs_, 0.0);
            osm.apply(time, dt, vOld, vNew);

            // increment
            solutionTrajectory_.push_back(vNew);
            vOld = vNew;
        }
        solveflag_ = true;
    }

    void writeVTK(std::string filename = "transportdgsolution")
    {
        // assertion for solver
        // warning: the time scaling should NEVER be subject to change in the pTree during compute
        if (!solveflag_)
            DUNE_THROW(Dune::Exception, "Run the solver first!");

        // vtk setup
        using SVTKWriter = Dune::SubsamplingVTKWriter<GV>;
        Dune::RefinementIntervals subsampling(pTree_.get<double>("visualization.subsamplingDG"));
        auto stationaryVTKWriter = std::make_shared<SVTKWriter>(gv_, subsampling);
        using VTKWriter = Dune::VTKSequenceWriter<GV>;
        VTKWriter vtkwriter(stationaryVTKWriter, "transportdgsolution", "", "");
        DGVector solutionCoefficients(dggfs_, 0.0);
        DGDGF solution(dggfs_, solutionCoefficients);
        using VTKGridAdapter = Dune::PDELab::VTKGridFunctionAdapter<DGDGF>;
        vtkwriter.addVertexData(std::make_shared<VTKGridAdapter>(solution, "uh"));

        // time stepping
        double time = pTree_.get<double>("time.time");
        double dt = pTree_.get<double>("time.dt");
        double T = pTree_.get<double>("time.T");
        for (int timeStep = 0 ; timeStep < solutionTrajectory_.size() ; ++timeStep )
        {
            solutionCoefficients = solutionTrajectory_[timeStep];
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
    std::vector<DGVector> solutionTrajectory_;
    bool solveflag_;
};

#endif // DUNE_DARCYFLOW_TRANSPORTSOLVER_HH