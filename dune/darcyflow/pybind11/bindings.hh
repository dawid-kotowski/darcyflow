// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_PYBIND11_BINDINGS_HH
#define DUNE_DARCYFLOW_PYBIND11_BINDINGS_HH

#include <array>
#include <map>
#include <memory>
#include <sstream>
#include <string>

#include <dune/darcyflow.hh>

namespace py = pybind11;

namespace Dune
{
  namespace DarcyFlow
  {
    namespace Python
    {

      /**
       * \brief Convert any Handle to a String
       */
      std::string toString(py::handle handle)
      {
        std::string type = handle.get_type().attr("__name__").cast<std::string>();
        if (type == "str")
          return handle.cast<std::string>();
        else if (type == "int")
          return std::to_string(handle.cast<int>());
        else if (type == "float") {
          std::stringstream sstr;
          sstr << handle.cast<double>();
          return sstr.str();
        } else if (type == "bool")
          return std::to_string(handle.cast<bool>());
        else if ((type == "list") || (type == "tuple")) {
          std::stringstream str;
          unsigned int i = 0;
          for (auto it = handle.begin(); it != handle.end(); ++it, ++i) {
            if (i > 0)
              str << " ";
            str << toString(*it);
          }
          return str.str();
        }
        DUNE_THROW(Dune::Exception, "type \"" << type << "\" not supported");
      }

      /**
       * \brief Convert a Python Dict to a String Map
       */
      std::map<std::string, std::string> toStringMap(py::dict dict)
      {
        std::map<std::string, std::string> map;
        for (const auto& item : dict) {
          std::string type = item.second.get_type().attr("__name__").cast<std::string>();
          std::string key = toString(item.first);
          if (type == "dict") {
            auto sub = toStringMap(item.second.cast<py::dict>());
            for (const auto& k : sub) {
              map[key + "." + k.first] = k.second;
            }
          } else {
            try {
              map[key] = toString(item.second);
            } catch (Dune::Exception&) {
              // ignore the entry. will be triggered for numpy arrays, which we do not want to be
              // converted to string
            }
          }
        }
        return map;
      }

      /**
        \brief Converts a python dictionary to a Dune::ParameterTree.
               Copied from dune-hypercut
      */
      Dune::ParameterTree toParameterTree(py::dict dict)
      {
        Dune::ParameterTree tree;
        auto map = toStringMap(dict);
        for (const auto& k : map)
        {
          tree[k.first] = k.second;
        }
        return tree;
      }
      
      /**
       * \brief Operator Type for wrapping up abstract grid operators
       */
      template<typename Impl>
      class Operator {
      public:
        using DomainType = typename Impl::domain_type;
        using RangeType = typename Impl::range_type;

        Operator(const Impl& op) : op_(op),
          dimDomain(op.getmat().N()), dimRange(op.getmat().M())
        {}

        void apply(const DomainType& x, RangeType& y) const {
          op_.apply(x,y);
        }

      protected:
        const Impl& op_;
      public:
        std::size_t dimDomain;
        std::size_t dimRange;
      };

      /**
       * \brief FullSolver initializing Grid and ParameterSolver with a 
       *        Python dictionary
       * 
       * \param config Python Dictionary containing non-parametric problem specification
       */
      template<typename RF, int dim>
      class FullSolver
      {

      private:
        using Grid = Dune::YaspGrid<dim>;
        using GV = typename Grid::LeafGridView;
        
        using Solver = ParameterSolver<GV>;
        using SolutionType = typename Solver::SolutionType;
        using WrappedVectorType = typename SolutionType::value_type;
        using NativeVectorType = Dune::PDELab::Backend::Native<WrappedVectorType>;
        using PythonSolutionType = std::vector<NativeVectorType>;
        
        using ParameterType = std::array<RF, ParameterParser::parameterSize>;
        
      private:
        using MatrixBackendType = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
        using L2Type = Dune::PDELab::L2;
        using L2GridOperator = Dune::PDELab::GridOperator<
          typename DGTraits<GV>::GFS, typename DGTraits<GV>::GFS,
          L2Type, MatrixBackendType, RF, RF, RF,
          typename DGTraits<GV>::ConstraintsType, typename DGTraits<GV>::ConstraintsType>;
        using L2MatrixType = Dune::PDELab::Backend::Native<typename L2GridOperator::Jacobian>;

      public:
        FullSolver(py::dict config)
          : pTree_(toParameterTree(config))
        {
          // grid setup
          Dune::FieldVector<double, dim> domain(1.0);
          std::array<int, dim> domainDims;
          domainDims[0] = pTree_.get<int>("grid.yasp_x");
          domainDims[1] = pTree_.get<int>("grid.yasp_y");

          grid_ = std::make_unique<Grid>(domain, domainDims);
          gv_ = std::make_unique<GV>(grid_->leafGridView());

          // set up solver
          solver_ = std::make_unique<Solver>(*gv_, pTree_);

          // assemble L2-product on the test space
          using GFSDomain = typename L2GridOperator::Domain;
          GFSDomain zero(solver_->getGfs(), 0.0);
          MatrixBackendType mbe(1<<(dim+1));
          L2Type mass;
          auto l2massgo = std::make_shared<L2GridOperator>(
            solver_->getGfs(), solver_->getConstraints(), 
            solver_->getGfs(), solver_->getConstraints(),
            mass, mbe);
          typename L2GridOperator::Jacobian l2mat(*l2massgo, 0.0);
          l2massgo->jacobian(zero, l2mat);
          gramMatrix_ = Dune::PDELab::Backend::native(l2mat);

          // set public members
          dimSource = solver_->getGfs().globalSize();
          dimRange = dimSource;
        }

        FullSolver(FullSolver&& other) = delete;

        PythonSolutionType solve(const ParameterType& mu)
        {
          SolutionType solution = solver_->template solve<ParameterType>(mu);
          PythonSolutionType nativeSolution;
          nativeSolution.reserve(solution.size());
          for (auto& state : solution)
            nativeSolution.push_back(Dune::PDELab::Backend::native(state));
          return nativeSolution;
        }

        const L2MatrixType& getL2MassMatrix() const
        {
          return gramMatrix_;
        }

        void visualize(const PythonSolutionType& solution,
          const std::string filename="solution")
        {
          SolutionType wrappedSolution;
          wrappedSolution.reserve(solution.size());
          for (const auto& state : solution)
          {
            WrappedVectorType wrappedState(solver_->getGfs());
            Dune::PDELab::Backend::native(wrappedState) = state;
            wrappedSolution.push_back(std::move(wrappedState));
          }
          solver_->visualize(wrappedSolution, filename);
        }

      protected:
        Dune::ParameterTree pTree_;
        std::unique_ptr<Grid> grid_;
        std::unique_ptr<GV> gv_;
        std::unique_ptr<Solver> solver_;
        L2MatrixType gramMatrix_;

      public:
        std::size_t dimSource;
        std::size_t dimRange;
      };
    }
  }
}

/**
  \brief Registers the FullSolver class into the pybind11 module m
*/
void registerSolverIntoModule(py::module m)
{
  constexpr int dim = 2;
  using T = Dune::DarcyFlow::Python::FullSolver<double, dim>;

  const std::string clsName = "DarcyFlowSolver";
  auto cls = py::class_<T, std::shared_ptr<T>>(m, clsName.c_str());

  cls.def(py::init<py::dict>(), py::arg("config"));
  cls.def_readonly("dim_source", &T::dimSource);
  cls.def_readonly("dim_range", &T::dimRange);

  cls.def("solve", &T::solve, "solve for a given parameter set", py::arg("mu"));
  cls.def("visualize", &T::visualize, "write graphic into vtk",
    py::arg("solution"),
    py::arg("filename") = "solution");
  cls.def("getL2MassMatrix", &T::getL2MassMatrix,
    py::return_value_policy::reference_internal);
}

#endif  // DUNE_DARCYFLOW_PYBIND11_BINDINGS_HH
