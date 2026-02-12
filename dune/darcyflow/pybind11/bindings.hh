// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_PYBIND11_BINDINGS_HH
#define DUNE_DARCYFLOW_PYBIND11_BINDINGS_HH

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
            } catch (Dune::Exception& ex) {
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
        
        using ParameterType = std::array<RF, ParameterParser::parameterSize>;

      public:
        FullSolver(py::dict config)
          : pTree_(toParameterTree(config))
        {
          // grid setup
          Dune::FieldVector<double, dim> domain({1.0, 1.0});
          std::array<int, dim> domainDims = {
            domainDims[0] = pTree_.get<int>("grid.yasp_x"),
            domainDims[1] = pTree_.get<int>("grid.yasp_y")
          };

          grid_ = std::make_unique<Grid>(domain, domainDims);
          const auto gv = grid_->leafGridView();

          // set up solver
          solver_ = std::make_unique<Solver>(gv, pTree_);

          // set public members
          dimSource = solver_->getGfs().globalSize();
          dimRange = dimSource;
        }

        FullSolver(FullSolver&& other) = delete;

        SolutionType solve(const ParameterType& mu)
        {
          return solver_->template solve<ParameterType>(mu);
        }

        void visualize(const SolutionType& solution,
          const std::string filename="solution")
        {
          solver_->visualize(solution, filename);
        }

      protected:
        Dune::ParameterTree pTree_;
        std::unique_ptr<Solver> solver_;
        std::unique_ptr<Grid> grid_;

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
}

#endif  // DUNE_DARCYFLOW_PYBIND11_BINDINGS_HH
