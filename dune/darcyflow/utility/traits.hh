// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_UTILITY_TRAITS_HH
#define DUNE_DARCYFLOW_UTILITY_TRAITS_HH

/**
 * \brief Traits Struct for Raviart Thomas Function space of the Darcy Problem
 *        for darcy velocity
 */
template<typename GV>
struct DarcyTraits
{
  static const int dim = GV::dimension;

  using DomainField = typename GV::Grid::ctype;
  using RangeField = double;

  using VBE = Dune::PDELab::ISTL::VectorBackend<>;

  static const int order_rt = 0;
  using FEM = Dune::PDELab::LocalFiniteElementMap<GV,DF,RF, order_rt>;
  using ConstraintsType = Dune::PDELab::RT0Constraints;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,RT0CON,VBE>;
  using VectorType = Dune::PDELab::Backend::Vector<GridFunctionSpace, DF>;
  using DiscreteGridFunction = Dune::PDELab::DiscreteGridFunctionPiola<GridFunctionSpace, CoefficientVector>;
};

/**
 * \brief Traits Struct for Discontinuous Galerkin space of the Darcy Problem
 *        for pressure
 */
template<typename GV>
struct DGTraits
{
  static const int dim = GV::dimension;

  using DomainField = typename GV::Grid::ctype;
  using RangeField = double;

  using VBE = Dune::PDELab::ISTL::VectorBackend<>;

  static const int order_dg = 0;
  using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<DF,RF,order_dg, dim, Dune::PDELab::QkDGBasisPolynomial::lagrange>;
  using ConstraintsType = Dune::PDELab::NoConstraints;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,ConstraintsType,VBE>;
  using VectorType = Dune::PDELab::Backend::Vector<GridFunctionSpace, DF>;
  using DiscreteGridFunction = Dune::PDELab::DiscreteGridFunctionPiola<GridFunctionSpace, CoefficientVector>;
};

#endif // DUNE_DARCYFLOW_UTILITY_TRAITS_HH