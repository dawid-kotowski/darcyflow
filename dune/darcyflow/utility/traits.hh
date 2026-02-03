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
  using FEM = Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV,DomainField,RangeField, order_rt>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,Dune::PDELab::RT0Constraints,VBE>;
  using ConstraintsType = typename GFS::template ConstraintsContainer<RangeField>::Type;
  using VectorType = Dune::PDELab::Backend::Vector<GFS, RangeField>;
  using DiscreteGridFunction = Dune::PDELab::DiscreteGridFunctionPiola<GFS, VectorType>;
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
  using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<DomainField,RangeField,order_dg, dim, Dune::PDELab::QkDGBasisPolynomial::lagrange>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,Dune::PDELab::NoConstraints,VBE>;
  using ConstraintsType = typename GFS::template ConstraintsContainer<RangeField>::Type;
  using VectorType = Dune::PDELab::Backend::Vector<GFS, DomainField>;
  using DiscreteGridFunction = Dune::PDELab::DiscreteGridFunction<GFS, VectorType>;
};

#endif // DUNE_DARCYFLOW_UTILITY_TRAITS_HH