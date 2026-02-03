// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_TRANSPORT_TRANSPORTPROBLEM_HH
#define DUNE_DARCYFLOW_TRANSPORT_TRANSPORTPROBLEM_HH

#include <cmath>

#include <dune/grid/utility/hierarchicsearch.hh>
#include <dune/pdelab.hh>

/**
 * \brief BoundaryProfiles for dirichlet conditioning of the Catalysator example
 */
namespace BoundaryProfiles {
  using RF = double;

  RF C2bump(const RF& x) {
    return (x>=0.25)*(x<=0.75)*(256*pow(x,4) - 512*pow(x,3) + 352*pow(x,2) - 96*x + 9);
  }

  RF sineBump(const RF& x, const unsigned int n=3) {
    return pow(sin(n*M_PI*x), 2);
  }

  RF L2bump(const RF& x, const unsigned int n=3) {
    using std::abs;
    const RF tol = 1e-10;
    return abs(n*x - int(n*x) - 0.5) <= 0.25 + tol;
  }
}

/**
 * \brief TransportProblem base for any stationary convection diffusion reaction
 *        equation. Contains standard initialization for the Catalysator example
 * 
 * \param pTree Dune::ParameterTree containing static parametrization information
 */
template <typename GV, typename RF>
class TransportProblem 
  : public Dune::PDELab::ConvectionDiffusionModelProblem<GV, RF>
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<RF>
{
public:
  using Base = Dune::PDELab::ConvectionDiffusionModelProblem<GV,RF>;
  using Traits = typename Base::Traits;

  TransportProblem(Dune::ParameterTree pTree) :
    Base(), pTree_(pTree),
    openingHeight_(pTree_.get<RF>("problem.openingHeight")),
    coatingHeight_(pTree_.get<RF>("problem.coatingHeight")),
    halfReactionBlockHeight_(0.5 - openingHeight_ - coatingHeight_)
  {}

  // no diffusion
  template<typename Element, typename X>
  auto A (const Element& el, const X& x) const
  {
    return typename Traits::PermTensorType(0.0);
  }

  // Boundary condition type
  template<typename Element, typename X>
  auto bctype(const Element& el, const X& x) const
  {
    auto global = el.geometry().global(x);
    const double tol = 1e-10;
    if (global[1] > 1-tol)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else if (global[1] < tol)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::None;
  }

  // Poiseuille profile
  template<typename Element, typename X>
  auto b (const Element& el, const X& x) const
  {
    const auto global = el.geometry().global(x);
    using std::abs;
    const auto r = abs(global[0] - 0.5);
    const double R = 0.5;
    const double eta = pTree_.template get<double>("problem.eta");
    using std::pow;
    return typename Traits::RangeType({0.0, -(pow(R,2) - pow(r,2))/(4*eta)});
  }

  // reaction coefficient
  template<typename Element, typename X>
  auto c (const Element& el, const X& x) const
  {
    const auto& global = el.geometry().center();
    using std::abs;
    auto d = abs(global[1]-0.5);
    return typename Traits::RangeType({RF(d <= halfReactionBlockHeight_),
      RF(d > halfReactionBlockHeight_ and d<=(halfReactionBlockHeight_ + coatingHeight_))});
  }

  // Dirichlet condition
  template<typename Element, typename X>
  RF g (const Element& el, const X& x) const
  {
    auto global = el.geometry().global(x);
    const int nBumps = pTree_.get<int>("problem.nInflowBumps");
    if (pTree_.get<bool>("problem.discontinuousInflow")){
      return BoundaryProfiles::L2bump(global[0], nBumps)
        +    BoundaryProfiles::L2bump(global[1], nBumps);
    }
    else {
      return BoundaryProfiles::sineBump(global[0], nBumps)
        +    BoundaryProfiles::sineBump(global[1], nBumps);
    }
  }

private:
  Dune::ParameterTree pTree_;
  const RF openingHeight_;
  const RF coatingHeight_;;
  const RF halfReactionBlockHeight_;
};

#endif // DUNE_DARCYFLOW_TRANSPORT_TRANSPORTPROBLEM_HH
