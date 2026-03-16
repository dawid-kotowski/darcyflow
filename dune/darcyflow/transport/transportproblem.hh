// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_TRANSPORT_TRANSPORTPROBLEM_HH
#define DUNE_DARCYFLOW_TRANSPORT_TRANSPORTPROBLEM_HH

#include <cmath>

#include <dune/grid/utility/hierarchicsearch.hh>
#include <dune/pdelab.hh>


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

  TransportProblem(Dune::ParameterTree& pTree) :
    Base(), pTree_(pTree),
    coatingReaction_(pTree_.get<RF>("problem.parametric.coatingReaction")),
    minReaction_(pTree_.get<RF>("problem.parametric.minReaction")),
    openingHeight_(pTree_.get<RF>("problem.non-parametric.openingHeight")),
    coatingHeight_(pTree_.get<RF>("problem.parametric.coatingHeight")),
    halfReactionBlockHeight_(0.5 - openingHeight_ - coatingHeight_)
  {}

  void update()
  {
    coatingReaction_ = pTree_.get<RF>("problem.parametric.coatingReaction");
    minReaction_ = pTree_.get<RF>("problem.parametric.minReaction");
    coatingHeight_ = pTree_.get<RF>("problem.parametric.coatingHeight");
    halfReactionBlockHeight_ = RF(0.5 - openingHeight_ - coatingHeight_);
  }

  template<typename X>
  bool isInflowBoundary(const X& global) const
  {
    return (global[0] < tol_) and (global[1] > 1 - openingHeight_ - tol_);
  }

  template<typename X>
  bool isOutflowBoundary(const X& global) const
  {
    return (global[0] > 1 - tol_) and (global[1] < openingHeight_ + tol_);
  }

  template<typename X>
  bool isOnOpening(const X& global) const
  {
    return isInflowBoundary(global) or isOutflowBoundary(global);
  }

  template<typename Element>
  auto elementBounds(const Element& el) const
  {
    const auto& geometry = el.geometry();
    auto lower = geometry.corner(0);
    auto upper = lower;
    for (int i = 1; i < geometry.corners(); ++i)
    {
      const auto corner = geometry.corner(i);
      for (int j = 0; j < Traits::dimDomain; ++j)
      {
        lower[j] = std::min(lower[j], corner[j]);
        upper[j] = std::max(upper[j], corner[j]);
      }
    }
    return std::make_pair(lower, upper);
  }

  template<typename Element>
  bool isInflowElement(const Element& el) const
  {
    const auto [lower, upper] = elementBounds(el);
    return (lower[0] < tol_) and (upper[1] > 1.0 - openingHeight_ - tol_);
  }

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
    if (isInflowBoundary(global))
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else if (isOutflowBoundary(global))
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
    const auto identity = typename Traits::RangeType({RF(1.0), RF(1.0)});
    using std::abs;
    if (abs(global[1] - 0.5) < halfReactionBlockHeight_ + tol_)
      return minReaction_ * identity;
    else if (abs(global[1] - 0.5) < halfReactionBlockHeight_ + coatingHeight_ + tol_)
      return coatingReaction_ * identity;
    else
      return identity * 0.0;
  }

  // Dirichlet condition
  template<typename Element, typename X>
  RF g (const Element& el, const X& x) const
  {
    auto global = el.geometry().global(x);
    const auto [lower, upper] = elementBounds(el);
    if (not isInflowElement(el))
      return 0.0;

    const RF center = 1.0 - 0.5 * openingHeight_;
    const RF radius = 0.5 * openingHeight_;
    const RF sample = std::min<RF>(std::max<RF>(global[1], std::max<RF>(lower[1], 1.0 - openingHeight_)),
                                   std::min<RF>(upper[1], 1.0));
    const RF r = std::abs(sample - center);
    if (r >= radius)
      return 0.0;

    const RF s = r / radius;
    return 0.5 * (1.0 + std::cos(M_PI * s));
  }

private:
  const double tol_ = 1e-10;
  Dune::ParameterTree& pTree_;
  RF minReaction_;
  RF coatingReaction_;
  RF openingHeight_;
  RF coatingHeight_;;
  RF halfReactionBlockHeight_;
};

#endif // DUNE_DARCYFLOW_TRANSPORT_TRANSPORTPROBLEM_HH
