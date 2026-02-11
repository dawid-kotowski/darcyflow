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
    openingHeight_(pTree_.get<RF>("problem.parametric.openingHeight")),
    coatingHeight_(pTree_.get<RF>("problem.parametric.coatingHeight")),
    halfReactionBlockHeight_(0.5 - openingHeight_ - coatingHeight_)
  {}

  void update()
  {
    openingHeight_ = pTree_.get<RF>("problem.parametric.openingHeight");
    coatingHeight_ = pTree_.get<RF>("problem.parametric.coatingHeight");
    halfReactionBlockHeight_ = RF(0.5 - openingHeight_ - coatingHeight_);
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
    if (global[1] > 1-tol_)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else if (global[1] < tol_)
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

    // TODO: write this into some more parametric interface
    
    // Smooth compact bump centered at the top-left corner.
    const RF x0 = 0.0;
    const RF y0 = 1.0;
    const RF dx = global[0] - x0;
    const RF dy = global[1] - y0;
    const RF r = std::sqrt(dx*dx + dy*dy);

    // Use openingHeight as a natural length scale for the bump radius.
    const RF r0 = 0.5 * openingHeight_;
    if (r >= r0)
      return 0.0;

    const RF s = r / r0;
    return 0.5 * (1.0 + std::cos(M_PI * s));
  }

private:
  const double tol_ = 1e-10;
  Dune::ParameterTree& pTree_;
  RF openingHeight_;
  RF coatingHeight_;;
  RF halfReactionBlockHeight_;
};

#endif // DUNE_DARCYFLOW_TRANSPORT_TRANSPORTPROBLEM_HH
