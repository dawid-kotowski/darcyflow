#ifndef DUNE_ULTRAWEAK_DARCY_VELOCITY_DARCY_VELOCITY_ADAPTER_HH
#define DUNE_ULTRAWEAK_DARCY_VELOCITY_DARCY_VELOCITY_ADAPTER_HH

#include <dune/pdelab.hh>

#include <dune/darcyflow/transport/transportproblem.hh>

/**
 \brief TransportProblemAdapter for exposing darcy velocity from
        given DiscreteGridFunction

  \param transportProblem TransportProblem base containing non-velocity information
  \param velocityDgf DiscreteGridFunction of the velocity information
*/
template<typename DGFType, typename GV, typename RF>
class TransportProblemAdapter : public TransportProblem<GV, RF>
{
public:
  using Base = TransportProblem<GV, RF>;
  using Traits = typename Base::Traits;

  TransportProblemAdapter(const Base& transportProblem, const DGFType& velocityDgf)
    : Base(transportProblem), transportProblem_(transportProblem), velocityDgf_(velocityDgf) {}

  template<typename Element, typename X>
  auto bctype(const Element& el, const X& x) const
  {
    const auto& insidePos = el.geometryInInside().global(x);

    const auto& vel = b(el.inside(), insidePos);
    const auto val = vel.dot(el.unitOuterNormal(x));
    const double tol = 1e-10;
    if (val > tol)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;
    else if (val < -tol)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::None;
  }

  template<typename Element, typename X>
  auto b (const Element& el, const X& x) const
  {
    // warning: refined grids are not implemented!
    typename Traits::RangeType ret(0.0);
    velocityDgf_.evaluate(el, x, ret);
    return ret;
  }

private:
  const Base& transportProblem_;
  const DGFType& velocityDgf_;
};

#endif  // DUNE_ULTRAWEAK_DARCY_VELOCITY_DARCY_VELOCITY_ADAPTER_HH
