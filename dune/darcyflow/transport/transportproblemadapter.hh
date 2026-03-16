// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

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
    : Base(transportProblem), velocityDgf_(&velocityDgf) {}

  void update()
  {
    Base::update();
  }

  void updateDgf(const DGFType& dgf)
  {
    velocityDgf_ = &dgf;
  }

  template<typename Element, typename X>
  auto bctype(const Element& el, const X& x) const
  {
    const auto& insidePos = el.geometryInInside().global(x);
    return Base::bctype(el.inside(), insidePos);
  }

  template<typename Element, typename X>
  auto b (const Element& el, const X& x) const
  {
    // warning: refined grids are not implemented!
    typename Traits::RangeType ret(0.0);
    velocityDgf_->evaluate(el, x, ret);
    return ret;
  }

private:
  const DGFType* velocityDgf_;
};

#endif  // DUNE_ULTRAWEAK_DARCY_VELOCITY_DARCY_VELOCITY_ADAPTER_HH
