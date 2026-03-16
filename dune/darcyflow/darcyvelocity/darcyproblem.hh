// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_DARCYVELOCITY_DARCYPROBLEM_HH
#define DUNE_DARCYFLOW_DARCYVELOCITY_DARCYPROBLEM_HH

#include <dune/pdelab.hh>

/**
   Pressure gradient from top left to bottom right with a horizontal
   low permeability block in the middle

   \brief Catalysator Example for Darcy Velocity

   \param pTree Dune::ParameterTree IF for problem specification
*/
template <typename GV, typename RF>
class DarcyProblem : public Dune::PDELab::ConvectionDiffusionModelProblem<GV, RF>
{
public:
  using Base = Dune::PDELab::ConvectionDiffusionModelProblem<GV, RF>;
  using Traits = typename Base::Traits;

  DarcyProblem(Dune::ParameterTree& pTree) : Base(),
    pTree_(pTree), I_(0.0),
    minPerm_(pTree_.get<RF>("problem.parametric.minPermeability")),
    coatingPerm_(pTree_.get<RF>("problem.parametric.coatingPermeability")),
    openingHeight_(pTree_.get<RF>("problem.non-parametric.openingHeight")),
    coatingHeight_(pTree_.get<RF>("problem.parametric.coatingHeight")),
    inflowAngle_(pTree_.get<RF>("problem.parametric.inflowAngle")),
    inflowVelocity_(pTree_.template get<RF>("problem.inflowVelocity")),
    halfReactionBlockHeight_(0.5 - openingHeight_ - coatingHeight_)
  {
    // precompute unity tensor
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      I_[i][i] = 1.0;
  }

  void update()
  {
    minPerm_ = pTree_.get<RF>("problem.parametric.minPermeability");
    coatingPerm_ = pTree_.get<RF>("problem.parametric.coatingPermeability");
    coatingHeight_ = pTree_.get<RF>("problem.parametric.coatingHeight");
    inflowAngle_ = pTree_.get<RF>("problem.parametric.inflowAngle");
    halfReactionBlockHeight_ = RF(0.5 - openingHeight_ - coatingHeight_);
  }

  // Boundary condition type
  template<typename Element, typename X>
  auto bctype(const Element& el, const X& x) const
  {
    auto global = el.geometry().global(x);
    if ((global[0] > 1-tol_) and (global[1] < openingHeight_ + tol_)){
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    }
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }

  // Permeability tensor (isotropic)
  template<typename Element, typename X>
  auto A (const Element& el, const X& x) const
  {
    const auto& global = el.geometry().center();

    using std::abs;
    if (abs(global[1] - 0.5) < halfReactionBlockHeight_ + tol_)
      return minPerm_ * I_;
    else if (abs(global[1] - 0.5) < halfReactionBlockHeight_ + coatingHeight_ + tol_)
      return coatingPerm_ * I_;
    else
      return I_;
  }

  // Dirichlet condition
  template<typename Element, typename X>
  auto g (const Element& el, const X& x) const
  {
    const auto global = el.geometry().global(x);
    const auto gval = 1.0 - global[0];
    return -gval;
  }

  //! Neumann boundary condition
  template<typename Intersection, typename X>
  auto j (const Intersection& is, const X& x) const
  {
    const auto global = is.geometry().global(x);

    using Vec = Dune::FieldVector<double,2>;
    const double eta = pTree_.template get<double>("problem.eta");
    using std::cos;
    using std::sin;
    const Vec direction({cos(inflowAngle_), sin(inflowAngle_)});

    // Poiseuille profile
    auto profile = [eta, this, direction](const auto& r){
      using std::pow;
      return direction * (inflowVelocity_ * (pow(0.5*openingHeight_,2) - pow(r,2)) / (4*eta));
    };

    using std::abs;
    Vec ret(0.0);
    if ((global[0] < tol_ and global[1] > 1-openingHeight_ - tol_)) {
      const auto& radius = abs(global[1] - (1-0.5*openingHeight_));
      ret = profile(radius);
    }
    else if ((global[0] > 1-tol_ and global[1] < openingHeight_ + tol_)) {
      const auto& radius = abs(global[1] - (0.5*openingHeight_));
      ret = profile(radius);
    }

    return ret;
  }

private:
  const double tol_ = 1e-10;
  typename Dune::ParameterTree& pTree_;
  typename Traits::PermTensorType I_;
  RF minPerm_;
  RF coatingPerm_;
  RF openingHeight_;
  RF coatingHeight_;
  RF inflowAngle_;
  RF inflowVelocity_;
  RF halfReactionBlockHeight_;
};

#endif // DUNE_DARCYFLOW_DARCYVELOCITY_DARCYPROBLEM_HH
