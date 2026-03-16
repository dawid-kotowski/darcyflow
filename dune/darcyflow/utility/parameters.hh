// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_UTILITY_PARAMETERS_HH
#define DUNE_DARCYFLOW_UTILITY_PARAMETERS_HH

#include <dune/darcyflow.hh>

class ParameterParser
{

public:
  static const int parameterSize = 4;

public:
  ParameterParser(Dune::ParameterTree& pTree)
    : pTree_(pTree)
  {}

  template<typename ParamType>
  ParamType shiftScale(ParamType mu, const double a, const double b)
  {
    return static_cast<ParamType>(a) * mu + static_cast<ParamType>(b);
  }

  template<typename Param>
  void parse(const Param& mu)
  {
    pTree_["problem.parametric.coatingHeight"] = std::to_string(
      shiftScale(mu[0], 0.3, 0.0)); // [0, 1] -> [0, 0.3]
    pTree_["problem.parametric.inflowAngle"] = std::to_string(
      shiftScale(mu[1], M_PI, M_PI * -0.5)); // [0, 1] -> [-pi/2, pi/2]
    pTree_["problem.parametric.minReaction"] = std::to_string(
      shiftScale(mu[2], 0.3, 0.3)); // [0, 1] -> [0.3, 0.6]
    pTree_["problem.parametric.coatingReaction"] = std::to_string(
      shiftScale(mu[3], 0.2, 0.1)); // [0, 1] -> [0.1, 0.3]
  }

  Dune::ParameterTree getParameterTree() const
  {
    return pTree_;
  }

protected:
  Dune::ParameterTree pTree_;
};

#endif  // DUNE_DARCYFLOW_UTILITY_PARAMETERS_HH
