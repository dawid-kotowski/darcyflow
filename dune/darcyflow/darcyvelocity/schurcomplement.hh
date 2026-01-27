// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_SCHURCOMPLEMENT_HH
#define DUNE_DARCYFLOW_SCHURCOMPLEMENT_HH

template<typename... T>
class Schurcomplement;

template<typename Y, typename X, typename BMatrixType>
class Schurcomplement<Y,X,BMatrixType> : public Dune::LinearOperator<Y,Y> {
public:
  using field_type = typename BMatrixType::field_type;

  Schurcomplement(Dune::InverseOperator<X,X>& aInvOp_, const BMatrixType& bMat_) : aInvOp(aInvOp_),
                  bMat(bMat_) {}

  void apply(const Y& x, Y& y) const {
    X by(bMat.N());
    X aby(bMat.N());  // initialize
    by = 0.0;
    aby = 0.0;

    bMat.mv(x, by);
    Dune::InverseOperatorResult res;
    aInvOp.apply(aby, by, res);
    bMat.mtv(aby, y);
  }

  void applyscaleadd(field_type alpha, const Y& x, Y& y) const {
    X by(bMat.N());
    X aby(bMat.N());  // initialize
    by = 0.0;
    aby = 0.0;

    bMat.mv(x, by);
    Dune::InverseOperatorResult res;
    aInvOp.apply(aby, by, res);
    bMat.usmtv(alpha, aby, y);
  }

  Dune::SolverCategory::Category category() const {
    return Dune::SolverCategory::sequential;
  }

private:
  Dune::InverseOperator<X,X>& aInvOp;
  const BMatrixType bMat; // TODO: replace by LinearOperator? how to transpose?
};

template<typename Y, typename X, typename BMatrixType, typename CMatrixType>
class Schurcomplement<Y,X,BMatrixType,CMatrixType> : public Dune::LinearOperator<Y,Y> {
public:
  using field_type = typename BMatrixType::field_type;

  Schurcomplement(Dune::InverseOperator<X,X>& aInvOp_, const BMatrixType& bMat_, const CMatrixType& cMat_) : aInvOp(aInvOp_),
                  bMat(bMat_), cMat(cMat_) {}

  void apply(const Y& x, Y& y) const {
    X by(bMat.N());
    X aby(bMat.N());  // initialize
    by = 0.0;
    aby = 0.0;

    bMat.mv(x, by);
    Dune::InverseOperatorResult res;
    aInvOp.apply(aby, by, res);
    bMat.mtv(aby, y);
    cMat.umv(x, y);
  }

  void applyscaleadd(field_type alpha, const Y& x, Y& y) const {
    X by(bMat.N());
    X aby(bMat.N());  // initialize
    by = 0.0;
    aby = 0.0;

    bMat.mv(x, by);
    Dune::InverseOperatorResult res;
    aInvOp.apply(aby, by, res);
    bMat.usmtv(alpha, aby, y);
    cMat.usmv(alpha, x, y);
  }

  Dune::SolverCategory::Category category() const {
    return Dune::SolverCategory::sequential;
  }

private:
  Dune::InverseOperator<X,X>& aInvOp;
  const BMatrixType bMat;
  const CMatrixType cMat;
};

#endif // DUNE_DARCYFLOW_SCHURCOMPLEMENT_HH
