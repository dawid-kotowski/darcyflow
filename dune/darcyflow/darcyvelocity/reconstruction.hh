// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_RECONSTRUCTION_HH
#define DUNE_DARCYFLOW_RECONSTRUCTION_HH

template<typename Y, typename X, typename MatrixType>
class Reconstruction : public Dune::LinearOperator<Y,X> {
public:
  using DomainType = Y;
  using RangeType = X;
  using field_type = typename Y::field_type;

  Reconstruction(Dune::InverseOperator<X,X>& invOp, const MatrixType& mat,
                 const int pm) : invOp_(invOp), mat_(mat), pm_(pm)
                 {
                    X tmp(mat.N());
                    tmp = 0.0;
                    rhsX_(tmp);
                 }

  Reconstruction(Dune::InverseOperator<X,X>& invOp, const MatrixType& mat,
                 const X& rhsX, const int pm) : invOp_(invOp), mat_(mat),
                                                   rhsX_(rhsX), pm_(pm) {}

  void apply(const Y& y, X& x) const {
    X temp(rhsX_);
    mat_.usmv(pm_, y, temp);

    Dune::InverseOperatorResult res;
    invOp_.apply(x, temp, res);
  }

  void applyscaleadd(field_type alpha, const Y& y, X& x) const {
    X temp1(rhsX_);
    mat_.usmv(pm_, y, temp1);

    X temp2(mat_.N());
    temp2 = 0.0;
    Dune::InverseOperatorResult res;
    invOp_.apply(temp2, temp1, res);
    x.axpy(alpha,temp2);
  }

  std::shared_ptr<X> getReconstructionVector() const {
    return std::make_shared<X>(rhsX_.size());
  }

  Dune::SolverCategory::Category category() const {
    return Dune::SolverCategory::sequential;
  }

private:
  Dune::InverseOperator<X,X>& invOp_;
  const MatrixType mat_;
  const X rhsX_;
  const int pm_;
};

#endif  // DUNE_DARCYFLOW_RECONSTRUCTION_HH
