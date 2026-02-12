// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_MASSLUMPING_HH
#define DUNE_DARCYFLOW_MASSLUMPING_HH


// TODO: do we want the A-Matrix dependency here?
template<typename AMatrixType, typename BMatrixType>
struct MassLumpedMatrixType {
  using type = typename Dune::TransposedMatMultMatResult<BMatrixType,BMatrixType>::type;
};

// TODO: what about a lower right block?
template<typename AMatrixType, typename BMatrixType>
auto getMassLumpedMatrix(const AMatrixType& aMat, const BMatrixType& bMat,
                         typename MassLumpedMatrixType<AMatrixType, BMatrixType>::type& lumpedMat) {
  using field_type = typename AMatrixType::field_type;

  BMatrixType scaledBMat(bMat);

  // perform mass lumping
  for (auto rowIt = aMat.begin(); rowIt != aMat.end(); ++rowIt) {
    field_type rowSum = 0.0;
    for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
      using std::abs;
      // TODO: allow other norms here
      rowSum += abs(*colIt); // l1-Norm of row
    }
    scaledBMat[rowIt.index()] /= rowSum;
  }

  Dune::transposeMatMultMat(lumpedMat, bMat, scaledBMat);
}

#endif // DUNE_DARCYFLOW_MASSLUMPING_HH
