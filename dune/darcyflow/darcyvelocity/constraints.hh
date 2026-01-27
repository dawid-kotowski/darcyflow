// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DARCYFLOW_CONSTRAINTS_HH
#define DUNE_DARCYFLOW_CONSTRAINTS_HH

// ATTENTION: this is some hardcoded stuff that will likely not work in other contexts
// In particular, mat may at the moment not have block structure
// symmetrize the matrix mat by eliminating the columns for Dirichlet dofs
template<typename Mat, typename CC>
void eliminateColumns(Mat& mat, const CC& cc) {
  using BlockType = typename Mat::row_type::member_type;

  for (const auto& dof : cc) {
    auto idx = dof.first[0]; // TODO: not working for true multiindices!
    for (auto row=mat.begin(); row!=mat.end(); row++) {
      if(row.index() != idx) {
        auto it = row->find(idx);
        if(!it.equals(row->end()))
          *it = BlockType(0.0);
      }
    }
  }
}

#endif  // DUNE_DARCYFLOW_CONSTRAINTS_HH
