#pragma once
#include "HeuristicSearch.h"

// function prototypes

namespace ntl {
double norm(NTL::Vec<NTL::ZZ> &v);
}

namespace fplll {
double mat_norm(mat_zz &M, int n);

void gen_identity(mat_zz &M, int n, int p);
void init_matrix(mat_zz &M, NTL::Vec<NTL::ZZ> &a, int p, int n);
} // namespace fplll
