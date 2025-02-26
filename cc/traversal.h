#pragma once
#include "HeuristicSearch.h"

// function prototypes
namespace traversal {
// increment patterns
void cube(NTL::Vec<NTL::ZZ> &a, NTL::ZZ b, long n);
void simplex(NTL::Vec<NTL::ZZ> &a, NTL::ZZ b, long n);
void diagonal(NTL::Vec<NTL::ZZ> &a, double r, long n);
void rdiagonal(NTL::Vec<NTL::ZZ> &a, NTL::ZZ p, double r, long n);
void cdiagonal(NTL::Vec<NTL::ZZ> &a, double r, NTL::Vec<NTL::ZZ> b1, NTL::Vec<NTL::ZZ> b2, long n);
} // namespace traversal

NTL::RR q(NTL::ZZ p, long n);
NTL::RR minkowski(long n);