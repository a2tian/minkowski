#pragma once
#include "HeuristicSearch.h"

// function prototypes
namespace traversal {
// increment patterns
void cube(NTL::Vec<NTL::ZZ> &a, int b, int n);
void simplex(NTL::Vec<NTL::ZZ> &a, int b, int n);
void diagonal(NTL::Vec<NTL::ZZ> &a, double r, int n);
void rdiagonal(NTL::Vec<NTL::ZZ> &a, int p, double r, int n);
void cdiagonal(NTL::Vec<NTL::ZZ> &a, double r, NTL::Vec<int> b1, NTL::Vec<int> b2, int n);
} // namespace traversal

double q(int p, int N);
double minkowski(int N);