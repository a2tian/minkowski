#pragma once
#include "HeuristicSearch.h"

// function prototypes
namespace traversal {
// increment patterns
void cube(NTL::Vec<NTL::ZZ> &a, int b, int n);
void simplex(NTL::Vec<NTL::ZZ> &a, int b, int n);
void diagonal(NTL::Vec<NTL::ZZ> &a, double r, int n);
void rdiagonal(NTL::Vec<NTL::ZZ> &a, int p, int N, double r);
void center(NTL::Vec<NTL::ZZ> &a, int p, NTL::Vec<NTL::ZZ> &center, NTL::ZZ &r, int N);
}

double q(int p, int N);
double minkowski(int N);