#ifndef __DP_HPP__
#define __DP_HPP__

#include "mtxutl.hpp"
#include "poa.hpp"
#include <algorithm>
#include <mutex>
#include <cstring>
#include <iostream>
#include <future>
#include <set>
#include <vector>
#include <tuple>
#include <cstdio>

#define DP_LONG_LEN 0

typedef int dp_t;

typedef std::vector <dp_t> dp_v_t;
typedef std::set <dp_v_t> dp_vs_t;
typedef std::set <dp_t> dp_s_t;

void print_vec(dp_t *v, size_t len);

char *strrev_(char *str);
dp_t* AllocateDPVec(size_t v);
dp_t readDPdata(FILE *F);
void printDPdata(dp_t x, FILE *F);
void FreeDPVec(dp_t *v);
void linear_dp_hirschberg(const char* seq1, const char* seq2, dp_t* CC, dp_t* DD, const int line1, const int len2, dp_t& M, dp_t& X, dp_t& O, dp_t& E, dp_t& tg);
void rev_linear_dp_hirschberg(const char* seq1, const char* seq2, dp_t* CC, dp_t* DD, const int line1, const int len1, const int len2, dp_t& M, dp_t& X, dp_t& O, dp_t& E, dp_t& tg);

void linear_poa_dp_hirschberg(const char* seq1, const Graph* g2, dp_t* CC, dp_t* DD, const int line1, dp_t& M, dp_t& X, dp_t& O, dp_t& E, dp_t& tg);
void rev_linear_poa_dp_hirschberg(const char* seq1, const Graph *g2, dp_t* CC, dp_t* DD, const int line1, const int len1, dp_t& M, dp_t& X, dp_t& O, dp_t& E, dp_t& tg);

Alignment graph_to_sequence(const char *seq1, const Graph *g2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, int *score);

#endif
