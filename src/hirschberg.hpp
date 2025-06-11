#ifndef __HIRSCHBERG__
#define __HIRSCHBERG__

#include "dp.hpp"
#include "mtset.hpp"
#include "../include/dpthreadpool/BS_thread_pool_light.hpp"
#include <queue>
#include "poa.hpp"

using namespace BS;

void hirschberg_multi_init(ConcurrentSet_hirschberg_status_t *&S, int len, thread_pool_light *&pool1, thread_pool_light *&pool2, int dp_threads, std::mutex *&mtx1, std::mutex *&mtx2, bool multidp);
void hirschberg_multi_start(char *seq1, Graph* G, int len1, dp_t O, dp_t E, dp_t M, dp_t X, ConcurrentSet_hirschberg_status_t *&status_set, bool multidp, int dp_threads,
                            thread_pool_light *pool1, thread_pool_light *pool2, std::mutex *mtx1, std::mutex *mtx2, int threads);
void hirschberg_multi_free(ConcurrentSet_hirschberg_status_t *&S, thread_pool_light *pool1, thread_pool_light *pool2, std::mutex *mtx1, std::mutex *mtx2);

#endif