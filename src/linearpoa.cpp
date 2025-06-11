#include "linearpoa.h"
#include "poa.hpp"
#include "dp.hpp"
#include "seqio.hpp"

#include <cassert>
#include <cstring>

DLL void poa_generate_debug(char **seqs, int seq_num, dp_t O, dp_t E, dp_t M, dp_t X, int threads, int dp_threads, short multidp, char ***outseqs)
{
    std::unique_ptr<AlignmentEngine> alignment_engine = AlignmentEngine::Create(M, X, O, E);
    size_t max_len = 0;
    size_t *seq_lens = new size_t[seq_num + 1];
    for(int i = 0; i < seq_num; ++ i)
    {
        seq_lens[i] = strlen(seqs[i]);
        max_len = (std::max)(max_len, seq_lens[i]);
    }
    Graph G{}, G2{};
    Alignment align_res, align_res2;
    int32_t score, score_linear, score_raw;
    G.AddAlignment(align_res, seqs[0], seq_lens[0]);
    G2.AddAlignment(align_res2, seqs[0], seq_lens[0]);
    for(int i = 1; i < seq_num; ++ i)
    {
        // align_res = alignment_engine->linearspace_Align(seqs[i], seq_lens[i], G, &score, threads, dp_threads, multidp ? true : false);
        align_res = graph_to_sequence(seqs[i], &G, M, X, O, E, &score);
        align_res2 = alignment_engine->linearspace_Align(seqs[i], seq_lens[i], G2, &score_linear, threads, dp_threads, multidp ? true : false);
        score_linear = G.CalculateAlignmentScore(align_res2, seqs[i], seq_lens[i], O, E, M, X);
        score_raw    = G.CalculateAlignmentScore(align_res,  seqs[i], seq_lens[i], O, E, M, X);
        G.AddAlignment(align_res, seqs[i], seq_lens[i]);
        G2.AddAlignment(align_res2, seqs[i], seq_lens[i]);
        if (score_raw < score_linear)
        {
            fprintf(stderr, "score = %d (func calc = %d), score_linear = %d\n", score, score_raw, score_linear);
#if 1
            G.PrintDot("./1.dot");
            G2.PrintDot("./2.dot");
            auto tmp_res = G.GenerateMultipleSequenceAlignment(false);
            for (int j = 0; j <= i; ++j)
                fprintf(stderr, "%s\n", tmp_res[j].c_str());
            fprintf(stderr, "\n");
#endif
        }
        else fprintf(stderr, "same/better\n");
    }
    // process string to char*
    auto msa_res = G.GenerateMultipleSequenceAlignment(false);
    *outseqs = (char**)malloc(sizeof(char*) * (seq_num + 1));
    if (*outseqs == NULL) { fprintf(stderr, "Error: Can not alloc enough space for save results.\n"); return; }
    for(int i = 0; i < seq_num; ++ i)
    {
        char* correct_str = strdup(msa_res[i].c_str());
        gap_insert(seqs[i], correct_str, seq_lens[i], strlen(correct_str));
        *(*outseqs + i) = correct_str;
        msa_res[i].clear();
    }
    *(*outseqs + seq_num) = NULL;
}

DLL void poa_generate_with_consensus(char **seqs, int seq_num, dp_t O, dp_t E, dp_t M, dp_t X, int threads, int dp_threads, short multidp, short linear_method, short gen_msa, short gen_cons, char ***outseqs, char **cons, int min_cov)
{
    std::unique_ptr<AlignmentEngine> alignment_engine = AlignmentEngine::Create(M, X, O, E);
    size_t max_len = 0;
    size_t *seq_lens = new size_t[seq_num + 1];
    for(int i = 0; i < seq_num; ++ i)
    {
        seq_lens[i] = strlen(seqs[i]);
        max_len = (std::max)(max_len, seq_lens[i]);
    }
    std::vector <int> seq_ids(seq_num); std::iota(seq_ids.begin(), seq_ids.end(), 0);
    std::sort(seq_ids.begin(), seq_ids.end(), [&](int &a, int &b){ return seq_lens[a] > seq_lens[b]; });
    Graph G{};
    Alignment align_res;
    int32_t score;
    G.AddAlignment(align_res, seqs[seq_ids[0]], seq_lens[seq_ids[0]]);
    for(int i = 1; i < seq_num; ++ i)
    {
        if (linear_method) align_res = alignment_engine->linearspace_Align(seqs[seq_ids[i]], seq_lens[seq_ids[i]], G, &score, threads, dp_threads, multidp ? true : false);
        else align_res = graph_to_sequence(seqs[seq_ids[i]], &G, M, X, O, E, &score);
        G.AddAlignment(align_res, seqs[seq_ids[i]], seq_lens[seq_ids[i]]);
        if (i % 10 == 0 || i == seq_num - 1) fprintf(stderr, "Step [ %d / %d ]\n", i, seq_num - 1);
    }
    if(gen_msa)
    {
        // process string to char*
        auto msa_res = G.GenerateMultipleSequenceAlignment(false);
        *outseqs = (char**)malloc(sizeof(char*) * (seq_num + 1));
        if (*outseqs == NULL) { fprintf(stderr, "Error: Can not alloc enough space for save results.\n"); return; }
        char *correct_str;
        for(int i = 0; i < seq_num; ++ i)
        {
            correct_str = strdup(msa_res[i].c_str());
            gap_insert(seqs[seq_ids[i]], correct_str, seq_lens[seq_ids[i]], strlen(correct_str));
            *(*outseqs + seq_ids[i]) = correct_str;
            msa_res[i].clear();
        }
        *(*outseqs + seq_num) = NULL;
    }
    if(gen_cons)
    {
        auto cons_res = G.GenerateConsensus(min_cov);
        *cons = strdup(cons_res.c_str());
    }
    fprintf(stderr, "Done!\n");
}
