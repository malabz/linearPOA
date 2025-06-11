#include "russians.hpp"
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <string>
#include <vector>
#include <algorithm>

char *seq1, *seq2; // needed in dp.cpp

int main()
{
    dp_t M = 0, X = 6, O = 5, E = 3;
    sequence_t seq_type;
    char* table = nullptr; int table_size;
#if 1
    seq1 = (char*)"ggaugaccagccacacuggaacugagacacgguccagacuccuacgggaggcagcaguggggaauauugcacaaugggcgcaagccugaugcagccaugccgcguguaugaagaaggccu";
    seq2 = (char*)"aauguuggguuaagucccgcaacgagcgcaacccuuauccuuuguugccagcgguccggcaacgucggaagaccaaagagggggaccuucgggccucucgccaucggaugugccc";
    seq_type = DNA;
#elif 1
    seq1 = (char*)"SGNAKIGHPAPSFKATAVMPDGQFKDISLSDYKGKYVVFFFYPLDFTFVCPTEIIAFSDRAEEFKKLNCQVIGASVDSHFSHLAWINTPKKQGGLGPMNIPLVSDPKRTIAQDYGVLKADEGISFRGLFIIDDKGILRQITINDLPVGRSVDEILRLVQAFQFTDKHGEVCPA";
    seq2 = (char*)"LLLGDVAPNFEANTTVGRIRFHDFLGDSWGILFSHPRDFTPVTTELGRAAKLAPEFAKRNVKLIALSIDSVEDHLAWSKDINAYNSEEPTEKLPFPIIDDRNRELAILLGMLDPAEKDEKGMPVTARVVFVFGPDKKLKLSILYPATTGRNFDEILRVVISLQLTAEKRVATPVDWKDGDSVMVLPTIPEEEAKKLFPKGVFTKELPSGKKYLRYTPQP";
    seq_type = Protein;
#elif 0
    seq1 = (char*)"ATATCT";
    seq2 = (char*)"ATGTCTTC";
    seq_type = DNA;
#else
    seq1 = (char*)"AGTAC";
    seq2 = (char*)"AAG";
    seq_type = DNA;
#endif
    size_t len1 = strlen(seq1), len2 = strlen(seq2);
    table = russian_table_init(seq_type, table, table_size);
//#define BLOCK_SIMPLE
#ifdef BLOCK_SIMPLE
    F_simple_t F;
    dp_t *final_C = AllocateDPVec(len2 + 1);
    linear_block_simple_main_dp(seq1, len1, seq2, len2, M, X, O, 3, &F, final_C, table, table_size, false);
    for(size_t i = 0; i <= len2; ++ i) fprintf(stderr, "%lld ", final_C[i]); fputc('\n', stderr);
    dp_vs_t CC;
    simple_dp_russian(seq1, seq2, len1, len2, M, X, O, CC);
#else
    F_affine_t F;
    dp_t *final_C = AllocateDPVec(len2 + 1), *final_D = AllocateDPVec(len2 + 1);
    linear_block_affine_main_dp(seq1, len1, seq2, len2, M, X, O, E, O, 3, &F, final_C, final_D, table, table_size, false);
    fprintf(stderr, "C = "); for(size_t i = 0; i <= len2; ++ i) fprintf(stderr, "%lld ", final_C[i]); fputc('\n', stderr);
    fprintf(stderr, "D = "); for(size_t i = 0; i <= len2; ++ i) fprintf(stderr, "%lld ", final_D[i]); fputc('\n', stderr);
    dp_vs_t CC, DD, II, CDI;
    affine_dp_russian(seq1, seq2, len1, len2, M, X, O, E, CC, DD, II, CDI);
#endif

    return 0;
}
