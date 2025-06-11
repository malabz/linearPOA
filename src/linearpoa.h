#ifndef __HIRSCHBERG_H__
#define __HIRSCHBERG_H__
#ifdef __cplusplus
extern "C"
{
#endif

#if (defined(_WIN32) || defined(_WIN64))
#ifdef RUSSIAN_LIB
#define DLL __declspec(dllexport)
#else
#define DLL __declspec(dllimport)
#endif // #ifdef RUSSIAN_LIB
#else
#define DLL
#endif

typedef int dp_t;

#ifndef sequence_t
typedef enum { DNA, Protein } sequence_t;
#endif

DLL void poa_generate_debug(char **seqs, int seq_num, dp_t O, dp_t E, dp_t M, dp_t X, int threads, int dp_threads, short multidp, char ***outseqs);
DLL void poa_generate_with_consensus(char **seqs, int seq_num, dp_t O, dp_t E, dp_t M, dp_t X, int threads, int dp_threads, short multidp, short linear_method, short gen_msa, short gen_cons, char ***outseqs, char **cons, int min_cov = -1);
#ifdef __cplusplus
}
#endif
#endif