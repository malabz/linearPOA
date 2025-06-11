#include "../include/kseq/kseq.h"
#include "dp.hpp"
#include "seqio.hpp"
#include "linearpoa.h"
#include <cstdio>
#include <string>

#if (defined(__linux__) || defined(__APPLE__))
#include <unistd.h>
#include <getopt.h>
#else
#include <io.h>
#include <process.h>
#include "../include/getopt9/include/getopt.h"
#endif

#define ___DEBUG_ 0

char* c_input = NULL, * c_output = NULL;
int O, E, M, X, threads;
char **seqs, **names;
short align_mode;
size_t len1 = 0, len2 = 0;
bool linear;

KSEQ_INIT(int, read)

void version()
{
    fprintf(stderr, "Version 1.0.0.0\n");
    exit(0);
}

void usage(char name[])
{
    fprintf(stderr, "POA alignment speedups based on parallel Hirschberg method with affine gap penalty\n");
    fprintf(stderr, "\n\tUsage: %s [options]\n\n", name);
    fprintf(stderr, "DP bits = %zu\n", sizeof(dp_t));
    fprintf(stderr, "Available options:\n");
    fprintf(stderr, "\t--in        FILE      sequence file name (Required)\n");
    fprintf(stderr, "\t--out       FILE      output file name (Required)\n");
    fprintf(stderr, "\t--threads   N         use N threads (N >= 1, default: 1)\n");
    fprintf(stderr, "\t--open      O         gap open penalty (default: 3)\n");
    fprintf(stderr, "\t--ext       E         gap extension penalty (default: 1)\n");
    fprintf(stderr, "\t--match     M         match score (default: 0)\n");
    fprintf(stderr, "\t--mismatch  X         mismatch score (default: 2)\n");
    fprintf(stderr, "\t--nolinear            do not use linear method (default: disabled)\n");
    fprintf(stderr, "\t--genmode   M         generate mode, 1: generate MSA, 2: generate consensus, 3: generate MSA+consensus (default: 1)\n");
    fprintf(stderr, "\t--help                print help message\n");
    fprintf(stderr, "\t--version             show program version\n");
    fprintf(stderr, "Example:\n\t%s --in seq.fasta --out seq_out.fasta\n", name);
}

int read_seqs()
{
    FILE* f_pointer = fopen(c_input, "r");
    if (f_pointer == NULL) { fprintf(stderr, "Error: file %s cannot open. Program will exit\n", c_input); exit(1); }
    kseq_t* file_t_ = kseq_init(fileno(f_pointer));
    int seq_cnt = 0, max_seqs = 2;
    seqs = (char**)malloc(max_seqs * sizeof(char*));
    names = (char**)malloc(max_seqs * sizeof(char*));
    while (kseq_read(file_t_) >= 0)
    {
        seq_cnt ++;
        if(max_seqs < seq_cnt)
        {
            max_seqs <<= 1;
            char **new_item = (char**)malloc(max_seqs * sizeof(char*));
            if(new_item == NULL) { fprintf(stderr, "Error: can not allocate enough memory. Program will exit.\n"); exit(1); }
            memcpy(new_item, seqs, (max_seqs >> 1) * sizeof(char*));
            free(seqs); seqs = new_item;
            new_item = (char**)malloc(max_seqs * sizeof(char*));
            if(new_item == NULL) { fprintf(stderr, "Error: can not allocate enough memory. Program will exit.\n"); exit(1); }
            memcpy(new_item, names, (max_seqs >> 1) * sizeof(char*));
            free(names); names = new_item;
        }
        seqs[seq_cnt - 1] = AllocateCharVec(file_t_ -> seq.l + 1);
        if (seqs[seq_cnt - 1] == NULL)
        {
            fprintf(stderr, "Error: Can not save sequence 2. Program will exit.\n");
            exit(1);
        }
        strncpy(seqs[seq_cnt - 1], file_t_->seq.s, file_t_->seq.l);
        seqs[seq_cnt - 1][file_t_->seq.l] = 0;
        names[seq_cnt - 1] = AllocateCharVec(file_t_ -> name.l + file_t_ -> comment.l + 1);
        strncpy(names[seq_cnt - 1], file_t_ -> name.s, file_t_ -> name.l);
        strncpy(names[seq_cnt - 1] + file_t_ -> name.l, file_t_ -> comment.s, file_t_ -> comment.l);
    }
    kseq_destroy(file_t_);
    fclose(f_pointer);
    return seq_cnt;
}

int main(int argc, char** argv)
{
    threads = 1; O = 6; E = 2; M = 0; X = 4; align_mode = 1; linear = true;
#if ( defined(_WIN32) || defined(_WIN64) )
    std::string consts[] = { "in", "out", "help", "version", "threads", "open", "ext", "match", "mismatch", "nolinear", "genmode" };
#endif
    int c;

    while (1)
    {
        int option_index = 0;
#if ( defined(_WIN32) || defined(_WIN64) )
        static struct option long_options[] =
        {
            {(char*)consts[0].c_str(),  required_argument, 0, 'i'},
            {(char*)consts[1].c_str(),  required_argument, 0, 'o'},
            {(char*)consts[2].c_str(),  no_argument,       0, 'h'},
            {(char*)consts[3].c_str(),  no_argument,       0, 'v'},
            {(char*)consts[4].c_str(),  required_argument, 0, 't'},
            {(char*)consts[5].c_str(),  required_argument, 0, 'O'},
            {(char*)consts[6].c_str(),  required_argument, 0, 'E'},
            {(char*)consts[7].c_str(),  required_argument, 0, 'M'},
            {(char*)consts[8].c_str(),  required_argument, 0, 'X'},
            {(char*)consts[9].c_str(),  no_argument,       0, 'N'},
            {(char*)consts[10].c_str(), required_argument, 0, 'g'},
            {0,                         0,                 0,  0 }
        };
#else
        static struct option long_options[] =
        {
            {"in",       required_argument, 0, 'i'},
            {"out",      required_argument, 0, 'o'},
            {"help",     no_argument,       0, 'h'},
            {"version",  no_argument,       0, 'v'},
            {"threads",  required_argument, 0, 't'},
            {"open",     required_argument, 0, 'O'},
            {"ext",      required_argument, 0, 'E'},
            {"match",    required_argument, 0, 'M'},
            {"mismatch", required_argument, 0, 'X'},
            {"nolinear", no_argument,       0, 'N'},
            {"genmode",  required_argument, 0, 'g'},
            {0,          0,                 0,  0 }
        };
#endif
        c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1) break;
        switch (c)
        {
        case 0:
            fprintf(stderr, "Not supported: option %s", long_options[option_index].name);
            if (optarg) fprintf(stderr, " with arg %s", optarg);
            fprintf(stderr, "\n");
            break;
        case 'h':
            usage(argv[0]);
            version();
            break;
        case 'i':
            c_input = argv[optind - 1];
            break;
        case 'o':
            c_output = argv[optind - 1];
            break;
        case 't':
            threads = atoi(argv[optind - 1]);
            break;
        case 'O':
            O = atoi(argv[optind - 1]);
            break;
        case 'E':
            E = atoi(argv[optind - 1]);
            break;
        case 'M':
            M = atoi(argv[optind - 1]);
            break;
        case 'X':
            X = atoi(argv[optind - 1]);
            break;
        case 'N':
            linear = false;
            break;
        case 'g':
            align_mode = atoi(argv[optind - 1]);
            if(!(1 <= align_mode && align_mode <= 3)) { fprintf(stderr, "Error: align mode is invalid. Please check your arugment and run again.\n"); exit(2); }
            break;
        case 'v':
            version();
        }
    }
    if(c_input == NULL || c_output == NULL)
    {
        fprintf(stderr, "Error: not specified input or output file. Check your arugments and try again.\n");
        exit(1);
    }
    int seq_cnt = read_seqs();
    if(seq_cnt <= 1) { fprintf(stderr, "Warning: detected less than 1 sequence(s). Program will exit.\n"); exit(0); }
    fprintf(stderr, "Total: %d strings\n", seq_cnt);
    char **out_seqs = nullptr, *cons = nullptr;
#if ___DEBUG_
    poa_generate_debug(seqs, seq_cnt, O, E, M, X, threads, 1, false, &out_seqs);
#else
    bool calc_msa = (align_mode & 1), calc_cons = ((align_mode >> 1) & 1);
    FILE *final_file = fopen(c_output, "wb");
    if(final_file == NULL) { fprintf(stderr, "Error: can not write to final file. Program will exit.\n"); exit(1); }
    poa_generate_with_consensus(seqs, seq_cnt, O, E, M, X, threads, 1, false, linear, calc_msa, calc_cons, &out_seqs, &cons);
    if(calc_msa)
    {
        for(int i = 0; i < seq_cnt; ++ i) fprintf(final_file, ">%s\n%s\n", names[i], out_seqs[i]);
    }
    if(calc_cons) fprintf(final_file, ">Cons\n%s\n", cons);
    fclose(final_file);
#endif
	return 0;
}
