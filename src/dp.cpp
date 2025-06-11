#include "dp.hpp"
#define DP_DEBUG 0
#define DP_PRINT_MARTIX 0
#define DP_DEBUG_CONDITION 0

#define MATCH_SCORE(a, b, M, X) ((a) == (b) ? (M) : (X))

const dp_t max_score = std::numeric_limits<dp_t>::max() - 1024;

inline dp_t** AllocateDPMtx(size_t a, size_t b)
{
#if DP_LONG_LEN == 2
    return AllocateInt128Mtx(a, b);
#elif DP_LONG_LEN == 1
    return AllocateLongLongMtx(a, b);
#else
    return AllocateIntMtx(a, b);
#endif
}

inline void FreeDPMtx(dp_t **m)
{
#if DP_LONG_LEN == 2
    FreeInt128Mtx(m);
#elif DP_LONG_LEN == 1
    FreeLongLongMtx(m);
#else
    FreeIntMtx(m);
#endif
}

dp_t* AllocateDPVec(size_t v)
{
#if DP_LONG_LEN == 2
    return AllocateInt128Vec(v);
#elif DP_LONG_LEN == 1
    return AllocateLongLongVec(v);
#else
    return AllocateIntVec(v);
#endif
}

dp_t readDPdata(FILE *F)
{
    dp_t x;
#if DP_LONG_LEN == 2
    x = 0;
    dp_t f = 1;
    char c = fgetc(F);
    while(! isdigit(c)){ if(c == '-') f = -1;  c = fgetc(F); }
    while(isdigit(c))  { x = x * 10 + c - '0'; c = fgetc(F); }
#elif DP_LONG_LEN == 1
    fscanf(F, "%lld", &x);
#else
    fscanf(F, "%d", &x);
#endif
    return x;
}

void printDPdata(dp_t x, FILE *F)
{
#if DP_LONG_LEN == 2
    fprintf(F, "%s", x.str().c_str());
#elif DP_LONG_LEN == 1
    fprintf(F, "%lld", x);
#else
    fprintf(F, "%d", x);
#endif
}

void FreeDPVec(dp_t *v)
{
#if DP_LONG_LEN == 2
    FreeInt128Vec(v);
#elif DP_LONG_LEN == 1
    FreeLongLongVec(v);
#else
    FreeIntVec(v);
#endif

}

void print_matrix(dp_t **mat, size_t len1, size_t len2)
{
    for(size_t i = 0; i <= len1; ++ i, std::cerr << "\n")
        for(size_t j = 0; j <= len2; ++ j, std::cerr << ' ')
            std::cerr << mat[i][j];
    std::cerr << std::endl;
}

void print_vec(dp_t *v, size_t len)
{
    for(size_t i = 0; i <= len; ++ i, std::cerr << ' ') std::cerr << v[i];
    std::cerr << std::endl;
}

char *strrev_(char *str)
{
    // code from https://stackoverflow.com/a/8534275 - it works!
    char *p1, *p2;
    if (! str || ! *str)
        return str;
    for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
    {
        *p1 ^= *p2;
        *p2 ^= *p1;
        *p1 ^= *p2;
    }
    return str;
}

void linear_dp_hirschberg(const char* seq1, const char* seq2, dp_t *CC, dp_t *DD, const int line1, const int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t &tg)
{
    dp_t t, s, e, c;
    CC[0] = 0;
    t = O;
    for(int j = 1; j <= len2; ++ j)
    {
        t = t + E;
        CC[j] = t;
        DD[j] = t + O;
    }
#if DP_DEBUG
    fprintf(stderr, "CC[0] = "); print_vec(CC, len2);
    fprintf(stderr, "DD[0] = "); print_vec(DD, len2);
#endif
    t = tg; // [*]
    for(int i = 1; i <= line1; ++ i)
    {
        s = CC[0];
        t = t + E;
        c = t;
        CC[0] = c; DD[0] = CC[0];
        e = t + O;
        for(int j = 1; j <= len2; ++ j)
        {
            e = std::min(e, c + O) + E;
            DD[j] = std::min(DD[j], CC[j] + O) + E;
            c = std::min(std::min(DD[j], e), s + MATCH_SCORE(seq1[i - 1], seq2[j - 1], M, X));
            s = CC[j];
            CC[j] = c;
        }
#if DP_DEBUG
        fprintf(stderr, "CC[%zu] = ", i); print_vec(CC, len2);
        fprintf(stderr, "DD[%zu] = ", i); print_vec(DD, len2);
#endif
    }
}

void rev_linear_dp_hirschberg(const char* seq1, const char* seq2, dp_t *CC, dp_t *DD, const int line1, const int len1, const int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t &tg)
{
    dp_t t, s, e, c;
    CC[0] = 0;
    t = O;
    for(int j = 1; j <= len2; ++ j)
    {
        t = t + E;
        CC[j] = t;
        DD[j] = t + O;
    }
#if DP_DEBUG
    fprintf(stderr, "CC[0] = "); print_vec(CC, len2);
    fprintf(stderr, "DD[0] = "); print_vec(DD, len2);
#endif
    t = tg; // [*]
    for(int i = 1; i <= line1; ++ i)
    {
        s = CC[0];
        t = t + E;
        c = t;
        CC[0] = c; DD[0] = CC[0];
        e = t + O;
        for(int j = 1; j <= len2; ++ j)
        {
            e = std::min(e, c + O) + E;
            DD[j] = std::min(DD[j], CC[j] + O) + E;
            c = std::min(std::min(DD[j], e), s + MATCH_SCORE(seq1[len1 - i], seq2[len2 - j], M, X));
            s = CC[j];
            CC[j] = c;
        }
#if DP_DEBUG
        fprintf(stderr, "CC[%zu] = ", i); print_vec(CC, len2);
        fprintf(stderr, "DD[%zu] = ", i); print_vec(DD, len2);
#endif
    }
}

void linear_poa_dp_hirschberg(const char* seq1, const Graph* g2, dp_t* CC, dp_t* DD, const int line1, dp_t& M, dp_t& X, dp_t& O, dp_t& E, dp_t& tg)
{
    auto const &graphorder = g2->rank_to_node();
    int len2 = graphorder.size();
    std::vector<std::uint32_t> node_to_rank(len2 + 1, 0);
    for(int i = 0; i < len2; ++ i) node_to_rank[graphorder[i]->id] = i;
    dp_t *II = AllocateDPVec(len2 + 1), *CC_pre = AllocateDPVec(len2 + 1);
    dp_t t, c;
    // init
    CC[0] = 0;
    t = O;
    for(int j = 1; j <= len2; ++ j)
    {
        CC[j] = max_score;
        for(const auto &item: graphorder[j - 1]->inedges)
        {
            if(! item->tail->code)
            {
                CC[j] = std::min(CC[j], O + E);
                // fprintf(stderr, "ID: %d, pre: %d\n", j - 1, 0);
                continue;
            }
            auto pre_id = node_to_rank[item->tail->id] + 1;
            CC[j] = std::min(CC[j], CC[pre_id] + E);
        }
        if(CC[j] == max_score) CC[j] = O + E;
        DD[j] = CC[j] + O;
    }
#if DP_DEBUG
    fprintf(stderr, "CC[0] = "); print_vec(CC, len2);
    fprintf(stderr, "DD[0] = "); print_vec(DD, len2);
#endif
    t = tg; // [*]
    for(int i = 1; i <= line1; ++ i)
    {
        CC_pre[0] = CC[0];
        CC[0] = c = t = t + E;
        II[0] = t + O;
        for(int j = 1; j <= len2; ++ j)
        {
            DD[j] = std::min(DD[j], CC[j] + O) + E;
            c = II[j] = max_score;
            const auto &this_node = graphorder[j - 1];
            for(const auto &item: this_node->inedges)
            {
                auto pre_id = item->tail->code ? node_to_rank[item->tail->id] + 1 : 0;
                II[j] = std::min(std::min(II[pre_id], CC[pre_id] + O) + E, II[j]);
                c = std::min(c, CC_pre[pre_id] + MATCH_SCORE(seq1[i - 1], g2->decoder(this_node->code), M, X));
#if DP_DEBUG
                fprintf(stderr, "seq1[%d] = %c, g2[%d] = %c\n", i - 1, seq1[i - 1], item->tail->id, g2->decoder(this_node->code));
#endif
            }
            c = std::min(std::min(DD[j], II[j]), c);
            CC_pre[j] = CC[j];
            CC[j] = c;
        }
#if DP_DEBUG
        fprintf(stderr, "CC[%d] = ", i); print_vec(CC, len2);
        fprintf(stderr, "DD[%d] = ", i); print_vec(DD, len2);
#endif
    }
    DD[0] = CC[0];
    FreeDPVec(II);
    FreeDPVec(CC_pre);
}

void rev_linear_poa_dp_hirschberg(const char* seq1, const Graph *g2, dp_t* CC, dp_t* DD, const int line1, const int len1, dp_t& M, dp_t& X, dp_t& O, dp_t& E, dp_t& tg)
{
    auto const &graphorder = g2->rank_to_node();
    int len2 = graphorder.size();
    std::vector<std::uint32_t> node_to_rank(len2 + 1, 0);
    for(int i = 0; i < len2; ++ i) node_to_rank[graphorder[i]->id] = len2 - i - 1;
    dp_t *II = AllocateDPVec(len2 + 1), *CC_pre = AllocateDPVec(len2 + 1);
    dp_t t;
    // init
    CC[0] = 0;
    t = O;
    for(int j = 1; j <= len2; ++ j)
    {
        CC[j] = max_score;
        for(const auto &item: graphorder[len2 - j]->outedges)
        {
            if(! item->head->code)
            {
                CC[j] = std::min(CC[j], O + E);
                // fprintf(stderr, "ID: %d, pre: %d\n", j - 1, 0);
                continue;
            }
            auto pre_id = node_to_rank[item->head->id] + 1;
            CC[j] = std::min(CC[j], CC[pre_id] + E);
        }
        if(CC[j] == max_score) CC[j] = O + E;
        DD[j] = CC[j] + O;
    }
    t = tg; // [*]
    for(int i = 1; i <= line1; ++ i)
    {
        CC_pre[0] = CC[0];
        CC[0] = t = t + E;
        II[0] = t + O;
        for(int j = 1; j <= len2; ++ j)
        {
            DD[j] = std::min(DD[j], CC[j] + O) + E;
            CC_pre[j] = CC[j];
            CC[j] = II[j] = max_score;
            const auto &this_node = graphorder[len2 - j];
            for(const auto &item: this_node->outedges)
            {
                auto pre_id = item->head->code ? node_to_rank[item->head->id] + 1 : 0;
                II[j] = std::min(std::min(II[pre_id], CC[pre_id] + O) + E, II[j]);
                CC[j] = std::min(std::min(CC[j], II[j]), CC_pre[pre_id] + MATCH_SCORE(seq1[len1 - i], g2->decoder(this_node->code), M, X));
            }
            CC[j] = std::min(DD[j], CC[j]);
        }
    }
    DD[0] = CC[0];
#if DP_DEBUG
    fprintf(stderr, "CC[len2] = "); print_vec(CC, len2);
    fprintf(stderr, "DD[len2] = "); print_vec(DD, len2);
#endif
    FreeDPVec(II);
    FreeDPVec(CC_pre);
}

void print_matrix(dp_t** m, int len1, int len2)
{
    for (int i = 0; i <= len1; ++i, fputc('\n', stderr))
        for (int j = 0; j <= len2; ++j, fputc(' ', stderr))
            fprintf(stderr, "%d", m[i][j]);
};

Alignment graph_to_sequence(const char *seq1, const Graph *g2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, int *score)
{
    auto const &graphorder = g2->rank_to_node();
    size_t len1 = strlen(seq1), len2 = graphorder.size();
    dp_t **C = AllocateDPMtx(len1 + 1, len2 + 1),
         **D = AllocateDPMtx(len1 + 1, len2 + 1),
         **I = AllocateDPMtx(len1 + 1, len2 + 1);
    if(len2 == 0) return Alignment();
    // init
    C[0][0] = 0;
    dp_t t = O, min_score = max_score;
    int32_t min_start_i, min_start_j;
    std::vector<std::uint32_t> node_to_rank(len2 + 1, 0);
    for(int i = 0; i < len2; ++ i) node_to_rank[graphorder[i]->id] = i;
    for(int j = 1; j <= len2; ++ j)
    {
        C[0][j] = max_score;
        for(const auto &item: graphorder[j - 1]->inedges)
        {
            if(! item->tail->code)
            {
                C[0][j] = std::min(C[0][j], O + E);
                // fprintf(stderr, "ID: %d, pre: %d\n", j - 1, 0);
                continue;
            }
            auto pre_id = node_to_rank[item->tail->id] + 1;
            C[0][j] = std::min(C[0][j], C[0][pre_id] + E);
            // fprintf(stderr, "ID: %d, pre: %d\n", j - 1, pre_id);
        }
        D[0][j] = C[0][j] + O;
    }
    t = O; // [*]
    for(int i = 1; i <= len1; ++ i)
    {
        C[i][0] = t = t + E;
        I[i][0] = t + O;
        for(int j = 1; j <= len2; ++ j)
        {
            // single mode:
            // I[i][j] = std::min(I[i][j - 1], C[i][j - 1] + O) + E;
            // C[i][j] = std::min(std::min(D[i][j], I[i][j]), C[i - 1][j - 1] + MATCH_SCORE(seq1[i - 1], g2->decoder(graphorder[j - 1]->code), M, X));
            D[i][j] = std::min(D[i - 1][j], C[i - 1][j] + O) + E;

            C[i][j] = D[i][j];
            I[i][j] = max_score;
            const auto &this_node = graphorder[j - 1];
            for(const auto &item: this_node->inedges)
            {
                // if pre node is start node, link to row 0
                auto pre_id = item->tail->code ? node_to_rank[item->tail->id] + 1 : 0;
                I[i][j] = std::min(std::min(I[i][pre_id], C[i][pre_id] + O) + E, I[i][j]);
                C[i][j] = std::min(C[i][j], C[i - 1][pre_id] + MATCH_SCORE(seq1[i - 1], g2->decoder(this_node->code), M, X));
            }
            C[i][j] = std::min(I[i][j], C[i][j]);
            if(i == len1) // graph is in end, update the last row
            {
                bool has_end = false;
                for(const auto &item: this_node->outedges)
                    if(! item->head->code)
                    {
                        has_end = true;
                        break;
                    }
                if(has_end && C[i][j] < min_score)
                {
                    min_score = C[i][j];
                    *score = min_score;
                    min_start_i = i;
                    min_start_j = j;
                    // fprintf(stderr, "Info: found backtrack start point (%d, %d), node_id = %d\n", i, j, this_node->id);
                }
            }
        }
    }
    int32_t i = min_start_i, j = min_start_j;
    auto nw_condition = [&i, &j] () -> bool
    {
        return (i == 0 && j == 0) ? false : true;
    };
    Alignment res;
    typedef enum {affine_matrix_M, affine_matrix_D, affine_matrix_I} matrix_type;
    matrix_type now_type = affine_matrix_M;
    while(nw_condition())
    {
        int32_t pred_j;
        const auto &this_node = graphorder[j - 1];
        bool ismatch = false;
        switch(now_type)
        {
            case affine_matrix_M:
                // first determine Match
                if(this_node->inedges.empty())
                {
                    // must not be match because null node align with sequence
                    ;
                }
                else
                {
                    for(const auto &item: this_node->inedges)
                    {
                        pred_j = item->tail->code ? node_to_rank[item->tail->id] + 1 : 0;
                        if(C[i - 1][pred_j] + M == C[i][j] && seq1[i - 1] == g2->decoder(this_node->code))
                        {
                            ismatch = true;
                            break;
                        }
                    }
                }
                if(ismatch)
                {
                    res.emplace_back(this_node->id, i - 1); // first graph, then sequence
#if DP_DEBUG_CONDITION
                    fprintf(stderr, "Condition %s: %d %d, will jump to (%d, %d)\n", "Match", i, j, i - 1, pred_j);
#endif
                    i = i - 1;
                    j = pred_j;
                }
                else
                {
                    if(C[i][j] == D[i][j])
                    {
                        now_type = affine_matrix_D;
                    }
                    else if(C[i][j] == I[i][j])
                    {
                        now_type = affine_matrix_I;
                    }
                    else // must be mismatch
                    {
                        bool mismatch = false;
                        if(this_node->inedges.empty())
                        {
                            if(C[i][j] == C[i - 1][0] + X) // must be different
                            {
                                mismatch = true;
                                pred_j = 0;
                            }
                        }
                        else
                        {
                            for(const auto &item: this_node->inedges)
                            {
                                pred_j = item->tail->code ? node_to_rank[item->tail->id] + 1 : 0;
                                if(C[i - 1][pred_j] + X == C[i][j] && seq1[i - 1] != g2->decoder(this_node->code))
                                {
                                    mismatch = true;
                                    break;
                                }
                            }
                        }
                        assert(mismatch);
                        res.emplace_back(this_node->id, i - 1); // first graph, then sequence
#if DP_DEBUG_CONDITION
                        fprintf(stderr, "Condition %s: %d %d, will jump to (%d, %d)\n", "Mismatch", i, j, i - 1, pred_j);
#endif
                        i = i - 1;
                        j = pred_j;
                    }
                }
                break;
            case affine_matrix_D:
                res.emplace_back(-1, i - 1);
#if DP_DEBUG_CONDITION
                fprintf(stderr, "Condition Delete 2: %d %d, will jump to (%d, %d)\n", i + 1, j, i, j);
#endif
                if(D[i - 1][j] + E != D[i][j]) now_type = affine_matrix_M;
                i--;
                break;
            case affine_matrix_I:
                res.emplace_back(this_node->id, -1);
#if DP_DEBUG_CONDITION
                fprintf(stderr, "Condition %s: %d %d, will jump to (%d, %d)\n", "Insert 2", i, j, i, pred_j);
#endif
                bool continue_INSERT = false;
                // try to find predecessor                            
                if(this_node->inedges.empty())
                {
                    pred_j = 0;
                    if(I[i][pred_j] + E == I[i][j]) continue_INSERT = true;
                }
                else
                {
                    pred_j = -1;
                    for(const auto &item: this_node->inedges)
                    {
                        pred_j = item->tail->code ? node_to_rank[item->tail->id] + 1 : 0;
                        if(I[i][pred_j] + E == I[i][j])
                        {
                            continue_INSERT = true;
                            break;
                        }
                    }
                }
                if(!continue_INSERT) now_type = affine_matrix_M;
                j = pred_j;
                break;
        }
        
        if(i == 0) // align before nodes to gap
        {
            while(j != 0)
            {
                const auto &this_node = graphorder[j - 1];
                std::uint32_t pred_j;
                for(const auto &it: this_node->inedges)
                {
                    if(it->tail->code)
                    {
                        pred_j = node_to_rank[it->tail->id] + 1;
                        if(C[i][j] == C[i][pred_j] + E)
                        {
                            j = pred_j;
                            break;
                        }
                    }
                    else
                    {
                        pred_j = 0;
                        if(C[i][j] == C[i][0] + O + E)
                        {
                            j = pred_j;
                            break;
                        }
                    }
                }
                res.emplace_back(this_node->id, -1);
            }
        }
        if(j == 0) // align
        {
            while(i != 0) res.emplace_back(-1, --i);
        }
    }
    std::reverse(res.begin(), res.end());
    *score = min_score;
#if DP_PRINT_MARTIX
    fprintf(stderr, "I\n");
    print_matrix(I, len1, len2);
    fprintf(stderr, "D\n");
    print_matrix(D, len1, len2);
    fprintf(stderr, "C\n");
    print_matrix(C, len1, len2);
#endif
    FreeDPMtx(C);
    FreeDPMtx(D);
    FreeDPMtx(I);
    return res;
}
