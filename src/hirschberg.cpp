#include "hirschberg.hpp"
#include "../include/threadpool/include/scheduler.hpp"
#include "linearpoa.h"
#define HIRSCHBERG_DEBUG 0
#define HIRSCHBERG_SET_DEBUG 0
#define HIRSCHBERG_CONQUER_DEBUG 0
#define HIRSCHBERG_DEBUG_PRINT 0
#define HIRSCHBERG_ASSERT 0

#define MATCH_SCORE(a, b, M, X) ((a) == (b) ? (M) : (X))
#define HIRSCHBERG_GAP(N, O, E) ((N) <= 0 ? 0 : (O) + (E) * (N))

constexpr std::int32_t kInfinity = std::numeric_limits<std::int32_t>::max() - 1024;

inline dp_t hirschberg_gap(dp_t &N, dp_t &O, dp_t &E)
{
    return N <= 0 ? 0 : O + E * N;
}

using namespace staccato;

class HirschbergTask
{
public:
    HirschbergTask(char *seq1, int b1, int l1, Graph *G, dp_t tb, dp_t te, dp_t &O, dp_t &E, dp_t &M, dp_t &X, bool multidp, int dp_threads,
                   thread_pool_light *dp_pool_fwd, thread_pool_light *dp_pool_rev, std::mutex *mtx_fwd, std::mutex *mtx_rev,
                   ConcurrentSet_hirschberg_status_t *&status_set, bool first) :
                   seq1(seq1), b1(b1), l1(l1), G(G), tb(tb), te(te), O(O), E(E), M(M), X(X), multidp(multidp), dp_threads(dp_threads),
                   dp_pool_fwd(dp_pool_fwd), dp_pool_rev(dp_pool_rev), mtx_fwd(mtx_fwd), mtx_rev(mtx_rev),
                   status_set(status_set), first(first) {}
    char *seq1;
    Graph* G;
    int b1, l1, dp_threads;
    dp_t tb, te, O, E, M, X;
    thread_pool_light *dp_pool_fwd, *dp_pool_rev;
    std::mutex *mtx_fwd, *mtx_rev;
    ConcurrentSet_hirschberg_status_t *status_set;
    bool multidp, first;
};

class linear_dp_Task
{
public:
    linear_dp_Task(dp_t *CC, dp_t *DD, char *seq1, int line1, int b1, int l1, Graph *G, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t tg,
                   bool rev, bool multidp, int &dp_threads, thread_pool_light *pool, std::mutex *mtx) :
                   CC(CC), DD(DD), seq1(seq1), line1(line1), b1(b1), l1(l1), G(G), M(M), X(X), O(O), E(E), tg(tg),
                   rev(rev), multidp(multidp), dp_threads(dp_threads), pool(pool), mtx(mtx) {}
    dp_t *CC, *DD;
    Graph* G;
    char *seq1;
    int line1, b1, l1, dp_threads;
    dp_t M, X, O, E, tg;
    bool rev, multidp;
    std::mutex *mtx;
    thread_pool_light *pool;
};

union task_info_t
{
    linear_dp_Task dpTask;
    HirschbergTask HirTask;
    task_info_t(linear_dp_Task dpTask):  dpTask(dpTask)   {}
    task_info_t(HirschbergTask HirTask): HirTask(HirTask) {}
};

class TaskRunner : public staccato::task<TaskRunner>
{
public:
    TaskRunner(linear_dp_Task dpTask) : Info(dpTask)  { type_ = 1; }
    TaskRunner(HirschbergTask HirTask): Info(HirTask) { type_ = 2; }
    void execute();
    void execute_linear_dp();
    void execute_hirschberg();
private:
    task_info_t Info;
    int type_;
};

void hirschberg_multi_init(ConcurrentSet_hirschberg_status_t *&S, int len, thread_pool_light *&pool1, thread_pool_light *&pool2, int dp_threads, std::mutex *&mtx1, std::mutex *&mtx2, bool multidp)
{
    S = new ConcurrentSet_hirschberg_status_t(len <= 100007 ? len : 100007); // TODO: less len to find error on sorting
    if(multidp)
    {
        pool1 = new thread_pool_light(dp_threads);
        pool2 = new thread_pool_light(dp_threads);
        mtx1 = new std::mutex;
        mtx2 = new std::mutex;
    }
    else pool1 = pool2 = nullptr, mtx1 = mtx2 = nullptr;
}

void hirschberg_multi_free(ConcurrentSet_hirschberg_status_t *&S, thread_pool_light *pool1, thread_pool_light *pool2, std::mutex *mtx1, std::mutex *mtx2)
{
    delete S;
    delete pool1;
    delete pool2;
    delete mtx1;
    delete mtx2;
}

// status: first is seq id, second is graph id
inline void del_status(int b1, int l1, ConcurrentSet_hirschberg_status_t &S)
{
#if HIRSCHBERG_SET_DEBUG
    fprintf(stderr, "%s seq1 = %\n", "del  ", b1, l1);
#endif
    for(int i = 0; i < l1; ++ i) S.insert(b1 + i, -1);
}

inline void ins_status(int b2, ConcurrentSet_hirschberg_status_t &S)
{
#if HIRSCHBERG_SET_DEBUG
    fprintf(stderr, "%s g2 = %d\n", "ins  ", b2);
#endif
    S.insert(-1, b2);
}

inline void match_status(int b1, int b2, ConcurrentSet_hirschberg_status_t &S)
{
#if HIRSCHBERG_SET_DEBUG
    fprintf(stderr, "%s seq1 = %d g2 = %d\n", "match", b1, b2);
#endif
    S.insert(b1, b2);
}

void TaskRunner::execute_linear_dp()
{
#define F(X) this->Info.dpTask.X
    dp_t *CC = F(CC), *DD = F(DD);
    char *seq1 = F(seq1);
    Graph *G = F(G);
    int line1 = F(line1), b1 = F(b1), l1 = F(l1), len2 = G->rank_to_node().size();
    dp_t M = F(M), X = F(X), O = F(O), E = F(E), tg = F(tg);
    bool rev = F(rev), multidp = F(multidp);
    int dp_threads = F(dp_threads);
    thread_pool_light *pool = F(pool);
    std::mutex *mtx = F(mtx);
#undef F
    if(rev)
    {
        if (multidp && dp_threads > 1 && (line1 >= 1000 && len2 >= 1000))
            rev_linear_poa_dp_hirschberg(seq1 + b1, G, CC, DD, line1, l1, M, X, O, E, tg); // TODO: may need to provide multithread version
        else rev_linear_poa_dp_hirschberg(seq1 + b1, G, CC, DD, line1, l1, M, X, O, E, tg);
    }
    else
    {
        if(multidp && dp_threads > 1 && (line1 >= 1000 && len2 >= 1000))
            linear_poa_dp_hirschberg(seq1 + b1, G, CC, DD, line1, M, X, O, E, tg); // TODO: may need to provide multithread version
        else linear_poa_dp_hirschberg(seq1 + b1, G, CC, DD, line1, M, X, O, E, tg);
    }
}

void TaskRunner::execute_hirschberg()
{
#define F(X) this->Info.HirTask.X
    auto b1 = F(b1), l1 = F(l1);
    auto G = F(G);
    auto tb = F(tb), te = F(te), O = F(O), E = F(E), M = F(M), X = F(X);
    auto status_set = F(status_set);
    auto seq1 = F(seq1);
    auto multidp = F(multidp), first = F(first);
    auto dp_threads = F(dp_threads);
    auto *dp_pool_fwd = F(dp_pool_fwd), *dp_pool_rev = F(dp_pool_rev);
    auto *mtx_fwd = F(mtx_fwd), *mtx_rev = F(mtx_rev);
#undef F
    /* if Graph B is empty.... */
    int l2 = G->rank_to_node().size();
    if(l2 <= 0)
    {
        /* if sequence A is not empty.... */
        if(l1 > 0)
        {
            int id = 0;
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "del l1: b1 = %zu, l1 = %zu, b2 = %zu (final)\n", b1, l1, G->nodes()[0].get()->raw_id);
#endif
            // get the pre node id for saving node id
            del_status(b1, l1, *status_set);
        }
        return;
    }
    else if(l1 <= 1)
    {
        /* if sequence A is empty.... */
        if(l1 <= 0)
        {
            /* insert residues from all nodes in Graph B */
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "add l2: b2 = %zu, l2 = %zu (final)\n", G->rank_to_node()[0]->raw_id, l2);
#endif
            // maybe do nothing?
            assert(l2 > 1);
            ins_status(G->rank_to_node()[0]->raw_id, *status_set);
            return;
        }
        /* if sequence A has just one residue.... */
        // maybe only need to simply sequence-to-graph alignment??
        // with only one sequence align to graph
        const auto &graphorder = G->rank_to_node(), &revgraphorder = G->rank_to_node_rev();
        l2 = graphorder.size() + 1; // calculate length in calling graphorder
        std::vector<std::uint32_t> node_to_rank(graphorder.size() + 1, 0), rev_node_to_rank(graphorder.size() + 1, 0);
        for(int i = 0; i < graphorder.size(); ++ i) node_to_rank[graphorder[i]->id] = i;
        for(int i = 0; i < revgraphorder.size(); ++ i) rev_node_to_rank[revgraphorder[i]->id] = i;
        int midj = 0;
        dp_t min_score = kInfinity, cur;
        for(auto &item: graphorder)
        {
            l2 = l2 > node_to_rank[item->id] + rev_node_to_rank[item->id] ? node_to_rank[item->id] + rev_node_to_rank[item->id] : l2;
            cur = HIRSCHBERG_GAP(node_to_rank[item->id], tb, E) + MATCH_SCORE(seq1[b1 + l1 - 1], G->decoder(item->code), M, X) + HIRSCHBERG_GAP(rev_node_to_rank[item->id], te, E);
            if(cur < min_score)
            {
                min_score = cur;
                midj = item->raw_id;
            }
        }
        // calculate the border condition
        cur = (tb + E) + HIRSCHBERG_GAP(l2, te, E);
        if(min_score > cur) min_score = cur, midj = -1;
        cur = (te + E) + HIRSCHBERG_GAP(l2, tb, E);
        if(min_score > cur) min_score = cur, midj = -1;

        if(midj == -1) // graph has not aligned with sequence
        {
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "del 1: b1 = %zu\n", b1);
#endif
            del_status(b1, 1, *status_set);
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "add l2: b2 = %zu, l2 = %zu\n", G->rank_to_node()[0]->raw_id, l2);
#endif
            for(auto &item: graphorder) ins_status(item->raw_id, *status_set);
        }
        else // aligned a node
        {
            match_status(b1 + l1 - 1, midj, *status_set);
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "match 1: seq1 = %d, g2 = %d (in)\n", l1 + b1, midj);
#endif
            for(auto &item: graphorder)
            {
                if(item->raw_id == midj) continue;
#if HIRSCHBERG_DEBUG
                fprintf(stderr, "add 1: g2 = %d\n", item->raw_id);
#endif
                ins_status(item->raw_id, *status_set);
            }
        }
        // graph clean by pre calling
        return;
    }
    else
    {
        int imid = l1 >> 1;
        std::uint32_t jmid, jmid_nxt;
        dp_t *CC = AllocateDPVec(l2 + 1), *DD = AllocateDPVec(l2 + 1), *RR = AllocateDPVec(l2 + 1), *SS = AllocateDPVec(l2 + 1);
        linear_dp_Task dp_pre(CC, DD, seq1, imid,      b1, l1, G, M, X, O, E, tb, false, multidp, dp_threads, dp_pool_fwd, mtx_fwd),
                       dp_rev(RR, SS, seq1, l1 - imid, b1, l1, G, M, X, O, E, te, true,  multidp, dp_threads, dp_pool_rev, mtx_rev);
        spawn(new(child()) TaskRunner(dp_pre));
        spawn(new(child()) TaskRunner(dp_rev));
        wait();

        bool type, has_multiple_starts = false; dp_t min_val, cur_val;

        const auto &graphorder = G->rank_to_node();
        std::vector<std::uint32_t> node_to_rank(l2 + 1, 0);
        for(int i = 0; i < l2; ++ i) node_to_rank[graphorder[i]->id] = i;

        /* find midj, such that CC[j] + RR[len2 - next(j)] or DD[j] + SS[len2 - next(j)] - gap is the max */
        min_val = std::numeric_limits<dp_t>::max() - 1024; jmid = 0; jmid_nxt = 0; type = false;
        // determine the head nodes
        for(const auto &item: G->nodes()[0].get()->outedges)
        {
            // type 1
            const auto nxt = item->head->code ? node_to_rank[item->head->id] : 0;
            cur_val = CC[0] + RR[l2 - nxt];
            if(cur_val <= min_val)
            {
                if(cur_val < min_val || (CC[0] != DD[0] && RR[l2 - nxt] == SS[l2 - nxt]))
                {
                    min_val = cur_val;
                    jmid = std::numeric_limits<std::uint32_t>::max(); // contain all first nodes
                    jmid_nxt = graphorder[nxt]->id;
                    type = false;
                }
            }
            // type 2
            cur_val = DD[0] + SS[l2 - nxt] - O;
            if(cur_val < min_val)
            {
                min_val = cur_val;
                jmid = std::numeric_limits<std::uint32_t>::max();
                jmid_nxt = graphorder[nxt]->id;
                type = true;
            }
        }

        for(int j = 1; j <= l2; ++ j)
        {
            const auto &this_node = graphorder[j - 1];
            for(const auto &item: graphorder[j - 1]->outedges)
            {
                unsigned int nxt;
                if(! item->head->code) // end node
                    nxt = l2;
                else
                    nxt = node_to_rank[item->head->id];
                // type 1
                cur_val = CC[j] + RR[l2 - nxt];
                if(cur_val <= min_val)
                {
                    if((cur_val < min_val) || (CC[j] != DD[j] && RR[l2 - nxt] == SS[l2 - nxt])) // consider type 1 more than type 2
                    {
                        min_val = cur_val;
                        jmid = this_node->id;
                        jmid_nxt = item->head->code ? item->head->id : std::numeric_limits<std::uint32_t>::max();
                        type = false;
                    }
                }
                // type 2
                cur_val = DD[j] + SS[l2 - nxt] - O;
                if(cur_val < min_val)
                {
                    min_val = cur_val;
                    jmid = this_node->id;
                    jmid_nxt = item->head->code ? item->head->id : std::numeric_limits<std::uint32_t>::max();
                    type = true;
                }
            }
        }
#if HIRSCHBERG_CONQUER_DEBUG
        for (size_t i = 0; i < l1; ++i) fputc(seq1[i + b1], stderr); fputc('\n', stderr);
        fprintf(stderr, "O = %d, E = %d, M = %d, X = %d\n", O, E, M, X);
        fprintf(stderr, "tb = "); printDPdata(tb, stderr); fputc('\n', stderr);
        fprintf(stderr, "te = "); printDPdata(te, stderr); fputc('\n', stderr);
        fprintf(stderr, "CC    DD    RR    SS\n");
        for(size_t i = 0; i <= l2; ++ i)
            fprintf(stderr, "%-5d %-5d %-5d %-5d\n", CC[i], DD[i], RR[i], SS[i]);
        fprintf(stderr, "min_val = %d, split point = %d, %d, type = %d\n", min_val, jmid, jmid_nxt, type ? 2 : 1);
#endif
        FreeDPVec(CC); FreeDPVec(DD); FreeDPVec(RR); FreeDPVec(SS);
        /* Conquer recursively around midpoint */
        Graph G_left, G_right;
        std::vector <const Graph::Node*> tmp_graph;
        if(! type)
        {
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "This function will choose type 1 with (%zu, %zu) => (%d, %zu) and (%zu, %zu) => (%zu, %zu); max = %d, (i*, j*) = (%d, %d)\n",
                            b1, imid, G->rank_to_node()[0]->raw_id, jmid, b1 + imid, l1 - imid, jmid == std::numeric_limits<std::uint32_t>::max() ? G->nodes()[0].get()->raw_id : G->rank_to_node()[jmid]->raw_id, l2 - jmid, min_val, imid, jmid);
#endif
            G_left = G->Subgraph_left(jmid);
            G_right = G->Subgraph_right(jmid_nxt, jmid);
            if(!first) G->Clear();
            HirschbergTask left_task (seq1, b1,        imid,      &G_left, tb,  O, O, E, M, X, multidp, dp_threads, dp_pool_fwd, dp_pool_rev, mtx_fwd, mtx_rev, status_set, false),
                           right_task(seq1, b1 + imid, l1 - imid, &G_right, O, te, O, E, M, X, multidp, dp_threads, dp_pool_fwd, dp_pool_rev, mtx_fwd, mtx_rev, status_set, false);
            spawn(new(child()) TaskRunner(left_task));
            spawn(new(child()) TaskRunner(right_task));
            wait();
            // Implicit Destruction for G_left and G_right (garbage collect)
        }
        else
        {
            del_status(imid + b1 - 1, 2, *status_set);
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "This function will choose type 2 with (%zu, %zu) => (%zu, %zu), with 2 deletions, and (%zu, %zu) => (%zu, %zu); max = %d, (i*, j*) = (%zu, %zu)\n",
                            b1, imid - 1, G->rank_to_node()[0]->raw_id, jmid, b1 + imid + 1, l1 - imid - 1, jmid == std::numeric_limits<std::uint32_t>::max() ? G->nodes()[0].get()->raw_id : G->rank_to_node()[jmid]->raw_id, l2 - jmid, min_val, imid, jmid);
            fprintf(stderr, "del 2: midi = %zu\n", imid + b1);
#endif
#undef HIRSCHBERG_DEBUG
            G_left = G->Subgraph_left(jmid);
            G_right = G->Subgraph_right(jmid_nxt, jmid);
            if(!first) G->Clear();
            HirschbergTask left_task (seq1, b1,            imid - 1,      &G_left, tb,  0, O, E, M, X, multidp, dp_threads, dp_pool_fwd, dp_pool_rev, mtx_fwd, mtx_rev, status_set, false),
                           right_task(seq1, b1 + imid + 1, l1 - imid - 1, &G_right, 0, te, O, E, M, X, multidp, dp_threads, dp_pool_fwd, dp_pool_rev, mtx_fwd, mtx_rev, status_set, false);
            spawn(new(child()) TaskRunner(left_task));
            spawn(new(child()) TaskRunner(right_task));
            wait();
            // Implicit Destruction for G_left and G_right (garbage collect)
            }
    }
}

void TaskRunner::execute()
{
    switch (this->type_)
    {
    case 1:
        this->execute_linear_dp();
        break;
    case 2:
        this->execute_hirschberg();
        break;
    default:
        break;
    }
}

void hirschberg_multi_start(char *seq1, Graph* G, int len1, dp_t O, dp_t E, dp_t M, dp_t X,
                            ConcurrentSet_hirschberg_status_t *&status_set, bool multidp, int dp_threads,
                            thread_pool_light *pool1, thread_pool_light *pool2, std::mutex *mtx1, std::mutex *mtx2, int threads)
{
    scheduler<TaskRunner> ____(4, threads);
    ____.spawn(new(____.root()) TaskRunner(HirschbergTask(seq1, 0, len1, G, O, O, O, E, M, X, multidp, dp_threads, pool1, pool2, mtx1, mtx2, status_set, true)));
    ____.wait();
}

// TODO: support type Alignment, and return Alignment
void write_final_result(int len1, const char *seq1, Graph *g2, ConcurrentSet_hirschberg_status_t *S, Alignment *res, int &score)
{
    // if not priority_queue, sort the data
    std::vector <hirschberg_status_t*> v;
    S->getFinalSnapshotAndFree(v);
    // the status has only two parts
    const auto &graphorder = g2->rank_to_node();
    std::vector<std::uint32_t> node_to_rank(graphorder.size() + 1, 0);
    for(int i = 0; i < graphorder.size(); ++ i) node_to_rank[graphorder[i]->id] = i;

    Alignment aligned_pair, unaligned_pair;
    aligned_pair.reserve(graphorder.size());
    unaligned_pair.reserve(graphorder.size());
    for(const auto& p : v)
    {
        if(p->seq1 == -1 || p->g2 == -1) unaligned_pair.emplace_back(p->g2, p->seq1);
        else aligned_pair.emplace_back(p->g2, p->seq1);
    }
    // align_pair_t: first g2, second seq1
    std::sort(aligned_pair.begin(), aligned_pair.end(), [&](const align_pair_t& a, const align_pair_t& b)
        {
            if (a.first != b.first) return node_to_rank[a.first] < node_to_rank[b.first];
            fprintf(stderr, "Warning: Graph node must be different in aligned_pair!!\n");
            return a.second < b.second; // may be useless because the order must fit the first condition
        });
    std::sort(unaligned_pair.begin(), unaligned_pair.end(), [&](const auto &a, const auto &b) 
             {
                 auto get_key = [&](const auto& p) -> int32_t {
                        if (p.first != -1) return node_to_rank[p.first];
                        return p.second;
                    };
                 auto a_key = get_key(a);
                 auto b_key = get_key(b);
                 // 主排序：键值升序
                 if (a_key != b_key) return a_key < b_key;
                 // 次排序：优先非-1在second位置
                 return (a.second != -1) > (b.second != -1);
             });
    for (const auto& p : unaligned_pair)
    {
        auto insert_key = (p.first != -1) ? p.first : p.second;
        bool use_first  = (p.first != -1);
        auto it = std::lower_bound(aligned_pair.begin(), aligned_pair.end(), insert_key,
                                   [use_first, &node_to_rank](const auto &elem, auto key)
                                   {
                                       return use_first ? (node_to_rank[elem.first] < node_to_rank[key]) : elem.second < key;
                                   });
        aligned_pair.insert(it, p);
    }
    std::copy(aligned_pair.begin(), aligned_pair.end(), std::inserter(*res, res->begin()));
}

Alignment AlignmentEngine::linearspace_Align(
    const char* seq1, std::uint32_t len1, Graph &G,
    std::int32_t *score,
    int threads, int dp_threads, bool multidp)
{
    const int len2 = G.nodes().size();
    if (len1 > std::numeric_limits<int32_t>::max())
    {
        throw std::invalid_argument("[linearPOA::linearspace_Align] Error: too large sequence!");
    }
    if (! len1 || ! len2) return Alignment();

    if (WorstCaseAlignmentScore(len1, len2) * -1 > kInfinity)
    {
        throw std::invalid_argument("[linearPOA::linearspace_Align] error: possiblly overflow, we throw this exception for not aligning sequences.");
    }
    ConcurrentSet_hirschberg_status_t* s;
    thread_pool_light *pool1, *pool2;
    std::mutex *mtx1, *mtx2;
    hirschberg_multi_init(s, len1 + len2, pool1, pool2, dp_threads, mtx1, mtx2, multidp ? true : false);
    hirschberg_multi_start((char*)seq1, &G, len1, g_, e_, m_, n_, s, multidp, dp_threads, pool1, pool2, mtx1, mtx2, threads);
    Alignment res;
    write_final_result(len1, seq1, &G, s, &res, *score);
    hirschberg_multi_free(s, pool1, pool2, mtx1, mtx2);
    return res;
}
