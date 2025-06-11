#include "poa.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <fstream>
#include <stack>
#include <stdexcept>
#include <unordered_set>
#include <iostream>
#include <cassert>
#include <climits>

constexpr std::int32_t kNegativeInfinity = std::numeric_limits<std::int32_t>::min() + 1024;

// AlignmentEngine start

AlignmentEngine::AlignmentEngine(std::int8_t m, std::int8_t n, std::int8_t g, std::int8_t e) : m_(m), n_(n), g_(g), e_(e)
{
    this->pimpl_ = std::make_unique<Implementation>();
}

void AlignmentEngine::Prealloc(std::uint32_t max_sequence_len, std::uint8_t alphabet_size)
{
    if (max_sequence_len > std::numeric_limits<int32_t>::max())
    {
        throw std::invalid_argument("[spoa::SisdAlignmentEngine::Prealloc] error: too large sequence!");
    }
    try
    {
        Realloc(static_cast<std::uint64_t>(max_sequence_len) + 1,
            static_cast<std::uint64_t>(max_sequence_len) * alphabet_size + alphabet_size,  // NOLINT
            alphabet_size);
    }
    catch (std::bad_alloc& ba)
    {
        throw std::invalid_argument("[spoa::SisdAlignmentEngine::Prealloc] error: insufficient memory!");
    }
}

void AlignmentEngine::Realloc(std::uint64_t matrix_width, std::uint64_t matrix_height, std::uint8_t num_codes)
{
    if (pimpl_->node_id_to_rank.size() < matrix_height - 1)
    {
        pimpl_->node_id_to_rank.resize(matrix_height - 1, 0);
    }
    if (pimpl_->sequence_profile.size() < num_codes * matrix_width)
    {
        pimpl_->sequence_profile.resize(num_codes * matrix_width, 0);
    }
    // Allocate Affine matrix
    if (pimpl_->M.size() < 3 * matrix_height * matrix_width)
    {
        pimpl_->M.resize(3 * matrix_width * matrix_height, 0);
        pimpl_->H = pimpl_->M.data();
        pimpl_->F = pimpl_->H + matrix_width * matrix_height;
        pimpl_->E = pimpl_->F + matrix_width * matrix_height;
    }
}

void AlignmentEngine::Initialize(const char* sequence, std::uint32_t sequence_len, const Graph& graph) noexcept
{
    std::uint32_t matrix_width = sequence_len + 1;
    std::uint32_t matrix_height = graph.nodes().size() + 1;

    for (std::uint32_t i = 0; i < graph.num_codes(); ++i)
    {
        char c = graph.decoder(i);
        pimpl_->sequence_profile[i * matrix_width] = 0;
        for (std::uint32_t j = 0; j < sequence_len; ++j)
            pimpl_->sequence_profile[i * matrix_width + (j + 1)] = (c == sequence[j] ? m_ : n_);
    }

    const auto& rank_to_node = graph.rank_to_node();
    for (std::uint32_t i = 0; i < rank_to_node.size(); ++i) {
        pimpl_->node_id_to_rank[rank_to_node[i]->id] = i;
    }

    // initialize secondary matrices
    pimpl_->F[0] = 0;
    pimpl_->E[0] = 0;
    for (std::uint32_t j = 1; j < matrix_width; ++j)
    {
        pimpl_->F[j] = kNegativeInfinity;
        pimpl_->E[j] = g_ + (j - 1) * e_;
    }
    for (std::uint32_t i = 1; i < matrix_height; ++i)
    {
        const auto& edges = rank_to_node[i - 1]->inedges;
        std::int32_t penalty = edges.empty() ? g_ - e_ : kNegativeInfinity;
        for (const auto& it : edges)
        {
            std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
            penalty = std::max(penalty, pimpl_->F[pred_i * matrix_width]);
        }
        pimpl_->F[i * matrix_width] = penalty + e_;
        pimpl_->E[i * matrix_width] = kNegativeInfinity;
    }
    pimpl_->H[0] = 0;

    // initialize primary matrix
    for (std::uint32_t j = 1; j < matrix_width; ++j)
        pimpl_->H[j] = pimpl_->E[j];
    for (std::uint32_t i = 1; i < matrix_height; ++i)
        pimpl_->H[i * matrix_width] = pimpl_->F[i * matrix_width];
}

Alignment AlignmentEngine::Align(const char* sequence, std::uint32_t sequence_len, const Graph& graph, std::int32_t* score)
{
    if (sequence_len > std::numeric_limits<int32_t>::max())
    {
        throw std::invalid_argument("[spoa::SisdAlignmentEngine::Align] error: too large sequence!");
    }

    if (graph.nodes().empty() || sequence_len == 0) return Alignment();

    if (WorstCaseAlignmentScore(sequence_len, graph.nodes().size()) < kNegativeInfinity)
    {
        throw std::invalid_argument("[spoa::SisdAlignmentEngine::Align] error: possible overflow!");
    }

    try
    {
        Realloc(sequence_len + 1, graph.nodes().size() + 1, graph.num_codes());
    }
    catch (std::bad_alloc& ba)
    {
        throw std::invalid_argument("[spoa::SisdAlignmentEngine::Align] error: insufficient memory!");
    }
    Initialize(sequence, sequence_len, graph);
    std::uint64_t matrix_width = sequence_len + 1;
    const auto& rank_to_node = graph.rank_to_node();

    std::int32_t max_score = kNegativeInfinity;
    std::uint32_t max_i = 0;
    std::uint32_t max_j = 0;
    auto update_max_score = [&max_score, &max_i, &max_j]
                            (std::int32_t* H_row, std::uint32_t i, std::uint32_t j) -> void
    {
        if (max_score < H_row[j])
        {
            max_score = H_row[j];
            max_i = i;
            max_j = j;
        }
        return;
    };

    // alignment
    for (const auto& it : rank_to_node)
    {
        const auto& char_profile = &(pimpl_->sequence_profile[it->code * matrix_width]);

        std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
        std::uint32_t pred_i = it->inedges.empty() ? 0 :
            pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

        std::int32_t* H_row = &(pimpl_->H[i * matrix_width]);
        std::int32_t* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

        std::int32_t* F_row = &(pimpl_->F[i * matrix_width]);
        std::int32_t* F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

        // update F and H
        for (std::uint64_t j = 1; j < matrix_width; ++j)
        {
            F_row[j] = std::max(H_pred_row[j] + g_, F_pred_row[j] + e_);
            H_row[j] = H_pred_row[j - 1] + char_profile[j];
        }
        // check other predeccessors
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p)
        {
            pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

            H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
            F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

            for (std::uint64_t j = 1; j < matrix_width; ++j)
            {
                F_row[j] = std::max(F_row[j], std::max(H_pred_row[j] + g_, F_pred_row[j] + e_));
                H_row[j] = std::max(H_row[j], H_pred_row[j - 1] + char_profile[j]);
            }
        }

        // update E and H
        std::int32_t* E_row = &(pimpl_->E[i * matrix_width]);
        for (std::uint64_t j = 1; j < matrix_width; ++j)
        {
            E_row[j] = std::max(H_row[j - 1] + g_, E_row[j - 1] + e_);
            H_row[j] = std::max(H_row[j], std::max(F_row[j], E_row[j]));

            if (it->outedges.empty() && j == matrix_width - 1)
                update_max_score(H_row, i, j);
        }
    }

    if (max_i == 0 && max_j == 0) return Alignment();

    if (score) *score = max_score;

    // backtrack
    Alignment alignment;
    std::uint32_t i = max_i;
    std::uint32_t j = max_j;

    auto nw_condition = [&i, &j] () -> bool
    {
        return (i == 0 && j == 0) ? false : true;
    };

    std::uint32_t prev_i = 0;
    std::uint32_t prev_j = 0;

    while (nw_condition())
    {
        auto H_ij = pimpl_->H[i * matrix_width + j];
        bool predecessor_found = false, extend_left = false, extend_up = false;

        if (i != 0 && j != 0)
        {
            const auto& it = rank_to_node[i - 1];
            std::int32_t match_cost = pimpl_->sequence_profile[it->code * matrix_width + j];

            std::uint32_t pred_i = it->inedges.empty() ? 0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

            if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost)
            {
                prev_i = pred_i;
                prev_j = j - 1;
                predecessor_found = true;
            }
            else
            {
                for (std::uint32_t p = 1; p < it->inedges.size(); ++p)
                {
                    pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

                    if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost)
                    {
                        prev_i = pred_i;
                        prev_j = j - 1;
                        predecessor_found = true;
                        break;
                    }
                }
            }
        }

        if (!predecessor_found && i != 0)
        {
            const auto& it = rank_to_node[i - 1];

            std::uint32_t pred_i = it->inedges.empty() ? 0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

            if ((extend_up = H_ij == pimpl_->F[pred_i * matrix_width + j] + e_) ||
                             H_ij == pimpl_->H[pred_i * matrix_width + j] + g_)
            {
                prev_i = pred_i;
                prev_j = j;
                predecessor_found = true;
            }
            else
            {
                for (std::uint32_t p = 1; p < it->inedges.size(); ++p)
                {
                    pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

                    if ((extend_up = H_ij == pimpl_->F[pred_i * matrix_width + j] + e_) ||
                                     H_ij == pimpl_->H[pred_i * matrix_width + j] + g_)
                    {
                        prev_i = pred_i;
                        prev_j = j;
                        predecessor_found = true;
                        break;
                    }
                }
            }
        }

        if (!predecessor_found && j != 0)
        {
            if ((extend_left = H_ij == pimpl_->E[i * matrix_width + j - 1] + e_) ||
                               H_ij == pimpl_->H[i * matrix_width + j - 1] + g_)
            {
                prev_i = i;
                prev_j = j - 1;
                predecessor_found = true;
            }
        }

        alignment.emplace_back(i == prev_i ? -1 : rank_to_node[i - 1]->id, j == prev_j ? -1 : j - 1);

        i = prev_i;
        j = prev_j;

        if (extend_left)
        {
            while (true)
            {
                alignment.emplace_back(-1, j - 1);
                --j;
                if (pimpl_->E[i * matrix_width + j] + e_ !=
                    pimpl_->E[i * matrix_width + j + 1]) break;
            }
        }
        else if (extend_up)
        {
            while (true)
            {
                bool stop = false;
                prev_i = 0;
                for (const auto& it : rank_to_node[i - 1]->inedges)
                {
                    std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;

                    if ((stop = pimpl_->F[i * matrix_width + j] == pimpl_->H[pred_i * matrix_width + j] + g_) ||
                                pimpl_->F[i * matrix_width + j] == pimpl_->F[pred_i * matrix_width + j] + e_)
                    {
                        prev_i = pred_i;
                        break;
                    }
                }

                alignment.emplace_back(rank_to_node[i - 1]->id, -1);
                i = prev_i;
                if (stop || i == 0) break;
            }
        }
    }

    std::reverse(alignment.begin(), alignment.end());
    return alignment;
}

std::int64_t AlignmentEngine::WorstCaseAlignmentScore(std::int64_t i, std::int64_t j) const
{
    auto gap_score = [&] (std::int64_t len) -> std::int64_t
    {
        return len == 0 ? 0 : g_ + (len - 1) * e_;
    };
    return std::min(-1 * (m_ * std::min(i, j) + gap_score(std::abs(i - j))), gap_score(i) + gap_score(j));
}

std::unique_ptr<AlignmentEngine> AlignmentEngine::Create(std::int8_t m, std::int8_t n, std::int8_t g, std::int8_t e)
{
    return std::unique_ptr<AlignmentEngine>(new AlignmentEngine(m, n, g, e));
}

// AlignmentEngine end

// Graph start

Graph::Node::Node(std::uint32_t id, std::uint32_t code)
    : id(id), code(code), raw_id(id), inedges(), outedges(), aligned_nodes(), first_node(false), last_node(false) {}

Graph::Node::Node(std::uint32_t id, std::uint32_t code, std::uint32_t raw_id)
    : id(id), code(code), raw_id(raw_id), inedges(), outedges(), aligned_nodes(), first_node(false), last_node(false) {}


Graph::Node* Graph::Node::Successor(std::uint32_t label) const
{
    for (const auto& it : outedges)
    {
        auto jt = std::find(it->labels.begin(), it->labels.end(), label);
        if (jt != it->labels.end()) return it->head;
    }
    return nullptr;
}

std::uint32_t Graph::Node::Coverage() const
{
    std::unordered_set<std::uint32_t> labels;
    for (const auto& it : inedges)
        std::copy(it->labels.begin(), it->labels.end(), std::inserter(labels, labels.end()));
    for (const auto& it : outedges)
        std::copy(it->labels.begin(), it->labels.end(), std::inserter(labels, labels.end()));
    return labels.size();
}

Graph::Edge::Edge(Node* tail, Node* head, std::uint32_t label, std::uint32_t weight)
    : tail(tail), head(head), labels(1, label), weight(weight) {}

void Graph::Edge::AddSequence(std::uint32_t label, std::uint32_t w)
{
    labels.emplace_back(label);
    weight += w;
}

Graph::Graph() : num_codes_(0), coder_(256, -1), decoder_(256, -1),
                 sequences_(), nodes_(), edges_(), rank_to_node_(), consensus_()
{
    // Add Virtual start and end node
    auto start_node = AddNode(0); // start node
    auto end_node   = AddNode(0); // end node
    coder_[0] = num_codes_;
    decoder_[num_codes_++] = 0;
    AddEdge(start_node, end_node, 1); // connect start and end edge, but ignored this edge
}

Graph::Node* Graph::AddNode(std::uint32_t code)
{
    nodes_.emplace_back(new Node(nodes_.size() - 2, code));
    return nodes_.back().get();
}

Graph::Node* Graph::AddNode(std::uint32_t code, std::uint32_t raw_id) {
  nodes_.emplace_back(new Node(nodes_.size() - 2, code, raw_id));
  return nodes_.back().get();
}

void Graph::AddEdge(Node* tail, Node* head, std::uint32_t weight)
{
    for (const auto& it : tail->outedges)
        if (it->head == head)
        {
            it->AddSequence(sequences_.size(), weight);
            return;
        }

    edges_.emplace_back(new Edge(tail, head, sequences_.size(), weight));
    tail->outedges.emplace_back(edges_.back().get());
    head->inedges.emplace_back(edges_.back().get());
}

Graph::Node* Graph::AddSequence(const char* sequence, const std::vector<std::uint32_t>& weights,
                                std::uint32_t begin, std::uint32_t end, bool connect_start, bool connect_end)
{
    if (begin == end) return nullptr;
    Node* prev = connect_start ? nodes_[0].get() : nullptr;
    for (std::uint32_t i = begin; i < end; ++i)
    {
        auto curr = AddNode(coder_[sequence[i]]);
        if (connect_start)
        {
            connect_start = false;
            AddEdge(prev, curr, 1);
            curr->first_node = true;
        }
        else if (prev)
        {
            // both nodes contribute to the weight
            AddEdge(prev, curr, weights[i - 1] + weights[i]);
        }
        prev = curr;
    }
    if(connect_end) AddEdge(prev, nodes_[1].get(), 1), prev->last_node = true;
    return nodes_[nodes_.size() - (end - begin)].get();
}

void Graph::AddAlignment(const Alignment& alignment, const std::string& sequence, std::uint32_t weight)
{
    AddAlignment(alignment, sequence.c_str(), sequence.size(), weight);
}

void Graph::AddAlignment(const Alignment& alignment, const char* sequence, std::uint32_t sequence_len, std::uint32_t weight)
{
    std::vector<std::uint32_t> weights(sequence_len, weight);
    AddAlignment(alignment, sequence, sequence_len, weights);
}

void Graph::AddAlignment(const Alignment& alignment, const std::string& sequence, const std::string& quality)
{
    AddAlignment(alignment, sequence.c_str(), sequence.size(), quality.c_str(), quality.size());
}

void Graph::AddAlignment(const Alignment& alignment, const char* sequence, std::uint32_t sequence_len, const char* quality, std::uint32_t quality_len)
{
    std::vector<std::uint32_t> weights;
    for (std::uint32_t i = 0; i < quality_len; ++i)
    {
        weights.emplace_back(quality[i] - 33);  // Phred quality
    }
    AddAlignment(alignment, sequence, sequence_len, weights);
}

void Graph::AddAlignment(const Alignment& alignment, const std::string& sequence, const std::vector<std::uint32_t>& weights)
{
    AddAlignment(alignment, sequence.c_str(), sequence.size(), weights);
}

void Graph::AddAlignment(const Alignment& alignment, const char* sequence, std::uint32_t sequence_len, const std::vector<std::uint32_t>& weights)
{
    if (sequence_len == 0) return;
    if (sequence_len != weights.size())
    {
        throw std::invalid_argument(
            "[spoa::Graph::AddAlignment] error: "
            "sequence and weights are of unequal size!");
    }

    for (std::uint32_t i = 0; i < sequence_len; ++i)
    {
        if (coder_[sequence[i]] == -1)
        {
            coder_[sequence[i]] = num_codes_;
            decoder_[num_codes_++] = sequence[i];
        }
    }

    if (alignment.empty())
    {
        sequences_.emplace_back(AddSequence(sequence, weights, 0, sequence_len, true, true));
        TopologicalSort();
        return;
    }

    std::vector<std::uint32_t> valid;
    for (const auto& it : alignment)
    {
        if (it.second != -1)
        {
            if (it.second < 0 || it.second >= static_cast<std::int32_t>(sequence_len))
            {
                throw std::invalid_argument("[spoa::Graph::AddAlignment] error: invalid alignment");
            }
            valid.emplace_back(it.second);
        }
    }
    if (valid.empty()) throw std::invalid_argument("[spoa::Graph::AddAlignment] error: missing sequence in alignment");

    // add unaligned bases
    Node* begin = AddSequence(sequence, weights, 0, valid.front(), true, false);
    Node* prev = begin ? nodes_.back().get() : nullptr;
    Node* last = AddSequence(sequence, weights, valid.back() + 1, sequence_len, false, true);

    // add aligned bases
    for (const auto& it : alignment)
    {
        if (it.second == -1) continue;

        std::uint32_t code = coder_[sequence[it.second]];
        Node* curr = nullptr;
        if (it.first == -1) curr = AddNode(code);
        else
        {
            auto jt = nodes_[it.first + 2].get();
            if (jt->code == code) curr = jt;
            else
            {
                for (const auto& kt : jt->aligned_nodes)
                {
                    if (kt->code == code)
                    {
                        curr = kt;
                        break;
                    }
                }
                if (!curr)
                {
                    curr = AddNode(code);
                    for (const auto& kt : jt->aligned_nodes)
                    {
                        kt->aligned_nodes.emplace_back(curr);
                        curr->aligned_nodes.emplace_back(kt);
                    }
                    jt->aligned_nodes.emplace_back(curr);
                    curr->aligned_nodes.emplace_back(jt);
                }
            }
        }
        if (!begin) begin = curr, curr->first_node = true;
        if (prev)
        {
            // both nodes contribute to weight
            AddEdge(prev, curr, weights[it.second - 1] + weights[it.second]);
        }
        else
        {
            // link to start node
            AddEdge(nodes_[0].get(), curr, 1 + weights[it.second]);
        }
        prev = curr;
    }
    if (last) AddEdge(prev, last, weights[valid.back()] + weights[valid.back() + 1]);
    else AddEdge(prev, nodes_[1].get(), weights[valid.back()] + 1), prev->last_node = true;
    sequences_.emplace_back(begin);

    TopologicalSort();
}

void Graph::TopologicalSort()
{
    rank_to_node_.clear();

    std::vector<std::uint8_t> marks(nodes_.size(), 0);
    std::vector<bool> ignored(nodes_.size(), 0);

    std::stack<Node*> stack;
    for (const auto& it : nodes_)
    {
        if(it->code == 0) continue; // ignore virtual nodes
        if (marks[it->id] != 0) continue;
        stack.push(it.get());
        while (!stack.empty())
        {
            auto curr = stack.top();
            bool is_valid = true;
            if (marks[curr->id] != 2)
            {
                for (const auto& jt : curr->inedges)
                    if (jt->tail->code && marks[jt->tail->id] != 2)
                    {
                        stack.push(jt->tail);
                        is_valid = false;
                    }
                if (!ignored[curr->id])
                {
                    for (const auto& jt : curr->aligned_nodes)
                        if (marks[jt->id] != 2)
                        {
                            stack.push(jt);
                            ignored[jt->id] = true;
                            is_valid = false;
                        }
                }
                assert((is_valid || marks[curr->id] != 1) && "Graph is not a DAG");
                if (is_valid)
                {
                    marks[curr->id] = 2;
                    if (!ignored[curr->id])
                    {
                        rank_to_node_.emplace_back(curr);
                        for (const auto& jt : curr->aligned_nodes) rank_to_node_.emplace_back(jt);
                    }
                }
                else
                {
                    marks[curr->id] = 1;
                }
            }

            if (is_valid) stack.pop();
        }
    }
    assert(IsTopologicallySorted() && "Graph is not topologically sorted");
}

void Graph::rev_TopologicalSort()
{
    rank_to_node_rev_.clear();

    std::vector<std::uint8_t> marks(nodes_.size(), 0);
    std::vector<bool> ignored(nodes_.size(), 0);

    std::stack<Node*> stack;
    for (const auto& it : nodes_)
    {
        if(it->code == 0) continue; // ignore virtual nodes
        if (marks[it->id] != 0) continue;
        stack.push(it.get());
        while (!stack.empty())
        {
            auto curr = stack.top();
            bool is_valid = true;
            if (marks[curr->id] != 2)
            {
                for (const auto& jt : curr->outedges)
                    if (jt->head->code && marks[jt->head->id] != 2)
                    {
                        stack.push(jt->head);
                        is_valid = false;
                    }
                if (!ignored[curr->id])
                {
                    for (const auto& jt : curr->aligned_nodes)
                        if (marks[jt->id] != 2)
                        {
                            stack.push(jt);
                            ignored[jt->id] = true;
                            is_valid = false;
                        }
                }
                assert((is_valid || marks[curr->id] != 1) && "Graph is not a DAG");
                if (is_valid)
                {
                    marks[curr->id] = 2;
                    if (!ignored[curr->id])
                    {
                        rank_to_node_rev_.emplace_back(curr);
                        for (const auto& jt : curr->aligned_nodes) rank_to_node_rev_.emplace_back(jt);
                    }
                }
                else
                {
                    marks[curr->id] = 1;
                }
            }

            if (is_valid) stack.pop();
        }
    }
}

bool Graph::IsTopologicallySorted() const
{
    // virtual nodes are not sorted
    assert(nodes_.size() - 2 == rank_to_node_.size() && "Topological sort not called ");

    std::vector<bool> visited(nodes_.size(), 0);
    for (const auto& it : rank_to_node_)
    {
        for (const auto& jt : it->inedges)
            if (jt->tail->code && !visited[jt->tail->id]) return false;
        visited[it->id] = 1;
    }

    return true;
}

std::vector<std::uint32_t> Graph::InitializeMultipleSequenceAlignment(std::uint32_t* row_size) const
{
    std::vector<std::uint32_t> dst(nodes_.size());
    std::uint32_t j = 0;
    for (std::uint32_t i = 0; i < rank_to_node_.size(); ++i, ++j)
    {
        auto it = rank_to_node_[i];
        dst[it->id] = j;
        for (const auto& jt : it->aligned_nodes)
        {
            dst[jt->id] = j;
            ++i;
        }
    }
    if (row_size) *row_size = j;
    return dst;
}

std::vector<std::string> Graph::GenerateMultipleSequenceAlignment(bool include_consensus)
{
    std::uint32_t row_size = 0;
    auto node_id_to_column = InitializeMultipleSequenceAlignment(&row_size);

    std::vector<std::string> dst;
    for (std::uint32_t i = 0; i < sequences_.size(); ++i)
    {
        std::string row(row_size, '-');
        auto it = sequences_[i];
        while (true)
        {
            row[node_id_to_column[it->id]] = decoder_[it->code];
            it = it->Successor(i);
            if (it->id + 1 == 0 || !it) break;
        }
        dst.emplace_back(row);
    }
    if (include_consensus)
    {
        TraverseHeaviestBundle();
        std::string row(row_size, '-');
        for (const auto& it : consensus_) row[node_id_to_column[it->id]] = decoder_[it->code];
        dst.emplace_back(row);
    }

    return dst;
}

std::string Graph::GenerateConsensus()
{
    TraverseHeaviestBundle();
    std::string dst{};
    for (const auto& it : consensus_)
        dst += decoder_[it->code];
    return dst;
}

std::string Graph::GenerateConsensus(std::int32_t min_coverage)
{
    TraverseHeaviestBundle();
    std::string dst{};
    for (const auto& it : consensus_)
        if (static_cast<std::int32_t>(it->Coverage()) >= min_coverage)
            dst += decoder_[it->code];
    return dst;
}

std::string Graph::GenerateConsensus(std::int32_t min_coverage, std::vector<std::uint32_t> *summary)
{
    if (!summary)
        throw std::invalid_argument("[spoa::Graph::GenerateConsensus] error: invalid ptr to summary");
    summary->clear();

    auto dst = GenerateConsensus(min_coverage);
    for (const auto &it : consensus_)
    {
        if (static_cast<std::int32_t>(it->Coverage()) >= min_coverage)
        {
            summary->emplace_back(0);
            summary->back() += it->Coverage();
            for (const auto &jt : it->aligned_nodes) summary->back() += jt->Coverage();
        }
    }

    return dst;
}

std::string Graph::GenerateConsensus(std::vector<std::uint32_t>* summary, bool verbose)
{
    if (!summary)
        throw std::invalid_argument("[spoa::Graph::GenerateConsensus] error: invalid ptr to summary");

    auto dst = GenerateConsensus();

    summary->clear();
    if (!verbose)
    {
        for (const auto& it : consensus_)
        {
            summary->emplace_back(0);
            summary->back() += it->Coverage();
            for (const auto& jt : it->aligned_nodes) summary->back() += jt->Coverage();
        }
    }
    else
    {
        summary->resize((num_codes_ + 1) * consensus_.size(), 0);
        auto node_id_to_column = InitializeMultipleSequenceAlignment();

        for (std::uint32_t i = 0; i < sequences_.size(); ++i)
        {
            Node* it = sequences_[i];
            std::uint32_t c = 0, p, column = node_id_to_column[it->id];
            bool is_gap = false;
            while (true)
            {
                for (; c < consensus_.size(); ++c)
                {
                    if (node_id_to_column[consensus_[c]->id] < column) continue;
                    else
                    {
                        if (node_id_to_column[consensus_[c]->id] == column)
                        {
                            if (is_gap)
                            {
                                for (std::uint32_t j = p + 1; j < c; ++j)
                                    ++(*summary)[num_codes_ * consensus_.size() + j];
                            }
                            is_gap = true;
                            p = c;
                            ++(*summary)[it->code * consensus_.size() + c];
                        }
                        break;
                    }
                }
                if (c == consensus_.size() || !(it = it->Successor(i))) break;
                column = node_id_to_column[it->id];
            }
        }
    }

    return dst;
}

void Graph::TraverseHeaviestBundle()
{
    if (rank_to_node_.empty()) return;

    std::vector<Node*> predecessors(nodes_.size(), nullptr);
    std::vector<std::int64_t> scores(nodes_.size(), -1);
    Node* max = nullptr, *max_pre = nullptr;

    for (const auto& it : rank_to_node_)
    {
        for (const auto& jt : it->inedges)
            if ((scores[it->id] < jt->weight) ||
                ((jt->tail->id != UINT32_MAX && jt->tail->id != UINT32_MAX - 1) // may contain loop???
                    && scores[it->id] == jt->weight && scores[predecessors[it->id]->id] <= scores[jt->tail->id]))
            {
                if(! jt->tail->code) continue;
                scores[it->id] = jt->weight;
                predecessors[it->id] = jt->tail;
            }
        if (predecessors[it->id]) scores[it->id] += scores[predecessors[it->id]->id];
        if (!max || scores[max->id] < scores[it->id]) max = it;
    }

    if (!max->outedges.empty())
    {
        std::vector<std::uint32_t> node_id_to_rank(nodes_.size(), 0);
        for (std::uint32_t i = 0; i < rank_to_node_.size(); ++i)
            node_id_to_rank[rank_to_node_[i]->id] = i;
        while (max && !max->outedges.empty())
        {
            max_pre = max;
            max = BranchCompletion(node_id_to_rank[max->id], &scores, &predecessors);
        }
        max = max_pre;
    }

    // traceback
    consensus_.clear();
    while (predecessors[max->id])
    {
        consensus_.emplace_back(max);
        max = predecessors[max->id];
    }
    consensus_.emplace_back(max);

    std::reverse(consensus_.begin(), consensus_.end());
}

Graph::Node* Graph::BranchCompletion(std::uint32_t rank, std::vector<std::int64_t>* scores, std::vector<Node*>* predecessors)
{
    auto start = rank_to_node_[rank];
    for (const auto& it : start->outedges)
    {
        if (!it->head->code) continue;
        for (const auto& jt : it->head->inedges)
            if (jt->tail != start) (*scores)[jt->tail->id] = -1;
    }

    Node* max = nullptr;
    for (std::uint32_t i = rank + 1; i < rank_to_node_.size(); ++i)
    {
        auto it = rank_to_node_[i];
        (*scores)[it->id] = -1;
        (*predecessors)[it->id] = nullptr;

        for (const auto& jt : it->inedges)
        {
            if (! jt->tail->code) continue;
            if ((*scores)[jt->tail->id] == -1) continue;
            if (((*scores)[it->id] < jt->weight) ||
                ((*scores)[it->id] == jt->weight && (*scores)[(*predecessors)[it->id]->id] <= (*scores)[jt->tail->id]))
            {
                (*scores)[it->id] = jt->weight;
                (*predecessors)[it->id] = jt->tail;
            }
        }
        if ((*predecessors)[it->id]) (*scores)[it->id] += (*scores)[(*predecessors)[it->id]->id];
        if (!max || (*scores)[max->id] < (*scores)[it->id]) max = it;
    }

    return max;
}

std::vector<bool> Graph::ExtractSubgraph(const Node* begin, const Node* end) const
{
    std::vector<bool> dst(nodes_.size(), false);
    std::stack<const Node*> stack;
    stack.push(begin);
    // std::cerr << "begin.id = " << begin->id << ", end.id = " << end->id << std::endl;

    while (!stack.empty())
    {
        auto curr = stack.top();
        stack.pop();
        if(curr->id == end->id) { if(curr->code) dst[curr->id] = true; continue; }

        if (!dst[curr->id])
        {
            bool same_step = false;
            // std::cerr << "curr = " << curr->id << " with " << curr->aligned_nodes.size() << " nodes." << std::endl;
            for (const auto& it : curr->aligned_nodes)
            {
                if( it->id == end->id )
                {
                    same_step = true;
                    break;
                }
            }
            if(same_step)
            {
                // std::cerr << "Has same step: " << curr->id << std::endl;
                continue;
            }
            for (const auto& it : curr->inedges) stack.push(it->tail);
            dst[curr->id] = true;
        }
    }

    return dst;
}

Graph Graph::Subgraph(std::uint32_t begin, std::uint32_t end, std::vector<const Node*>* subgraph_to_graph) const
{
    if (!subgraph_to_graph)
        throw std::invalid_argument("[spoa::Graph::Subgraph] error: invalid ptr to subgraph_to_graph");

    auto is_in_subgraph = ExtractSubgraph(nodes_[end].get(), nodes_[begin].get());

    // init subgraph
    Graph subgraph{};
    subgraph.num_codes_ = num_codes_;
    subgraph.coder_ = coder_;
    subgraph.decoder_ = decoder_;
    // subgraph.sequences_ = TODO(rvaser) maybe add sequences

    // create a map from subgraph nodes to graph nodes and vice versa
    subgraph_to_graph->clear();
    subgraph_to_graph->resize(nodes_.size(), nullptr);

    std::vector<Node*> graph_to_subgraph(nodes_.size(), nullptr);

    for (const auto& it : nodes_)
    {
        if (!is_in_subgraph[it->id]) continue;
        subgraph.AddNode(it->code, it->raw_id);
        graph_to_subgraph[it->id] = subgraph.nodes_.back().get();
        (*subgraph_to_graph)[subgraph.nodes_.back()->id] = it.get();
    }

    // connect nodes
    for (const auto& it : nodes_)
    {
        if (!is_in_subgraph[it->id]) continue;
        auto jt = graph_to_subgraph[it->id];
        for (const auto& kt : it->inedges)
        {
            if (graph_to_subgraph[kt->tail->id])
                subgraph.AddEdge(graph_to_subgraph[kt->tail->id], jt, kt->weight);
        }
        for (const auto& kt : it->aligned_nodes)
        {
            if (graph_to_subgraph[kt->id])
                jt->aligned_nodes.emplace_back(graph_to_subgraph[kt->id]);
        }
    }

    subgraph.TopologicalSort();

    return subgraph;
}

std::vector<bool> Graph::ExtractSubgraphLeft(const Graph::Node* begin) const
{
    std::vector<bool> dst(nodes_.size(), false);
    std::stack<const Graph::Node*> stack;
    stack.push(begin);
    // std::cerr << "begin.id = " << begin->id << ", end.id = " << end->id << std::endl;

    while (!stack.empty())
    {
        auto curr = stack.top();
        stack.pop();
        if(! curr->code)
        {
            // Graph comes to end
            continue;
        }

        if (!dst[curr->id])
        {
            // std::cerr << "curr = " << curr->id << " with " << curr->aligned_nodes.size() << " nodes." << std::endl;
            for (const auto& it : curr->inedges) stack.push(it->tail);
            dst[curr->id] = true;
        }
    }

    return dst;
}

std::vector<bool> Graph::ExtractSubgraphRight(const Graph::Node* begin) const
{
    std::vector<bool> dst(nodes_.size(), false);
    std::stack<const Graph::Node*> stack;
    stack.push(begin);
    // std::cerr << "begin.id = " << begin->id << ", end.id = " << end->id << std::endl;

    while (!stack.empty())
    {
        auto curr = stack.top();
        stack.pop();
        if(! curr->code)
        {
            // Graph comes to end
            continue;
        }

        if (!dst[curr->id])
        {
            // std::cerr << "curr = " << curr->id << " with " << curr->aligned_nodes.size() << " nodes." << std::endl;
            for (const auto& it : curr->outedges) stack.push(it->head);
            dst[curr->id] = true;
        }
    }

    return dst;
}

Graph Graph::Subgraph_left(std::uint32_t end) const
{
    if(end == std::numeric_limits<std::uint32_t>::max())
    {
        auto G = Graph();
        G.nodes()[0].get()->raw_id = nodes_[0].get()->raw_id; // copy raw id for alignment
        return G;
    }
    auto is_in_subgraph = ExtractSubgraphLeft(nodes_[end + 2].get());
    nodes_[end + 2]->last_node = true; // this is the end node

    // init subgraph
    Graph subgraph{};
    subgraph.nodes()[0].get()->raw_id = nodes_[0].get()->raw_id; // copy raw id for alignment
    subgraph.num_codes_ = num_codes_;
    subgraph.coder_ = coder_;
    subgraph.decoder_ = decoder_;

    // create a map from subgraph nodes to graph nodes and vice versa
    std::vector<const Node*> subgraph_to_graph;
    subgraph_to_graph.resize(nodes_.size(), nullptr);

    std::vector<Node*> graph_to_subgraph(nodes_.size(), nullptr);

    for (const auto& it : nodes_)
    {
        if (!it->code || !is_in_subgraph[it->id]) continue;
        subgraph.AddNode(it->code, it->raw_id);
        graph_to_subgraph[it->id] = subgraph.nodes_.back().get();
        subgraph_to_graph[subgraph.nodes_.back()->id] = it.get();
    }

    // connect nodes
    Node *raw_graph_start_node = subgraph.nodes().at(0).get(), *raw_graph_end_node = subgraph.nodes().at(1).get();
    for (const auto& it : nodes_)
    {
        if (!it->code || !is_in_subgraph[it->id]) continue;
        auto jt = graph_to_subgraph[it->id];
        for (const auto& kt : it->inedges)
        {
            if (kt->tail->code && graph_to_subgraph[kt->tail->id])
                subgraph.AddEdge(graph_to_subgraph[kt->tail->id], jt, kt->weight);
        }
        for (const auto& kt : it->aligned_nodes)
        {
            if (graph_to_subgraph[kt->id])
                jt->aligned_nodes.emplace_back(graph_to_subgraph[kt->id]);
        }
        // connect start && end virtual node
        if (it->first_node)
        {
            subgraph.AddEdge(raw_graph_start_node, jt, 1); // TODO: Sub alignment may no need to add weight
            jt->first_node = true;
        }
        if (it->last_node)
        {
            subgraph.AddEdge(jt, raw_graph_end_node, 1); // TODO: same as before
            jt->last_node = true;
        }
    }

    subgraph.TopologicalSort();

    return subgraph;
}

Graph Graph::Subgraph_right(std::uint32_t begin, std::uint32_t prev_id) const
{
    if(begin == std::numeric_limits<std::uint32_t>::max())
    {
        auto G = Graph();
        G.nodes()[0].get()->raw_id = (prev_id == std::numeric_limits<std::uint32_t>::max() ? nodes_[0].get()->raw_id : nodes_[prev_id + 2].get()->raw_id); // copy raw id for alignment
        return G;
    }
    auto is_in_subgraph = ExtractSubgraphRight(nodes_[begin + 2].get());
    nodes_[begin + 2]->first_node = true; // this is the first node

    // init subgraph
    Graph subgraph{};
    subgraph.nodes()[0].get()->raw_id = (prev_id == std::numeric_limits<std::uint32_t>::max() ? nodes_[0].get()->raw_id : nodes_[prev_id + 2].get()->raw_id); // copy raw id for alignment
    subgraph.num_codes_ = num_codes_;
    subgraph.coder_ = coder_;
    subgraph.decoder_ = decoder_;

    // create a map from subgraph nodes to graph nodes and vice versa
    std::vector<const Node*> subgraph_to_graph;
    subgraph_to_graph.resize(nodes_.size(), nullptr);

    std::vector<Node*> graph_to_subgraph(nodes_.size(), nullptr);

    for (const auto& it : nodes_)
    {
        if (!it->code || !is_in_subgraph[it->id]) continue;
        subgraph.AddNode(it->code, it->raw_id);
        graph_to_subgraph[it->id] = subgraph.nodes_.back().get();
        subgraph_to_graph[subgraph.nodes_.back()->id] = it.get();
    }

    // connect nodes
    Node *raw_graph_start_node = subgraph.nodes().at(0).get(), *raw_graph_end_node = subgraph.nodes().at(1).get();
    for (const auto& it : nodes_)
    {
        if (!it->code || !is_in_subgraph[it->id]) continue;
        auto jt = graph_to_subgraph[it->id];
        for (const auto& kt : it->inedges)
        {
            if (kt->tail->code && graph_to_subgraph[kt->tail->id])
                subgraph.AddEdge(graph_to_subgraph[kt->tail->id], jt, kt->weight);
        }
        for (const auto& kt : it->aligned_nodes)
        {
            if (graph_to_subgraph[kt->id])
                jt->aligned_nodes.emplace_back(graph_to_subgraph[kt->id]);
        }
        // connect start && end virtual node
        if (it->first_node)
        {
            subgraph.AddEdge(raw_graph_start_node, jt, 1); // TODO: Sub alignment may no need to add weight
            jt->first_node = true;
        }
        if (it->last_node)
        {
            subgraph.AddEdge(jt, raw_graph_end_node, 1); // TODO: same as before
            jt->last_node = true;
        }

    }

    subgraph.TopologicalSort();

    return subgraph;

}

std::uint32_t Graph::depth() const
{
    // TODO: modify it
    return rank_to_node_.size();
}

void Graph::UpdateAlignment(const std::vector<const Node*>& subgraph_to_graph, Alignment* alignment) const
{
    for (auto& it : *alignment)
    {
        if (it.first != -1)
            it.first = subgraph_to_graph[it.first]->id;
    }
}

#define MATCH_SCORE(a, b, M, X) ((a) == (b) ? (M) : (X))

int Graph::CalculateAlignmentScore(const Alignment& alignment, const char* sequence, size_t str_len, int gap_open, int gap_extend, int match, int mismatch) const
{
    // TODO: real go on graph instead of reading the alignment result
    int score = 0;
    bool in_gap_graph = false;  // graph deletion
    bool in_gap_seq = false;    // sequence insertion
    auto graph_size = nodes_.size();
    int seq_place = 0; // seq id
    Node* node_point = nodes_[0].get(), *new_node_point; // node id

    while(seq_place != str_len || node_point->outedges.empty())
    {
        if(node_point->outedges.empty())
        {
            // must be in seq_place
            if (!in_gap_seq)
            {
                score += gap_open + gap_extend;
                in_gap_seq = true;
            }
            else
            {
                score += gap_extend;
            }
            in_gap_graph = false;
            seq_place ++;
        }
        else
        {
            if(seq_place == str_len)
            {
                if (!in_gap_graph)
                {
                    score += gap_open + gap_extend;
                    in_gap_graph = true;
                }
                else
                {
                    score += gap_extend;
                }
                in_gap_seq = false;
                new_node_point = node_point;
                for(auto & edge : node_point->outedges)
                {
                    if (edge->head->id == UINT_MAX) continue;
                    if(std::find_if(alignment.begin(), alignment.end(), [&](const align_pair_t &x) {return x.first == edge->head->id;}) != alignment.end())
                    {
                        new_node_point = nodes_[edge->head->id + 2].get();
                        break;
                    }
                }
                node_point = new_node_point;
            }
            else
            {
                bool found = false;
                new_node_point = node_point;
                for(auto &edge: node_point->outedges)
                {
                    if (edge->head->id == UINT_MAX) continue;
                    if(std::find_if(alignment.begin(), alignment.end(), [&](const align_pair_t &x) {return x.first == edge->head->id && x.second == seq_place;}) != alignment.end()) // match/mismatch
                    {
                        char node_char = decoder_.at(node_point->code);  // may throw std::out_of_range in nodes_
                        char seq_char = sequence[seq_place];
                        score += MATCH_SCORE(node_char, seq_char, match, mismatch);
                        new_node_point = nodes_[edge->head->id + 2].get();
                        seq_place  ++;
                        // reset gap condition
                        in_gap_graph = false;
                        in_gap_seq = false;
                        found = true;
                        break;
                    }
                    else if(std::find_if(alignment.begin(), alignment.end(), [&](const align_pair_t& x) {return x.first == edge->head->id; }) != alignment.end()) // gap aligned with graph node
                    {
                        if (!in_gap_graph)
                        {
                            score += gap_open + gap_extend;
                            in_gap_graph = true;
                        }
                        else
                        {
                            score += gap_extend;
                        }
                        in_gap_seq = false;
                        new_node_point = nodes_[edge->head->id + 2].get();
                        found = true;
                        break;
                    }
                }
                node_point = new_node_point;
                if(! found) // gap aligned with sequence
                {
                    if (!in_gap_seq)
                    {
                        score += gap_open + gap_extend;
                        in_gap_seq = true;
                    }
                    else
                    {
                        score += gap_extend;
                    }
                    in_gap_graph = false;
                    seq_place ++;
                }
            }
        }
    }
    return score;

    for (const auto& step : alignment)
    {
        const auto& graph_idx = step.first;
        const auto& seq_pos = step.second;

        if (graph_idx == -1 && seq_pos == -1)
        {
            throw std::invalid_argument("Invalid alignment step: both indices are -1");
        }

        // Case 1: match/mismatch
        if (graph_idx != -1 && seq_pos != -1)
        {
            if (graph_idx < 0 || static_cast<std::size_t>(graph_idx) >= graph_size)
            {
                throw std::out_of_range("Graph node index out of range: " + std::to_string(graph_idx));
            }
            if (seq_pos < 0 || static_cast<std::size_t>(seq_pos) >= str_len)
            {
                throw std::out_of_range("Sequence position out of range: " + std::to_string(seq_pos));
            }

            try
            {
                char node_char = decoder_.at(rank_to_node_[graph_idx]->code);  // may throw std::out_of_range in nodes_
                char seq_char = sequence[seq_pos];
                score += MATCH_SCORE(node_char, seq_char, match, mismatch);
            }
            catch (const std::out_of_range&)
            {
                throw std::out_of_range("[CalcualteAlignmentScore] Error: Graph Node id = " + std::to_string(graph_idx) + "is out of range");
            }

            // reset gap condition
            in_gap_graph = false;
            in_gap_seq = false;
        }
        // Case 2: deletion
        else if (seq_pos == -1)
        {
            if (graph_idx < 0 || static_cast<std::size_t>(graph_idx) >= graph_size)
            {
                throw std::out_of_range("[CalcualteAlignmentScore] Error: Invalid graph index in deletion: " + std::to_string(graph_idx));
            }

            if (!in_gap_graph)
            {
                score += gap_open + gap_extend;
                in_gap_graph = true;
            }
            else
            {
                score += gap_extend;
            }
            in_gap_seq = false;
        }
        // Case 3: insertion
        else if (graph_idx == -1)
        {
            if (seq_pos < 0 || static_cast<std::size_t>(seq_pos) >= str_len)
            {
                throw std::out_of_range("Invalid sequence position in insertion: " + std::to_string(seq_pos));
            }

            if (!in_gap_seq)
            {
                score += gap_open + gap_extend;
                in_gap_seq = true;
            }
            else
            {
                score += gap_extend;
            }
            in_gap_graph = false;
        }
    }

    return score;
}

#undef MATCH_SCORE

void Graph::PrintDot(const std::string& path) const
{
    if (path.empty()) return;
    std::ofstream os(path);

    std::vector<std::int32_t> consensus_rank(nodes_.size(), -1);
    std::int32_t rank = 0;
    for (const auto& it : consensus_)
        consensus_rank[it->id] = rank++;

    os << "digraph " << sequences_.size() << " {\n"
       << "  graph [rankdir = LR]" << std::endl;
    for (const auto& it : nodes_)
    {
        if (it->id == UINT_MAX -1 || it->id == UINT_MAX) continue;
        os << "  " << it->id << "[label = \"" << it->raw_id << " - "
           << static_cast<char>(decoder_[it->code]) << "\"";
        if (consensus_rank[it->id] != -1) os << ", style = filled, fillcolor = goldenrod1";
        os << "]" << std::endl;

        for (const auto& jt : it->outedges)
        {
            if (jt->head->id == UINT_MAX - 1 || jt->head->id == UINT_MAX) continue;
            os << "  " << it->id << " -> " << jt->head->id
               << " [label = \"" << jt->weight << "\"";
            if (consensus_rank[it->id] + 1 == consensus_rank[jt->head->id]) os << ", color = goldenrod1";
            os << "]" << std::endl;
        }
        for (const auto& jt : it->aligned_nodes)
        {
            if (jt->id > it->id)
                os << "  " << it->id << " -> " << jt->id << " [style = dotted, arrowhead = none]" << std::endl;
        }
    }
    os << "}" << std::endl;

    os.close();
}

void Graph::Clear()
{
    num_codes_ = 0;
    std::fill(coder_.begin(), coder_.end(), -1);
    std::fill(decoder_.begin(), decoder_.end(), -1);
    sequences_.clear();
    nodes_.clear();
    edges_.clear();
    rank_to_node_.clear();
    consensus_.clear();
}

// Graph end
