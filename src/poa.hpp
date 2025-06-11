// based on spoa
// Copyright (c) 2020 Robert Vaser


#ifndef __POA_DEFINE__
#define __POA_DEFINE__

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

typedef std::pair<std::int32_t, std::int32_t> align_pair_t;
using Alignment = std::vector<align_pair_t>;

class Graph
{
public:
    Graph();

    Graph(const Graph&) = delete;
    Graph& operator=(const Graph&) = delete;

    Graph(Graph&&) = default;
    Graph& operator=(Graph&&) = default;

    ~Graph() = default;

    struct Node;
    struct Edge;

    struct Node
    {
    public:
        Node(std::uint32_t id, std::uint32_t code);
        Node(std::uint32_t id, std::uint32_t code, std::uint32_t raw_code);

        Node(const Node&) = default;
        Node& operator=(const Node&) = default;

        Node(Node&&) = default;
        Node& operator=(Node&&) = default;

        Node* Successor(std::uint32_t label) const;

        std::uint32_t Coverage() const;

        std::uint32_t id;
        std::uint32_t code;
        std::uint32_t raw_id;
        std::vector<Edge*> inedges;
        std::vector<Edge*> outedges;
        std::vector<Node*> aligned_nodes;
        bool first_node, last_node;
    };
    struct Edge
    {
    public:
        Edge(Node* tail, Node* head, std::uint32_t label, std::uint32_t weight);

        Edge(const Edge&) = default;
        Edge& operator=(const Edge&) = default;

        Edge(Edge&&) = default;
        Edge& operator=(Edge&&) = default;

        void AddSequence(std::uint32_t label, std::uint32_t weight = 1);

        Node* tail;
        Node* head;
        std::vector<std::uint32_t> labels;
        std::int64_t weight;
    };

    const std::vector<std::unique_ptr<Node>>& nodes() const { return nodes_; }

    const std::vector<std::unique_ptr<Edge>>& edges() const { return edges_; }

    const std::vector<Node*>& rank_to_node() const { return rank_to_node_; }

    const std::vector<Node*>& rank_to_node_rev() { rev_TopologicalSort(); return rank_to_node_rev_; }

    const std::vector<Node*>& sequences() const { return sequences_; }

    std::uint32_t num_codes() const { return num_codes_; }

    std::uint8_t coder(std::uint8_t c) const { return coder_[c]; }

    std::uint8_t decoder(std::uint8_t code) const { return decoder_[code]; }

    std::uint32_t depth() const;

    const std::vector<Node*>& consensus() const { return consensus_; }

    void AddAlignment(
        const Alignment& alignment,
        const std::string& sequence,
        std::uint32_t weight = 1);

    void AddAlignment(
        const Alignment& alignment,
        const std::string& sequence,
        const std::vector<std::uint32_t>& weights);

    void AddAlignment(
        const Alignment& alignment,
        const std::string& sequence,
        const std::string& quality);

    void AddAlignment(
        const Alignment& alignment,
        const char* sequence, std::uint32_t sequence_len,
        std::uint32_t weight = 1);

    void AddAlignment(
        const Alignment& alignment,
        const char* sequence, std::uint32_t sequence_len,
        const std::vector<std::uint32_t>& weights);

    void AddAlignment(
        const Alignment& alignment,
        const char* sequence, std::uint32_t sequence_len,
        const char* quality, std::uint32_t quality_len);

    std::vector<std::string> GenerateMultipleSequenceAlignment(
        bool include_consensus = false);

    std::string GenerateConsensus();

    std::string GenerateConsensus(std::int32_t min_coverage);

    std::string GenerateConsensus(std::int32_t min_coverage, std::vector<std::uint32_t> *summary);

    std::string GenerateConsensus(std::vector<std::uint32_t>* summary, bool verbose = false);

    Graph Subgraph(std::uint32_t begin, std::uint32_t end, std::vector<const Node*>* subgraph_to_graph) const;
    Graph Subgraph_left(std::uint32_t end) const;
    Graph Subgraph_right(std::uint32_t begin, std::uint32_t prev_id) const;

    void UpdateAlignment(const std::vector<const Node*>& subgraph_to_graph, Alignment* alignment) const;

    int CalculateAlignmentScore(const Alignment& alignment, const char* sequence, size_t str_len, int gap_open, int gap_extend, int match, int mismatch) const;

    // print with Graphviz
    void PrintDot(const std::string& path) const;

    void Clear();

private:
    Node* AddNode(std::uint32_t code);
    Node* AddNode(std::uint32_t code, std::uint32_t raw_id);


    void AddEdge(Node* tail, Node* head, std::uint32_t weight);

    Node* AddSequence(
        const char* sequence,
        const std::vector<std::uint32_t>& weights,
        std::uint32_t begin,
        std::uint32_t end,
        bool connect_start,
        bool connect_end);

    void TopologicalSort();
    void rev_TopologicalSort();

    bool IsTopologicallySorted() const;

    void TraverseHeaviestBundle();

    Node* BranchCompletion(std::uint32_t rank, std::vector<std::int64_t>* scores, std::vector<Node*>* predecessors);

    std::vector<bool> ExtractSubgraph(const Node* begin, const Node* end) const;
    std::vector<bool> ExtractSubgraphLeft(const Graph::Node* begin) const;
    std::vector<bool> ExtractSubgraphRight(const Graph::Node* begin) const;

    std::vector<std::uint32_t> InitializeMultipleSequenceAlignment(std::uint32_t* row_size = nullptr) const;

    std::uint32_t num_codes_;
    std::vector<std::int32_t> coder_;
    std::vector<std::int32_t> decoder_;
    std::vector<Node*> sequences_;
    std::vector<std::unique_ptr<Node>> nodes_;
    std::vector<std::unique_ptr<Edge>> edges_;
    std::vector<Node*> rank_to_node_;
    std::vector<Node*> rank_to_node_rev_; // calculate broader condition
    std::vector<Node*> consensus_;
};

class AlignmentEngine
{
public:
    ~AlignmentEngine() = default;

    static std::unique_ptr<AlignmentEngine> Create(
        std::int8_t m,   // match
        std::int8_t n,   // mismatch
        std::int8_t g,   // gap open
        std::int8_t e);  // gap extend


    void Prealloc(
        std::uint32_t max_sequence_len,
        std::uint8_t alphabet_size);


    Alignment Align(
        const char* sequence, std::uint32_t sequence_len,
        const Graph& graph,
        std::int32_t* score = nullptr);
    Alignment linearspace_Align(
        const char* seq, std::uint32_t len, Graph &G, 
        std::int32_t *score, 
        int threads, int dp_threads, bool multidp);

protected:
    AlignmentEngine(
        std::int8_t m,
        std::int8_t n,
        std::int8_t g,
        std::int8_t e);

    std::int64_t WorstCaseAlignmentScore(
        std::int64_t sequence_len,
        std::int64_t graph_len) const;

    void Realloc(
      std::uint64_t matrix_width,
      std::uint64_t matrix_height,
      std::uint8_t num_codes);

    void Initialize(
        const char* sequence, std::uint32_t sequence_len,
        const Graph& graph) noexcept;

    struct Implementation
    {
        std::vector<std::uint32_t> node_id_to_rank;
        std::vector<std::int32_t> sequence_profile;
        std::vector<std::int32_t> M;
        std::int32_t* H;
        std::int32_t* F;
        std::int32_t* E;

        Implementation()
            : node_id_to_rank(),
            sequence_profile(),
            M(),
            H(nullptr),
            F(nullptr),
            E(nullptr) {}
    };
    std::unique_ptr<Implementation> pimpl_;
    std::int8_t m_;
    std::int8_t n_;
    std::int8_t g_;
    std::int8_t e_;
};

#endif