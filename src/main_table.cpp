#include <atomic>
#include <thread>
#include <shared_mutex>
#include <mutex>
#include <vector>
#include <cassert>
#include <algorithm>

#include "../include/murmurhash/murmurhash3.h"
#include "../include/boost/include/boost/multiprecision/cpp_int.hpp"

using namespace boost::multiprecision;

class ConcurrentSet_int_t
{
private:
    struct Node
    {
        int value;
        std::atomic<Node*> next;

        Node(int val) : value(val), next(nullptr) {}
    };

    std::vector<std::atomic<Node*>> *table;
    int capacity;
    std::atomic<size_t> set_size;  // 添加一个原子计数器
    std::shared_mutex* mutexes;  // 读写锁数组

public:
    ConcurrentSet_int_t(int cap = 967) : capacity(cap)
    {
        table = new std::vector<std::atomic<Node*>>(capacity);
        for (int i = 0; i < capacity; ++i) {
            std::atomic_init(&(*table)[i], nullptr);
        }
        mutexes = new std::shared_mutex[capacity];  // 初始化读写锁数组
        set_size.store(0);
    }

    ~ConcurrentSet_int_t() {
        delete table;      // 在析构函数中释放内存
        delete[] mutexes;  // 删除读写锁数组
    }

    ConcurrentSet_int_t(const ConcurrentSet_int_t&) = delete;
    ConcurrentSet_int_t(ConcurrentSet_int_t &&) = default;

    void insert(int value)
    {
        int index = (value > 0 ? value : -value) % capacity;

        Node* newNode = new Node(value);
        newNode->next.store(nullptr);
        std::lock_guard<std::shared_mutex> lock(mutexes[index]);  // 获取写锁

        Node* curr = (*table)[index];
        Node* prev = nullptr;

        while (curr != nullptr) {
            if (curr->value == value) {
                delete newNode;
                return;
            }
            prev = curr;
            curr = curr->next.load();
        }

        newNode->next.store((*table)[index]);
        (*table)[index] = newNode;
        set_size.fetch_add(1);

        return;

    }
    bool getFinalSnapshotAndFree(std::vector<int> &v)
    {
        std::vector<std::shared_lock<std::shared_mutex>> locks(mutexes, mutexes + capacity);  // 获取所有哈希表的读锁
        Node *curr, *pre;
        v.reserve(set_size.load());
        for(int i = 0; i < capacity; ++ i)
        {
            curr = (*table)[i];
            while(curr != nullptr)
            {
                v.emplace_back(std::move(curr -> value));
                pre = curr;
                curr = curr->next.load();
                delete pre;
            }
        }
        std::sort(v.begin(), v.end());
        return true;
    }

    size_t getSize()
    {
        std::vector<std::shared_lock<std::shared_mutex>> locks(mutexes, mutexes + capacity);  // 获取所有哈希表的读锁
        return set_size.load();
    }
};

void get_vec(int128_t x, std::vector <int> &ans, int one_size, int all_size)
{
    ans.clear();
    int mod;
    while(x)
    {
        mod = int(x % all_size);
        ans.emplace_back(mod - one_size);
        x /= all_size;
    }
}

int hash(const std::vector<int>& K)
{
    int ans;
    MurmurHash3_x86_32(K.data(), K.size() * sizeof(int), 0, &ans);
    return ans;
}


void block_process(ConcurrentSet_int_t &s, int128_t begin_, int128_t end_, int max_len, int one_size, int all_size)
{
    // fprintf(stderr, "%s %s\n", begin_.str().c_str(), end_.str().c_str());
    std::vector <int> tmp;
    for(int128_t i = begin_; i < end_; ++ i)
    {
        get_vec(i, tmp, one_size, all_size);
        if(tmp.size() > max_len) return;
        else
        {
            while(tmp.size() < max_len) tmp.emplace_back(-one_size);
            // for(auto &j: tmp) fprintf(stderr, "%d ", j); fprintf(stderr, "\n");
            // fprintf(stderr, "-> %d\n", hash(tmp));
            s.insert(hash(tmp));
        }
    }
}


int main(int argc, char **argv)
{
    int vector_size = atoi(argv[1]), max_size = atoi(argv[2]), threads = atoi(argv[3]);
    if(! vector_size || ! max_size)
    {
        fprintf(stderr, "Error: maybe lose arugments.\n");
        fprintf(stderr, "Arguments List: %s vector_size max_size [threads]\n", argv[0]);
    }
    if(threads < 0) threads = 1;
    fprintf(stderr, "vector size = %d, max range = %d, use %d threads\n", vector_size, max_size, threads);
    int128_t now_max = 2 * max_size + 1, max_range = 1;
    for(int i = 0; i < vector_size; ++ i) max_range *= now_max;
    fprintf(stderr, "has %s items\n", max_range.str().c_str());
    ConcurrentSet_int_t s;

    std::vector <std::thread> threads_;
    int128_t begin_, end_;

    for(int i = 0; i < threads; ++ i)
    {
        begin_ = max_range / threads * i;
        end_ = max_range / threads * (i + 1);
        threads_.emplace_back(block_process, std::ref(s), begin_, end_, vector_size, max_size, 2 * max_size + 1);
    }
    for(auto &entry: threads_) entry.join();
    block_process(s, max_range / threads * threads, max_range, vector_size, max_size, 2 * max_size + 1);

    fprintf(stderr, "has %zu items\n", s.getSize());
    std::vector <int> final_v;

    // s.getFinalSnapshotAndFree(final_v);
    // for(auto &i: final_v) fprintf(stderr, "%d ", i); fprintf(stderr, "\n");

    return 0;
}