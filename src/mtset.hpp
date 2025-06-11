#ifndef __MUTLITHREAD_SET__
#define __MUTLITHREAD_SET__

#include <atomic>
#include <thread>
#include <shared_mutex>
#include <mutex>
#include <vector>
#include <cassert>
#include <algorithm>
#include <cstddef>

typedef int32_t align_status_t;

typedef struct hirschberg_status_t
{
    align_status_t seq1, g2; // first is graph id, second is sequence id
    bool operator < (const hirschberg_status_t& b) const // TODO: modify it
    {
        return (seq1 == (align_status_t)(-1) || g2 == (align_status_t)(-1)) ? true : false;
    }
    bool operator == (const hirschberg_status_t& b) const
    {
        return seq1 == b.seq1 && g2 == b.g2;
    }
    bool operator > (const hirschberg_status_t& b) const
    {
        return !(*this == b || *this < b);
    }
    hirschberg_status_t() {}
    hirschberg_status_t(int a, int b) : seq1(a), g2(b) {}
} hirschberg_status_t;

class ConcurrentSet_hirschberg_status_t
{
private:
    struct Node
    {
        hirschberg_status_t *value;
        std::atomic<Node*> next;

        Node(hirschberg_status_t *val) : value(val), next(nullptr) {}
    };

    std::vector<std::atomic<Node*>> *table;
    int capacity;
    std::atomic<size_t> set_size;  // 添加一个原子计数器
    std::mutex* mutexes;  // 读写锁数组

public:
    ConcurrentSet_hirschberg_status_t(int cap = 23) : capacity(cap)
    {
        table = new std::vector<std::atomic<Node*>>(capacity);
        for (int i = 0; i < capacity; ++i) {
            std::atomic_init(&(*table)[i], nullptr);
        }
        mutexes = new std::mutex[capacity];  // 初始化读写锁数组
        set_size.store(0);
    }

    ~ConcurrentSet_hirschberg_status_t() {
        delete table;      // 在析构函数中释放内存
        delete[] mutexes;  // 删除读写锁数组
    }

    void insert(int a, int b);
    bool getFinalSnapshotAndFree(std::vector<hirschberg_status_t*> &v);
    
    size_t getSize();
};


#endif