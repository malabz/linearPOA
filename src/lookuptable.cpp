#include <iostream>
#include <unordered_set>
#include <atomic>
#include <thread>
#include <string>
#include <vector>
#include <shared_mutex>
#include <mutex>
#include <algorithm>
#include <execution>

class ConcurrentSet {
private:
    struct Node {
        size_t value;
        std::atomic<Node*> next;

        Node(size_t val) : value(val), next(nullptr) {}
    };

    std::vector<std::atomic<Node*>> *table;
    int capacity;
    std::atomic<size_t> size;  // 添加一个原子计数器
    std::shared_mutex* mutexes;  // 读写锁数组

public:
    ConcurrentSet(int cap = 29) : capacity(cap), size(0) {
        table = new std::vector<std::atomic<Node*>>(capacity);
        for (int i = 0; i < capacity; ++i) {
            std::atomic_init(&(*table)[i], nullptr);
        }
        mutexes = new std::shared_mutex[capacity];  // 初始化读写锁数组
    }

    ~ConcurrentSet() {
        delete table;  // 在析构函数中释放内存
        delete[] mutexes;  // 删除读写锁数组
    }

    ConcurrentSet(const ConcurrentSet&) = delete;
    ConcurrentSet(ConcurrentSet &&) = default;

    void insert(size_t value) {
        int index = hash(value);

        Node* newNode = new Node(value);
        newNode->next.store(nullptr);

        std::unique_lock<std::shared_mutex> lock(mutexes[index], std::defer_lock);  // 获取写锁
        while (true) {
            lock.lock();  // 加锁
            Node* curr = (*table)[index];
            Node* prev = nullptr;

            while (curr != nullptr) {
                if (curr->value == value) {
                    lock.unlock();  // 解锁
                    delete newNode;
                    return;
                }
                prev = curr;
                curr = curr->next.load();
            }

            newNode->next.store((*table)[index]);
            (*table)[index] = newNode;
            size.fetch_add(1);

            lock.unlock();  // 解锁
            return;
        }
    }

    void erase(size_t value) {
        int index = hash(value);

        std::unique_lock<std::shared_mutex> lock(mutexes[index], std::defer_lock);  // 获取写锁
        lock.lock();  // 加锁

        Node* curr = (*table)[index];
        Node* prev = nullptr;

        while (curr != nullptr) {
            if (curr->value == value) {
                Node* next = curr->next.load();
                if (prev == nullptr) {
                    (*table)[index] = next;
                } else {
                    prev->next.store(next);
                }
                delete curr;
                size.fetch_sub(1);

                lock.unlock();  // 解锁
                return;
            }
            prev = curr;
            curr = curr->next.load();
        }

        lock.unlock();  // 解锁
    }

    bool contains(size_t value) {
        int index = hash(value);

        std::shared_lock<std::shared_mutex> lock(mutexes[index], std::defer_lock);  // 获取读锁
        lock.lock();  // 加锁

        Node* curr = (*table)[index];

        while (curr != nullptr) {
            if (curr->value == value) {
                lock.unlock();  // 解锁
                return true;
            }
            curr = curr->next.load();
        }

        lock.unlock();  // 解锁
        return false;
    }

    size_t getSize() const {
        std::vector<std::shared_lock<std::shared_mutex>> locks(mutexes, mutexes + capacity);  // 获取所有哈希表的读锁
        return size.load();
    }

    bool getFinalSnapshotAndFree(std::vector<size_t> &v)
    {
        std::vector<std::shared_lock<std::shared_mutex>> locks(mutexes, mutexes + capacity);  // 获取所有哈希表的读锁
        Node *curr, *pre;
        v.reserve(size.load());
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
        std::sort(std::execution::par, v.begin(), v.end());
        return true;
}

private:
    int hash(size_t value) {
        return value % capacity;
    }
};


size_t compress(std::string &s, int sz)
{
    std::vector <int> v(sz, 0);
    int count = 1;
    size_t ans = 0;
    for(auto &i: s)
    {
        if(! v[i - '0'])
        {
            v[i - '0'] = count ++;
        }
        ans = ans * 10 + v[i - '0'];
    }
    return ans;
}

std::string to_string_X(size_t val, int X)
{
    int now;
    std::string s;
    while(val)
    {
        now = int(val % X);
        s += (now + '0');
        val /= X;
    }
    return s;
}

void block_process(ConcurrentSet &s, size_t begin_, size_t end_, size_t max_len, size_t one_size)
{
    std::string tmp;
    int real_one_size = int(one_size);
    for(size_t i = begin_; i < end_; ++ i)
    {
        tmp = to_string_X(i, real_one_size);
        if(tmp.size() <= max_len) tmp.insert(0, (max_len - tmp.size()), '0');
        else return;
        s.insert(compress(tmp, 4));
    }
}

inline char get_bit(size_t raw_id, size_t bit)
{
    while(bit)
    {
        raw_id /= 10;
        -- bit;
    }
    return raw_id % 10 + '0';
}

int main() {
    ConcurrentSet s;
    std::vector <std::thread> threads;
    const size_t str_len = 6, one_size = 4, thread_id = 12;
    size_t all_size = 1;
    for(int i = 1; i <= str_len; ++ i) all_size *= one_size;
    for(size_t i = 0, begin_, end_; i < thread_id; ++ i) 
    {
        begin_ = all_size / thread_id * i;
        end_ = all_size / thread_id * (i + 1);
        threads.emplace_back(block_process, std::ref(s), begin_, end_, str_len, one_size);
    }
    for(auto &entry: threads) entry.join();
    std::vector <size_t> v;
    s.getFinalSnapshotAndFree(v);
    for(auto &item: v)
    {
        std::cout << item << '\n';
    }
    std::cout << std::flush;
    return 0;
}
