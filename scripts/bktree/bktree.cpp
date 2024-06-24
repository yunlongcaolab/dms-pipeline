#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <queue>
#include <cstdint>
#include <string>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <random>

const uint64_t mask1 = 0x5555555555555555ull; // 64bit 010101...
const uint64_t mask2 = 0xAAAAAAAAAAAAAAAAull; // 64bit 101010...

inline int hamming_distance_bit(uint64_t x, uint64_t y) {
    uint64_t _xor = x^y;
    return __builtin_popcountll((_xor&mask1)|((_xor&mask2)>>1));
}

const std::string chars = "ATCG";
std::unordered_map<char, uint64_t> char2int = {
    {'A', 0}, {'T', 1}, {'C', 2}, {'G', 3}
};

std::string decode(uint64_t x, int len=26) {
    std::string s = "";
    int highest = len*2-2;

    if (x > (0xFFFFFFFFFFFFFFFFull >> (64-2*len))) {
        return "NotFound";
    }

    for (int i = highest; i >= 0; i -= 2) {
         s += chars[(x & (3ll << i)) >> i];
    }
    return s;
}

uint64_t encode(std::string x) {
    uint64_t res = 0;
    for (int i = 0; i < x.size(); i++) {
        res |= (char2int[x[x.size()-i-1]] << 2*i);
    }
    return res;
}

// A: 00, T: 01, C: 10, G:11

struct node {
    int max_dist;
    int* has_child;
    int n_child;

    node** child;
    uint64_t s;

    node(uint64_t s, int max_dist) {
        this->max_dist = max_dist;
        this->n_child = 0;
        this->has_child = new int[max_dist+1];

        this->child = new node*[max_dist+1];
        for (int i = 0; i <= max_dist; i++) {
            this->child[i] = nullptr;
        }
        this->s = s;
    }
    ~node() {
        for (int i = 0; i <= this->max_dist; i++) {
            delete this->child[i];
        }
    }
};

class BKTree {
    int max_dist;
    node* root;

    public:
    BKTree(int max_dist) {
        this->max_dist = max_dist;
        this->root = nullptr;
    }
    
    void insert(std::string str) {
        uint64_t s = encode(str);
        node* cur = this->root;
        if (cur == nullptr) {
            this->root = new node(s, this->max_dist);
            return;
        }
        int _dist = hamming_distance_bit(s, cur->s);
        while (cur->child[_dist] != nullptr) {
            cur = cur->child[_dist];
            _dist = hamming_distance_bit(cur->s, s);
        }
        cur->has_child[cur->n_child++] = _dist;
        cur->child[_dist] = new node(s, this->max_dist);
    }

    std::pair<std::pair<int, std::string>, std::pair<int, std::string>> find2(std::string str, int required_dist = 114514) {
        std::queue<node*> q;
        uint64_t s = encode(str);
        node* minNode = nullptr;
        node* secNode = nullptr;
        int _min_dist = std::min(required_dist, this->max_dist)+1;
        int _sec_dist = std::min(required_dist, this->max_dist)+1;
        int cnt = 1;
        q.push(this->root);
        
        while (!q.empty()) {
            node* cur = q.front();
            q.pop();
            int _dist = hamming_distance_bit(cur->s, s);
            if (_dist <= _min_dist) {
                _sec_dist = _min_dist;
                _min_dist = _dist;
                secNode = minNode;
                minNode = cur;
            }
            else if (_dist < _sec_dist) {
                _sec_dist = _dist;
                secNode = cur;
            }
            for (int i = 0; i < cur->n_child; i++) {
                if (cur->has_child[i] >= _dist - _sec_dist && cur->has_child[i] <= _dist + _sec_dist) {
                    q.push(cur->child[cur->has_child[i]]);
                    cnt += 1;
                }
            }
        }
        // std::cout << cnt << std::endl;
        
        std::string s1 = decode(minNode == nullptr ? 0xFFFFFFFFFFFFFFFFull : minNode->s, this->max_dist);
        std::string s2 = decode(secNode == nullptr ? 0xFFFFFFFFFFFFFFFFull : secNode->s, this->max_dist);
        
        return std::make_pair(std::make_pair(_min_dist, s1), std::make_pair(_sec_dist, s2));
    }

    std::pair<int, std::string> find_nearest(std::string str, int required_dist = 114514) {
        std::queue<node*> q;
        uint64_t s = encode(str);
        node* minNode = nullptr;
        int _min_dist = std::min(required_dist, this->max_dist)+1;
        q.push(this->root);
        while (!q.empty()) {
            node* cur = q.front();
            q.pop();
            int _dist = hamming_distance_bit(cur->s, s);
            if (_dist < _min_dist) {
                _min_dist = _dist;
                minNode = cur;
            }
            for (int i = 0; i < cur->n_child; i++) {
                if (cur->has_child[i] >= _dist - _min_dist && cur->has_child[i] <= _dist + _min_dist) {
                    q.push(cur->child[cur->has_child[i]]);
                }
            }
        }
        std::string s1 = decode(minNode == nullptr ? 0xFFFFFFFFFFFFFFFFull : minNode->s, this->max_dist);
        return std::make_pair(_min_dist, s1);
    }
};

namespace py=pybind11;

PYBIND11_MODULE(bktree, m) {
    py::class_<BKTree>(m, "BKTree")
        .def(py::init<int>())
        .def("insert", &BKTree::insert)
        .def("find2", &BKTree::find2)
        .def("find_nearest", &BKTree::find_nearest);
}
