#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <queue>
#include <cstdint>
#include <string>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <random>
#include <map>
#include <vector>
#include <algorithm>

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
    for (unsigned long i = 0; i < x.size(); i++) {
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

        while (this->root && !q.empty()) {
            node* cur = q.front();
            q.pop();
            int _dist = hamming_distance_bit(cur->s, s);
            if (_dist < _min_dist) {
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
                if (cur->has_child[i] > _dist - _sec_dist && cur->has_child[i] < _dist + _sec_dist) {
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
        while (this->root && !q.empty()) {
            node* cur = q.front();
            q.pop();
            int _dist = hamming_distance_bit(cur->s, s);
            if (_dist < _min_dist) {
                _min_dist = _dist;
                minNode = cur;
            }
            for (int i = 0; i < cur->n_child; i++) {
                if (cur->has_child[i] > _dist - _min_dist && cur->has_child[i] < _dist + _min_dist) {
                    q.push(cur->child[cur->has_child[i]]);
                }
            }
        }
        std::string s1 = decode(minNode == nullptr ? 0xFFFFFFFFFFFFFFFFull : minNode->s, this->max_dist);
        return std::make_pair(_min_dist, s1);
    }
};

int minhash(std::string x, std::map<std::pair<int, char>, int> &permutation) {
  int hash = INT_MAX;
  for (unsigned long i = 0; i < x.size(); i++) {
    hash = std::min(hash, permutation[std::make_pair(i, x[i])]);
  }
  return hash;
}

class PBKTree {
  const int num_tree;
  const int num_hash;
  const int max_dist;

  std::map<std::pair<int, char>, int> **permutations;
  std::map<std::vector<int>, BKTree*> *trees;

public:
  PBKTree(int num_tree, int num_hash, int max_dist)
      : num_tree(num_tree), num_hash(num_hash), max_dist(max_dist) {
    std::vector<std::pair<int, char>> items;
    for (int i = 0; i < max_dist; i++) {
      items.push_back(std::make_pair(i, 'A'));
      items.push_back(std::make_pair(i, 'T'));
      items.push_back(std::make_pair(i, 'C'));
      items.push_back(std::make_pair(i, 'G'));
    }

    permutations = new std::map<std::pair<int, char>, int> *[num_tree];
    for (int i = 0; i < num_tree; i++) {
      permutations[i] = new std::map<std::pair<int, char>, int>[num_hash];
      for (int j = 0; j < num_hash; j++) {
        std::random_shuffle(items.begin(), items.end());
        for (unsigned long k = 0; k < items.size(); k++) {
          permutations[i][j][items[k]] = k;
        }
      }
    }

    trees = new std::map<std::vector<int>, BKTree*>[num_tree];
  }

  ~PBKTree() {
    for (int i = 0; i < num_tree; i++) {
      delete[] permutations[i];
    }
    delete[] permutations;

                for (int i = 0; i < num_tree; i++) {
                        for (auto [hashes, ptr] : trees[i]) {
                                delete ptr;
                        }
                }
    delete[] trees;
  }

  void insert(std::string str) {
                std::vector<int> hashes[num_tree];
    for (int i = 0; i < num_tree; i++) {
      for (int j = 0; j < num_hash; j++) {
        hashes[i].push_back(minhash(str, permutations[i][j]));
      }
    }

    for (int i = 0; i < num_tree; i++) {
                        if (trees[i][hashes[i]] == nullptr) {
                                trees[i][hashes[i]] = new BKTree(max_dist);
                        }
      trees[i][hashes[i]]->insert(str);
    }
  }

  std::pair<int, std::string> find_nearest(std::string str,
                                           int required_dist = INT_MAX) {
                std::vector<int> hashes[num_tree];
    for (int i = 0; i < num_tree; i++) {
      for (int j = 0; j < num_hash; j++) {
        hashes[i].push_back(minhash(str, permutations[i][j]));
      }
    }

    int min_dist = INT_MAX;
    std::string min_str;

    for (int i = 0; i < num_tree; i++) {
                        if (trees[i][hashes[i]] == nullptr) {
                                trees[i][hashes[i]] = new BKTree(max_dist);
                        }
      auto [_min_dist, _min_str] =
          trees[i][hashes[i]]->find_nearest(str, required_dist);
      if (_min_dist < min_dist) {
        min_dist = _min_dist;
        min_str = _min_str;
      }
    }
    return std::make_pair(min_dist, min_str);
  }

  std::pair<std::pair<int, std::string>, std::pair<int, std::string>>
  find2(std::string str, int required_dist = INT_MAX) {
                std::vector<int> hashes[num_tree];
    for (int i = 0; i < num_tree; i++) {
      for (int j = 0; j < num_hash; j++) {
        hashes[i].push_back(minhash(str, permutations[i][j]));
      }
    }

    int min_dist = INT_MAX;
    int sec_dist = INT_MAX;
    std::string min_str;
    std::string sec_str;

    for (int i = 0; i < num_tree; i++) {
                        if (trees[i][hashes[i]] == nullptr) {
                                trees[i][hashes[i]] = new BKTree(max_dist);
                        }
      auto [_min_pair, _sec_pair] = trees[i][hashes[i]]->find2(str, required_dist);
      auto [_min_dist, _min_str] = _min_pair;
      auto [_sec_dist, _sec_str] = _sec_pair;

      if (_min_dist < min_dist) {
        sec_dist = min_dist;
        sec_str = min_str;
        min_dist = _min_dist;
        min_str = _min_str;
      } else if (_min_dist < sec_dist && _min_str != min_str) {
        sec_dist = _min_dist;
        sec_str = _min_str;
      }

      if (_sec_dist < sec_dist) {
        sec_dist = _sec_dist;
        sec_str = _sec_str;
      }
    }
    return std::make_pair(std::make_pair(min_dist, min_str),
                          std::make_pair(sec_dist, sec_str));
  }
};

namespace py = pybind11;

PYBIND11_MODULE(bktree, m) {
    py::class_<BKTree>(m, "BKTree")
        .def(py::init<int>())
        .def("insert", &BKTree::insert)
        .def("find2", &BKTree::find2)
        .def("find_nearest", &BKTree::find_nearest);

  py::class_<PBKTree>(m, "PBKTree")
      .def(py::init<int, int, int>())
      .def("insert", &PBKTree::insert)
      .def("find2", &PBKTree::find2)
      .def("find_nearest", &PBKTree::find_nearest);
}
