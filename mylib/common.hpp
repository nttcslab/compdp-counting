#ifndef SIMFRONTIER_COMMON_HPP
#define SIMFRONTIER_COMMON_HPP

#include <cmath>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <utility>

using addr_t = long long int;
using Bint = long long int;

// FNV hash constant
constexpr uint64_t FNV_OFFSET_BASIS_64 = 14695981039346656037ULL;
constexpr uint64_t FNV_PRIME_64 = 1099511628211ULL;

// hash function of pair of ints
class HashPI{
public:
  size_t operator()(const std::pair<int, int>& x) const {
    uint64_t h = FNV_OFFSET_BASIS_64;
    h = FNV_PRIME_64* h ^ x.first;
    h = FNV_PRIME_64* h ^ x.second;
    return h;
  }
};

class HashPaddr{
public:
  size_t operator()(const std::pair<addr_t, addr_t>& x) const {
    uint64_t h = FNV_OFFSET_BASIS_64;
    h = FNV_PRIME_64* h ^ x.first;
    h = FNV_PRIME_64* h ^ x.second;
    return h;
  }
};

#endif // SIMFRONTIER_COMMON_HPP