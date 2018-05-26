
#ifndef LS_UTIL_H
#define LS_UTIL_H

#include <string>
#include <sstream>
#include <bitset>

std::string to_string(char function, uint64_t i, uint64_t m) {
    std::stringstream ss;
    ss << function << '_' << i;
    ss << '(';
    if (m > 0) ss << 'a' << '_' << 1;
    if (m > 2) ss << ',' << " ... ";
    if (m > 1) ss << ", " << 'a' << '_' << m;
    ss << ')';
    return ss.str();
}

std::string to_string(char family, char function, uint64_t n, uint64_t m) {
    std::stringstream ss;
    ss << family << " = ";
    ss << '{';
    if (n > 0) ss << to_string(function, 1, m);
    if (n > 2) ss << ',' << " ... ";
    if (n > 1) ss << ", " << to_string(function, n, m);
    ss << '}';
    return ss.str();
}

inline uint64_t extract(uint64_t family, uint64_t n, uint64_t k, uint64_t i) {
    return (family >> k * (n - i - 1)) & ((1ul << k) - 1);
}

inline bool _is_const(uint64_t f, uint64_t k) {
    return f == 0 || f == (1ul << k) - 1;
}

inline bool _is_equal(uint64_t f, uint64_t k, uint64_t a, uint64_t b) {
    return !((f >> (k - 1 - a) & 1ul) ^ (f >> (k - 1 - b) & 1ul));
}

inline bool _is_dummy(uint64_t f, uint64_t n, uint64_t k, uint64_t i) {
    const uint64_t lc = i;
    const uint64_t rc = n - i - 1;
    const uint64_t lm = 1ul << lc;
    const uint64_t rm = 1ul << rc;
    if (i == 0) {
        for (uint64_t r = 0; r < rm; r++) {
            uint64_t a = 0u << rc | r;
            uint64_t b = 1u << rc | r;
            if (!_is_equal(f, k, a, b)) return false;
        }
        return true;
    }
    if (i == k - 1) {
        for (uint64_t l = 0; l < lm; l++) {
            uint64_t a = l << 1u | 0u;
            uint64_t b = l << 1u | 1u;
            if (!_is_equal(f, k, a, b)) return false;
        }
        return true;
    }
    for (uint64_t l = 0; l < lm; l++) {
        for (uint64_t r = 0; r < rm; r++) {
            uint64_t a = (l << 1u | 0u) << rc | r;
            uint64_t b = (l << 1u | 1u) << rc | r;
            if (!_is_equal(f, k, a, b)) return false;
        }
    }
    return true;
}

inline void to_vector(uint64_t family, uint64_t k, uint64_t n, std::vector <uint64_t> &v) {
    v.resize(n);
    for (auto i = 0; i < n; i++) {
        v[i] = extract(family, n, k, i);
    }
}

inline std::vector <uint64_t> to_vector(uint64_t family, uint64_t k, uint64_t n) {
    std::vector <uint64_t> v(n);
    to_vector(family, k, n, v);
    return v;
}

inline uint64_t ith(uint64_t x, uint64_t n, uint64_t i) {
    return (x >> (n - 1 - i)) & 1ul;
}

#endif //LS_UTIL_H
