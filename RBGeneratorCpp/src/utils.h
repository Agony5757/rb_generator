#pragma once
#include <complex>
#include <vector>
#include <iostream>
#include <array>
#include <set>
#include <iomanip>
#include <tuple>
#include <map>
#include <random>
#include <chrono>
#include <omp.h>
#include <sstream>
#include <algorithm>

constexpr double SQRT2 = 1.0 / 1.414213562373095145474621858739;

namespace std {
    inline bool operator<(const ::std::complex<double>& x, const ::std::complex<double>& y)
    {
        return std::forward_as_tuple(x.real(), x.imag()) < std::forward_as_tuple(y.real(), y.imag());
    }
}

constexpr double EPSILON = 1e-8; 

template<typename T>
inline bool is_close(T a, T b) {
    return std::abs(a - b) < EPSILON;
}

template<typename T>
inline bool is_close_to_zero(T a) {
    return std::abs(a) < EPSILON;
}

#define MAKE_POSSIBLE_VALUES(m, value) if (is_close(m.real(), value)) m.real(value); if (is_close(m.imag(), value)) m.imag(value)
#define MAKE_POSSIBLE_VALUES_ZERO(m) if (is_close_to_zero(m.real())) m.real(0); if (is_close_to_zero(m.imag())) m.imag(0)

inline void clarify(std::complex<double>& m)
{
    MAKE_POSSIBLE_VALUES_ZERO(m);
    MAKE_POSSIBLE_VALUES(m, SQRT2);
    MAKE_POSSIBLE_VALUES(m, -SQRT2);
    MAKE_POSSIBLE_VALUES(m, 0.5);
    MAKE_POSSIBLE_VALUES(m, -0.5);
    MAKE_POSSIBLE_VALUES(m, 1.0);
    MAKE_POSSIBLE_VALUES(m, -1.0);
}

template<typename T>
void print_group(const std::map<T, std::string>& groupset)
{
    for (auto&& [elem, generate_str] : groupset)
    {
        elem.print();
        std::cout << "(generator: " << generate_str << ")\n";
    }
}

template<typename T>
std::string vec2str(const std::vector<T>& vec)
{
    std::stringstream ss;
    ss << "[";
    for (int i=0;i<vec.size();++i)
    {
        ss << vec[i];
        if (i != vec.size() - 1)
            ss << ", ";
    }
    ss << "]";
    return ss.str();
}