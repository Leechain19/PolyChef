//
// Created by AnthonyZhang on 2025/1/11.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_HASH_H
#define ATOM_SEARCH_CPP_HASH_H

#include <array>
#include <vector>
#include <Eigen/Dense>

namespace hashing {
    template <typename T, std::size_t N>
    struct array_hash {
        std::size_t operator()(const std::array<T, N> &arr) const {
            std::size_t seed = 0;
            for (const auto& elem : arr) {
                seed ^= std::hash<T>{}(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    template<typename T>
    struct [[maybe_unused]] vector_hash {
        std::size_t operator()(const std::vector<T> &arr) const {
            std::size_t seed = 0;
            for (const auto& elem : arr) {
                seed ^= std::hash<T>{}(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    template<typename T>
    struct pair_hash {
        std::size_t operator()(const std::pair<T, T>& pr) const {
            std::size_t seed = 0;
            seed ^= (std::hash<T>{}(pr.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
            seed ^= (std::hash<T>{}(pr.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
            return seed;
        }
    };

    template <typename T>
    struct eigen_vector_hash {
        std::size_t operator()(const Eigen::Vector<T, 3>& vec) const {
            std::size_t h1 = std::hash<T>()(vec(0));
            std::size_t h2 = std::hash<T>()(vec(1));
            std::size_t h3 = std::hash<T>()(vec(2));
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };

    template <typename T>
    struct eigen_vector_equal {
        bool operator()(const Eigen::Vector<T, 3>& a, const Eigen::Vector<T, 3>& b) const {
            return a.isApprox(b, 1e-3);
        }
    };
}

#endif //ATOM_SEARCH_CPP_HASH_H



