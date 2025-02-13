//
// Created by AnthonyZhang on 2025/1/7.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_GRID_H
#define ATOM_SEARCH_CPP_GRID_H

#include "hash.h"
#include "atom.h"
#include "graph.h"
#include <omp.h>
#include <Eigen/Dense>
#include <unordered_map>
#include <limits>
#include <memory>
#include <shared_mutex>
#include <mutex>

namespace grid {
    constexpr float inf = std::numeric_limits<float>::max();
}

struct GridCell {
    std::vector<std::shared_ptr<Atom>> atoms;
//    std::shared_mutex cell_rw_lock;
};

// A structure based on hashmap to maintain the positions of points for collision detection
class Grid {
private:
    float _interval;
    int _size;
    std::unordered_map<std::array<int, 3>, GridCell, hashing::array_hash<int, 3>> mp;
//    std::shared_mutex global_rw_lock;
public:
    explicit Grid(float interval = 5.0);
    virtual ~Grid();
    [[nodiscard]] int _getpos(float x) const;
    [[nodiscard]] float interval() const;
    [[nodiscard]] int size() const;
    std::pair<std::shared_ptr<Atom>, float> search_nn(const Position &point);
    void add(const Position &point, const std::string &symbol, bool including_hydrogen = true);
    void add(const std::shared_ptr<Atom>& atom_ptr, bool including_hydrogen = true);
    void add_mol(const Graph& g, bool including_hydrogen = true);
    bool isCollision(const Position& point, float threshold);
};

#endif //ATOM_SEARCH_CPP_GRID_H