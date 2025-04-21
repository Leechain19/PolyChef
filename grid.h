//
// Created by AnthonyZhang on 2025/1/7.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_GRID_H
#define ATOM_SEARCH_CPP_GRID_H

#include "hash.h"
#include "atom.h"
#include "utility.h"
#include "exception.h"
#include <omp.h>
#include <Eigen/Dense>
#include <unordered_map>
#include <limits>
#include <memory>
#include <shared_mutex>
#include <mutex>
#include <fstream>

namespace grid {
    constexpr float inf = std::numeric_limits<float>::max();
}

struct GridCell {
    std::vector<std::shared_ptr<Atom>> atoms;
};

// A structure based on hashmap to maintain the positions of points for collision detection
class Grid {
private:
    float _interval;
    int _size;
    std::unordered_map<std::array<int, 3>, GridCell, hashing::array_hash<int, 3>> mp;
public:
    explicit Grid(float interval = 4.0);
    virtual ~Grid();
    [[nodiscard]] int _getpos(float x) const;
    [[nodiscard]] float interval() const;
    [[nodiscard]] int size() const;
    std::pair<std::shared_ptr<Atom>, float> search_nn(const Position &point, float collision_threshold = 2.5);
    void add(const Position &point, const std::string &symbol, bool including_hydrogen = true);
    void add(std::shared_ptr<Atom> atom_ptr, bool including_hydrogen = true);
    void erase(const std::shared_ptr<Atom>& atom_ptr);

    template<typename T>
    void add_mol(const T& g, bool including_hydrogen = true);

    template<typename T>
    void add_mol(const std::shared_ptr<T>& g_ptr, bool including_hydrogen = true);

    template<typename T>
    void erase_mol(const T& g);

    template<typename T>
    void erase_mol(const std::shared_ptr<T>& g_ptr);

    bool isCollision(const Position& point, float threshold);

    void addCollisionMol2(const std::string& s);

    void addCollisionMol2(const std::vector<std::string>& vec);
};


template<typename T>
void Grid::add_mol(const T& g, bool including_hydrogen) {
    for (const auto& ptr : g.getAtomVec()) {
        add(ptr, including_hydrogen);
    }
}

template<typename T>
void Grid::add_mol(const std::shared_ptr<T> &g_ptr, bool including_hydrogen) {
    for (const auto& ptr : g_ptr->getAtomVec()) {
        add(ptr, including_hydrogen);
    }
}

template<typename T>
void Grid::erase_mol(const T& g) {
    for (auto it = g.getAtomVec().rbegin(); it != g.getAtomVec().rend(); it ++) {
        erase(*it);
    }
}

template<typename T>
void Grid::erase_mol(const std::shared_ptr<T> &g_ptr) {
    for (auto it = g_ptr->getAtomVec().rbegin(); it != g_ptr->getAtomVec().rend(); it ++) {
        erase(*it);
    }
}


#endif //ATOM_SEARCH_CPP_GRID_H