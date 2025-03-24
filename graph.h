//
// Created by AnthonyZhang on 2025/1/12.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_GRAPH_H
#define ATOM_SEARCH_CPP_GRAPH_H

#include "atom.h"
#include "exception.h"
#include "hash.h"
#include "mathfunc.h"
#include "chemtable.h"
#include "grid.h"
#include "loadlib.h"
#include "progresscpp/ProgressBar.hpp"
#include <utility>
#include <vector>
#include <numeric>
#include <memory>
#include <unordered_set>
#include <set>
#include <random>
#include <Eigen/Dense>
#include <cassert>
#include <queue>
#include <fstream>


class [[maybe_unused]] DSU {
private:
    int n;
    std::vector<int> f, siz;
public:
    [[maybe_unused]] explicit DSU(int n);
    ~DSU() = default;

    int leader(int x);
    bool same(int x, int y);
    bool merge(int x, int y);
    int size(int x);
};

class PolyManager {
private:
    std::shared_ptr<Poly> front_ptr, back_ptr;
public:
    explicit PolyManager();

    explicit PolyManager(std::shared_ptr<Poly> ptr);

    PolyManager(std::shared_ptr<Poly> ptr1, std::shared_ptr<Poly> ptr2);
    ~PolyManager() = default;

    PolyManager(const PolyManager& other);
    PolyManager(PolyManager&& other) noexcept;

    PolyManager& operator=(const PolyManager& other);
    PolyManager& operator=(PolyManager&& other) noexcept;

    std::shared_ptr<Poly>& operator[](int index);
    std::shared_ptr<Poly>& operator()(int index);

    [[nodiscard]] const std::shared_ptr<Poly>& front() const ;
    [[nodiscard]] const std::shared_ptr<Poly>& back() const ;

    void pop_front();
    void pop_back();
    void pop(int index);

    void push_back(std::shared_ptr<Poly> ptr);
    void push_front(std::shared_ptr<Poly> ptr);

    [[nodiscard]] int size() const;

    void translation(float dx, float dy, float dz);
    void randomize();
    void rotate(const Eigen::Matrix3f& rod, const Position& ver);
};

class Graph {
private:
    int n;
    std::vector<std::shared_ptr<Atom>> vertices;
    std::vector<std::vector<std::shared_ptr<Edge>>> edges;
    std::vector<int> monos;
    PolyManager polys;
    bool is_period{};
    std::unordered_set<std::pair<int, int>, hashing::pair_hash<int>> ring_edges;
    std::vector<int> is_on_main_chain;
    std::vector<int> is_ar;
    std::vector<int> mono_types;
    std::vector<std::vector<int>> mono2ver;

    void _addMono2verTable(int mono_idx, int atom_idx);

    void _modifyMonomer(int mono_idx, int atom_idx);

    static Position _point_rotate_scale_translation(const Position& pos, const Position& ver, const Eigen::Matrix3f& rod, float scale_ratio, const Vector& dvec);

    static Position _point_rotate_scale(const Position& pos, const Position& ver, const Eigen::Matrix3f& rod, float scale_ratio);

    static Position _point_rotate(const Position& pos, const Position& ver, const Eigen::Matrix3f& rod);

    static Position _point_scale(const Position& pos, const Position& ver, float scale_ratio);

    void _bfs_process(int root, int fa, const std::function<Position(const Position&, const Position&)>& f);

    void _make_end(const std::string& end_symbol = "H", int poly_index = -1, bool poly_delete = true);

public:
    Graph();
    Graph(int n, std::vector<std::shared_ptr<Atom>> atoms, std::vector<int> monos, std::vector<int> is_ar, int mono_type = 0);
    virtual ~Graph() = default;

    Graph(const Graph& other);
    Graph(Graph&& other) noexcept;

    Graph& operator=(const Graph& o);
    Graph& operator=(Graph&& other) noexcept ;

    void addAtom(std::shared_ptr<Atom> atom_ptr, int mono, int mono_type = 0, bool ar = false);

    void addAtom(int index, std::shared_ptr<Atom> atom_ptr, int mono, int mono_type = 0, bool ar = false);

    void setMonoType(int index, int mono_type);

    void setMonoTypeAll(int mono_type);

    void addEdge(int from, int to, const std::string& type, bool onring = false);

    bool checkOnRing(int from, int to);

    bool checkOnMainChain(int idx);

    const std::vector<int>& getOnMainChain() const;

    int size() const;

    bool is_periodic() const;

    void addRingEdge(int from, int to);

    std::shared_ptr<Atom> getAtom(int index);

    Position getAtomPosition(int index);

    const std::vector<std::shared_ptr<Atom>>& getAtomVec() const;

    int getMonomer(int index) const;

    int getMonomerType(int index) const;

    const std::vector<int>& getMonomerVec() const;

    const std::vector<int>& getMono2verIndex(int index) const;

    const std::vector<std::vector<int>>& getMono2verTable() const;

    bool isAr(int index) const;

    const std::vector<int>& getArVec() const;

    const std::vector<int>& getMonoTypeVec() const ;

    float atomDistance(int idx1, int idx2);

    const std::vector<std::shared_ptr<Edge>>& getEdge(int index) const;

    int getPolysSize() const;

    int getEdgesSize() const;

    const std::vector<std::vector<std::shared_ptr<Edge>>>& getEdgesVec() const;

    void addPoly(float x, float y, float z, int neigh);

    std::shared_ptr<Poly>& getPoly(int index);

    const std::shared_ptr<Poly>& polyFront() const;

    const std::shared_ptr<Poly>& polyBack() const;

    Position getPolyPosition(int index);

    const std::unordered_set<std::pair<int, int>, hashing::pair_hash<int>>& getRingEdges() const ;

    float polyDistanceWithNeigh(int index);

    void translation(float dx, float dy, float dz);

    void translation(const Vector & vec);

    void polyRandomize();

    void polyPopBack();

    void polyPopFront();

    void rotate(const Eigen::Matrix3f& rod, const Position& ver);

    void startRandomDirection();

    Vector getRotateAxis(int poly_index);

    void calMainChain();

    bool isOnMainChainIdx(int index) const;

    void bfsRotate(int root, int fa, const Eigen::Matrix3f& rod);

    void bfsScale(int root, int fa, float scale_ratio);

    void bfsRotateScale(int root, int fa, const Eigen::Matrix3f& rod, float scale_ratio);

    void bfsRotateScaleTranslation(int root, int fa, const Eigen::Matrix3f& rod, float scale_ratio, Vector dvec);

    void attachPoint(const Position& target_poly_point, const Vector& target_direction);

    void attract(const std::shared_ptr<Graph>& g, int poly_index = -1);

    void connect(const std::shared_ptr<Graph>& g, int poly_index = -1, bool poly_delete = true);

    void makeStart(std::shared_ptr<Graph>& g_ptr);

    void makePeriodicBoundaryCondition(const std::array<float, 3> &box_size);

    void makeEnd(const std::string &symbol);

    friend std::ostream& operator<<(std::ostream &os, const Graph &g);
};

#endif //ATOM_SEARCH_CPP_GRAPH_H
