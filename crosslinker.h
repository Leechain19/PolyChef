//
// Created by AnthonyZhang on 2025/2/17.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_CROSSLINKER_H
#define ATOM_SEARCH_CPP_CROSSLINKER_H

#include <vector>
#include <memory>
#include "atom.h"
#include "exception.h"
#include "graph.h"
#include "curve.h"
#include "spreading.h"

class CrossLinker {
public:
    CrossLinker() = default;
    CrossLinker(int n, std::vector<std::shared_ptr<Atom>> atoms, std::vector<int> is_ar, int mono_type);
    ~CrossLinker() = default;

    CrossLinker(const CrossLinker& other);
    CrossLinker(CrossLinker&& other) noexcept ;

    CrossLinker& operator=(const CrossLinker &other);
    CrossLinker& operator=(CrossLinker &&other) noexcept ;

    [[nodiscard]] int size() const ;

    [[nodiscard]] int getPolysSize() const;

    [[nodiscard]] int getEdgesSum() const ;

    [[nodiscard]] std::shared_ptr<Atom> getAtom(int index) const ;

    [[nodiscard]] Position getAtomPosition(int index) const ;

    [[nodiscard]] const std::vector<std::shared_ptr<Atom>>& getAtomVec() const;

    [[nodiscard]] std::shared_ptr<Poly> getPoly(int index) const ;

    [[nodiscard]] Position getPolyPosition(int index) const ;

    [[nodiscard]] const std::vector<std::shared_ptr<Poly>>& getPolyVec() const;

    void addAtom(std::shared_ptr<Atom> atom_ptr, bool ar = false);

    void addAtom(int index, std::shared_ptr<Atom> atom_ptr, bool ar = false);

    void setMonoType(int new_mono_type);

    [[nodiscard]] int getMonomerType() const ;

    void addEdge(int from, int to, Bond_type type);

    void addPoly(float x, float y, float z, int neigh);

    [[nodiscard]] bool isAr(int index) const;

    [[nodiscard]] const std::vector<int>& getArVec() const;

    [[nodiscard]] const std::vector<std::shared_ptr<Edge>>& getEdge(int index) const;

    [[nodiscard]] const std::vector<std::vector<std::shared_ptr<Edge>>>& getEdgeVec() const;

    void makeEnd(int poly_index, const std::string& end_symbol = "H");

    void makePolyUsedFlag(int index);

    bool checkPolyUsedFlag(int index);

    friend std::ostream& operator<<(std::ostream& os, const CrossLinker& cl);

private:
    int n{};
    std::vector<std::shared_ptr<Atom>> vertices;
    std::vector<std::vector<std::shared_ptr<Edge>>> edges;
    std::vector<std::shared_ptr<Poly>> polys;
    std::vector<int> poly_seen;
    std::vector<int> is_ar;
    int mono_type{};
};

class CrosslinkingSystem {
public:
    CrosslinkingSystem() = delete;
    CrosslinkingSystem(std::vector<std::shared_ptr<CrossLinker>> crosslinkers, std::vector<std::array<int, 4>> crosslinker_network,
                       std::vector<std::vector<Position>> point_lists, std::shared_ptr<Grid> tree, float window_distance);

    CrosslinkingSystem(const CrosslinkingSystem& other) = delete;
    CrosslinkingSystem& operator=(const CrosslinkingSystem& other) = delete;

    ~CrosslinkingSystem() = default;

    [[nodiscard]] int getAtomSize() const ;

    [[nodiscard]] int getEdgeSize() const ;

    [[nodiscard]] int getCrosslinkerNumber() const;

    [[nodiscard]] int getCrosslinkerNetworkNumber() const ;

    [[nodiscard]] std::shared_ptr<CrossLinker> getCrosslinkGraph(int index) const;

    [[nodiscard]] const std::vector<std::shared_ptr<CrossLinker>> &getCrosslinkGraphVec() const;

    [[nodiscard]] std::shared_ptr<Graph> getChainGraph(int index) const;

    [[nodiscard]] const std::vector<std::shared_ptr<Graph>> &getChainGraphVec() const;

    [[nodiscard]] const std::array<int, 4>& getCrosslinkerNetworkInIdx(int index) const ;

    [[nodiscard]] const std::vector<std::array<int, 4>>& getCrosslinkerNetwork() const ;

    void spreadingChain(int chain_index, const std::vector<std::shared_ptr<Graph>> &sequence,
                        int degree_polymerization, const std::string& pool_choice, float para_A,
                        float para_B, bool random_polymerization, int optimize_size);

    void calcChainGraphs(const std::vector<std::vector<std::shared_ptr<Graph>>>& sequences, int degree_polymerization,
                         const std::string& pool_choice, float para_A, float para_B, bool random_polymerization, int optimize_size);

    void makeEnd(const std::string& end_system);

    friend std::ostream& operator<<(std::ostream& os, const CrosslinkingSystem& cls);

private:
    std::vector<std::shared_ptr<CrossLinker>> crosslinkers_;
    std::vector<std::array<int, 4>> crosslinker_network_;
    std::vector<std::vector<Position>> point_lists_;
    std::shared_ptr<Grid> tree_ptr_;
    std::vector<std::shared_ptr<Graph>> chain_graphs_;
    float window_distance_ = 5.0f;
};

#endif //ATOM_SEARCH_CPP_CROSSLINKER_H
