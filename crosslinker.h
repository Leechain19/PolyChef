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

class CrossLinker {
public:
    CrossLinker() = default;
    CrossLinker(int n, std::vector<std::shared_ptr<Atom>> atoms, std::vector<int> is_ar);
    ~CrossLinker() = default;

    CrossLinker(const CrossLinker& other);
    CrossLinker(CrossLinker&& other) noexcept ;

    CrossLinker& operator=(const CrossLinker &other);
    CrossLinker& operator=(CrossLinker &&other) noexcept ;

    [[nodiscard]] int size() const ;

    [[nodiscard]] int getPolysSize() const;

    [[nodiscard]] int getEdgesSize() const ;

    [[nodiscard]] std::shared_ptr<Atom> getAtom(int index) const ;

    [[nodiscard]] Position getAtomPosition(int index) const ;

    [[nodiscard]] const std::vector<std::shared_ptr<Atom>>& getAtomVec() const;

    [[nodiscard]] std::shared_ptr<Poly> getPoly(int index) const ;

    [[nodiscard]] Position getPolyPosition(int index) const ;

    [[nodiscard]] const std::vector<std::shared_ptr<Poly>>& getPolyVec() const;

    void addAtom(std::shared_ptr<Atom> atom_ptr, int mono, bool ar = false);

    void addAtom(int index, std::shared_ptr<Atom> atom_ptr, int mono, bool ar = false);

    void addEdge(int from, int to, const std::string& type);

    [[nodiscard]] bool isAr(int index) const;

    [[nodiscard]] const std::vector<int>& getArVec() const;

    [[nodiscard]] const std::vector<std::shared_ptr<Edge>>& getEdge(int index) const;

    [[nodiscard]] const std::vector<std::vector<std::shared_ptr<Edge>>>& getEdgeVec() const;

    friend std::ostream& operator<<(std::ostream& os, const CrossLinker& cl);

private:
    int n{};
    std::vector<std::shared_ptr<Atom>> vertices;
    std::vector<std::vector<std::shared_ptr<Edge>>> edges;
    std::vector<std::shared_ptr<Poly>> polys;
    std::vector<int> is_ar;
};


class CrosslinkingSystem {
public:
    CrosslinkingSystem(std::vector<std::shared_ptr<CrossLinker>> crosslink_graph, std::vector<std::array<int, 4>> cross_network);
    CrosslinkingSystem(const CrosslinkingSystem& other) = delete;
    CrosslinkingSystem& operator=(const CrosslinkingSystem& other) = delete;

    ~CrosslinkingSystem() = default;

    [[nodiscard]] int getCrosslinkSize() const;

    [[nodiscard]] std::shared_ptr<CrossLinker> getCrosslinkGraph(int index) const;

    [[nodiscard]] const std::vector<std::shared_ptr<CrossLinker>> &getCrosslinkGraphVec() const;

    [[nodiscard]] std::shared_ptr<Graph> getChainGraph(int index) const;

    [[nodiscard]] const std::vector<std::shared_ptr<Graph>> &getChainGraphVec() const;

    [[nodiscard]] const std::array<int, 4>& getEdge(int index) const ;

    [[nodiscard]] const std::vector<std::array<int, 4>>& getEdgeVec() const ;

    void addEdge(const std::array<int, 4>& e);

    bool calcChainGraphs(const std::array<int, 4>& e, const CustomFunctionLoader& custfunc, const std::vector<std::shared_ptr<Graph>>& sequence,
                         int degree_polymerization, bool random_polymerization = false, int optimize_size = 1);



private:
    int n_;
    std::vector<std::shared_ptr<CrossLinker>> crosslink_graphs_;
    std::vector<std::array<int, 4>> edges_;
    std::shared_ptr<Grid> tree_ptr_;
    std::vector<std::shared_ptr<Graph>> chain_graphs_;
    std::vector<std::vector<int>> poly_seen_;
};

#endif //ATOM_SEARCH_CPP_CROSSLINKER_H
