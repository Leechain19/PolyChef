//
// Created by AnthonyZhang on 2025/1/18.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_CURVE_H
#define ATOM_SEARCH_CPP_CURVE_H

#include "atom.h"
#include "Eigen/Dense"
#include "grid.h"
#include <memory>
#include <utility>
#include <vector>
#include <queue>
#include <cmath>
#include <omp.h>
#include <Eigen/Dense>
#include <functional>
#include <string>

struct Node {
    Position pos;
    std::shared_ptr<Node> fa;
    Node(Position pos, std::shared_ptr<Node> fa);
    virtual ~Node() = default;
    [[nodiscard]] std::string getPositionToString() const ;
};

struct AStarNode : public Node {
    float cost, h;
    AStarNode(Position pos, float cost, float h, std::shared_ptr<Node> fa = nullptr);
    [[nodiscard]] float f_cost() const ;
    ~AStarNode() override = default;
};

namespace AStar {
    float heuristic(const Position& start, const Position &goal);

    std::vector<Position> getNeighbors(const Position &position, float step);

    std::vector<Position> AStarSearch(const Position& start, const Position& goal, const std::shared_ptr<Grid> &tree, float step, float collisionThreshold = 3.5f);
}

namespace RRT {

}

#endif //ATOM_SEARCH_CPP_CURVE_H
