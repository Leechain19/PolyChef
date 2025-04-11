//
// Created by AnthonyZhang on 2025/1/18.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_CURVE_H
#define ATOM_SEARCH_CPP_CURVE_H

#include "atom.h"
#include "Eigen/Dense"
#include "grid.h"
#include "loadlib.h"
#include "kdtree.hpp"
#include "progresscpp/ProgressBar.hpp"
#include <memory>
#include <utility>
#include <vector>
#include <queue>
#include <cmath>
#include <omp.h>
#include <Eigen/Dense>
#include <functional>
#include <string>

std::vector<Position> calCRS(const std::vector<Position>& path, float distance_threshold = 0.8f);
std::vector<Position> getScatterFromCSV(const std::string& path, bool header = true, float distance_upper_bound = 0.3f);

using myKDTree = kdtree::KDTree<int, 3>;

struct Node {
    Position pos;
    std::shared_ptr<Node> fa;
    Node(Position pos, std::shared_ptr<Node> fa);
    virtual ~Node() = default;

    [[maybe_unused]] [[nodiscard]] std::string getPositionToString() const ;
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

// RRT 命名空间
namespace RRT {
    // 双向 RRT 规划
    std::vector<Position> bidirectionalRRT(const Position &start, const Position &goal, const std::shared_ptr<Grid>& grid_tree, const CustomFunctionLoader& custfunc,
                                           float step_size, float goal_tolerance, float to_end_possibility, float collisionThreshold = 3.5, int max_trial = 10000);

}


#endif //ATOM_SEARCH_CPP_CURVE_H
