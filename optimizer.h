//
// Created by AnthonyZhang on 2025/1/12.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_OPTIMIZER_H
#define ATOM_SEARCH_CPP_OPTIMIZER_H

#include "de/DifferentialEvolution.h"
#include "grid.h"
#include "atom.h"
#include <utility>
#include <ctime>

struct Pointer {
    int left, right;
    explicit Pointer(int left = 0, int right = 0);
};

class Optimizer {
private:
    float LJ_weight;
    std::shared_ptr<Grid> tree;
    const std::vector<Position>& target_points;
    const Pointer& pointer;
    int cal_len;

public:
    std::vector<Position> atoms_list;
    Position root_position{};
    Vector K{};

    Optimizer(float LJ_weight, std::shared_ptr<Grid> tree, const std::vector<Position>& target_points, const Pointer& pointer, int cal_len);

    std::pair<float, float> objective_fcn_pair(float angle);
    float objective_fcn(float angle);
};

class TargetFunction : public de::IOptimizable {
private:
    std::shared_ptr<Optimizer> opti_ptr;
public:
    explicit TargetFunction(std::shared_ptr<Optimizer> ptr);
    ~TargetFunction() override = default;
    [[nodiscard]] double EvaluateCost(const Eigen::VectorXd& inputs) const override;
    [[nodiscard]] unsigned int NumberOfParameters() const override;
    [[nodiscard]] std::vector<de::IOptimizable::Constraints> GetConstraints() const override;
};

class TestFunction : public de::IOptimizable {
public:
    TestFunction() = default;
    ~TestFunction() override = default;
    [[nodiscard]] double EvaluateCost(const Eigen::VectorXd& inputs) const override;
    [[nodiscard]] unsigned int NumberOfParameters() const override;
    [[nodiscard]] std::vector<de::IOptimizable::Constraints> GetConstraints() const override;
};

namespace optimizer {
    constexpr int inf = 1000000000;
    constexpr double costThreshold = 0.01;
    constexpr int convergenceThreshold = 20;
    constexpr double convergenceEpsilon = 0.01;
    float optimize(std::shared_ptr<Optimizer>& opti_opt, std::vector<Position> atom_list, Position root_position_, Vector K_, bool verbose = false);
    float testOptimize(bool verbose = true);
    bool EarlyTerminationFunction(const de::DifferentialEvolution& diffevo);
}
#endif //ATOM_SEARCH_CPP_OPTIMIZER_H
