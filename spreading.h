//
// Created by AnthonyZhang on 2025/1/18.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_SPREADING_H
#define ATOM_SEARCH_CPP_SPREADING_H

#include "graph.h"
#include "optimizer.h"
#include "grid.h"
#include "mathfunc.h"
#include "chemtable.h"
#include "progresscpp/ProgressBar.hpp"
#include <queue>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <chrono>
#include <utility>

enum class Spreading_stat {
    EARLY_STOP,
    FULL_STOP,
    BAD_SIGNAL
};

class MolGenerator {
private:
    std::vector<std::shared_ptr<Graph>> sequence;
    bool is_random;
    std::mt19937 gen;
    std::uniform_int_distribution<int> dist;
protected:
    int cur;
    int len;
public:
    MolGenerator(std::vector<std::shared_ptr<Graph>> sequence, bool is_random, unsigned int seed = 0x9A3F7B1C);
    virtual ~MolGenerator() = default;
    [[nodiscard]] std::shared_ptr<Graph> getSequence(int index) const;
    virtual std::shared_ptr<Graph> getNext();
};


class FastaGenerator : public MolGenerator {
private:
    const std::string fasta;
public:
    [[maybe_unused]] FastaGenerator(std::vector<std::shared_ptr<Graph>> sequence, std::string fasta);
    ~FastaGenerator() override = default;
    std::shared_ptr<Graph> getNext() override;
};

void curveSpreading(const std::vector<Position>& target_points, std::shared_ptr<Graph> g, std::shared_ptr<Grid> tree, const std::vector<std::shared_ptr<Graph>>& sequence,
                    int degree_of_polymerization, const std::string& pool_choice, float para_A, float para_B, float window_distance, int optimize_atom_number, bool random_polymerization, int optimize_size,
                    bool verbose, double bad_cost, const std::unique_ptr<std::vector<std::pair<double, double>>>& loss_vector_ptr = nullptr);



#endif //ATOM_SEARCH_CPP_SPREADING_H
