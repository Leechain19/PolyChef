//
// Created by AnthonyZhang on 2025/1/18.
//

#include "spreading.h"
#include <list>
#include <cmath>
#include <iostream>
#include <chrono>

MolGenerator::MolGenerator(std::vector<std::shared_ptr<Graph>> sequence, bool is_random, unsigned int seed) : sequence(std::move(sequence)), is_random(is_random) {
    this->len = (int)this->sequence.size();
    this->cur = 0;

    if (is_random) {
        if (seed == 0x9A3F7B1C) {
            seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
        }
        gen.seed(seed);
        dist.param(std::uniform_int_distribution<int>::param_type(0, this->len - 1));
    }
}

std::shared_ptr<Graph> MolGenerator::getNext() {
    auto copied_ptr = std::make_shared<Graph>(*this->getSequence(cur));
    if (is_random) {
        cur = dist(gen);
    }
    else {
        cur += 1;
        if (cur >= len) cur -= len;
    }
    return copied_ptr;
}

std::shared_ptr<Graph> MolGenerator::getSequence(int index) const {
    return this->sequence[index];
}

[[maybe_unused]] FastaGenerator::FastaGenerator(std::vector<std::shared_ptr<Graph>> sequence, std::string fasta) : MolGenerator(std::move(sequence), false, 0), fasta(std::move(fasta)) {
    this->len = (int)fasta.size();
}

std::shared_ptr<Graph> FastaGenerator::getNext() {
    auto copied_ptr = std::make_shared<Graph>(*this->getSequence(fasta[cur] - 'A'));
    cur += 1;
    if (cur >= len)
        cur -= len;
    return copied_ptr;
}


void curveSpreading(const std::vector<Position>& target_points, std::shared_ptr<Graph> g, std::shared_ptr<Grid> tree, const std::vector<std::shared_ptr<Graph>>& sequence,
                    int degree_of_polymerization, const std::string& pool_choice, float para_A, float para_B, float window_distance, int optimize_atom_number,
                    bool random_polymerization, int optimize_size, bool verbose, double bad_cost, const std::unique_ptr<std::vector<std::pair<double, double>>>& loss_vector_ptr) {
    assert(optimize_atom_number > 2);
    Pointer pointer(0, 0);

    auto mol_genrator_ptr = std::make_shared<MolGenerator>(sequence, random_polymerization, 114514);
    auto start_mol = mol_genrator_ptr->getNext();

    auto start_target_direction = atom::positionMinusPosition(target_points[0], target_points[1]);
    start_mol->attachPoint(target_points[0], start_target_direction);

    g->makeStart(start_mol);
    std::list<int> bone_line;
    bone_line.push_back(g->polyFront()->getNeigh());

    std::unordered_map<int, bool> seen;
    seen[g->polyFront()->getNeigh()] = true;

    // 标记下一个该加入grid的原子
    std::unordered_set<int> already_add_tree_set;
    auto tree_index = g->polyFront()->getNeigh();

    // 定义Optimize的后序原子数

    auto dfs = [&seen, &g, &bone_line](auto&& self, int u) -> bool {
        // target: -1到poly
        if (seen.count(u)) {
            return seen[u];
        }
        seen[u] = false;
        for (const auto& e : g->getEdge(u)) {
            int atom_index = e->getTo();
            if (g->checkOnMainChain(atom_index) && self(self, atom_index)) {
                seen[u] = true;
                bone_line.push_back(u);
                return true;
            }
        }
        return false;
    };

    auto add_tree = [&g, &tree, &already_add_tree_set](int cur, int aim) -> int {
        while (cur != aim) {
            tree->add(g->getAtom(cur), true);
            already_add_tree_set.insert(cur);
            auto nxt = aim;
            for (const auto& e : g->getEdge(cur)) {
                int x = e->getTo();
                // 判断是不是在主链上
                if (g->checkOnMainChain(x)) {
                    if (already_add_tree_set.count(x)) continue;
                    nxt = x;
                }
                else {
                    tree->add(g->getAtom(x), true);
                }
            }
            cur = nxt;
        }
        return cur;
    };

    auto opt_ptr = std::make_shared<Optimizer>(1.0f, tree, target_points, pointer, optimize_size, pool_choice, para_A, para_B);

    Spreading_stat stat = Spreading_stat::FULL_STOP;

    auto optimize_process = [&bone_line, &g, &pointer, &target_points, &window_distance, &opt_ptr, &add_tree, &verbose,
                             &stat, &bad_cost, &loss_vector_ptr](int tree_index, int size) -> int {
        while ((int)bone_line.size() >= size && stat != Spreading_stat::BAD_SIGNAL) {
            int u = bone_line.front();
            bone_line.pop_front();
            int v = bone_line.front();

            auto position_u = g->getAtomPosition(u);
            auto position_v = g->getAtomPosition(v);

            std::vector<Position> atom_list;
            for (auto x : bone_line) {
                if (x == v) continue;
                atom_list.emplace_back(g->getAtomPosition(x));
            }

            while (pointer.right < (int)target_points.size() && atom::positionDistance(target_points[pointer.right], position_v) < window_distance) {
                pointer.right += 1;
            }

            float nearest_distance = std::numeric_limits<float>::max();
            int nearest_index = -1;

            for (int i = pointer.left; i < pointer.right; i ++) {
                float cur_distance = atom::positionDistance(target_points[i], position_v);
                if (cur_distance < nearest_distance) {
                    nearest_distance = cur_distance;
                    nearest_index = i;
                }
            }

            while (pointer.left < nearest_index) {
                pointer.left += 1;
            }

            if (!g->checkOnRing(u, v)) {
                auto K = atom::positionMinusPosition(position_v, position_u);
                auto theta = optimizer::optimize(opt_ptr, atom_list, position_u, K, verbose);
                auto [val_cost, lj_cost] = opt_ptr->objective_fcn_pair(theta);
                if (val_cost + lj_cost > bad_cost) {
                    stat = Spreading_stat::BAD_SIGNAL;
                    break;
                }
                auto R = rodrigues(K, theta);
                g->bfsRotate(v, u, R);
                if (loss_vector_ptr) {
                    loss_vector_ptr->emplace_back(val_cost, lj_cost);
                }
            }
            tree_index = add_tree(tree_index, *std::next(bone_line.begin()));
        }
        return tree_index;
    };

    dfs(dfs, g->polyBack()->getNeigh());

//    auto start = std::chrono::steady_clock::now();
//    std::vector<long long> time_memo;

    progresscpp::ProgressBar progress_bar(degree_of_polymerization-1, 70);
    for (int epoch = 1; epoch < degree_of_polymerization; epoch ++) {
        auto cur_mol = mol_genrator_ptr->getNext();
        g->attract(cur_mol);
        g->connect(cur_mol, -1, true);

        dfs(dfs, g->polyBack()->getNeigh());
        optimize_process(tree_index, optimize_atom_number);
        ++ progress_bar;
        progress_bar.display();

        if (pointer.left + 1 == pointer.right && pointer.right == (int)target_points.size() &&
        atom::positionMinusPosition(target_points.back(), g->getAtomPosition(g->polyBack()->getNeigh())).dot(atom::positionMinusPosition(g->getPolyPosition(1), g->getAtomPosition(g->polyBack()->getNeigh()))) <= 0) {
            stat = Spreading_stat::EARLY_STOP;
            break;
        }
        if (stat == Spreading_stat::BAD_SIGNAL) {
            break;
        }
//        if (epoch == 100 || epoch == 500 || epoch == 1000 || epoch == 2000 || epoch % 3000 == 0) {
//            auto end = std::chrono::steady_clock::now();
//            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//            time_memo.emplace_back(duration.count());
//            std::cout << "Mono: " << epoch << " Time Cost: " << duration.count() << " ms" << std::endl;
//        }
    }
//    std::cout << "[";
//    for (int i = 0; i < time_memo.size(); i ++) {
//        std::cout << time_memo[i];
//        if (i < (int)time_memo.size() - 1) std::cout << ", ";
//    }
//    std::cout << "]";

    progress_bar.done();
    if (stat == Spreading_stat::EARLY_STOP)
        std::cout << "Early stop!" << std::endl;
    else if (stat == Spreading_stat::FULL_STOP)
        std::cout << "Full stop!" << std::endl;
    else
        std::cout << "Warning: Bad signal!" << std::endl;

    optimize_process(tree_index, 3);
    assert(bone_line.size() == 2);
}