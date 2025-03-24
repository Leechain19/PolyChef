//
// Created by AnthonyZhang on 2025/1/23.
//
#include "optimizer.h"

#include <utility>

Pointer::Pointer(int left, int right) : left(left), right(right) {}

bool optimizer::EarlyTerminationFunction(const de::DifferentialEvolution& diffevo) {
    // 检查是否达到成本阈值
    if (diffevo.GetBestCost() <= costThreshold) {
//        std::cout << "Termination condition met: Best cost is below threshold." << std::endl;
        return true;
    }

    // 检查是否收敛
    static std::vector<double> previousCosts;
    previousCosts.push_back(diffevo.GetBestCost());

    // 限制记录的成本数量
    if (previousCosts.size() > convergenceThreshold) {
        previousCosts.erase(previousCosts.begin());
    }

    // 检查最近几次迭代的最佳成本是否有显著变化
    if (previousCosts.size() == convergenceThreshold) {
        double maxDiff = 0.0;
        for (size_t i = 1; i < previousCosts.size(); ++i) {
            maxDiff = std::max(maxDiff, std::abs(previousCosts[i] - previousCosts[i - 1]));
        }

        if (maxDiff < convergenceEpsilon) {
//            std::cout << "Termination condition met: Convergence detected." << std::endl;
            return true;
        }
    }
    return false;
}

Optimizer::Optimizer(float LJ_weight, std::shared_ptr<Grid> tree, const std::vector<Position>& target_points, const Pointer& pointer, int cal_len) :
LJ_weight(LJ_weight), tree(std::move(tree)), target_points(target_points), pointer(pointer), cal_len(cal_len) {}

std::pair<float, float> Optimizer::objective_fcn_pair(float angle) {
    auto R = rodrigues(K, angle);
    int sz = std::min((int)atoms_list.size(), cal_len);
    if (!sz) {
        return {0.0f, 0.0f};
    }

    std::vector<Position> cur_atoms;
    int cur = 0;
    for (const auto& p : atoms_list) {
        auto vec = atom::positionMinusPosition(p, root_position);
        auto cur_position = root_position + R * vec;
        cur_atoms.emplace_back(cur_position);
        if (++ cur >= sz) break;
    }

    float lj = 0.0f;
    float val = 0.0f;
    auto sz_mul = 1.0f / (float)sz;

    for (const auto& position : cur_atoms) {
        float sum = 0.0;
        for (int i = pointer.left; i < pointer.right; i ++) {
            sum += atom::positionDistance(position, target_points[i]);
        }
        sum /= static_cast<float>(pointer.right - pointer.left);
        val += sum;

        auto [node, min_dist] = tree->search_nn(position);
        if (!node) {
            min_dist = 1e6f;
        }
        const float sigma = 3.5f;
        lj += LJ_weight * static_cast<float>(qmi(sigma / min_dist, 6));
    }
    return {val * sz_mul, lj * sz_mul};
}

float Optimizer::objective_fcn(float angle) {
    auto [val, lj] = objective_fcn_pair(angle);
    return val + lj;
}

TargetFunction::TargetFunction(std::shared_ptr<Optimizer> ptr) : opti_ptr(std::move(ptr)) {}

unsigned int TargetFunction::NumberOfParameters() const {
    return 1;
}

std::vector<de::IOptimizable::Constraints> TargetFunction::GetConstraints() const {
    std::vector<de::IOptimizable::Constraints> constr(NumberOfParameters());
    for (auto& c : constr) {
        c = Constraints(-M_PI, M_PI, true);
    }
    return constr;
}

double TargetFunction::EvaluateCost(const Eigen::VectorXd& inputs) const {
    assert((int)inputs.size() == 1);
    auto x = static_cast<float>(inputs[0]);
    return opti_ptr->objective_fcn(x);
}

float optimizer::optimize(std::shared_ptr<Optimizer>& opti_opt, std::vector<Position> atom_list, Position root_position_, Vector K_, bool verbose) {
    opti_opt->atoms_list = std::move(atom_list);
    opti_opt->root_position = std::move(root_position_);
    opti_opt->K = std::move(K_);

    TargetFunction cost(opti_opt);
    de::DifferentialEvolution diffevo(cost, 32, EarlyTerminationFunction, 1);
    bool ok = diffevo.Optimize(1000, false); // turns + verbose

    auto best_xs = diffevo.GetBestAgent();
    double best_cost = diffevo.GetBestCost();

    if (verbose) {
        std::cout << "Best Solution: ";
        for (double val : best_xs) {
            std::cout << val << " ";
        }
        std::cout << '\n';
        std::cout << "Best Cost: " << best_cost << '\n';
        std::cout << (ok ? "Status: OK" : "Status: ") << std::endl;
    }
    return static_cast<float>(best_xs[0]);
}

double TestFunction::EvaluateCost(const Eigen::VectorXd& inputs) const {
    assert((int)inputs.size() == 1);
    auto x = static_cast<float>(inputs[0]);
    return x * x;
}

unsigned int TestFunction::NumberOfParameters() const {
    return 1;
}

std::vector<de::IOptimizable::Constraints> TestFunction::GetConstraints() const {
    std::vector<de::IOptimizable::Constraints> constr(NumberOfParameters());
    for (auto& c : constr) {
        c = Constraints(-M_PI, M_PI, true);
    }
    return constr;
}

float optimizer::testOptimize(bool verbose) {
    /*
     *  这个是Test!!!!
     */
    TestFunction cost;
    de::DifferentialEvolution diffevo(cost, 32, EarlyTerminationFunction, 1);
    bool ok = diffevo.Optimize(1000, verbose);

    auto best_xs = diffevo.GetBestAgent();
    double best_cost = diffevo.GetBestCost();

    if (verbose) {
        std::cout << "Best Solution: ";
        for (double val : best_xs) {
            std::cout << val << " ";
        }
        std::cout << '\n';
        std::cout << "Best Cost: " << best_cost << '\n';
        std::cout << (ok ? "Status: OK" : "Status: ") << std::endl;
    }

    return static_cast<float>(best_xs[0]);
}
