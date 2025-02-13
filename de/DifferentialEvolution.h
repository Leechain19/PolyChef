/**
 * \file DifferentialEvolution.h
 * \author Milos Stojanovic Stojke (milsto)
 * \modified by AnthonyZhang
 *
 * Implementation of Differential evolution algorithm.
 */

#pragma once

#include <iostream>
#include <utility>
#include <vector>
#include <cassert>
#include <random>
#include <iomanip>
#include <utility>
#include <memory>
#include <limits>
#include <functional>
#include <omp.h>
#include <Eigen/Dense>

namespace de {
    class IOptimizable {
    public:
        struct Constraints {
            double lower;
            double upper;
            bool isConstrained;

            explicit Constraints(double lower = 0.0, double upper = 1.0, bool isConstrained = false) : lower(lower), upper(upper),
                isConstrained(isConstrained) {}

            [[nodiscard]] bool Check(double candidate) const {
                return !isConstrained || (candidate >= lower && candidate <= upper);
            }
        };

        IOptimizable() = default;
        [[nodiscard]] virtual double EvaluateCost(const Eigen::VectorXd& inputs) const = 0;
        [[nodiscard]] virtual unsigned int NumberOfParameters() const = 0;
        [[nodiscard]] virtual std::vector<Constraints> GetConstraints() const = 0;
        virtual ~IOptimizable() = default;
    };

    class DifferentialEvolution {
    public:
        /**
         * Construct Differential Evolution optimizer
         *
         * \param costFunction Cost function to minimize
         * \param populationSize Number of agents in each optimization step
         * \param randomSeed Set random seed to a fix value to have repeatable (non stochastic) experiments
         * \param shouldCheckConstraints Should constraints bee checked on for each new candidate.
         * This check check may be turned off to increase performance if the cost function is defined
         * and has no local minimum outside of the constraints.
         * \param callback Optional callback to be called after each optimization iteration has finished.
         * Optimization iteration is defined as processing of single population with SelectionAndCorssing method.
         */

        DifferentialEvolution(const IOptimizable& costFunction, unsigned int populationSize, int max_trial) : DifferentialEvolution(costFunction, populationSize, 123, true, nullptr, nullptr, max_trial) {}

        DifferentialEvolution(const IOptimizable& costFunction, unsigned int populationSize, std::function<bool(const DifferentialEvolution&)> terminationCondition, int max_trial) : DifferentialEvolution(costFunction, populationSize, 123, true, nullptr, std::move(terminationCondition), max_trial) {}

        DifferentialEvolution(const IOptimizable& costFunction, unsigned int populationSize, int randomSeed, bool shouldCheckConstraints, std::function<void(const DifferentialEvolution&)> callback,
                                std::function<bool(const DifferentialEvolution&)> terminationCondition, int max_trial) :
            m_cost(costFunction), m_populationSize(populationSize), m_F(0.8), m_CR(0.9), m_bestAgentIndex(0), rng(randomSeed), max_trial(max_trial),
            m_minCost(std::numeric_limits<double>::infinity()), m_shouldCheckConstraints(shouldCheckConstraints),
            m_callback(std::move(callback)), m_terminationCondition(std::move(terminationCondition)) {
            assert(m_populationSize >= 4);
            m_numberOfParameters = m_cost.NumberOfParameters();
            m_population = Eigen::MatrixXd::Zero(m_populationSize, m_numberOfParameters);
            m_minCostPerAgent.resize(m_populationSize);
            m_constraints = costFunction.GetConstraints();
        }

        void InitPopulation() {
            #pragma omp parallel for
            for (unsigned int i = 0; i < m_populationSize; ++i) {
                thread_local std::mt19937 thread_rng(rng());
                for (unsigned int j = 0; j < m_numberOfParameters; ++j) {
                    double lower = m_constraints[j].isConstrained ? m_constraints[j].lower :  -std::numeric_limits<double>::infinity();
                    double upper = m_constraints[j].isConstrained ? m_constraints[j].upper :  std::numeric_limits<double>::infinity();
                    std::uniform_real_distribution<double> distribution(lower, upper);
                    m_population(i, j) = distribution(thread_rng);
                }
                m_minCostPerAgent[i] = m_cost.EvaluateCost(m_population.row(i));
                #pragma omp critical
                {
                    if (m_minCostPerAgent[i] < m_minCost) {
                        m_minCost = m_minCostPerAgent[i];
                        m_bestAgentIndex = i;
                    }
                }
            }
        }

        void SelectionAndCrossing() {
            Eigen::MatrixXd newPopulation = m_population;

            #pragma omp parallel for
            for (int x = 0; x < m_populationSize; ++x) {
                for (int trial = 0; trial < max_trial; trial ++) {
                    bool updated = ProcessIndividual(x, newPopulation);
                    if (updated) {
                        #pragma omp critical
                        {
                            if (m_minCostPerAgent[x] < m_minCost) {
                                m_minCost = m_minCostPerAgent[x];
                                m_bestAgentIndex = x;
                            }
                        }
                        break;
                    }
                }
            }
            m_population = newPopulation;
        }

        [[nodiscard]] Eigen::VectorXd GetBestAgent() const {
            return m_population.row(m_bestAgentIndex).eval();
        }

        [[nodiscard]] double GetBestCost() const {
            return m_minCost;
        }

        [[nodiscard]] std::vector<std::pair<Eigen::VectorXd, double>> GetPopulationWithCosts() const {
            std::vector<std::pair<Eigen::VectorXd, double>> result;
            for (int i = 0; i < m_populationSize; ++i) {
                result.emplace_back(m_population.row(i).eval(), m_minCostPerAgent[i]);
            }
            return result;
        }

        [[maybe_unused]] void PrintPopulation() const {
            for (int i = 0; i < m_populationSize; ++i) {
                for (int j = 0; j < m_numberOfParameters; ++j) {
                    std::cout << m_population(i, j) << " ";
                }
                std::cout << std::endl;
            }
        }

        bool Optimize(int iterations, bool verbose = true) {
            InitPopulation();

            // Optimization loop
            for (int i = 0; i < iterations; i ++) {
                // Optimization step
                SelectionAndCrossing();

                if (verbose) {
                    std::cout << std::fixed << std::setprecision(5);
                    std::cout << "Current minimal cost: " << m_minCost << "\t\t";
                    std::cout << "Best agent: ";
                    for (int j = 0; j < m_numberOfParameters; j ++) {
                        std::cout<< m_population(m_bestAgentIndex, j) << " ";
                    }
                    std::cout << std::endl;
                }

                if (m_callback) {
                    m_callback(*this);
                }

                if (m_terminationCondition && m_terminationCondition(*this)) {
                    if (verbose) {
                        std::cout << "Terminated due to positive evaluation of the termination condition." << std::endl;
                    }
                    return true;
                }
            }

            if (verbose) {
                std::cout << "Terminated due to exceeding total number of generations." << std::endl;
            }
            return false;
        }

    private:

        bool ProcessIndividual(int x, Eigen::MatrixXd& newPopulation) {
            // 线程局部随机数生成器
            thread_local std::mt19937 thread_rng(rng());
            std::uniform_int_distribution<int> distribution(0, (int)m_populationSize - 1);
            std::uniform_real_distribution<double> real_distribution(0.0, 1.0);

            // 随机选择 a, b, c
            int a, b, c;
            do {
                a = distribution(thread_rng);
                b = distribution(thread_rng);
                c = distribution(thread_rng);
            } while (a == x || b == x || c == x || a == b || a == c || b == c);

            // 形成中间解 z
            Eigen::VectorXd z = m_population.row(a) + m_F * (m_population.row(b) - m_population.row(c));
            Eigen::VectorXd newX = m_population.row(x);
            int R = distribution(thread_rng);

            // 交叉操作
            for (int i = 0; i < m_numberOfParameters; ++i) {
                if (real_distribution(thread_rng) < m_CR || i == R) {
                    newX(i) = z(i);
                }
            }

            // 检查约束条件
            if (m_shouldCheckConstraints && !CheckConstraints(newX)) {
                return false;  // 不满足约束条件，不更新
            }

            // 计算新成本
            double newCost = m_cost.EvaluateCost(newX);
            if (newCost < m_minCostPerAgent[x]) {
                newPopulation.row(x).noalias() = newX;
                m_minCostPerAgent[x] = newCost;
                return true;  // 更新成功
            }
            return false;  // 未更新
        }

        bool CheckConstraints(const Eigen::VectorXd& agent) {
            for (int i = 0; i < agent.size(); i ++) {
                if (!m_constraints[i].Check(agent[i])) {
                    return false;
                }
            }
            return true;
        }

        const IOptimizable& m_cost;
        unsigned int m_populationSize{};
        double m_F{};
        double m_CR{};
        unsigned int m_numberOfParameters{};
        bool m_shouldCheckConstraints{};

        std::function<void(const DifferentialEvolution&)> m_callback;
        std::function<bool(const DifferentialEvolution&)> m_terminationCondition;

        std::mt19937 rng;
        Eigen::MatrixXd m_population;

        std::vector<double> m_minCostPerAgent;
        std::vector<IOptimizable::Constraints> m_constraints;

        int m_bestAgentIndex;
        double m_minCost;
        int max_trial;
    };
}