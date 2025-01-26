//
// Created by AnthonyZhang on 2025/1/18.
//

#include "curve.h"
#include "hash.h"

Node::Node(Position pos, std::shared_ptr<Node> fa) : pos(std::move(pos)), fa(std::move(fa)) {}

std::string Node::getPositionToString() const {
    std::stringstream ss;
    ss << pos.format(Eigen::IOFormat(3, 0, ",", ",", "", "", "(", ")"));
    return ss.str();
}

AStarNode::AStarNode(Position pos, float cost, float h, std::shared_ptr<Node> fa) : Node(std::move(pos), std::move(fa)), cost(cost), h(h) {}

float AStarNode::f_cost() const {
    return cost + h;
}

float AStar::heuristic(const Position& start, const Position& goal) {
    return (start - goal).norm();
}

std::vector<Position> AStar::getNeighbors(const Position &position, float step) {
    std::vector<Position> neighbors;

    #pragma omp parallel for default(none), shared(neighbors, position, step)
    for (int dx = -1; dx <= 1; dx ++) {
        for (int dy = -1; dy <= 1; dy ++) {
            for (int dz = -1; dz <= 1; dz ++) {
                if (dx == 0 && dy == 0 && dz == 0) continue;  // 跳过自身
                Position neighbor = position + Vector((float)dx, (float)dy, (float)dz).normalized() * step;
                #pragma omp critical
                neighbors.emplace_back(neighbor);
            }
        }
    }
    return neighbors;
}

std::vector<Position> AStar::AStarSearch(const Position& start, const Position& goal, const std::shared_ptr<Grid> &tree, float step, float collisionThreshold) {

    auto cmp = [](const std::shared_ptr<AStarNode>& ptr1, const std::shared_ptr<AStarNode>& ptr2) -> bool {
        return ptr1->f_cost() > ptr2->f_cost();
    };

    // 优先队列，存储 shared_ptr<AStarNode>
    std::priority_queue<std::shared_ptr<AStarNode>, std::vector<std::shared_ptr<AStarNode>>, decltype(cmp)> openSet(cmp);

    // 哈希表，用于记录所有已访问的节点
    std::unordered_map<Position, std::shared_ptr<AStarNode>, hashing::eigen_vector_hash<float>, hashing::eigen_vector_equal<float>> allNodes;

    // 创建起始节点
    auto startNode = std::make_shared<AStarNode>(start, 0, heuristic(start, goal));
    openSet.push(startNode);
    allNodes[startNode->pos] = startNode;

    while (!openSet.empty()) {
        auto current = openSet.top();
        openSet.pop();

        // 检查是否到达目标
        if (current->pos.isApprox(goal, 1e-3)) {
            // 找到目标，回溯路径
            std::vector<Position> path;
            auto node = current;
            while (node) {
                path.push_back(node->pos);
                node = std::dynamic_pointer_cast<AStarNode>(node->fa);
            }
            std::reverse(path.begin(), path.end());
            return path;
        }

        // 生成邻居节点
        #pragma omp parallel for default(none), shared(start, goal, current, step, tree, allNodes, openSet, collisionThreshold)
        for (const auto& neighborPos : getNeighbors(current->pos, step)) {
            if (tree->isCollision(neighborPos, collisionThreshold)) continue;

            float next_cost = current->cost + step;

            #pragma omp critical
            {
                if (allNodes.find(neighborPos) == allNodes.end() || next_cost < allNodes[neighborPos]->cost) {
                    auto neighborNode = std::make_shared<AStarNode>(neighborPos, next_cost, heuristic(neighborPos, goal), current);
                    openSet.push(neighborNode);
                    allNodes[neighborPos] = neighborNode;
                }
            }

        }
    }

    return {};  // 没有找到路径
}
