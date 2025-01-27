//
// Created by AnthonyZhang on 2025/1/18.
//

#include "curve.h"
#include <algorithm>
#include <utility>
#include "hash.h"

Node::Node(Position pos, std::shared_ptr<Node> fa) : pos(std::move(pos)), fa(std::move(fa)) {}

[[maybe_unused]] std::string Node::getPositionToString() const {
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

bool RRT::extendTree(const Position &random_position, std::unique_ptr<KDTree> &kdt, const std::shared_ptr<Grid> &grid_tree, std::vector<Point*>& point_ptrs,
                     std::vector<int>& fa, std::unordered_map<Point*, int>& Ptr2id, float step_size, float collisionThreshold) {

    Point random_point(random_position);
    Point* nearest_point = kdt->find_nearest(&random_point);
    Position nearest_position = nearest_point->to_vector3f();

    // 计算新点
    Vector direction = (random_position - nearest_position).normalized();
    Position new_position = nearest_position + direction * step_size;

    // 碰撞检测（这里假设没有障碍物）
    if (!grid_tree->isCollision(new_position, collisionThreshold)) {
        // 添加新节点
        auto* new_point = new Point(new_position);
        point_ptrs.push_back(new_point);
        fa.push_back(Ptr2id[nearest_point]);

        // 更新 KD-Tree
        kdt->insert(new_point);
        return true;
    }
    return false;
}

std::vector<Position> RRT::bidirectionalRRT(const Position &start, const Position &goal, const std::shared_ptr<Grid>& grid_tree, const Position& box_size,
                                            float step_size, float goal_tolerance, float to_end_possibility, float collisionThreshold, int max_trial) {

    auto kdt_start = std::make_unique<KDTree>(AABB(Point(0, 0, 0), Point(box_size)));
    auto kdt_goal = std::make_unique<KDTree>(AABB(Point(0, 0, 0), Point(box_size)));

    std::vector<Point*> point_ptrs_start = {new Point(start)};
    std::vector<Point*> point_ptrs_goal = {new Point(goal)};

    kdt_start->insert(point_ptrs_start[0]);
    kdt_goal->insert(point_ptrs_goal[0]);

    // 建立指针到下标的哈希映射
    // 建立父节点表
    std::unordered_map<Point*, int> Ptr2id_start = {{point_ptrs_start[0], 0}};
    std::unordered_map<Point*, int> Ptr2id_goal = {{point_ptrs_goal[0], 0}};

    std::vector<int> fa_start(1, -1);
    std::vector<int> fa_goal(1, -1);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis0(0.0, box_size.x());
    std::uniform_real_distribution<float> dis1(0.0, box_size.y());
    std::uniform_real_distribution<float> dis2(0.0, box_size.z());
    std::uniform_real_distribution<float> dis_possibility(0.0, 1.0);

    for (int i = 0; i < max_trial; ++i) { // 最大迭代次数
        auto random_position = Position(dis0(gen), dis1(gen), dis2(gen));
        if (dis_possibility(gen) < to_end_possibility) {
            random_position = goal;
        }

        if (extendTree(random_position, kdt_start, grid_tree, point_ptrs_start, fa_start, Ptr2id_start, step_size, collisionThreshold)) {
            // 检查是否连接到目标树
            auto* new_point = point_ptrs_start.back();
            auto* nearest_point = kdt_goal->find_nearest(new_point);

            if ((*new_point - *nearest_point).length() < goal_tolerance) {
                // 找到路径
                std::cout << "Path found!" << std::endl;
                std::vector<Position> path;

                int cur_idx = (int)point_ptrs_start.size() - 1;
                while (cur_idx >= 0) {
                    path.emplace_back(point_ptrs_start[cur_idx]->to_vector3f());
                    cur_idx = fa_start[cur_idx];
                }

                std::reverse(path.begin(), path.end());

                cur_idx = Ptr2id_goal[nearest_point];
                while (cur_idx >= 0) {
                    path.emplace_back(point_ptrs_goal[cur_idx]->to_vector3f());
                    cur_idx = fa_goal[cur_idx];
                }

                if (point_ptrs_goal[0]->to_vector3f() == path.front()) {
                    std::reverse(path.begin(), path.end());
                }

                // 析构
                for (auto& ptr : point_ptrs_start) {
                    delete ptr;
                }
                for (auto& ptr : point_ptrs_goal) {
                    delete ptr;
                }

                return path;
            }
        }

        // 交换
        std::swap(kdt_start, kdt_goal);
        std::swap(point_ptrs_start, point_ptrs_goal);
        std::swap(Ptr2id_start, Ptr2id_goal);
        std::swap(fa_start, fa_goal);
    }

    std::cout << "Path not found." << std::endl;

    // 析构
    for (auto& ptr : point_ptrs_start) {
        delete ptr;
    }
    for (auto& ptr : point_ptrs_goal) {
        delete ptr;
    }

    return {};
}

