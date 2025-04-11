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

    static constexpr std::array<std::array<int, 3>, 26> offsets = {{
        {-1, 0, 0},  {0, -1, 0},
        {0, 0, -1},   {1, 0, 0},   {0, 1, 0},
        {0, 0, 1},    {-1, -1, 0}, {-1, 0, -1},
        {-1, 0, 1},   {0, -1, 1},  {0, -1, -1},
        {-1, 1, 0},   {0, 1, 1},   {0, 1, -1},
        {1, -1, 0},   {1, 0, -1},  {1, 1, 0},
        {1, 0, 1},    {-1, -1, 1}, {-1, 1, -1},
        {-1, 1, 1},   {1, -1, -1}, {1, -1, 1},
        {-1, -1, -1}, {1, 1, -1},  {1, 1, 1}
    }};


    for (const auto& vec : offsets) {
        auto [dx, dy, dz] = vec;
        Position neighbor = position + Vector((float)dx, (float)dy, (float)dz).normalized() * step;
        neighbors.emplace_back(neighbor);
    }
    return neighbors;
}

std::vector<Position> AStar::AStarSearch(const Position& start, const Position& goal, const std::shared_ptr<Grid> &tree, float step, float collisionThreshold) {

//    auto cmp = [](const std::shared_ptr<AStarNode>& ptr1, const std::shared_ptr<AStarNode>& ptr2) -> bool {
//        return ptr1->f_cost() > ptr2->f_cost();
//    };
//
//    // 优先队列，存储 shared_ptr<AStarNode>
//    std::priority_queue<std::shared_ptr<AStarNode>, std::vector<std::shared_ptr<AStarNode>>, decltype(cmp)> openSet(cmp);
//
//    // 哈希表，用于记录所有已访问的节点
//    std::unordered_map<Position, std::shared_ptr<AStarNode>, hashing::eigen_vector_hash<float>, hashing::eigen_vector_equal<float>> allNodes;
//
//    // 创建起始节点
//    auto startNode = std::make_shared<AStarNode>(start, 0, heuristic(start, goal));
//    openSet.push(startNode);
//    allNodes[startNode->pos] = startNode;
//
//    while (!openSet.empty()) {
//        auto current = openSet.top();
//        openSet.pop();
//
//        // 检查是否到达目标
//        if (current->pos.isApprox(goal, 1e-3)) {
//            // 找到目标，回溯路径
//            std::vector<Position> path;
//            auto node = current;
//            while (node) {
//                path.push_back(node->pos);
//                node = std::dynamic_pointer_cast<AStarNode>(node->fa);
//            }
//            std::reverse(path.begin(), path.end());
//            return path;
//        }
//
//        // 生成邻居节点
//        for (const auto& neighborPos : getNeighbors(current->pos, step)) {
//            if (tree->isCollision(neighborPos, collisionThreshold)) continue;
//
//            float next_cost = current->cost + step;
//            if (allNodes.find(neighborPos) == allNodes.end() || next_cost < allNodes[neighborPos]->cost) {
//                auto neighborNode = std::make_shared<AStarNode>(neighborPos, next_cost, heuristic(neighborPos, goal), current);
//                openSet.push(neighborNode);
//                allNodes[neighborPos] = neighborNode;
//            }
//        }
//    }

    return {};  // 没有找到路径
}

std::vector<Position> RRT::bidirectionalRRT(const Position &start, const Position &goal, const std::shared_ptr<Grid>& grid_tree, const CustomFunctionLoader& custfunc,
                                            float step_size, float goal_tolerance, float to_end_possibility, float collisionThreshold, int max_trial) {

    auto kdt_start = std::make_unique<myKDTree>();
    auto kdt_goal = std::make_unique<myKDTree>();

    std::vector<Position> points_start = {start};
    std::vector<Position> points_goal = {goal};

    kdt_start->addPoint(start, 0);
    kdt_goal->addPoint(goal, 0);

    std::vector<int> fa_start(1, -1);
    std::vector<int> fa_goal(1, -1);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis0((float)custfunc.x_min(), (float)custfunc.x_max());
    std::uniform_real_distribution<float> dis1((float)custfunc.y_min(), (float)custfunc.y_max());
    std::uniform_real_distribution<float> dis2((float)custfunc.z_min(), (float)custfunc.z_max());
    std::uniform_real_distribution<float> dis_possibility(0.0, 1.0);

    progresscpp::ProgressBar progress_bar(max_trial, 70);

    for (int i = 0; i < max_trial; ++i) { // 最大迭代次数
        auto random_position = Position(dis0(gen), dis1(gen), dis2(gen));
        if (dis_possibility(gen) < to_end_possibility) {
            random_position = goal;
        }
        ++ progress_bar;
        progress_bar.display();

        auto cur_nearest_payload = kdt_start->searchKnn(random_position, 1).at(0);
        const Position& nearest_position = points_start.at(cur_nearest_payload.payload);
        Vector direction = (random_position - nearest_position).normalized();
        Position new_position = nearest_position + direction * step_size;

        bool flag = false;

        // 碰撞检测（这里假设没有障碍物）
        if (!grid_tree->isCollision(new_position, collisionThreshold)) {
            // 添加新节点
            points_start.push_back(new_position);
            fa_start.push_back(cur_nearest_payload.payload);

            // 更新 KD-Tree
            kdt_start->addPoint(points_start.back(), (int)points_start.size() - 1);
            flag = true;
        }

        if (flag) {
            // 检查是否连接到目标树
            const auto& new_point = points_start.back();
            auto nearest_payload = kdt_goal->searchKnn(new_point, 1).at(0);

            if (nearest_payload.distance < goal_tolerance) {
                // 找到路径
                std::cout << "Path found!" << std::endl;
                std::vector<Position> path;

                int cur_idx = (int)points_start.size() - 1;
                while (cur_idx >= 0) {
                    path.emplace_back(points_start.at(cur_idx));
                    cur_idx = fa_start.at(cur_idx);
                }
                std::reverse(path.begin(), path.end());

                cur_idx = nearest_payload.payload;
                while (cur_idx >= 0) {
                    path.emplace_back(points_goal.at(cur_idx));
                    cur_idx = fa_goal.at(cur_idx);
                }

                if (points_goal[0] == path.front()) {
                    std::reverse(path.begin(), path.end());
                }
                progress_bar.done();
                return path;
            }
        }

        // 交换
        std::swap(kdt_start, kdt_goal);
        std::swap(points_start, points_goal);
        std::swap(fa_start, fa_goal);
    }

    progress_bar.done();
    std::cout << "Path not found." << std::endl;
    return {};
}

std::vector<Position> calCRS(const std::vector<Position>& path, float distance_threshold) {
    const int n = (int)path.size();
    if (n == 0 || distance_threshold < 0.1) {
        throw exception::InvalidParameterException("Too small distance threshold");
    }
    Eigen::Matrix<float, 4, 4> CR_Matrix;
    CR_Matrix << 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,
    -3.0, -2.0, 3.0, -1.0,
    2.0, 1.0, -2.0, 1.0;

    Eigen::MatrixXf velocity(n, 3);
    for (int i = 0; i < n; i ++) {
        velocity.row(i) = 0.5 * (path[std::min(n-1, i+1)] - path[std::max(0, i-1)]);
    }

    std::vector<Position> all_points;

    for (int i = 0; i < n - 1; ++ i) {
        Eigen::Matrix<float, 4, 3> mtx;
        mtx.row(0) = path[i];
        mtx.row(1) = velocity.row(i);
        mtx.row(2) = path[i+1];
        mtx.row(3) = velocity.row(i+1);

        auto weight_matrix = CR_Matrix * mtx;
        float distance = atom::positionDistance(path[i], path[i+1]);
        int num_points = std::max(10, static_cast<int>(std::ceil(distance / distance_threshold)));

        for (int j = 0; j < num_points; j ++) {
            Eigen::Matrix<float, 1, 4> t;
            t(0) = 1.0f;
            for (int k = 1; k < 4; k ++) {
                t(k) = t(k-1) * (float)j / (float)num_points;
            }
            all_points.emplace_back((t * weight_matrix).transpose());
        }
    }
    return all_points;
}

std::vector<Position> getScatterFromCSV(const std::string& path, bool header, float distance_upper_bound) {
    std::vector<Position> points;
    std::ifstream file(path);
    std::string line;
    std::string token;

    // 跳过表头
    if (header)
        std::getline(file, line);

    int line_id = 1;
    while (std::getline(file, line)) {
        if (line.empty()) {
            line_id ++;
            continue;
        }
        std::istringstream ss(line);
        float x, y, z;
        bool valid = true;

        // 读取 x 坐标
        if (std::getline(ss, token, ',') &&!token.empty()) {
            try {
                x = static_cast<float>(std::stod(token));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid value for x: " << token << " in line: " << line_id << std::endl;
                valid = false;
            }
        } else {
            std::cerr << "Missing x value in line: " << line_id << std::endl;
            valid = false;
        }

        // 读取 y 坐标
        if (valid && std::getline(ss, token, ',') &&!token.empty()) {
            try {
                y = static_cast<float>(std::stod(token));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid value for y: " << token << " in line: " << line_id << std::endl;
                valid = false;
            }
        } else {
            std::cerr << "Missing y value in line: " << line_id << std::endl;
            valid = false;
        }

        // 读取 z 坐标
        if (valid && std::getline(ss, token) &&!token.empty()) {
            try {
                z = static_cast<float>(std::stod(token));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid value for z: " << token << " in line: " << line_id << std::endl;
                valid = false;
            }
        } else {
            std::cerr << "Missing z value in line: " << line_id << std::endl;
            valid = false;
        }

        // 如果所有值都有效，则添加到 points 向量中
        if (valid) {
            Position p(x, y, z);
            while (points.size() >= 2 && atom::positionDistance(p, points[points.size()-1]) < distance_upper_bound &&
                atom::positionDistance(p, points[points.size()-2]) < distance_upper_bound) {
                points.pop_back();
            }
            points.emplace_back(p);
        }
        line_id ++;
    }

    return calCRS(points, 0.8);
}