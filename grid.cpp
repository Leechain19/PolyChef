//
// Created by AnthonyZhang on 2025/1/7.
//
#include "grid.h"

Grid::Grid(float interval) : _interval(interval), _size(0) {}
Grid::~Grid() = default;

int Grid::_getpos(float x) const {
    return static_cast<int>(x / this->_interval);
}

float Grid::interval() const {
    return this->_interval;
}

int Grid::size() const {
    return this->_size;
}

//std::pair<std::shared_ptr<Atom>, float> Grid::search_nn(const Position &point) {
//    auto x = point[0], y = point[1], z = point[2];
//    auto xx = _getpos(x), yy = _getpos(y), zz = _getpos(z);
//    float min_dist = grid::inf;
//    std::shared_ptr<Atom> ret_ptr;
//
//    std::vector<std::array<int, 3>> neigh_coords;
//
//    for (int i = xx - 1; i <= xx + 1; i ++ ) {
//        for (int j = yy - 1; j <= yy + 1; j ++ ) {
//            for (int k = zz - 1; k <= zz + 1; k ++) {
//                neigh_coords.push_back({i, j, k});
//            }
//        }
//    }
//
//    for (const auto& coord : neigh_coords) {
//        auto it = mp.find(coord);
//        if (it == mp.end()) continue;
//        GridCell& cell = it->second;
//
//        for (const auto &ptr : cell.atoms) {
//            auto p = ptr->getPosition();
//            auto cur_dist = atom::positionDistance(p, point);
//            if (cur_dist < min_dist) {
//                min_dist = cur_dist;
//                ret_ptr = ptr;
//            }
//        }
//    }
//    return std::make_pair(ret_ptr, min_dist);
//}

std::pair<std::shared_ptr<Atom>, float> Grid::search_nn(const Position &point, float collision_threshold) {
    auto x = point[0], y = point[1], z = point[2];
    auto xx = _getpos(x), yy = _getpos(y), zz = _getpos(z);
    float min_dist_sq = std::numeric_limits<float>::max();
    std::shared_ptr<Atom> ret_ptr;

    // 预生成的27个偏移量
    static constexpr std::array<std::array<int, 3>, 27> offsets = {{
        {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1},
        {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1},
        {-1, 1, -1}, {-1, 1, 0}, {-1, 1, 1},
        {0, -1, -1}, {0, -1, 0}, {0, -1, 1},
        {0, 0, -1}, {0, 0, 0}, {0, 0, 1},
        {0, 1, -1}, {0, 1, 0}, {0, 1, 1},
        {1, -1, -1}, {1, -1, 0}, {1, -1, 1},
        {1, 0, -1}, {1, 0, 0}, {1, 0, 1},
        {1, 1, -1}, {1, 1, 0}, {1, 1, 1}
    }};

    float collision_threshold_sq = collision_threshold * collision_threshold;

    for (const auto& offset : offsets) {
        std::array<int, 3> coord{xx + offset[0], yy + offset[1], zz + offset[2]};
        auto it = mp.find(coord);
        if (it == mp.end()) continue;

        for (const auto& ptr : it->second.atoms) {
            auto p = ptr->getPosition();
            float dist_sq = atom::positionDistanceSquared(p, point);
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                ret_ptr = ptr;
                if (min_dist_sq < collision_threshold_sq) break; // 提前终止
            }
        }
        if (min_dist_sq < collision_threshold_sq) break;
    }

    return {ret_ptr, std::sqrt(min_dist_sq)};
}


void Grid::add(const Position &point, const std::string &symbol, bool including_hydrogen) {
    if (symbol.empty()) {
        throw exception::InvalidParameterException("The symbol string is empty");
    }
    if (!including_hydrogen && symbol[0] == 'H') return;
    int xx = _getpos(point.x()), yy = _getpos(point.y()), zz = _getpos(point.z());
    auto p = std::make_shared<Atom>(symbol, point);
    GridCell& cell = mp[{xx, yy, zz}];
    cell.atoms.emplace_back(p);
    _size += 1;
}

void Grid::add(const std::shared_ptr<Atom>& atom_ptr, bool including_hydrogen) {
    if (atom_ptr->getSymbol().empty()) {
        throw exception::InvalidParameterException("The symbol string is empty");
    }
    if (!including_hydrogen && atom_ptr->getSymbol()[0] == 'H') return;
    int xx = _getpos(atom_ptr->getx()), yy = _getpos(atom_ptr->gety()), zz = _getpos(atom_ptr->getz());
    GridCell& cell = mp[{xx, yy, zz}];
    cell.atoms.emplace_back(atom_ptr);
    _size += 1;
}

[[maybe_unused]] void Grid::add_mol(const Graph &g, bool including_hydrogen) {
    for (const auto& ptr : g.getAtomVec()) {
        add(ptr, including_hydrogen);
    }
}

bool Grid::isCollision(const Position &point, float threshold) {
    return search_nn(point).second < threshold;
}
