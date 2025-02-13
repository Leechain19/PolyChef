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

std::pair<std::shared_ptr<Atom>, float> Grid::search_nn(const Position &point) {
    auto x = point[0], y = point[1], z = point[2];
    auto xx = _getpos(x), yy = _getpos(y), zz = _getpos(z);
    float min_dist = grid::inf;
    std::shared_ptr<Atom> ret_ptr;

    std::vector<std::array<int, 3>> neigh_coords;

    for (int i = xx - 1; i <= xx + 1; i ++ ) {
        for (int j = yy - 1; j <= yy + 1; j ++ ) {
            for (int k = zz - 1; k <= zz + 1; k ++) {
                neigh_coords.push_back({i, j, k});
            }
        }
    }

    {
//        std::shared_lock global_read_lock(global_rw_lock);
        for (const auto& coord : neigh_coords) {
            auto it = mp.find(coord);
            if (it == mp.end()) continue;
            GridCell& cell = it->second;
//            std::shared_lock cell_read_lock(cell.cell_rw_lock);

            for (const auto &ptr : cell.atoms) {
                auto p = ptr->getPosition();
                auto cur_dist = atom::positionDistance(p, point);
                if (cur_dist < min_dist) {
                    min_dist = cur_dist;
                    ret_ptr = ptr;
                }
            }
        }
    }

    return std::make_pair(ret_ptr, min_dist);
}

void Grid::add(const Position &point, const std::string &symbol, bool including_hydrogen) {
    if (symbol.empty()) {
        throw exception::InvalidParameterException("The symbol string is empty");
    }
    if (!including_hydrogen && symbol[0] == 'H') return;
    int xx = _getpos(point.x()), yy = _getpos(point.y()), zz = _getpos(point.z());
    auto p = std::make_shared<Atom>(symbol, point);
    {
//        std::unique_lock global_write_lock(global_rw_lock);
        GridCell& cell = mp[{xx, yy, zz}];
//        std::unique_lock cell_write_lock(cell.cell_rw_lock);
        cell.atoms.emplace_back(p);
        _size += 1;
    }
}

void Grid::add(const std::shared_ptr<Atom>& atom_ptr, bool including_hydrogen) {
    if (atom_ptr->getSymbol().empty()) {
        throw exception::InvalidParameterException("The symbol string is empty");
    }
    if (!including_hydrogen && atom_ptr->getSymbol()[0] == 'H') return;
    int xx = _getpos(atom_ptr->getx()), yy = _getpos(atom_ptr->gety()), zz = _getpos(atom_ptr->getz());

    {
//        std::unique_lock global_write_lock(global_rw_lock);
        GridCell& cell = mp[{xx, yy, zz}];
//        std::unique_lock cell_write_lock(cell.cell_rw_lock);
        cell.atoms.emplace_back(atom_ptr);
        _size += 1;
    }
}

void Grid::add_mol(const Graph &g, bool including_hydrogen) {
    for (const auto& ptr : g.getAtomVec()) {
        add(ptr, including_hydrogen);
    }
}

bool Grid::isCollision(const Position &point, float threshold) {
    return search_nn(point).second < threshold;
}
