//
// Created by AnthonyZhang on 2025/1/21.
//

#include "graph.h"

// DSU
DSU::DSU(int n) : n(n), f(n), siz(n, 1) {
    std::iota(f.begin(), f.end(), 0);
}

int DSU::leader(int x) {
    while (x != f[x]) {
        x = f[x] = f[f[x]];
    }
    return x;
}

bool DSU::same(int x, int y) {
    return leader(x) == leader(y);
}

bool DSU::merge(int x, int y) {
    x = leader(x);
    y = leader(y);
    if (x == y) {
        return false;
    }
    siz[x] += siz[y];
    f[y] = x;
    return true;
}

int DSU::size(int x) {
    return siz[leader(x)];
}

// PolyManager
PolyManager::PolyManager() : front_ptr(nullptr), back_ptr(nullptr) {}

[[maybe_unused]] PolyManager::PolyManager(std::shared_ptr<Poly> ptr) : front_ptr(std::move(ptr)), back_ptr(nullptr) {}

[[maybe_unused]] PolyManager::PolyManager(std::shared_ptr<Poly> ptr1, std::shared_ptr<Poly> ptr2) : front_ptr(std::move(ptr1)), back_ptr(std::move(ptr2)) {}

PolyManager::PolyManager(const PolyManager& other) {
    front_ptr = (other.front_ptr ? std::make_shared<Poly>(*other.front_ptr) : nullptr);
    back_ptr = (other.back_ptr ? std::make_shared<Poly>(*other.back_ptr) : nullptr);
}

PolyManager::PolyManager(PolyManager&& other) noexcept : front_ptr(std::move(other.front_ptr)), back_ptr(std::move(other.back_ptr)) {
    other.front_ptr = nullptr;
    other.back_ptr = nullptr;
}

PolyManager& PolyManager::operator=(const PolyManager& other) {
    if (this != &other) {
        front_ptr = (other.front_ptr ? std::make_shared<Poly>(*other.front_ptr) : nullptr);
        back_ptr = (other.back_ptr ? std::make_shared<Poly>(*other.back_ptr) : nullptr);
    }
    return *this;
}

PolyManager& PolyManager::operator=(PolyManager&& other) noexcept {
    if (this != &other) {
        front_ptr = other.front_ptr;
        back_ptr = other.back_ptr;
        other.front_ptr = nullptr;
        other.back_ptr = nullptr;
    }
    return *this;
}

int PolyManager::size() const {
    return (int)(front_ptr != nullptr) + (int)(back_ptr != nullptr);
}

std::shared_ptr<Poly>& PolyManager::operator[](int index) {
    if (index >= size()) {
        throw exception::InvalidParameterException("Invalid index");
    }
    if (size() == 1) return front_ptr;
    if (index == 0)
        return front_ptr;
    return back_ptr;
}

std::shared_ptr<Poly>& PolyManager::operator()(int index) {
    if (index >= size()) {
        throw exception::InvalidParameterException("Invalid index");
    }
    if (size() == 1) return front_ptr;
    if (index == 0)
        return front_ptr;
    return back_ptr;
}

const std::shared_ptr<Poly>& PolyManager::front() const {
    if (size() == 0) {
        throw exception::InvalidParameterException("PolyManager(front) Empty PolyManager");
    }
    return front_ptr;
}

const std::shared_ptr<Poly>& PolyManager::back() const {
    if (size() == 0) {
        throw exception::InvalidParameterException("PolyManager(back) Empty PolyManager");
    }
    if (size() == 1) {
        return front_ptr;
    }
    return back_ptr;
}

void PolyManager::pop_front() {
    if (size() == 0) {
        throw exception::InvalidParameterException("PolyManager(pop_front) Empty PolyManager");
    }
    if (size() == 1) {
        front_ptr = nullptr;
        return;
    }
    front_ptr = std::move(back_ptr);
    back_ptr = nullptr;
}

void PolyManager::pop_back() {
    if (size() == 0) {
        throw exception::InvalidParameterException("PolyManager(pop_back) Empty PolyManager");
    }
    if (size() == 1) {
        front_ptr = nullptr;
        return;
    }
    back_ptr = nullptr;
}

void PolyManager::pop(int index) {
    if (index >= size()) {
        throw exception::InvalidParameterException("Invalid index");
    }
    if (size() == 1) {
        pop_front();
        return;
    }
    if (index == 0) {
        pop_front();
        return;
    }
    pop_back();
}

void PolyManager::push_back(std::shared_ptr<Poly> ptr) {
    if (size() >= 2) {
        throw exception::InvalidParameterException("PolyManager(push_back) Full PolyManager");
    }
    if (size() == 0) {
        front_ptr = std::move(ptr);
        return;
    }
    back_ptr = std::move(ptr);
}

void PolyManager::push_front(std::shared_ptr<Poly> ptr) {
    if (size() >= 2) {
        throw exception::InvalidParameterException("PolyManager(push_front) Full PolyManager");
    }
    if (size() == 1) {
        back_ptr = std::move(front_ptr);
    }
    front_ptr = std::move(ptr);
}

void PolyManager::translation(float dx, float dy, float dz) {
    if (size() == 0) return;
    if (size() == 2) {
        back_ptr->translation(dx, dy, dz);
    }
    front_ptr->translation(dx, dy, dz);
}

void PolyManager::randomize() {
    if (size() < 2) return;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);
    if (dis(gen)) {
        std::swap(front_ptr, back_ptr);
    }
}

void PolyManager::rotate(const Eigen::Matrix3f& rod, const Position& ver) {
    if (size() == 0) return;
    if (size() == 2) back_ptr->rotate(rod, ver);
    front_ptr->rotate(rod, ver);
}

void Graph::_addMono2verTable(int mono_idx, int atom_idx) {
    while (mono_idx >= (int)this->mono2ver.size()) {
        this->mono2ver.emplace_back();
    }
    this->mono2ver[mono_idx].push_back(atom_idx);
}

void Graph::_modifyMonomer(int mono_idx, int atom_idx) {
    auto old_mono_idx = this->getMonomer(atom_idx);
    monos[atom_idx] = mono_idx;
    mono2ver[old_mono_idx].erase(std::find(mono2ver[old_mono_idx].begin(), mono2ver[old_mono_idx].end(), atom_idx));
    _addMono2verTable(mono_idx, atom_idx);
}

Position Graph::_point_rotate_scale_translation(const Position& pos, const Position& ver, const Eigen::Matrix3f& rod, float scale_ratio, const Vector& dvec) {
    Vector vec = atom::positionMinusPosition(pos, ver);
    Position ans_p = ver + rod * vec * scale_ratio + dvec;
    return ans_p;
}

Position Graph::_point_rotate_scale(const Position& pos, const Position& ver, const Eigen::Matrix3f& rod, float scale_ratio) {
    Vector vec = atom::positionMinusPosition(pos, ver);
    Position ans_p = ver + rod * vec * scale_ratio;
    return ans_p;
}

Position Graph::_point_rotate(const Position& pos, const Position& ver, const Eigen::Matrix3f& rod) {
    Vector vec = atom::positionMinusPosition(pos, ver);
    Position ans_p = ver + rod * vec;
    return ans_p;
}

Position Graph::_point_scale(const Position& pos, const Position& ver, float scale_ratio) {
    Vector vec = atom::positionMinusPosition(pos, ver);
    Position ans_p = ver + vec * scale_ratio;
    return ans_p;
}
void Graph::_bfs_process(int root, int fa, const std::function<Position(const Position&, const Position&)>& f) {
    Position ver = getAtom(root)->getPosition();
    std::queue<int> q;
    q.emplace(root);

    std::unordered_set<int> visited;
    visited.insert(root);
    visited.insert(fa);

    while (!q.empty()) {
        auto cur = q.front();
        q.pop();
        for (const auto& ptr : getEdge(cur)) {
            auto nxt = ptr->getTo();
            if (visited.count(nxt)) continue;
            Position pos = getAtom(nxt)->getPosition();
            Position new_pos = f(pos, ver);
            vertices[nxt]->moveTo(new_pos);
            q.emplace(nxt);
            visited.insert(nxt);
        }
    }

    Position pos = polys.back()->getPosition();
    Position new_pos = f(pos, ver);
    polys.back()->moveTo(new_pos);
}

void Graph::_make_end(const std::string& end_symbol, int poly_index, bool poly_delete) {
    auto poly_ptr = polyBack();
    if (poly_delete) {
        polys.pop(poly_index);
    }
    // shorten C-C bond as C-H
    auto atom0_ptr = getAtom(poly_ptr->getNeigh());
    int mono_index = getMonomer(poly_ptr->getNeigh());

    Vector direct = atom::positionMinusPosition(poly_ptr->getPosition(), atom0_ptr->getPosition()).normalized() * BondTable::get_length(end_symbol, atom0_ptr->getSymbol(), "-");
    Position new_pos = atom0_ptr->getPosition() + direct;
    auto atom1_ptr = std::make_shared<Atom>(end_symbol, new_pos.x(), new_pos.y(), new_pos.z());
    this->addAtom(atom1_ptr, mono_index);
    int cur = n - 1;
    addEdge(cur, poly_ptr->getNeigh(), end_symbol);
    addEdge(poly_ptr->getNeigh(), cur, end_symbol);
}

Graph::Graph() : n(0), vertices(), edges(), monos(), polys(), is_period(false), ring_edges(), is_on_main_chain(), is_ar(), mono2ver() {}
Graph::Graph(int n, std::vector<std::shared_ptr<Atom>> atoms, std::vector<int> monos, std::vector<int> is_ar) : n(n), vertices(std::move(atoms)), edges(n), monos(std::move(monos)), polys(),
is_period(false), ring_edges(), is_on_main_chain(), is_ar(std::move(is_ar)), mono2ver() {
    for (int atom_idx = 0; atom_idx < n; atom_idx ++) {
        _addMono2verTable(0, atom_idx);
    }
}

Graph::Graph(const Graph& other) : n(other.n), edges(other.edges), monos(other.monos), polys(other.polys), is_period(other.is_period),
ring_edges(other.ring_edges), is_on_main_chain(other.is_on_main_chain), is_ar(other.is_ar), mono2ver(other.mono2ver) {
    for (const auto& ptr : other.vertices) {
        vertices.emplace_back(std::make_shared<Atom>(*ptr));
    }
}

Graph::Graph(Graph&& other) noexcept : n(other.n), vertices(std::move(other.vertices)), edges(std::move(other.edges)), monos(std::move(other.monos)), polys(std::move(other.polys)),
is_period(other.is_period), ring_edges(std::move(other.ring_edges)), is_on_main_chain(std::move(other.is_on_main_chain)), is_ar(std::move(other.is_ar)),
mono2ver(std::move(other.mono2ver)) {
    other.n = 0;
    other.is_period = false;
}

Graph& Graph::operator=(const Graph& other) {
    n = other.n;
    edges = other.edges;
    monos = other.monos;
    polys = other.polys;
    is_period = other.is_period;
    ring_edges = other.ring_edges;
    is_on_main_chain = other.is_on_main_chain;
    is_ar = other.is_ar;
    mono2ver = other.mono2ver;
    for (const auto& ptr : other.vertices) {
        vertices.emplace_back(std::make_shared<Atom>(*ptr));
    }

    return *this;
}

Graph& Graph::operator=(Graph&& other) noexcept {
    n = other.n;
    vertices = std::move(other.vertices);
    edges = std::move(other.edges);
    monos = std::move(other.monos);
    polys = std::move(other.polys);
    is_period = other.is_period;
    ring_edges = std::move(other.ring_edges);
    is_on_main_chain = std::move(other.is_on_main_chain);
    is_ar = std::move(other.is_ar);
    mono2ver = std::move(other.mono2ver);
    other.n = 0;
    other.is_period = false;

    return *this;
}

void Graph::addAtom(std::shared_ptr<Atom> atom_ptr, int mono, bool ar) {
    this->n += 1;
    this->vertices.push_back(std::move(atom_ptr));
    this->monos.push_back(mono);
    this->edges.emplace_back();
    this->is_ar.push_back(ar);
    this->_addMono2verTable(mono, n - 1);
}

void Graph::addAtom(int index, std::shared_ptr<Atom> atom_ptr, int mono, bool ar) {
    if (index >= n) {
        throw exception::InvalidParameterException("(Graph.addAtom) The index exceeds the array bound");
    }
    this->vertices[index] = std::move(atom_ptr);
    this->monos[index] = mono;
    this->is_ar[index] = ar;
    this->_modifyMonomer(mono, index);
}

void Graph::addEdge(int from, int to, const std::string& type, bool onring) {
    if (from < 0 || from >= n || to < 0 || to >= n) {
        throw exception::InvalidParameterException("(Graph.addEdge) Index exceeds bound");
    }
    auto ptr = std::make_shared<Edge>(to, type);
    edges[from].emplace_back(ptr);
    if (onring && from < to) {
        addRingEdge(from, to);
    }
}

bool Graph::checkOnRing(int from, int to) {
    if (from > to)
        std::swap(from, to);
    return this->ring_edges.count(std::make_pair(from, to));
}

bool Graph::checkOnMainChain(int idx) {
    return is_on_main_chain[idx];
}

const std::vector<int>& Graph::getOnMainChain() const {
    return is_on_main_chain;
}

int Graph::size() const {
    return n;
}

bool Graph::is_periodic() const {
    return is_period;
}

void Graph::addRingEdge(int from, int to) {
    if (from > to)
        std::swap(from, to);
    ring_edges.insert(std::make_pair(from, to));
}

std::shared_ptr<Atom> Graph::getAtom(int index) {
    if (index < 0 || index >= n) {
        throw exception::InvalidParameterException("(Graph.getAtom) Index exceeds bound");
    }
    return this->vertices[index];
}

Position Graph::getAtomPosition(int index) {
    auto p = getAtom(index);
    return p->getPosition();
}

const std::vector<std::shared_ptr<Atom>>& Graph::getAtomVec() const {
    return this->vertices;
}

int Graph::getMonomer(int index) const {
    if (index < 0 || index >= n) {
        throw exception::InvalidParameterException("(Graph.getMonomer) Index exceeds bound");
    }
    return monos[index];
}

const std::vector<int>& Graph::getMonomerVec() const {
    return monos;
}

const std::vector<int>& Graph::getMono2verIndex(int index) const {
    if (index < 0 || index >= (int)mono2ver.size()) {
        throw exception::InvalidParameterException("(Graph.getMono2verIndex) Index exceeds bound");
    }
    return mono2ver[index];
}

const std::vector<std::vector<int>>& Graph::getMono2verTable() const {
    return mono2ver;
}

bool Graph::isAr(int index) const {
    if (index < 0 || index >= n) {
        throw exception::InvalidParameterException("(Graph.isAr) Index exceeds bound");
    }
    return is_ar[index];
}

const std::vector<int>& Graph::getArVec() const {
    return is_ar;
}

float Graph::atomDistance(int idx1, int idx2) {
    auto p1 = getAtomPosition(idx1);
    auto p2 = getAtomPosition(idx2);
    return atom::positionDistance(p1, p2);
}

const std::vector<std::shared_ptr<Edge>>& Graph::getEdge(int index) const {
    if (index < 0 || index >= n) {
        throw exception::InvalidParameterException("(Graph.getEdge) Index exceeds bound");
    }
    return edges[index];
}

int Graph::getPolysSize() const {
    return polys.size();
}

int Graph::getEdgesSize() const {
    int size = 0;
    for (const auto& vec : this->getEdgesVec()) {
        size += (int)vec.size();
    }
    return size >> 1;
}

const std::vector<std::vector<std::shared_ptr<Edge>>>& Graph::getEdgesVec() const {
    return edges;
}

void Graph::addPoly(float x, float y, float z, int neigh) {
    this->polys.push_back(std::make_shared<Poly>(neigh, x, y, z));
}

std::shared_ptr<Poly>& Graph::getPoly(int index) {
    return polys(index);
}

const std::shared_ptr<Poly>& Graph::polyFront() const {
    return polys.front();
}

const std::shared_ptr<Poly>& Graph::polyBack() const {
    return polys.back();
}

Position Graph::getPolyPosition(int index) {
    return polys[index]->getPosition();
}

const std::unordered_set<std::pair<int, int>, hashing::pair_hash<int>>& Graph::getRingEdges() const {
    return ring_edges;
}

float Graph::polyDistanceWithNeigh(int index) {
    auto ptr = getPoly(index);
    int nei = ptr->getNeigh();
    auto poly_position = ptr->getPosition();
    auto atom_position = getAtom(nei)->getPosition();
    return atom::positionDistance(poly_position, atom_position);
}

void Graph::translation(float dx, float dy, float dz) {
    for (const auto& ptr : getAtomVec()) {
        ptr->translation(dx, dy, dz);
    }
    polys.translation(dx, dy, dz);
}

void Graph::translation(const Vector & vec) {
    translation(vec[0], vec[1], vec[2]);
}

void Graph::polyRandomize() {
    polys.randomize();
}

void Graph::polyPopBack() {
    polys.pop_back();
}

void Graph::polyPopFront() {
    polys.pop_front();
}

void Graph::rotate(const Eigen::Matrix3f& rod, const Position& ver) {
    for (const auto& ptr : getAtomVec()) {
        ptr->rotate(rod, ver);
    }
    polys.rotate(rod, ver);
}

void Graph::startRandomDirection() {
    auto K = Vector::Random(3).normalized();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 2.0);
    auto theta = static_cast<float>(dis(gen) * M_PI);
    Position ver{static_cast<float>(dis(gen)), static_cast<float>(dis(gen)), static_cast<float>(dis(gen))};
    auto rod = rodrigues(K, theta);
    rotate(rod, ver);
}

Vector Graph::getRotateAxis(int poly_index) {
    auto poly_ptr = polys[poly_index];
    int nei = poly_ptr->getNeigh();
    auto atom_ptr = getAtom(nei);

    Vector K;
    K << atom_ptr->getx() - poly_ptr->getx(), atom_ptr->gety() - poly_ptr->gety(), atom_ptr->getz() - poly_ptr->getz();
    return K;
}

void Graph::calMainChain() {
    assert(polys.size() == 2);
    std::vector<int> fa(n, -1);
    this->is_on_main_chain.resize(n, false);
    int start = polyFront()->getNeigh();
    int end = polyBack()->getNeigh();
    std::queue<int> q;
    q.emplace(start);
    fa[start] = start;
    is_on_main_chain[start] = true;

    while (!q.empty() && fa[end] == -1) {
        auto ver = q.front();
        q.pop();
        for (const auto& ptr : this->getEdge(ver)) {
            auto x = ptr->getTo();
            if (fa[x] >= 0) continue;
            fa[x] = ver;
            q.emplace(x);
        }
    }

    while (end != start) {
        is_on_main_chain[end] = true;
        end = fa[end];
    }

}

bool Graph::isOnMainChainIdx(int index) const {
    if ((int)is_on_main_chain.size() < n) {
        throw std::runtime_error("Graph.calMainChain() needs to be called before");
    }
    if (index >= n || index < 0) {
        throw exception::InvalidParameterException("Graph.isOnMainChainIdx Index exceeds bound");
    }
    return this->is_on_main_chain[index];
}

void Graph::bfsRotate(int root, int fa, const Eigen::Matrix3f& rod) {
    auto func = [&rod](const Position& pos, const Position& ver) {
        return _point_rotate(pos, ver, rod);
    };
    _bfs_process(root, fa, func);
}

void Graph::bfsScale(int root, int fa, float scale_ratio) {
    auto func = [&scale_ratio](const Position& pos, const Position& ver) {
        return _point_scale(pos, ver, scale_ratio);
    };
    _bfs_process(root, fa, func);
}

void Graph::bfsRotateScale(int root, int fa, const Eigen::Matrix3f& rod, float scale_ratio) {
    auto func = [&rod, &scale_ratio](const Position& pos, const Position& ver) {
        return _point_rotate_scale(pos, ver, rod, scale_ratio);
    };
    _bfs_process(root, fa, func);
}

void Graph::bfsRotateScaleTranslation(int root, int fa, const Eigen::Matrix3f& rod, float scale_ratio, float dx, float dy, float dz) {
    Vector dvec;
    dvec << dx, dy, dz;
    auto func = [&rod, &scale_ratio, &dvec](const Position& pos, const Position& ver) {
        return _point_rotate_scale_translation(pos, ver, rod, scale_ratio, dvec);
    };
    _bfs_process(root, fa, func);
}

void Graph::attachPoint(const Position& target_poly_point, const Vector& target_direction) {
    // 主动向目标点移动 target_direction是向外的 从atom指向poly的
    Vector vec = atom::positionMinusPosition(target_poly_point, polyFront()->getPosition());
    translation(vec);

    // rotate
    Position ver = getAtomPosition(polyFront()->getNeigh());
    vec = atom::positionMinusPosition(polyFront()->getPosition(), ver);
    auto R = rotateMatrix(vec, target_direction);
    rotate(R, ver);
}

void Graph::attract(const std::shared_ptr<Graph>& g, int poly_index) {
    // 吸引别的Graph过来
    auto poly1_ptr = getPoly(poly_index);
    auto g_poly_ptr = g->polyFront();
    auto g_head_atom_ptr = g->getAtom(g_poly_ptr->getNeigh());

    g->translation(poly1_ptr->getx() - g_head_atom_ptr->getx(), poly1_ptr->gety() - g_head_atom_ptr->gety(), poly1_ptr->getz() - g_head_atom_ptr->getz());

    // rotate
    Vector vec1 = g->getRotateAxis(0);
    vec1 = -vec1;
    Vector vec2 = this->getRotateAxis(poly_index);
    Position ver = this->getPolyPosition(poly_index);

    auto R = rotateMatrix(vec1, vec2);
    g->rotate(R, ver);
}

void Graph::connect(const std::shared_ptr<Graph>& g, int poly_index, bool poly_delete) {
    // please use self.attract before this function
    auto poly1_ptr = getPoly(poly_index);
    auto poly2_ptr = g->polyFront();

    // change index
    int old_size = this->size();
    int old_mono = (int)this->getMono2verTable().size();

    vertices.insert(vertices.end(), g->getAtomVec().begin(), g->getAtomVec().end());
    this->n = (int)vertices.size();
    monos.insert(monos.end(), g->getMonomerVec().begin(), g->getMonomerVec().end());
    for (int i = 0; i < g->size(); i ++) {
        monos[i+old_size] += old_mono;
        _addMono2verTable(monos[i+old_size], i+old_size);
    }

    is_ar.insert(is_ar.end(), g->getArVec().begin(), g->getArVec().end());
    is_on_main_chain.insert(is_on_main_chain.end(), g->getOnMainChain().begin(), g->getOnMainChain().end());

    for (const auto& vec : g->getEdgesVec()) {
        this->edges.emplace_back();
        for (const auto& ptr : vec) {
            int x = ptr->getTo() + old_size;
            this->edges.back().emplace_back(std::make_shared<Edge>(x, ptr->getType()));
        }
    }

    for (const auto& [from, to] : g->getRingEdges()) {
        ring_edges.insert(std::make_pair(from + old_size, to + old_size));
    }

    // add new bond
    addEdge(poly1_ptr->getNeigh(), poly2_ptr->getNeigh() + old_size, "1", false);
    addEdge(poly2_ptr->getNeigh() + old_size, poly1_ptr->getNeigh(), "1", false);

    if (poly_delete) {
        polys.pop(poly_index);
        this->polys.push_back(std::make_shared<Poly>(g->polyBack()->getNeigh() + old_size, g->polyBack()->getx(), g->polyBack()->gety(), g->polyBack()->getz()));
    }
}

void Graph::makeStart(std::shared_ptr<Graph>& g_ptr) {
    *this = std::move(*g_ptr);
}

void Graph::makePeriodicBoundaryCondition(const std::array<float, 3>& box_size) {
    if (is_period) return;
    for (const auto& ptr : getAtomVec()) {
        auto new_x = ptr->getx() - ptr->getx() / box_size[0] * box_size[0];
        auto new_y = ptr->gety() - ptr->gety() / box_size[1] * box_size[1];
        auto new_z = ptr->getz() - ptr->getz() / box_size[2] * box_size[2];
        ptr->moveTo(new_x, new_y, new_z);
    }
    is_period = true;
}

void Graph::makeEnd(const std::string &symbol) {
    while (polys.size() > 0) {
        this->_make_end(symbol, -1, true);
    }
}

std::ostream &operator<<(std::ostream &os, const Graph &g) {
    os << "Graph (Atom_num: " << g.size() << ", "
    << "Poly_num: " << g.getPolysSize() << ", "
    << "PolyFrontNeigh: " << g.polyFront()->getNeigh() << ' ' << "PolyBackNeigh: " << g.polyBack()->getNeigh() << std::endl
    << "Atoms: " << std::endl;
    for (const auto& ptr : g.getAtomVec()) {
        os << *ptr << std::endl;
    }
    os << "Polys: " << std::endl;
    if (g.polyFront()) {
        os << "Polyfront: " << *(g.polyFront()) << std::endl;
    }
    if (g.polyBack()) {
        os << "Polyback: " << *(g.polyBack()) << std::endl;
    }
    os << std::endl;
    return os;
}