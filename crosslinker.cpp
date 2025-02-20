//
// Created by AnthonyZhang on 2025/2/17.
//

#include "crosslinker.h"


CrossLinker::CrossLinker(int n, std::vector<std::shared_ptr<Atom>> atoms, std::vector<int> is_ar) :
n(n), vertices(std::move(atoms)), is_ar(std::move(is_ar)) {}

CrossLinker::CrossLinker(const CrossLinker &other) {
    n = other.n;
    vertices.resize(other.n);
    for (int i = 0; i < other.n; i ++) {
        vertices[i] = std::make_shared<Atom>(*other.vertices[i]);
    }
    polys.resize((int)other.polys.size());
    for (int i = 0; i < (int)other.polys.size(); i ++) {
        polys[i] = std::make_shared<Poly>(*other.polys[i]);
    }
    edges = other.edges;
    is_ar = other.is_ar;
}

CrossLinker::CrossLinker(CrossLinker &&other) noexcept {
    n = other.n;
    vertices = std::move(other.vertices);
    polys = std::move(other.polys);
    edges = std::move(other.edges);
    is_ar = std::move(other.is_ar);
    other.n = 0;
}

CrossLinker& CrossLinker::operator=(const CrossLinker &other) {
    n = other.n;
    vertices.resize(other.n);
    for (int i = 0; i < other.n; i ++) {
        vertices[i] = std::make_shared<Atom>(*other.vertices[i]);
    }
    polys.resize((int)other.polys.size());
    for (int i = 0; i < (int)other.polys.size(); i ++) {
        polys[i] = std::make_shared<Poly>(*other.polys[i]);
    }
    edges = other.edges;
    is_ar = other.is_ar;
    return *this;
}

CrossLinker& CrossLinker::operator=(CrossLinker &&other) noexcept {
    n = other.n;
    vertices = std::move(other.vertices);
    polys = std::move(other.polys);
    edges = std::move(other.edges);
    is_ar = std::move(other.is_ar);
    other.n = 0;
    return *this;
}

int CrossLinker::size() const {
    return n;
}

int CrossLinker::getPolysSize() const {
    return (int)polys.size();
}

int CrossLinker::getEdgesSize() const {
    int size = 0;
    for (const auto& vec : this->getEdgeVec()) {
        size += (int)vec.size();
    }
    return size >> 1;
}

std::shared_ptr<Atom> CrossLinker::getAtom(int index) const {
    if (index < 0 || index >= n) {
        throw exception::InvalidParameterException("(CrossLinker.getAtom) Index exceeds bound");
    }
    return this->vertices.at(index);
}

Position CrossLinker::getAtomPosition(int index) const {
    auto p = getAtom(index);
    return p->getPosition();
}

const std::vector<std::shared_ptr<Atom>>& CrossLinker::getAtomVec() const {
    return vertices;
}

std::shared_ptr<Poly> CrossLinker::getPoly(int index) const {
    if (index < 0 || index >= getPolysSize()) {
        throw exception::InvalidParameterException("(CrossLinker.getPoly) Index exceeds bound");
    }
    return this->polys.at(index);
}

Position CrossLinker::getPolyPosition(int index) const {
    auto p = getPoly(index);
    if (!p) {
        throw std::runtime_error("(getPolyPosition: Poly is null)");
    }
    return p->getPosition();
}

const std::vector<std::shared_ptr<Poly>>& CrossLinker::getPolyVec() const {
    return polys;
}

void CrossLinker::addAtom(std::shared_ptr<Atom> atom_ptr, bool ar) {
    this->n += 1;
    this->vertices.push_back(std::move(atom_ptr));
    this->edges.emplace_back();
    this->is_ar.push_back(ar);
}

void CrossLinker::addAtom(int index, std::shared_ptr<Atom> atom_ptr, bool ar) {
    if (index >= n) {
        throw exception::InvalidParameterException("(CrossLinker.addAtom) The index exceeds the array bound");
    }
    this->vertices[index] = std::move(atom_ptr);
    this->is_ar[index] = ar;
}

void CrossLinker::addEdge(int from, int to, const std::string &type) {
    if (from < 0 || from >= n || to < 0 || to >= n) {
        throw exception::InvalidParameterException("(CrossLinker.addEdge) Index exceeds bound");
    }
    auto ptr = std::make_shared<Edge>(to, type);
    edges[from].emplace_back(ptr);
}

void CrossLinker::addPoly(float x, float y, float z, int neigh) {
    this->polys.push_back(std::make_shared<Poly>(neigh, x, y, z));
}

bool CrossLinker::isAr(int index) const {
    if (index < 0 || index >= n) {
        throw exception::InvalidParameterException("(CrossLinker.isAr) Index exceeds bound");
    }
    return is_ar[index];
}

const std::vector<int>& CrossLinker::getArVec() const {
    return is_ar;
}

const std::vector<std::shared_ptr<Edge>>& CrossLinker::getEdge(int index) const {
    if (index < 0 || index >= n) {
        throw exception::InvalidParameterException("(CrossLinker.getEdge) Index exceeds bound");
    }
    return edges[index];
}

const std::vector<std::vector<std::shared_ptr<Edge>>>& CrossLinker::getEdgeVec() const {
    return edges;
}

std::ostream &operator<<(std::ostream &os, const CrossLinker &cl) {
    os << "CrossLinker: (Atom_num: " << cl.size() << ", "
    << "Poly_num: " << cl.getPolysSize() << std::endl;
    os << "Atoms: " << std::endl;
    for (const auto& ptr : cl.getAtomVec()) {
        if (!ptr) {
            os << "nullptr" << std::endl;
            continue;
        }
        os << *ptr << std::endl;
    }
    os << "Polys: " << std::endl;
    for (const auto& ptr : cl.getPolyVec()) {
        if (!ptr) {
            os << "nullptr" << std::endl;
            continue;
        }
        os << *ptr << std::endl;
    }
    return os;
}

void CrossLinker::makeEnd(int poly_index, const std::string& end_symbol) {
    auto poly_ptr = polys.at(poly_index);
    // shorten C-C bond as C-H
    auto atom0_ptr = getAtom(poly_ptr->getNeigh());

    Vector direct = atom::positionMinusPosition(poly_ptr->getPosition(), atom0_ptr->getPosition()).normalized() * BondTable::get_length(end_symbol, atom0_ptr->getSymbol(), "-");
    Position new_pos = atom0_ptr->getPosition() + direct;
    auto atom1_ptr = std::make_shared<Atom>(end_symbol, new_pos.x(), new_pos.y(), new_pos.z());
    this->addAtom(atom1_ptr);
    int cur = n - 1;
    addEdge(cur, poly_ptr->getNeigh(), end_symbol);
    addEdge(poly_ptr->getNeigh(), cur, end_symbol);
}

CrosslinkingSystem::CrosslinkingSystem(std::vector<std::shared_ptr<CrossLinker>> crosslinkers, std::vector<std::array<int, 4>> crosslinker_network,
                                       std::vector<std::vector<Position>> point_lists) :
                                       crosslinkers_(std::move(crosslinkers)), crosslinker_network_(std::move(crosslinker_network)),
                                       point_lists_(std::move(point_lists)) {
    assert((int)point_lists.size() == (int)crosslinker_network.size());
    tree_ptr_ = std::make_shared<Grid>();
    for (const auto& cl : crosslinkers_) {
        poly_seen_.emplace_back(std::vector<int>(cl->getPolysSize()));
        tree_ptr_->add_mol(cl, true);
    }
    chain_graphs_.resize((int)crosslinker_network.size(), nullptr);
}

int CrosslinkingSystem::getAtomSize() const {
    int sum = 0;
    for (const auto& cl : crosslinkers_) {
        sum += cl->size();
    }
    for (const auto& chain_ptr : chain_graphs_) {
        sum += chain_ptr->size();
    }
    return sum;
}

int CrosslinkingSystem::getEdgeSize() const {
    int sum = 0;
    for (const auto& cl : crosslinkers_) {
        sum += cl->getEdgesSize();
    }
    for (const auto& chain_ptr : chain_graphs_) {
        sum += chain_ptr->getEdgesSize();
    }
    sum >>= 1;
    sum += getCrosslinkerNetworkNumber();
    return sum;
}

int CrosslinkingSystem::getCrosslinkerNumber() const {
    return (int)crosslinkers_.size();
}

int CrosslinkingSystem::getCrosslinkerNetworkNumber() const {
    return (int)crosslinker_network_.size();
}

std::shared_ptr<CrossLinker> CrosslinkingSystem::getCrosslinkGraph(int index) const {
    return crosslinkers_.at(index);
}

const std::vector<std::shared_ptr<CrossLinker>>& CrosslinkingSystem::getCrosslinkGraphVec() const {
    return crosslinkers_;
}

std::shared_ptr<Graph> CrosslinkingSystem::getChainGraph(int index) const {
    if (index < 0 or index >= (int)chain_graphs_.size()) {
        throw exception::InvalidParameterException("CrosslinkingSystem: getChainGraph");
    }
    return chain_graphs_.at(index);
}

const std::vector<std::shared_ptr<Graph>>& CrosslinkingSystem::getChainGraphVec() const {
    return chain_graphs_;
}

const std::array<int, 4>& CrosslinkingSystem::getCrosslinkerNetworkInIdx(int index) const {
    if (index < 0 or index >= (int)crosslinker_network_.size()) {
        throw exception::InvalidParameterException("CrosslinkingSystem: getEdge");
    }
    return crosslinker_network_.at(index);
}

const std::vector<std::array<int, 4>>& CrosslinkingSystem::getCrosslinkerNetwork() const {
    return crosslinker_network_;
}

void CrosslinkingSystem::calcChainGraphs(int chain_index, const std::vector<std::shared_ptr<Graph>> &sequence,
                                    int degree_polymerization, bool random_polymerization, int optimize_size) {
    auto [cid1, cpid1, cid2, cpid2] = crosslinker_network_.at(chain_index);
    assert(!poly_seen_[cid1][cpid1] && !poly_seen_[cid2][cpid2]);

    auto& chain_ptr = chain_graphs_[chain_index];
    curveSpreading(point_lists_.at(chain_index), chain_ptr, tree_ptr_, sequence, degree_polymerization, 5.0f, 5, random_polymerization,
                   optimize_size, false);

    auto position1 = crosslinkers_.at(cid1)->getPolyPosition(cpid1);
    auto position2 = crosslinkers_.at(cid2)->getPolyPosition(cpid2);
    auto vec1 = chain_ptr->getAtomPosition(chain_ptr->polyBack()->getNeigh()) - chain_ptr->getAtomPosition(chain_ptr->polyBack()->getNeigh());
    auto vec2 = position2 - position1;

    auto rod = rotateMatrix(vec1, vec2);
    auto scale_ratio = vec2.norm() / vec1.norm();

    auto d_vec = position1 - chain_ptr->getAtomPosition(chain_ptr->polyBack()->getNeigh());
    auto dx = d_vec.x();
    auto dy = d_vec.y();
    auto dz = d_vec.z();

    tree_ptr_->erase_mol(chain_ptr);
    chain_ptr->bfsRotateScaleTranslation(chain_ptr->polyFront()->getNeigh(), -1, rod, scale_ratio, dx, dy, dz);
    tree_ptr_->add_mol(chain_ptr);

    poly_seen_[cid1][cpid1] = true;
    poly_seen_[cid2][cpid2] = true;
}

void CrosslinkingSystem::makeEnd(const std::string &end_system) {
    for (int i = 0; i < (int)crosslinkers_.size(); i ++) {
        const auto& cl = crosslinkers_[i];
        for (int j = 0; j < (int)cl->getPolysSize(); i ++) {
            if (!poly_seen_[i][j]) {
                poly_seen_[i][j] = true;
                cl->makeEnd(j, end_system);
            }
        }
    }
}

std::ostream &operator<<(std::ostream &os, const CrosslinkingSystem &cls) {
    os << "CrosslinkSystem:\n" <<
    "CrossLinker number: " << cls.crosslinkers_.size() << '\n' <<
    "Crosslinker_Network: " << cls.crosslinker_network_.size() << std::endl;
    return os;
}