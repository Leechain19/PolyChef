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
    return this->vertices[index];
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
    return this->polys[index];
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

void CrossLinker::addAtom(std::shared_ptr<Atom> atom_ptr, int mono, bool ar) {
    this->n += 1;
    this->vertices.push_back(std::move(atom_ptr));
    this->edges.emplace_back();
    this->is_ar.push_back(ar);
}

void CrossLinker::addAtom(int index, std::shared_ptr<Atom> atom_ptr, int mono, bool ar) {
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


CrosslinkingSystem::CrosslinkingSystem(std::vector<std::shared_ptr<CrossLinker>> crosslink_graphs, std::vector<std::array<int, 4>> cross_network) :
crosslink_graphs_(std::move(crosslink_graphs)), edges_(std::move(cross_network)) {
    n_ = (int)crosslink_graphs_.size();
    tree_ptr_ = std::make_shared<Grid>();
    for (const auto& cg : crosslink_graphs_) {
        poly_seen_.emplace_back(std::vector<int>(cg->getPolysSize()));
        tree_ptr_->add_mol(cg, true);
    }
}

int CrosslinkingSystem::getCrosslinkSize() const {
    return n_;
}

std::shared_ptr<CrossLinker> CrosslinkingSystem::getCrosslinkGraph(int index) const {
    if (index < 0 or index >= n_) {
        throw exception::InvalidParameterException("CrosslinkingSystem: getChainGraph");
    }
    return crosslink_graphs_[index];
}

const std::vector<std::shared_ptr<CrossLinker>>& CrosslinkingSystem::getCrosslinkGraphVec() const {
    return crosslink_graphs_;
}

std::shared_ptr<Graph> CrosslinkingSystem::getChainGraph(int index) const {
    if (index < 0 or index >= (int)chain_graphs_.size()) {
        throw exception::InvalidParameterException("CrosslinkingSystem: getChainGraph");
    }
    return chain_graphs_[index];
}

const std::vector<std::shared_ptr<Graph>>& CrosslinkingSystem::getChainGraphVec() const {
    return chain_graphs_;
}

const std::array<int, 4>& CrosslinkingSystem::getEdge(int index) const {
    if (index < 0 or index >= (int)edges_.size()) {
        throw exception::InvalidParameterException("CrosslinkingSystem: getEdge");
    }
    return edges_[index];
}

const std::vector<std::array<int, 4>>& CrosslinkingSystem::getEdgeVec() const {
    return edges_;
}

void CrosslinkingSystem::addEdge(const std::array<int, 4>& e) {
    edges_.emplace_back(e);
}

bool CrosslinkingSystem::calcChainGraphs(const std::array<int, 4>& e, const CustomFunctionLoader& custfunc, const std::vector<std::shared_ptr<Graph>>& sequence,
                                         int degree_polymerization, bool random_polymerization, int optimize_size) {
    if (sequence.empty()) {
        throw exception::InvalidParameterException("CorsslinkingSystem::calcChainGraphs: empty sequence");
    }
    return false;
}
