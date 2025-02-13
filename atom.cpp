//
// Created by AnthonyZhang on 2025/1/21.
//

#include "atom.h"

Atom::Atom() : symbol(std::string()), x(0.0), y(0.0), z(0.0) {}

Atom::Atom(std::string symbol, float x, float y, float z) : symbol(std::move(symbol)), x(x), y(y), z(z) {}

Atom::Atom(std::string symbol, const Position& pos) : Atom(std::move(symbol), pos.x(), pos.y(), pos.z()) {}

Atom::Atom(Atom&& o) noexcept : symbol(std::move(o.symbol)), x(o.x), y(o.y), z(o.z) {}

Atom& Atom::operator=(const Atom& other) {
    if (this != &other) {
        symbol = other.symbol;
        x = other.x;
        y = other.y;
        z = other.z;
    }
    return *this;
}

Atom& Atom::operator=(Atom&& other) noexcept {
    if (this != &other) {
        symbol = std::move(other.symbol);
        x = other.x;
        y = other.y;
        z = other.z;
        // 将被移动的对象置为有效但未定义的状态
        other.x = other.y = other.z = 0.0f;
        other.symbol = "";
    }
    return *this;
}

std::string Atom::getSymbol() const {
    return this->symbol;
}

float Atom::getx() const {
    return x;
}

float Atom::gety() const {
    return y;
}

float Atom::getz() const {
    return z;
}

void Atom::setSymbol(const std::string& new_symbol) {
    this->symbol = new_symbol;
}

void Atom::setx(float new_x) {
    this->x = new_x;
}

void Atom::sety(float new_y) {
    this->y = new_y;
}

void Atom::setz(float new_z) {
    this->z = new_z;
}

void Atom::set(float new_x, float new_y, float new_z) {
    this->x = new_x;
    this->y = new_y;
    this->z = new_z;
}

void Atom::translation(float dx, float dy, float dz) {
    this->x += dx;
    this->y += dy;
    this->z += dz;
}

void Atom::moveTo(float nx, float ny, float nz) {
    this->x = nx;
    this->y = ny;
    this->z = nz;
}

void Atom::moveTo(const Position& pos) {
    moveTo(pos.x(), pos.y(), pos.z());
}

void Atom::rotate(const Eigen::Matrix3f& rod, const Position& ver) {
    Vector vec;
    vec << x - ver[0], y - ver[1], z - ver[2];
    auto vec_rot = rod * vec;
    x = vec_rot[0] + ver[0];
    y = vec_rot[1] + ver[1];
    z = vec_rot[2] + ver[2];
}

Position Atom::getPosition() const {
    return {x, y, z};
}

std::ostream& operator<<(std::ostream& os, const Atom& atom) {
    os << "Atom(symbol=" << atom.symbol
    << ", x=" << atom.x
    << ", y=" << atom.y
    << ", z=" << atom.z << ")";
    return os;
}

Poly::Poly() : neigh(-1), x(0.0f), y(0.0f), z(0.0f) {}

Poly::Poly (int neigh, float x, float y, float z) : neigh(neigh), x(x), y(y), z(z) {}

Poly::Poly(Poly&& other) noexcept : neigh(other.neigh), x(other.x), y(other.y), z(other.z) {}

Poly& Poly::operator=(const Poly& other) {
    if (this != &other) {
        neigh = other.neigh;
        x = other.x;
        y = other.y;
        z = other.z;
    }
    return *this;
}

Poly& Poly::operator=(Poly&& other) noexcept {
    if (this != &other) {
        neigh = other.neigh;
        x = other.x;
        y = other.y;
        z = other.z;
        // 将被移动的对象置为有效但未定义的状态
        other.x = other.y = other.z = 0.0f;
        other.neigh = 0;
    }
    return *this;
}

int Poly::getNeigh() const {
    return this->neigh;
}

float Poly::getx() const {
    return x;
}

float Poly::gety() const {
    return y;
}

float Poly::getz() const {
    return z;
}

void Poly::setNeigh(int new_neigh) {
    this->neigh = new_neigh;
}

void Poly::setx(float new_x) {
    this->x = new_x;
}

void Poly::sety(float new_y) {
    this->y = new_y;
}

void Poly::setz(float new_z) {
    this->z = new_z;
}

void Poly::set(float new_x, float new_y, float new_z) {
    this->x = new_x;
    this->y = new_y;
    this->z = new_z;
}

void Poly::translation(float dx, float dy, float dz) {
    this->x += dx;
    this->y += dy;
    this->z += dz;
}

void Poly::moveTo(float nx, float ny, float nz) {
    this->x = nx;
    this->y = ny;
    this->z = nz;
}

void Poly::moveTo(const Position& pos) {
    moveTo(pos.x(), pos.y(), pos.z());
}

void Poly::rotate(const Eigen::Matrix3f& rod, const Position& ver) {
    Vector vec;
    vec << x - ver[0], y - ver[1], z - ver[2];
    auto vec_rot = rod * vec;
    x = vec_rot[0] + ver[0];
    y = vec_rot[1] + ver[1];
    z = vec_rot[2] + ver[2];
}

Position Poly::getPosition() const {
    return {x, y, z};
}

std::ostream& operator<<(std::ostream& os, const Poly& poly) {
    os << "Poly(neigh=" << poly.neigh
    << ", x=" << poly.x
    << ", y=" << poly.y
    << ", z=" << poly.z << ")";
    return os;
}

Edge::Edge(int to, std::string type) : to(to), type(std::move(type)) {}

Edge::Edge(Edge&& e)  noexcept {
    to = e.to;
    type = std::move(e.type);
}

Edge& Edge::operator=(const Edge& e) {
    to = e.to;
    type = e.type;
    return *this;
}

Edge& Edge::operator=(Edge&& e)  noexcept {
    to = e.to;
    e.to = -1;
    type = std::move(e.type);
    e.type.clear();
    return *this;
}

int Edge::getTo() const {
    return to;
}

void Edge::setTo(int new_to) {
    to = new_to;
}

const std::string& Edge::getType() {
    return type;
}

void Edge::setType(std::string new_type) {
    type = std::move(new_type);
}

std::ostream& operator<<(std::ostream& os, const Edge& e) {
    os << "Edge(to = " << e.to << " "
    << "type = " << e.type << ")";
    return os;
}

float atom::positionDistance(const Position& p1, const Position& p2) {
    return (p1 - p2).norm();
}

float atom::positionDistanceSquared(const Position& p1, const Position& p2) {
    return (p1 - p2).squaredNorm();
}

Vector atom::positionMinusPosition(const Position& p1, const Position& p2) {
    return p1 - p2;
}
