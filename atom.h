//
// Created by AnthonyZhang on 2025/1/11.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_ATOM_H
#define ATOM_SEARCH_CPP_ATOM_H

#include "Eigen/Dense"
#include <string>
#include <utility>
#include <iostream>
#include <array>
#include <cmath>

using Position = Eigen::Vector3f;
using Vector = Eigen::Vector3f;

class Atom {
private:
    std::string symbol;
    float x, y, z;
public:
    Atom();
    Atom(std::string symbol, float x, float y, float z);
    Atom(std::string symbol, const Position& pos);
    Atom(const Atom &o) = default;
    Atom(Atom&& o) noexcept ;
    ~Atom() = default;

    Atom& operator=(const Atom& other);

    Atom& operator=(Atom&& other) noexcept ;

    [[nodiscard]] std::string getSymbol() const ;

    [[nodiscard]] float getx() const ;

    [[nodiscard]] float gety() const ;

    [[nodiscard]] float getz() const ;

    void setSymbol(const std::string& new_symbol);

    void setx(float new_x);

    void sety(float new_y);

    void setz(float new_z);

    void set(float new_x, float new_y, float new_z);

    void translation(float dx, float dy, float dz);

    void moveTo(float nx, float ny, float nz);

    void moveTo(const Position& pos);

    void rotate(const Eigen::Matrix3f& rod, const Position& ver);

    [[nodiscard]] Position getPosition() const;

    friend std::ostream& operator<<(std::ostream& os, const Atom& atom);
};

class Poly {
private:
    int neigh;
    float x, y, z;

public:
    Poly();
    Poly (int neigh, float x, float y, float z);
    Poly(const Poly& other) = default;
    Poly(Poly&& other) noexcept ;
    ~Poly() = default;

    Poly& operator=(const Poly& other);

    Poly& operator=(Poly&& other) noexcept;

    [[nodiscard]] int getNeigh() const;

    [[nodiscard]] float getx() const;

    [[nodiscard]] float gety() const;

    [[nodiscard]] float getz() const;

    void setNeigh(int new_neigh);

    void setx(float new_x);

    void sety(float new_y);

    void setz(float new_z);

    void set(float new_x, float new_y, float new_z);

    void translation(float dx, float dy, float dz);

    void moveTo(float nx, float ny, float nz);

    void moveTo(const Position& pos);

    void rotate(const Eigen::Matrix3f& rod, const Position& ver);

    [[nodiscard]] Position getPosition() const;

    friend std::ostream& operator<<(std::ostream& os, const Poly& poly);
};


class Edge {

private:
    int to;
    std::string type;

public:
    Edge(int to, std::string type);
    ~Edge() = default;
    Edge(const Edge& e) = default;
    Edge(Edge&& e)  noexcept ;

    Edge& operator=(const Edge& e);

    Edge& operator=(Edge&& e)  noexcept;

    [[nodiscard]] int getTo() const;

    void setTo(int new_to);

    const std::string& getType();

    void setType(std::string new_type);

    friend std::ostream& operator<<(std::ostream& os, const Edge& e);
};

namespace atom {
    float positionDistance(const Position& p1, const Position& p2);
    float positionDistanceSquared(const Position& p1, const Position& p2);
    Vector positionMinusPosition(const Position& p1, const Position& p2);
}

#endif //ATOM_SEARCH_CPP_ATOM_H
