//
// Created by AnthonyZhang on 2025/1/11.
//
#pragma once

#ifndef ATOM_SEARCH_CPP_MATHFUNC_H
#define ATOM_SEARCH_CPP_MATHFUNC_H

#include <Eigen/Dense>
#include "atom.h"

namespace mathfunc {
    constexpr float eps = 1e-5;
}

Eigen::Matrix3f rodrigues(const Vector &K, float theta);

Eigen::Matrix3f rodriguesWithVector(const Vector &K, float theta, const Vector &vec1, const Vector &vec2);

float getAngleFromVector(const Vector &vec1, const Vector &vec2);

Eigen::Matrix3f rotateMatrix(const Vector &vec1, const Vector &vec2);

// 找到平面法向量
Vector getNormalVectorFromTriplePoints(const Vector &p1, const Vector &p2, const Vector &p3);

#endif //ATOM_SEARCH_CPP_MATHFUNC_H