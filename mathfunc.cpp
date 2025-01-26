//
// Created by AnthonyZhang on 2025/1/11.
//

#include "mathfunc.h"
#include <cmath>

Eigen::Matrix3f rodrigues(const Vector &K, float theta) {
    Vector norm_K = K.normalized();
    float k1 = norm_K(0), k2 = norm_K(1), k3 = norm_K(2);
    Eigen::RowVector3f norm_KT = norm_K.transpose();
    float cos_theta = std::cos(theta);
    float sin_theta = std::sin(theta);
    float one_minus_cos_theta = 1.0f - cos_theta;

    Eigen::Matrix3f R1, R2, R3;
    R1 << 1.0f, 0.0f, 0.0f,
          0.0f, 1.0f, 0.0f,
          0.0f, 0.0f, 1.0f;

    R2 = norm_K * norm_KT;

    R3 << 0.0f, -k3, k2,
          k3, 0.0f, -k1,
          -k2, k1, 0.0f;

    return R1 * cos_theta + R2 * one_minus_cos_theta + R3 * sin_theta;
}

Eigen::Matrix3f rodriguesWithVector(const Vector &K, float theta, const Vector &vec1, const Vector &vec2) {
    auto R = rodrigues(K, theta);
    auto norm_vec1 = vec1.normalized();
    auto norm_vec2 = vec2.normalized();
    auto norm_vec1_rot = R * norm_vec1;

    for (int i = 0; i < 3; i ++) {
        if (abs(norm_vec1_rot[i] - norm_vec2[i]) > mathfunc::eps) {
            return rodrigues(K, -theta);
        }
    }
    return R;
}

float getAngleFromVector(const Vector &vec1, const Vector &vec2) {
    auto norm_vec1 = vec1.normalized();
    auto norm_vec2 = vec2.normalized();
    //    zero angle check
    {
        auto diff = norm_vec1 - norm_vec2;
        if (abs(diff(0)) < mathfunc::eps && abs(diff(1)) < mathfunc::eps && abs(diff(2)) < mathfunc::eps) {
            return 0.0f;
        }
    }
    // 180 angle check
    {
        auto diff = norm_vec1 + norm_vec2;
        if (abs(diff(0)) < mathfunc::eps && abs(diff(1)) < mathfunc::eps && abs(diff(2)) < mathfunc::eps) {
            return M_PI;
        }
    }

    return static_cast<float>(acos(norm_vec1.dot(norm_vec2)));
}

Eigen::Matrix3f rotateMatrix(const Vector &vec1, const Vector &vec2) {
    auto theta = getAngleFromVector(vec1, vec2);
    auto norm_vec1 = vec1.normalized();
    auto norm_vec2 = vec2.normalized();

    // check theta nearly equals 0/pi
    Vector K;
    if (abs(theta) < mathfunc::eps || abs(theta - M_PI) < mathfunc::eps) {
        K << norm_vec1[1], -norm_vec1[0], 0.0;
    }
    else {
        K = norm_vec1.cross(norm_vec2);
    }

    auto R = rodrigues(K, theta);
    auto norm_vec1_rot = R * norm_vec1;

    for (int i = 0; i < 3; i ++) {
        if (abs(norm_vec1_rot[i] - norm_vec2[i]) > mathfunc::eps) {
            return rodrigues(K, -theta);
        }
    }
    return R;
}

// 找到平面法向量
Vector getNormalVectorFromTriplePoints(const Vector &p1, const Vector &p2, const Vector &p3) {
    auto vec1 = p2 - p1;
    auto vec2 = p3 - p1;
    return vec1.cross(vec2);
}