#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <utility>

class Point {
private:
    Eigen::Vector3f vec;
public:
    // 默认构造函数
    Point() : vec(0.0f, 0.0f, 0.0f) {}

    // 带参数的构造函数
    Point(float x, float y, float z) : vec(x, y, z) {}

    // 从 Eigen::Vector3f 构造
    explicit Point(const Eigen::Vector3f &vec) : vec(vec) {}

    // 转换为 Eigen::Vector3f
    [[nodiscard]] Eigen::Vector3f to_vector3f() const {
        return vec;
    }

    [[nodiscard]] float x() const {
        return vec.x();
    }

    [[nodiscard]] float y() const {
        return vec.y();
    }

    [[nodiscard]] float z() const {
        return vec.z();
    }

    void setx(float new_x) {
        vec(0) = new_x;
    }

    void sety(float new_y) {
        vec(1) = new_y;
    }

    void setz(float new_z) {
        vec(2) = new_z;
    }

    // 运算符重载
    Point operator+(const Point& other) const {
        return Point(vec + other.vec);
    }

    Point operator-(const Point& other) const {
        return Point(vec - other.vec);
    }

    Point operator*(float mul) const {
        return Point(vec * mul);
    }

    Point operator/(float divisor) const {
        return Point(vec / divisor);
    }

    bool operator==(const Point& other) const {
        return vec == other.vec;
    }

    // 其他功能
    [[nodiscard]] float length() const {
        return vec.norm();
    }

    [[nodiscard]] Point x_component() const {
        return {vec.x(), 0.0f, 0.0f};
    }

    [[nodiscard]] Point y_component() const {
        return {0.0f, vec.y(), 0.0f};
    }

    [[nodiscard]] Point z_component() const {
        return {0.0f, 0.0f, vec.z()};
    }

    // 输出运算符重载
    friend std::ostream& operator<<(std::ostream& os, const Point& p) {
        os << "(" << p.vec.x() << ", " << p.vec.y() << ", " << p.vec.z() << ")";
        return os;
    }
};

std::ostream &operator<<(std::ostream &os, const Point &p);
