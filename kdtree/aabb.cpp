#include "aabb.hpp"

#include <iostream>
#include <cmath>

AABB::AABB() : lower(), upper() {}

AABB::AABB(float min_x, float max_x, float min_y, float max_y, float min_z, float max_z) :
    lower(min_x, min_y, min_z), upper(max_x, max_y, max_z) {}

AABB::AABB(Point lower, Point upper) : lower(std::move(lower)), upper(std::move(upper)) {}

AABB::AABB(const std::vector<Point*>& points) {
    if (points.empty()) {
        lower = Point();
        upper = Point();
        return;
    }

    float min_x = std::numeric_limits<float>::max();
    float max_x = std::numeric_limits<float>::lowest();
    float min_y = std::numeric_limits<float>::max();
    float max_y = std::numeric_limits<float>::lowest();
    float min_z = std::numeric_limits<float>::max();
    float max_z = std::numeric_limits<float>::lowest();

    for (const auto& p : points) {
        if (!p) continue;
        min_x = p->x() < min_x ? p->x() : min_x;
        max_x = p->x() > max_x ? p->x() : max_x;
        min_y = p->y() < min_y ? p->y() : min_y;
        max_y = p->y() > max_y ? p->y() : max_y;
        min_z = p->z() < min_z ? p->z() : min_z;
        max_z = p->z() > max_z ? p->z() : max_z;
    }

    lower = Point(min_x, min_y, min_z);
    upper = Point(max_x, max_y, max_z);
}

AABB::AABB(const AABB& parent, const Point *split, int alignment, bool is_upper) : lower(parent.lower), upper(parent.upper) {
    if (!split) return;
    if (alignment == 0) {
        if (is_upper) lower.setx(split->x());
        else upper.setx(split->x());
        return;
    }
    if (alignment == 1) {
        if (is_upper) lower.sety(split->y());
        else upper.sety(split->y());
        return;
    }
    if (alignment == 2) {
        if (is_upper) lower.setz(split->z());
        else upper.setz(split->z());
    }
}

bool AABB::includes(const Point *p) const {
    if (!p) return false;
    return p->x() >= lower.x() && p->x() <= upper.x() &&
    p->y() >= lower.y() && p->y() <= upper.y() &&
    p->z() >= lower.z() && p->z() <= upper.z();
}

float AABB::distance_outside(const Point *p) const {
    if (!p) return std::numeric_limits<float>::max();
    float x_dist = p->x() >= lower.x() && p->x() <= upper.x() ? 0 : std::min(std::abs(p->x() - lower.x()), std::abs(p->x() - upper.x()));
    float y_dist = p->y() >= lower.y() && p->y() <= upper.y() ? 0 : std::min(std::abs(p->y() - lower.y()), std::abs(p->y() - upper.y()));
    float z_dist = p->z() >= lower.z() && p->z() <= upper.z() ? 0 : std::min(std::abs(p->z() - lower.z()), std::abs(p->z() - upper.z()));
    return std::sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist);
}

float AABB::distance_inside(const Point *p) const {
    if (!p) return std::numeric_limits<float>::max();
    float x_dist = std::min(std::abs(p->x() - lower.x()), std::abs(p->x() - upper.x()));
    float y_dist = std::min(std::abs(p->y() - lower.y()), std::abs(p->y() - upper.y()));
    float z_dist = std::min(std::abs(p->z() - lower.z()), std::abs(p->z() - upper.z()));
    return std::sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist);
}

std::ostream &operator<<(std::ostream &os, const AABB &aabb) {
    os << "AABB " << aabb.lower << "-" << aabb.upper;
    return os;
}
