#pragma once

#include "point.hpp"
#include <vector>
#include <memory>
#include <limits>

struct AABB {
public:
    Point lower;
    Point upper;
    AABB();
    AABB(float min_x, float max_x, float min_y, float max_y, float min_z, float max_z);
    explicit AABB(const std::vector<Point*>& points);
    AABB(Point lower, Point upper);
    AABB(const AABB& parent, const Point *split, int alignment, bool upper);

    [[nodiscard]] bool includes(const Point *p) const;
    [[nodiscard]] float distance_outside(const Point *p) const;
    [[nodiscard]] float distance_inside(const Point *p) const;
};

std::ostream &operator<<(std::ostream &os, const AABB &aabb);
