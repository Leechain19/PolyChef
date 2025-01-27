#include "point.hpp"
#include "aabb.hpp"
#include <memory>
#include <iostream>

class KDTree {
public:
    Point *point;
    int alignment = 0;
    std::unique_ptr<KDTree> upper_child; // 使用 unique_ptr 管理子树
    std::unique_ptr<KDTree> lower_child; // 使用 unique_ptr 管理子树
    AABB aabb;

    KDTree() : point(nullptr) {}
    explicit KDTree(int alignment) : point(nullptr), alignment(alignment) {}
    explicit KDTree(AABB aabb) : point(nullptr), aabb(std::move(aabb)) {}
    KDTree(Point *p, int alignment, AABB aabb) : point(p), alignment(alignment), aabb(std::move(aabb)) {}
    void insert(Point *p);
    Point *find_point(const Point *p) const;
    Point *find_nearest(const Point *p) const;
};

std::ostream &operator<<(std::ostream &os, const KDTree &tree);
