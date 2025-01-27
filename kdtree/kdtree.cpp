#include "kdtree.hpp"

void KDTree::insert(Point *p) {
    if (!point) {
        point = p;
        aabb = AABB(*p, *p);
        return;
    }
    if (*p == *point) return;

    bool is_upper = alignment == 0 ? (p->x() > point->x()) : (alignment == 1 ? (p->y() > point->y()) : (p->z() > point->z()));
    std::unique_ptr<KDTree>& correct_child = is_upper ? (upper_child) : (lower_child);
    if (!correct_child) {
        int child_alignment = (alignment + 1) % 3;
        correct_child = std::make_unique<KDTree>(p, child_alignment, AABB(aabb, point, child_alignment, is_upper));
    }
    else {
        correct_child->insert(p);
    }
}

Point *KDTree::find_point(const Point *p) const {
    if (!point) return nullptr;

    if (*p == *point) {
        return point;
    }

    bool is_upper = alignment == 0 ? (p->x() > point->x()) : (alignment == 1 ? (p->y() > point->y()) : (p->z() > point->z()));
    const std::unique_ptr<KDTree> &correct_child = is_upper ? upper_child : lower_child;
    return correct_child ? correct_child->find_point(p) : nullptr;
}

Point *KDTree::find_nearest(const Point *p) const {
    if (!point) return nullptr;

    // Set own point as closest
    Point *closest_p = point;
    float closest_distance = (*p - *point).length();

    const std::unique_ptr<KDTree> &primary_child = (alignment == 0 ? (p->x() > point->x()) : (alignment == 1 ? (p->y() > point->y()) : (p->z() > point->z()))) ? upper_child : lower_child;

    if (primary_child && primary_child->aabb.includes(p)) {
        Point *new_closest = primary_child->find_nearest(p);
        if (new_closest) {
            float new_distance = (*p - *new_closest).length();
            if (new_distance < closest_distance) {
                closest_p = new_closest;
                closest_distance = new_distance;
            }
        }
    }

    const std::unique_ptr<KDTree> &secondary_child = (primary_child == upper_child) ? lower_child : upper_child;
    if (secondary_child && secondary_child->aabb.distance_outside(p) < closest_distance) {
        Point *new_closest = secondary_child->find_nearest(p);
        if (new_closest) {
            float new_distance = (*p - *new_closest).length();
            if (new_distance < closest_distance) {
                closest_p = new_closest;
                closest_distance = new_distance;
            }
        }
    }
    return closest_p;
}

std::ostream &operator<<(std::ostream &os, const KDTree &tree) {

    os << "KDTree: (" << tree.aabb << ") with alignment " << tree.alignment;
    if (!tree.point) {
        os << " without point";
    }
    else {
        os << " and point " << *tree.point;
    }

    if (tree.lower_child) {
        os << "\nLower Child: " << *tree.lower_child;
    }
    if (tree.upper_child) {
        os << "\nUpper Child: " << *tree.upper_child;
    }

    return os;
}
