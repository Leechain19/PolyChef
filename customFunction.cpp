//
// Created by AnthonyZhang on 2025/2/13.
//

#include <Eigen/Dense>

extern "C" {
    bool inZ(const Eigen::Vector3f& point) {
        if (point[0] <= 40 && point[0] >= 10 && point[1] <= 45 && point[1] >= 5) {
            if (10<=point[1] && point[1] <=40 && 1.16*point[0]-point[1]<=0) return false;
            if (10<=point[1] && point[1] <=40 && 1.16*point[0]-point[1]>=10) return false;
            return true;
        }
        return false;
    }

    bool inJ(const Eigen::Vector3f& point) {
        if (point[0] <= 35 && point[0] >= 10 && point[1] <= 45 && point[1] >= 5) {
            if (point[1]<=15 && (point[0]-25)*(point[0]-25) + (point[1]-15)*(point[1]-15) <=25) return false;
            if (point[1]<=15 && (point[0]-25)*(point[0]-25) + (point[1]-15)*(point[1]-15) >=101) return false;
            if (point[0] <= 30 && point[1] >=15) return false;
            return true;
        }
        return false;
    }

    bool inU(const Eigen::Vector3f& point) {
        if (point[0] <= 40 && point[0] >= 10 && point[1] <= 45 && point[1] >= 5) {
            if (point[1]<=20 && (point[0]-25)*(point[0]-25) + (point[1]-20)*(point[1]-20) <=110) return false;
            if (point[1]<=20 && (point[0]-25)*(point[0]-25) + (point[1]-20)*(point[1]-20) >=235) return false;
            if (15<= point[0] && point[0] <= 35 && 20<=point[1] && point[1]<=45) return false;
            return true;
        }
        return false;
    }

    bool customFunction(const Eigen::Vector3f& point) {
        return inZ(point);
    }

    int xmin() { return 0;}

    int xmax() { return 50; }

    int ymin() { return 0; }

    int ymax() { return 50; }

    int zmin() { return 0; }

    int zmax() { return 20; }

}

// 编译命令: g++ -shared -fPIC -I/usr/local/include/eigen3 -O3 -march=native -o libcustomFunction.so customFunction.cpp