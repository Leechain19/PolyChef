# Atom-Search-cpp

### 介绍
cpp version for Atom-Search (PolyChef)

### 软件依赖

Cmake version 3.20+ | GCC 11.2 | C++17 | Eigen3 | Python 3.8 | RDKit 2020.09.01 | OpenMP

### 安装教程

TBD


### 使用说明

#### 编译动态库命令:
g++ -shared -fPIC -I/usr/local/include/eigen3 -O3 -march=native -o libcustomFunction.so customFunction.cpp

#### 运行指令
./poly run

./poly config [options]

### 参与贡献

1.  Fork 本仓库
2.  新建 Feat_xxx 分支
3.  提交代码
4.  新建 Pull Request
