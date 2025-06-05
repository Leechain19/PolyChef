# PolyChef

### Introduce

PolyChef is a toolkit for generating polymers based on reference curves. Use pSMILES and CSV files of reference curves as input. At present, it supports monomers, intrachain entanglement, interchain entanglement, dendritic polymers, hyperbranched polymers, cross-linked polymers, mechanically interlocked polymers, and semi interpenetrating polymers.

### Dependency

Cmake version 3.20+ | GCC 11.2 | C++17 | Eigen3 | Python 3.8 | RDKit 2020.09.01 | OpenMP

### Install

Dependent Libraries: Eigen / RDKit

It is recommended to install Eigen/RDKit/GCC/GXX/CMake in a conda virtual environment.

Use   `which g++`  ,   `which gcc`  , and   `which cmake`   to check if the compilation environment is correctly set up.

#### Make

```bash
mkdir build
cd build
cmake ..
make -j
```

### Example input file

```json
{
    // chain_curve_list 聚合物链曲线 为一个list
    "chain_curve_list": [
        "./3Dflower.csv"
    ],
    // 最大聚合度 如果不满足会提前停止
    "chain_max_polymerization_degree": 1000,
    // 每一条单链的组成 要求长度和chain_curve_list相等
    "chain_psmiles_list": [
        [
            "[*]CC[*]"
        ]
    ],
    // 聚合物 交联点的mol2 list
    "crosslinker_mol2_list": [
        "../inputs/newmol.mol2"
    ],
    // 交联曲线list
    "crosslinking_curves": [
        "../inputs/testCurve.csv"
    ],
    // 交联network 长度和交联曲线list相等
    "crosslinking_network": [
        [
            0,
            0,
            0,
            1
        ]
    ],

    // 备注信息 (optional)
    "file_info": "NO_INFO",
    // 优化原子长度 一般为 1-10
    "optimize_size": 1,
    // 输出文件夹
    "output_directory_path": "./",
    // 参数A (optional)
    "para_A": 1.0,
    // 参数B (optional)
    "para_B": 1.0,
    // 窗口大小 (optional)
    "window_distance": 5.0,
    // 聚合物类型 chain or crosslink
    "polymer_type": "chain",
    // 池化方案 ave or max or min
    "pool_choice": "ave",
    // 随机聚合选项 true or false
    "random_polymerization": false,
    // thread 默认是8
    "threads": 16,
    // 是否详细输出 (optional)
    "verbose": false,
    // 障碍物 (optional)
    "obstacle_list": []
}
```

#### Compile dynamic library command:

g++ -shared -fPIC -I/usr/local/include/eigen3 -O3 -march=native -o libcustomFunction.so customFunction.cpp

#### Run command

```bash
polychef run config.json
polychef config config.json [options]
```
