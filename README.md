# Atom-Search-cpp

### 介绍
cpp version for Atom-Search (PolyChef)

### 软件依赖

Cmake version 3.20+ | GCC 11.2 | C++17 | Eigen3 | Python 3.8 | RDKit 2020.09.01 | OpenMP

### 安装教程
依赖库: Eigen / RDKit

推荐在conda虚拟环境中安装Eigen / RDKit / GCC / GXX / CMake

先使用`which g++` `which gcc` `which cmake` 判断编译环境正确

#### 编译命令
```bash
mkdir build
cd build
cmake ..
make -j6
```

### 配置文件说明

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





#### 编译动态库命令:
g++ -shared -fPIC -I/usr/local/include/eigen3 -O3 -march=native -o libcustomFunction.so customFunction.cpp

#### 运行指令
```bash
polychef run config.json
polychef config config.json [options]
```

### 参与贡献

1.  Fork 本仓库
2.  新建 Feat_xxx 分支
3.  提交代码
4.  新建 Pull Request
