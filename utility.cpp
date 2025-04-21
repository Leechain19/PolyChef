//
// Created by AnthonyZhang on 2025/4/19.
//

#include "utility.h"

bool fileExists(const std::string& filename) {
    return std::filesystem::exists(filename);
}

bool startsWith(const std::string& s, const std::string& prefix) {
    if ((int)s.size() < (int)prefix.size()) return false;
    for (int i = 0; i < (int)prefix.size(); i ++) {
        if (s[i] != prefix[i]) return false;
    }
    return true;
}

bool endsWith(const std::string& s, const std::string& suffix) {
    if ((int)s.size() < (int)suffix.size()) return false;
    int n = (int)s.size(), m = (int)suffix.size();

    for (int i = 0; i < m; i ++) {
        if (s[n-i-1] != suffix[m-i-1]) return false;
    }
    return true;
}

std::string lower(const std::string& s) {
    std::string ret;
    for (auto c : s) {
        if (c >= 'A' and c <= 'Z') c = static_cast<char>(c + 'a' - 'A');
        ret.push_back(c);
    }
    return ret;
}

std::string upper(const std::string& s) {
    std::string ret;
    for (auto c : s) {
        if (c >= 'a' and c <= 'z') c = static_cast<char>(c + 'A' - 'a');
        ret.push_back(c);
    }
    return ret;
}

std::string getCurrentTimeAsString() {
    // 获取当前时间点
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);

    // 转换为tm结构体
    std::tm* ltm = std::localtime(&time_t_now);

    // 格式化为字符串
    char buffer[80];
    strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", ltm);
    return {buffer};
}

std::string connectStringVectorAsOne(const std::vector<std::string>& vec) {
    std::string ret;
    for (const auto& s : vec) {
        ret += s;
        ret.push_back('-');
    }
    if (!ret.empty()) ret.pop_back();
    return ret;
}

void createDirectory(const std::string& path) {
    try {
        if (fs::create_directory(path)) {
            std::cout << "Directory created successfully: " << path << std::endl;
        } else {
            std::cout << "Directory already exists or creation failed: " << path << std::endl;
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }
}

void InitializeOpenMP(int target_thread_num) {
    std::cout << "Threads: " << target_thread_num << std::endl;
    omp_set_num_threads(target_thread_num);
    //    eigen 不开启并行
    //    Eigen::setNbThreads(target_thread_num);
}

std::vector<std::string> my_split(const std::string &s, char delim) {
    std::vector<std::string> ret;
    std::string t;
    for (auto x : s) {
        if (x == delim) {
            if (!t.empty()) ret.push_back(t);
            t = "";
        }
        else t += x;
    }
    ret.push_back(t);
    return ret;
}

std::shared_ptr<Atom> getAtomFromMol2Line(const std::string& s) {
    auto v = my_split(s, ' ');
    assert(v.size() == 9);
    auto symbol = v[1];
    auto x = std::stof(v[2]), y = std::stof(v[3]), z = std::stof(v[4]);
    return std::make_shared<Atom>(symbol, x, y, z);
}

