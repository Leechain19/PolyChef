#include "chemtable.h"
#include "nlohmann/json.hpp"
#include "graph.h"
#include "grid.h"
#include "chemio.h"
#include "spreading.h"
#include "cxxopts.hpp"
#include "optimizer.h"
#include "exception.h"
#include <filesystem>
#include <omp.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <ctime>

#define STRINGIFY(x) (#x)

// debug
std::string to_string(const std::string& s) { return '"' + s + '"'; }
std::string to_string(const char *s) { return to_string((std::string) s); }
std::string to_string(bool b) { return (b ? "true" : "false"); }
template<typename A, typename B>
std::string to_string(std::pair<A, B> p) { return "(" + to_string(p.first) + ", " + to_string(p.second) + ")"; }
template<typename A>
std::string to_string(A v) { bool first = true; std::string res = "{"; for(const auto &x : v) { if(!first) { res += ", "; } first = false; res += to_string(x);} res += "}"; return res; }
void debug_out() { std::cout << std::endl; }
template<typename Head, typename... Tail> void debug_out(Head H, Tail... T) { std::cout << " " << to_string(H); debug_out(T...);}
#define dbg(x) std::cout << '[' << (#x) << ']' << (x) << std::endl
#define debug(...) std::cout << "[" << #__VA_ARGS__ << "]:", debug_out(__VA_ARGS__)
// debug finish

using json = nlohmann::json;
constexpr int inf = 1e9;
using Options = cxxopts::Options;

bool startsWith(const std::string& s, const std::string& prefix) {
    if ((int)s.size() < (int)prefix.size()) return false;
    for (int i = 0; i < (int)prefix.size(); i ++) {
        if (s[i] != prefix[i]) return false;
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

template<typename T>
void readFromJSON(const json& j, T& val, const std::string& name) {
    if (j.count(name)) {
        val = j[name].get<T>();
        return;
    }
    val = T{};
}

template<typename T>
void readVectorFromJSON(const json& j, std::vector<T>& vec, const std::string& name) {
    if (j.count(name) && j[name].is_array()) {
        for (const auto& x : j[name]) {
            vec.push_back(x.get<T>());
        }
    }
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

namespace fs = std::filesystem;
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

std::vector<Position> getScatterFromCSV(const std::string& path, bool header = true) {
    std::vector<Position> points;
    std::ifstream file(path);
    std::string line;
    std::string token;

    // 跳过表头
    if (header)
        std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        float x, y, z;
        std::getline(ss, token, ','); x = static_cast<float>(std::stod(token));
        std::getline(ss, token, ','); y = static_cast<float>(std::stod(token));
        std::getline(ss, token, ','); z = static_cast<float>(std::stod(token));
        points.emplace_back(x, y, z);
    }

    return points;
}

void buildSystemWithScatter(const std::string& path, int string_number, const std::vector<std::vector<Position>>& target_points_list, const std::vector<std::shared_ptr<Graph>>& sequence,
                            int degree_of_polymerization, bool random_polymerization = false, int optimize_size = 1, bool verbose = false, const std::string& file_info = "NO_INFO") {
    std::cout << "Build System..." << std::endl;
    auto tree = std::make_shared<Grid>(5.0);

    for (int i = 0; i < string_number; i ++) {
        const std::vector<Position>& target_points = target_points_list[i];
        std::cout << std::string(50, '=') << std::endl;
        std::cout << "String ID: " << i << std::endl;

        std::string file_path = path + "/string_" + std::to_string(i) + ".mol2";
        std::cout << "Path: " << file_path << std::endl;

        auto g_ptr = std::make_shared<Graph>();

        curveSpreading(target_points, g_ptr, tree, sequence, degree_of_polymerization, 5.0f, 5, random_polymerization, optimize_size, verbose);

        g_ptr->makeEnd("H");
        chemio::writeMol2File(file_path, g_ptr, file_info);
        std::cout << std::string(50, '=') << std::endl;
    }
}


void buildChainBasedCurve(const std::string& output_path, const std::vector<std::string>& smiles_list, const std::vector<std::string>& chain_curve_list, int degree_of_polymerization = 5000,
                          bool random_polymerization = true, int optimize_size = 1, bool verbose = false, const std::string& file_info = "NO_INFO") {
    std::vector<std::shared_ptr<Graph>> sequence;
    for (const auto& smiles : smiles_list) {
        sequence.emplace_back(chemio::buildGraphFromPSmiles(smiles));
//        std::cout << "monos: " << *(sequence.back()) << std::endl;
    }

    auto now_time = getCurrentTimeAsString();
    std::cout << "Time: " << now_time << std::endl;

    std::string system = connectStringVectorAsOne(smiles_list);
    std::cout << "System: " << system << std::endl;

    auto saving_dir = output_path;
    if (saving_dir.empty() or saving_dir.back() != '/') saving_dir.push_back('/');
    saving_dir += now_time;

    createDirectory(saving_dir);

    std::vector<std::vector<Position>> target_points_list;
    for (const auto& path : chain_curve_list) {
        target_points_list.emplace_back(getScatterFromCSV(path));
    }

    int string_number = (int)chain_curve_list.size();

    buildSystemWithScatter(saving_dir, string_number, target_points_list, sequence, degree_of_polymerization, random_polymerization, optimize_size, verbose, file_info);
    std::cout << "Finish!" << std::endl;
}

void InitializeOpenMP(int target_thread_num) {
    std::cout << "Threads: " << target_thread_num << std::endl;
//    omp_set_num_threads(target_thread_num);
    Eigen::setNbThreads(target_thread_num);
}

// 读取配置
bool read_config(json& config, const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Config file not found. Please run 'config' command first." << std::endl;
        return false;
    }
    file >> config;
    file.close();
    return true;
}

// 保存配置
bool save_config(const json& config, const std::string& path) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "Failed to save config." << std::endl;
        return false;
    }
    file << config.dump(4); // 保存为格式化的 JSON
    file.close();
    std::cout << "Success! Config saved to " << path << std::endl;
    return true;
}

void solve(const std::string& filename) {
    json j;
    if (!read_config(j, filename)) {
        return;
    }

    std::string polymer_type;
    readFromJSON(j, polymer_type, STRINGIFY(polymer_type));
    polymer_type = lower(polymer_type);

    unsigned int seed = 114514;
    readFromJSON(j, seed, STRINGIFY(seed));

    int threads = 1;
    readFromJSON(j, threads, STRINGIFY(threads));

    // OpenMP threads setting
    InitializeOpenMP(threads);

    std::vector<int> box_size;
    readVectorFromJSON(j, box_size, STRINGIFY(box_size));

    std::string output_directory_path;
    readFromJSON(j, output_directory_path, STRINGIFY(output_directory_path));

    std::vector<std::string> chain_curve_list;
    readVectorFromJSON(j, chain_curve_list, STRINGIFY(chain_curve_list));

    std::vector<std::string> chain_psmiles_list;
    readVectorFromJSON(j, chain_psmiles_list, STRINGIFY(chain_psmiles_list));

    bool random_polymerization = false;
    readFromJSON(j, random_polymerization, STRINGIFY(random_polymerization));

    int chain_max_polymerization_degree;
    readFromJSON(j, chain_max_polymerization_degree, STRINGIFY(chain_max_polymerization_degree));

    int optimize_size = 5000;
    readFromJSON(j, optimize_size, STRINGIFY(optimize_size));

    std::string file_info;
    readFromJSON(j, file_info, STRINGIFY(file_info));

    bool verbose = false;
    readFromJSON(j, verbose, STRINGIFY(verbose));

    if (startsWith(polymer_type, "cross")) {
        // type: crosslink polymer
        std::vector<std::string> crosslinks_psmiles_list;
        readVectorFromJSON(j, crosslinks_psmiles_list, STRINGIFY(crosslinks_psmiles_list));

        int crosslinks_numbers;
        readFromJSON(j, crosslinks_numbers, STRINGIFY(crosslinks_numbers));
    }

    else if (startsWith(polymer_type, "chain")) {
        buildChainBasedCurve(output_directory_path, chain_psmiles_list, chain_curve_list, chain_max_polymerization_degree, random_polymerization, optimize_size, verbose, file_info);
    }
}

void test() {
    std::string smi = "[*]CC[*]";
    auto mono_ptr = chemio::buildGraphFromPSmiles(smi);
    auto start_mol = std::make_shared<Graph>();

    auto cur_mol = std::make_shared<Graph>(*mono_ptr);
    cur_mol->translation(10,10,10);
}

void config(int argc, char* argv[], const std::string& config_filename) {
    Options options("polychef config", "Configure settings");
    options.add_options()
    ("t,threads", "Threads", cxxopts::value<int>())
    ("o,output", "Output directory path", cxxopts::value<std::string>())
    ("v,verbose", "Verbose", cxxopts::value<bool>())
    ("i,info", "File info", cxxopts::value<std::string>())
    ("h,help", "Show help");

    auto result = options.parse(argc - 1, argv + 1);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return;
    }

    json config;
    read_config(config, config_filename);

    if (result.count("threads")) {
        config["threads"] = result["threads"].as<int>();
        std::cout << "Threads number: " << config["threads"] << std::endl;
    }

    if (result.count("output")) {
        config["output_directory_path"] = result["output"].as<std::string>();
        std::cout << "Output directory path has been changed to: " << config["output_directory_path"] << std::endl;
    }

    if (result.count("verbose")) {
        config["verbose"] = result["verbose"].as<bool>();
        std::cout << "Verbose mode " << (config["verbose"] ? "enabled" : "disabled") << "." << std::endl;
    }

    if (result.count("info")) {
        config["file_info"] = result["info"].as<std::string>();
        std::cout << "File Info: " << config["file_info"] << std::endl;
    }

    if (save_config(config, config_filename)) {
        std::cout << "Configure settings Finish!" << std::endl;
    }
}

void get_struct() {

}

int main(int argc, char* argv[]) {
    std::string cwd = std::filesystem::current_path();
    std::cout << "Current working directory: " << cwd << std::endl;
    std::string config_filename = cwd + "/../config.json";

    if (argc < 2) {
        std::cerr << "Usage: polychef <command> [options]\n"
        << "Commands:\n"
        << "  run     Run a task\n"
        << ("  config  Configure settings in " + config_filename + "\n")
        << "  struct  Generate a 3D structure for PSmiles\n";
        return 1;
    }

    std::string command = argv[1];

    if (command == "run") {
        auto start = std::chrono::high_resolution_clock::now(); // start time
        solve(config_filename);
        auto end = std::chrono::high_resolution_clock::now(); // end time
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Elapsed time: " << duration.count() << " ms" << std::endl;
    }
    else if (command == "config") {
        config(argc, argv, config_filename);
    }
    else {
        // 抛出异常
        throw exception::InvalidParameterException("Unknown command: " + command);
    }

    return 0;
}