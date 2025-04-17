#include "chemtable.h"
#include "nlohmann/json.hpp"
#include "graph.h"
#include "curve.h"
#include "grid.h"
#include "atom.h"
#include "loadlib.h"
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
#include "curve.h"

#define STRINGIFY(x) (#x)

using json = nlohmann::json;
constexpr int inf = 1e9;
using Options = cxxopts::Options;

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

void InitializeOpenMP(int target_thread_num) {
    std::cout << "Threads: " << target_thread_num << std::endl;
    omp_set_num_threads(target_thread_num);
//    eigen 不开启并行
//    Eigen::setNbThreads(target_thread_num);
}

// 读取配置
bool read_config(json& config, const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Config file not found. Please run 'config' command first." << std::endl;
        return false;
    }
    if (file.fail()) {
        throw std::runtime_error("config file is corrupted");
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

    int threads = 1;
    readFromJSON(j, threads, STRINGIFY(threads));

    // OpenMP threads setting
    InitializeOpenMP(threads);

    std::string output_directory_path;
    readFromJSON(j, output_directory_path, STRINGIFY(output_directory_path));

    std::vector<std::string> chain_curve_list;
    readVectorFromJSON(j, chain_curve_list, STRINGIFY(chain_curve_list));

    std::vector<std::vector<std::string>> chain_psmiles_list;
    readVectorFromJSON(j, chain_psmiles_list, STRINGIFY(chain_psmiles_list));

    bool random_polymerization = false;
    readFromJSON(j, random_polymerization, STRINGIFY(random_polymerization));

    int chain_max_polymerization_degree = 5000;
    readFromJSON(j, chain_max_polymerization_degree, STRINGIFY(chain_max_polymerization_degree));

    int optimize_size = 1;
    readFromJSON(j, optimize_size, STRINGIFY(optimize_size));

    std::string file_info;
    readFromJSON(j, file_info, STRINGIFY(file_info));

    std::string pool_choice;
    readFromJSON(j, pool_choice, STRINGIFY(pool_choice));

    float para_A = 1.0f, para_B = 1.0f;
    readFromJSON(j, para_A, STRINGIFY(para_A));
    readFromJSON(j, para_B, STRINGIFY(para_B));

    bool verbose = false;
    readFromJSON(j, verbose, STRINGIFY(verbose));

    // ==================

    auto now_time = getCurrentTimeAsString();
    std::cout << "Time: " << now_time << std::endl;

    // ===================

    auto saving_dir = output_directory_path;
    if (saving_dir.size() > 1 && saving_dir.back() == '/') saving_dir.pop_back();

    // chain type
    if (startsWith(polymer_type, "chain")) {

        if (chain_psmiles_list.size() != chain_curve_list.size()) {
            std::string err_info = std::string("Input Error: ") + "Chain Curve list size: " +
                    std::to_string(chain_psmiles_list.size()) + " Chain psmiles list size: " + std::to_string(chain_psmiles_list.size());
            throw std::runtime_error(err_info);
        }

        std::cout << "Build System [Chains] ..." << std::endl;
        auto tree = std::make_shared<Grid>();
        int string_number = chain_curve_list.size();
        std::cout << "Path: " << saving_dir << std::endl;

        for (int i = 0; i < string_number; i ++) {
            std::vector<Position> target_points = getScatterFromCSV(chain_curve_list[i]);
            std::vector<std::shared_ptr<Graph>> sequence;
            PsmilesBuilder<Graph> builder;
            for (const auto& smiles: chain_psmiles_list[i]) {
                auto ptr = builder.build(smiles);
                int mono_type = builder.getMonomerType(smiles);
                ptr->setMonoTypeAll(mono_type);
                sequence.emplace_back(std::move(ptr));
            }

            std::cout << std::string(50, '=') << std::endl;
            std::cout << "String ID: " << i << std::endl;

            std::string mol2_file_name = saving_dir + "/string_" + std::to_string(i) + ".mol2";
            std::string adj_file_name = saving_dir + "/string_" + std::to_string(i) + "_adj";
            std::string loss_file_name = saving_dir + "/string_" + std::to_string(i) + "_loss";

            constexpr double bad_signal_cost = 10000.0;

            auto g_ptr = std::make_shared<Graph>();
            auto loss_vector_ptr = std::make_unique<std::vector<std::pair<double, double>>>();
            curveSpreading(target_points, g_ptr, tree, sequence, chain_max_polymerization_degree, pool_choice, para_A, para_B, 5.0f, 5, random_polymerization, optimize_size,
                           verbose, bad_signal_cost, loss_vector_ptr);

            g_ptr->makeEnd("H");
            chemio::writeMol2File(mol2_file_name, adj_file_name, g_ptr, file_info);
            chemio::writeLoss2File(loss_file_name, loss_vector_ptr);
            std::cout << std::string(50, '=') << std::endl;
        }
        std::cout << "Finish!" << std::endl;
        return;
    }

    // crosslink type
    if (startsWith(polymer_type, "cross")) {

        std::cout << "Build System [CrossLinks] ..." << std::endl;
        std::cout << "Path: " << saving_dir << std::endl;
        // type: crosslink polymer
        std::vector<std::string> crosslinker_mol2_list;
        readVectorFromJSON(j, crosslinker_mol2_list, STRINGIFY(crosslinker_mol2_list));

        int crosslinker_number = (int)crosslinker_mol2_list.size();

        std::vector<std::array<int, 4>> crosslinking_network;
        readVectorFromJSON(j, crosslinking_network, STRINGIFY(crosslinking_network));

        std::vector<std::string> crosslinking_curves;
        readVectorFromJSON(j, crosslinking_curves, STRINGIFY(crosslinking_curves));

        // check inputs
        if (chain_psmiles_list.size() != crosslinking_curves.size() || crosslinking_curves.size() != crosslinking_network.size()) {
            std::string err_info = std::string("Input Error: ") +
                                   " Chain psmiles list size: " + std::to_string(chain_psmiles_list.size()) +
                                   " Crosslinking Curves size: " + std::to_string(crosslinking_curves.size()) +
                                   " Crosslinking Network size: " + std::to_string(crosslinking_network.size());
            throw std::runtime_error(err_info);
        }

        // 查看文件是否存在
        for (const auto& file : crosslinker_mol2_list) {
            if (!fileExists(file)) {
                throw std::runtime_error("Unknown file: " + file);
            }
        }
        for (const auto& file : crosslinking_curves) {
            if (!fileExists(file)) {
                throw std::runtime_error("Unknown file: " + file);
            }
        }

        std::vector<std::shared_ptr<CrossLinker>> crosslinkers(crosslinker_number, nullptr);
        std::vector<std::vector<int>> seen(crosslinker_number);
        // 因为所有的交联点一定名称不同 设置Builder多此一举

        for (int i = 0; i < crosslinker_number; i ++) {
            int crosslinker_type = i + 1;
            crosslinkers[i] = chemio::buildCrossLinkerFromMol2(crosslinker_mol2_list[i]);
            seen[i].resize(crosslinkers[i]->getPolysSize(), false);
            crosslinkers[i]->setMonoType(crosslinker_type);
        }

        // 检查poly是否存在.
        auto checkPolyAndGetPosition = [&crosslinkers](int who, int poly_id) -> Position {
            int poly_size = crosslinkers.at(who)->getPolysSize();
            if (poly_id < 0 || poly_id >= poly_size) {
                throw std::runtime_error("Error wrong poly:" + std::to_string(who) + " Crosslinker_poly_size:" + std::to_string(poly_size) + " PolyID: " + std::to_string(poly_id));
            }
            return crosslinkers[who]->getPolyPosition(poly_id);
        };

        std::vector<std::vector<Position>> curve_points((int)crosslinking_curves.size());

        // 检查曲线是否正确
        for (int i = 0; i < (int)crosslinking_curves.size(); i ++ ) {
            auto& cv = curve_points[i];
            cv = getScatterFromCSV(crosslinking_curves[i]);

            const auto& net = crosslinking_network[i];
            if (seen[net[0]][net[1]]) {
                throw std::runtime_error("Repeated poly: " + std::to_string(net[0]) + "-" + std::to_string(net[1]));
            }
            if (seen[net[2]][net[3]]) {
                throw std::runtime_error("Repeated poly: " + std::to_string(net[2]) + "-" + std::to_string(net[3]));
            }
            seen[net[0]][net[1]] = true;
            seen[net[2]][net[3]] = true;

            auto poly_position1 = checkPolyAndGetPosition(net[0], net[1]);
            auto poly_position2 = checkPolyAndGetPosition(net[2], net[3]);

            if (
                atom::positionDistanceSquared(cv.front(), poly_position1) > atom::positionDistanceSquared(cv.back(), poly_position1) &&
                atom::positionDistanceSquared(cv.back(), poly_position2) > atom::positionDistanceSquared(cv.front(), poly_position2)
                ) {
                // 曲线反了
                std::reverse(cv.begin(), cv.end());
            }

            // 检查曲线是否合格
            if (atom::positionDistanceSquared(cv.front(), poly_position1) > 4.0f || atom::positionDistanceSquared(cv.back(), poly_position2) > 4.0f) {
                throw std::runtime_error("Unqualified curve: " + crosslinking_curves[i]);
            }
        }

        std::cout << "========== Curves Generation Finish =============" << std::endl;

        auto cls = std::make_unique<CrosslinkingSystem>(crosslinkers, crosslinking_network, curve_points);
        std::vector<std::vector<std::shared_ptr<Graph>>> sequences((int)chain_psmiles_list.size());
        PsmilesBuilder<Graph> builder;

        for (int i = 0; i < (int)chain_psmiles_list.size(); i ++) {
            const auto& vec = chain_psmiles_list.at(i);
            for (const auto& smiles : vec) {
                auto ptr = builder.build(smiles);
                auto mono_type = builder.getMonomerType(smiles);
                ptr->setMonoTypeAll(mono_type);
                sequences[i].emplace_back(std::move(ptr));
            }
        }

        cls->calcChainGraphs(sequences, chain_max_polymerization_degree, pool_choice, para_A, para_B, random_polymerization, optimize_size);
        cls->makeEnd(std::string("H"));

        std::string mol2_file_name = saving_dir + "/crosslink.mol2";
        std::string adj_file_name = saving_dir + "/crosslink_adj";

        chemio::writeMol2File(mol2_file_name, adj_file_name, cls, file_info);
        std::cout << "Finish!" << std::endl;

        int cnt_poly = 0, sum_poly = 0;
        for (const auto& vec : seen) {
            sum_poly += (int)vec.size();
            for (const auto& t : vec) {
                cnt_poly += t;
            }
        }
//        std::cout << "sum_poly: " << sum_poly << " cnt_poly: " << cnt_poly << std::endl;
        float degree_crosslinking = (sum_poly ? (float)cnt_poly / (float)sum_poly : 0.0f);
        std::cout << "Degree of crosslinking: " << std::fixed << std::setprecision(2) << degree_crosslinking * 100 << "%" << std::endl;
        return;
    }

    // Unknown type
    throw std::runtime_error("Unknown Polymer Type: " + polymer_type);
}

void config(int argc, char* argv[], const std::string& config_filename) {
    Options options("polychef config", "Configure settings");
    options.add_options()
    ("type", "Task type: Chain or Crosslinks", cxxopts::value<std::string>())
    ("i,input", "Input files list", cxxopts::value<std::vector<std::string>>())
    ("p,psmiles", "Chain PSMILES list", cxxopts::value<std::vector<std::vector<std::string>>>())
    ("r,random", "Random polymerization or not", cxxopts::value<bool>())
    ("m,maximum", "Maximum length of chain", cxxopts::value<int>())
    ("opti", "Optimize Size", cxxopts::value<int>())
    ("o,output", "Output directory path", cxxopts::value<std::string>())
    ("info", "File info", cxxopts::value<std::string>())
    ("pool", "pooling choice", cxxopts::value<std::string>())
    ("v,verbose", "Verbose", cxxopts::value<bool>())
    ("t,threads", "Threads", cxxopts::value<int>())
    ("para_A", "Parameter A", cxxopts::value<float>())
    ("para_B", "Parameter B", cxxopts::value<float>())
    ("h,help", "Show help");

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return;
    }

    if (!fileExists(config_filename)) {
        throw std::runtime_error("Configure File " + config_filename + " does not exist");
    }

    json config;
    read_config(config, config_filename);

    if (result.count("type")) {
        config["polymer_type"] = result["type"].as<std::string>();
        std::cout << "Task type: " << config["polymer_type"] << std::endl;
    }

    if (result.count("input")) {
        config["chain_curve_list"] = result["input"].as<std::vector<std::string>>();
        std::cout << "Input files list: [\n";
        for (const auto& ss : config["chain_curve_list"]) {
            if (!fileExists(ss)) {
                throw exception::InvalidParameterException("Input filename Error");
            }
            std::cout << "  " << ss << "\n";
        }
        std::cout << "  ]" << std::endl;
    }

    if (result.count("psmiles")) {
        config["chain_psmiles_list"] = result["psmiles"].as<std::vector<std::vector<std::string>>>();
        std::cout << "Chain PSMILES list: [\n";
        for (const auto& vec : config["chain_psmiles_list"]) {
            std::cout << "==> ";
            for (const auto& ss : vec) {
                std::cout << ss << " ";
            }
            std::cout << "\n";
        }
        std::cout << "  ]" << std::endl;
    }

    if (result.count("random")) {
        config["random_polymerization"] = result["random"].as<bool>();
        std::cout << "Random Polymerization: " << (config["random_polymerization"] ? "Yes" : "No") << std::endl;
    }

    if (result.count("maximum")) {
        config["chain_max_polymerization_degree"] = result["maximum"].as<int>();
        std::cout << "Chain_max_polymerization_degree: " << config["chain_max_polymerization_degree"] << std::endl;
    }

    if (result.count("opti")) {
        config["optimize_size"] = result["opti"].as<int>();
        std::cout << "Optimize Size: " << config["optimize_size"] << std::endl;
    }

    if (result.count("output")) {
        config["output_directory_path"] = result["output"].as<std::string>();
        std::cout << "Output directory path has been changed to: " << config["output_directory_path"] << std::endl;
    }

    if (result.count("pool")) {
        config["pool_choice"] = result["pool"].as<std::string>();
        std::cout << "Pooling choice has been changed to: " << config["pool_choice"] << std::endl;
    }

    if (result.count("info")) {
        config["file_info"] = result["info"].as<std::string>();
        std::cout << "File Info: " << config["file_info"] << std::endl;
    }

    if (result.count("verbose")) {
        config["verbose"] = result["verbose"].as<bool>();
        std::cout << "Verbose mode " << (config["verbose"] ? "enabled" : "disabled") << "." << std::endl;
    }

    if (result.count("threads")) {
        config["threads"] = result["threads"].as<int>();
        std::cout << "Threads number: " << config["threads"] << std::endl;
    }

    if (result.count("para_A")) {
        config["para_A"] = result["para_A"].as<float>();
        std::cout << "Parameter A: " << config["para_A"] << std::endl;
    }

    if (result.count("para_B")) {
        config["para_B"] = result["para_B"].as<float>();
        std::cout << "Parameter B: " << config["para_B"] << std::endl;
    }

    if (save_config(config, config_filename)) {
        std::cout << "Configure settings Finish!" << std::endl;
    }
}

void run(int argc, char* argv[], const std::string& config_filename) {
    if (!fileExists(config_filename)) {
        throw std::runtime_error("Configure File " + config_filename + " does not exist");
    }
    auto start = std::chrono::high_resolution_clock::now(); // start time
    solve(config_filename);
    auto end = std::chrono::high_resolution_clock::now(); // end time
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Elapsed time: " << duration.count() << " ms" << std::endl;
}

int main(int argc, char* argv[]) {
    std::string cwd = std::filesystem::current_path();
    std::cout << "Current working directory: " << cwd << std::endl;

    if (argc < 3) {
        std::cerr << "Usage: polychef <command> [options]\n"
        << "Commands:\n"
        << "  run     Run a task\n"
        << ("  config  Configure \n");
        return 1;
    }

    std::string command = argv[1];
    std::string config_filename = argv[2];

    if (command == "run") {
        run(argc, argv, config_filename);
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