//
// Created by AnthonyZhang on 2025/1/12.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_CHEMIO_H
#define ATOM_SEARCH_CPP_CHEMIO_H

#include "graph.h"
#include "crosslinker.h"
#include <memory>
#include <vector>
#include <string>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <Python.h>
#include <variant>
#include <any>
#include <filesystem>
#include <functional>
#include <type_traits>

class Index2CodePrinter {
public:
    explicit Index2CodePrinter(std::string header);
    static std::string index2code(int idx);
    std::string get(int mono_type);
private:
    std::vector<std::string> string_memo{};
    std::string header_;
};

namespace chemio {
    std::variant<long, double, std::string> convertPyObject(PyObject* obj);
    std::vector<std::variant<long, double, std::string>> extract_tuple(PyObject* tuple);
    std::vector<std::vector<std::variant<long, double, std::string>>> extract_list(PyObject* list);

    std::shared_ptr<Graph> PyInfoConvertToGraph(
            const std::vector<std::vector<std::variant<long, double, std::string>>> &atoms_vec,
            const std::vector<std::vector<std::variant<long, double, std::string>>> &edges_vec,
            const std::vector<std::vector<std::variant<long, double, std::string>>> &polys_vec,
            const std::vector<std::vector<std::variant<long, double, std::string>>> &ring_edges_vec
            );

    std::shared_ptr<CrossLinker> PyInfoConvertToCrossLinker(
            const std::vector<std::vector<std::variant<long, double, std::string>>> &atoms_vec,
            const std::vector<std::vector<std::variant<long, double, std::string>>> &edges_vec,
            const std::vector<std::vector<std::variant<long, double, std::string>>> &polys_vec,
            const std::vector<std::vector<std::variant<long, double, std::string>>> &ring_edges_vec
            );

    std::shared_ptr<Graph> buildGraphFromPSmiles(const std::string& psmiles);
    std::shared_ptr<Graph> buildGraphFromMol2(const std::string& path);

    std::shared_ptr<CrossLinker> buildCrossLinkerFromPSmiles(const std::string& psmiles);
    std::shared_ptr<CrossLinker> buildCrossLinkerFromMol2(const std::string& path);

    std::string getAtomType(const std::string& name, int bond_num, bool ar = false);
    void writeMol2File(const std::string& file_name, const std::string& adj_file_name, const std::shared_ptr<Graph>& g, const std::string& file_info);
    void writeMol2File(const std::string& file_name, const std::string& adj_file_name, const std::unique_ptr<CrosslinkingSystem>& cls, const std::string& file_info);
    void writeLoss2File(const std::string& loss_file_name, const std::unique_ptr<std::vector<std::pair<double, double>>>& loss_vec_ptr);
}

template <typename T>
struct TypeTraits;

template <>
struct TypeTraits<Graph> {
    static std::shared_ptr<Graph> buildFromPSmiles(const std::string& psmiles) {
        return chemio::buildGraphFromPSmiles(psmiles);
    }
    static std::shared_ptr<Graph> buildFromMol2(const std::string& path) {
        return chemio::buildGraphFromMol2(path);
    }
};

template <>
struct TypeTraits<CrossLinker> {
    static std::shared_ptr<CrossLinker> buildFromPSmiles(const std::string& psmiles) {
        return chemio::buildCrossLinkerFromPSmiles(psmiles);
    }
    static std::shared_ptr<CrossLinker> buildFromMol2(const std::string& path) {
        return chemio::buildCrossLinkerFromMol2(path);
    }
};

template<typename T>
class PsmilesBuilder {
public:
    PsmilesBuilder() = default;
    ~PsmilesBuilder() = default;
    PsmilesBuilder(const PsmilesBuilder<T>& o) = delete;
    PsmilesBuilder<T>& operator=(const PsmilesBuilder<T>& o) = delete;

    std::shared_ptr<T> build(const std::string& psmiles) {
        auto it = memo_.find(psmiles);
        if (it != memo_.end()) return (it->second).first;
        auto ret = TypeTraits<T>::buildFromPSmiles(psmiles);
        memo_[psmiles] = std::make_pair(ret, mono_type ++);
        return ret;
    }

    int getMonomerType(const std::string& psmiles) {
        return memo_[psmiles].second;
    }

    void print() {
        for (const auto& [ss, p] : memo_) {
            std::cout << "ss: " << ss << " index: " << p.second << std::endl;
        }
    }

private:
    std::unordered_map<std::string, std::pair<std::shared_ptr<T>, int>> memo_{};
    int mono_type = 1;
};


template<typename T>
class Mol2Builder {
public:
    Mol2Builder() = default;
    ~Mol2Builder() = default;
    Mol2Builder(const Mol2Builder<T>& o) = delete;
    Mol2Builder<T>& operator=(const Mol2Builder<T>& o) = delete;

    std::shared_ptr<T> build(const std::string& path) {
        auto it = memo_.find(path);
        if (it != memo_.end()) return (it->second).first;
        auto ret = TypeTraits<T>::buildFromMol2(path);
        memo_[path] = std::make_pair(ret, mono_type ++);
        return ret;
    }

    int getMonomerType(const std::string& psmiles) {
        return memo_[psmiles].second;
    }

    void print() {
        for (const auto& [ss, p] : memo_) {
            std::cout << "ss: " << ss << " index: " << p.second << std::endl;
        }
    }

private:
    std::unordered_map<std::string, std::pair<std::shared_ptr<T>, int>> memo_{};
    int mono_type = 1;
};

#endif //ATOM_SEARCH_CPP_CHEMIO_H
