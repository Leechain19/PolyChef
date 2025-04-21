//
// Created by AnthonyZhang on 2025/4/19.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_UTILITY_H
#define ATOM_SEARCH_CPP_UTILITY_H

#include "nlohmann/json.hpp"
#include "atom.h"
#include <memory>
#include <vector>
#include <cassert>
#include <string>

using json = nlohmann::json;

bool fileExists(const std::string& filename);

bool startsWith(const std::string& s, const std::string& prefix);

bool endsWith(const std::string& s, const std::string& suffix);

std::string lower(const std::string& s);

std::string upper(const std::string& s);

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

std::string getCurrentTimeAsString();

std::string connectStringVectorAsOne(const std::vector<std::string>& vec);

namespace fs = std::filesystem;
void createDirectory(const std::string& path);

void InitializeOpenMP(int target_thread_num);

std::vector<std::string> my_split(const std::string &s, char delim);

std::shared_ptr<Atom> getAtomFromMol2Line(const std::string& s);


#endif //ATOM_SEARCH_CPP_UTILITY_H
