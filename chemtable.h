//
// Created by AnthonyZhang on 2025/1/7.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_CHEMTABLE_H
#define ATOM_SEARCH_CPP_CHEMTABLE_H

#include <unordered_map>
#include <string>
#include <utility>
#include <vector>
#include <stdexcept>

namespace chemtable {
    std::string process(const std::string &s);
}

class BondTable {
private:
    static std::unordered_map<std::string, float> bond_table;
    static bool have_init;
public:
    static float get_length(const std::string &s1, const std::string &s2, const std::string &bond_type);
    static void init();
};

class AminoAcidTable {
private:
    static std::unordered_map<std::string, std::string> amino2psmiles;
    static bool have_init;
public:
    static void init();
    static std::string get_amino_psmiles(const std::string &symbol);
    static std::vector<std::pair<std::string, std::string>> get_table();
};

class AtomTable {
private:
    static std::unordered_map<std::string, float> atom_table;
    static bool have_init;
public:
    static void init();
    static float get_radius(const std::string &symbol);
};



#endif //ATOM_SEARCH_CPP_CHEMTABLE_H