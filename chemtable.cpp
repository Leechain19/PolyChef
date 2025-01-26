//
// Created by AnthonyZhang on 2025/1/7.
//

#include "chemtable.h"
#include "exception.h"

std::string chemtable::process(const std::string &s) {
    if (s.empty() || !isalpha(s[0])) {
        throw exception::IllegalStringException("Input string is empty or does not start with an alphabetic character");
    }
    std::string ret;
    for (const auto& c : s) {
        if (isdigit(c)) break;
        ret.push_back(c);
    }
    return ret;
}

bool BondTable::have_init = false;
std::unordered_map<std::string, float> BondTable::bond_table = {};

void BondTable::init() {
    if (have_init) {
        return;
    }
    have_init = true;
    bond_table = {
            {"Br+Br+-", 2.29}, {"C+B+-", 1.56}, {"C+Br+-", 1.94}, {"C+C+-", 1.54}, {"C+C+=", 1.34},
            {"C+C+≡", 1.2}, {"C+Cl+-", 1.77}, {"C+F+-", 1.38}, {"C+H+-", 1.09}, {"C+I+-", 2.14},
            {"C+N+-", 1.48}, {"C+N+=", 1.35}, {"C+N+≡", 1.16}, {"C+O+-", 1.43}, {"C+O+=", 1.2},
            {"C+P+-", 1.87}, {"C+S+-", 1.82}, {"C+S+=", 1.56}, {"C+Si+-", 1.86}, {"Cl+Cl+-", 1.99},
            {"Cs+I+-", 3.37}, {"F+F+-", 1.4}, {"H+H+-", 0.75}, {"H+Br+-", 1.42}, {"H+Cl+-", 1.27},
            {"H+F+-", 0.92}, {"H+I+-", 1.61}, {"I+I+-", 2.66}, {"K+Br+-", 2.82}, {"K+Cl+-", 2.67},
            {"N+H+-", 1.01}, {"N+N+-", 1.45}, {"N+N+=", 1.25}, {"N+N+≡", 1.1}, {"N+O+-", 1.46},
            {"N+O+=", 1.14}, {"O+H+-", 0.98}, {"O+O+-", 1.48}, {"O+O+=", 1.2}, {"P+Br+-", 2.2},
            {"P+Cl+-", 2.03}, {"P+O+-", 1.63}, {"P+O+=", 1.38}, {"S+H+-", 1.35}, {"S+O+=", 1.43},
            {"S+S+-", 2.07}, {"S+S+=", 1.89}, {"Se+Se+-", 2.32}, {"Se+H+-", 1.47}, {"Se+Se+=", 2.15}
    };
}

float BondTable::get_length(const std::string &s1, const std::string &s2, const std::string &bond_type) {
    init();
    auto ss1 = chemtable::process(s1);
    auto ss2 = chemtable::process(s2);
    {
        auto &&key = ss1 + "+" + ss2 + "+" + bond_type;
        if (bond_table.count(key)) {
            return bond_table[key];
        }
    }

    {
        auto &&key = ss2 + "+" + ss1 + "+" + bond_type;
        if (bond_table.count(key)) {
            return bond_table[key];
        }
    }
    throw exception::IllegalStringException("No such a chemical bond in BondTable");
}

bool AminoAcidTable::have_init = false;
std::unordered_map<std::string, std::string> AminoAcidTable::amino2psmiles = {};

void AminoAcidTable::init() {
    if (have_init) {
        return;
    }
    have_init = true;
    amino2psmiles = {
            {"A","C[C@H](N[*])C(=O)[*]"},
            {"R","N=C(N)NCCC[C@H](N[*])C(=O)[*]"},
            {"N","NC(=O)C[C@H](N[*])C(=O)[*]"},
            {"D","[*]N[C@@H](CC(=O)O)C(=O)[*]"},
            {"C","[*]N[C@@H](CS)C(=O)[*]"},
            {"Q","NC(=O)CC[C@H](N[*])C(=O)[*]"},
            {"E","[*]N[C@@H](CCC(=O)O)C(=O)[*]"},
            {"G","[*]NCC(=O)[*]"},
            {"H","[*]N[C@@H](Cc1c[nH]cn1)C(=O)[*]"},
            {"I","CC[C@H](C)[C@H](N[*])C(=O)[*]"},
            {"L","CC(C)C[C@H](N[*])C(=O)[*]"},
            {"K","NCCCC[C@H](N[*])C(=O)[*]"},
            {"M","CSCC[C@H](N[*])C(=O)[*]"},
            {"F","[*]N[C@@H](Cc1ccccc1)C(=O)[*]"},
            {"P","O=C([*])[C@@H]1CCCN1[*]"},
            {"S","[*]N[C@@H](CO)C(=O)[*]"},
            {"T","C[C@@H](O)[C@H](N[*])C(=O)[*]"},
            {"W","[*]N[C@@H](Cc1c[nH]c2ccccc12)C(=O)[*]"},
            {"Y","[*]N[C@@H](Cc1ccc(O)cc1)C(=O)[*]"},
            {"V","CC(C)[C@H](N[*])C(=O)[*]"}
    };
}

std::string AminoAcidTable::get_amino_psmiles(const std::string &symbol) {
    init();
    if ((int)symbol.size() != 1 || symbol[0] < 'A' || symbol[0] > 'Z') {
        throw exception::IllegalStringException("Illegal AminoAcid String");
    }
    if (!amino2psmiles.count(symbol)) {
        throw exception::IllegalStringException("No such an AminoAcid in AminoAcidTable");
    }
    return amino2psmiles[symbol];
}

std::vector<std::pair<std::string, std::string>> AminoAcidTable::get_table() {
    std::vector<std::pair<std::string, std::string>> ret(amino2psmiles.begin(), amino2psmiles.end());
    return ret;
}

bool AtomTable::have_init = false;
std::unordered_map<std::string, float> AtomTable::atom_table = {};

void AtomTable::init() {
    if (have_init) {
        return;
    }
    have_init = true;
    atom_table = {
            {"H", 1.1}, {"He", 1.4}, {"B", 2.13},
            {"C", 1.72}, {"N", 1.5}, {"O", 1.4},
            {"F", 1.35}, {"Ne", 1.54}, {"Al", 2.51}, {"Si", 2.1},
            {"P", 1.90}, {"S", 1.85}, {"Cl", 1.75}
    };
}

float AtomTable::get_radius(const std::string &symbol) {
    init();
    auto s = chemtable::process(symbol);
    if (!atom_table.count(s)) {
        throw exception::IllegalStringException("No such an Atom in AtomTable");
    }
    return atom_table[s];
}


