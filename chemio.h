//
// Created by AnthonyZhang on 2025/1/12.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_CHEMIO_H
#define ATOM_SEARCH_CPP_CHEMIO_H

#include "graph.h"
#include <memory>
#include <vector>
#include <string>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <Python.h>
#include <variant>
#include <any>

namespace chemio {
    std::variant<long, double, std::string> convertPyObject(PyObject* obj);
    std::vector<std::variant<long, double, std::string>> extract_tuple(PyObject* tuple);
    std::vector<std::vector<std::variant<long, double, std::string>>> extract_list(PyObject* list);

    std::shared_ptr<Graph> buildGraphFromPSmiles(const std::string& psmiles);
    std::string getAtomType(const std::string& name, int bond_num, bool ar = false);
    void writeMol2File(const std::string& path, const std::shared_ptr<Graph>& g, const std::string& file_info);
}

#endif //ATOM_SEARCH_CPP_CHEMIO_H
