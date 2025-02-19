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
#include <type_traits>

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

    template <class T>
    void PyInfoConvertToCppType(
            std::shared_ptr<T>& g,
            const std::vector<std::vector<std::variant<long, double, std::string>>> &atoms_vec,
            const std::vector<std::vector<std::variant<long, double, std::string>>> &edges_vec,
            const std::vector<std::vector<std::variant<long, double, std::string>>> &polys_vec,
            const std::vector<std::vector<std::variant<long, double, std::string>>> &ring_edges_vec ) {
        g = std::make_shared<T>();
        // atoms
        for (const auto& vec : atoms_vec) {
            std::string name = std::get<std::string>(vec[0]);
            float x = static_cast<float>(std::get<double>(vec[1]));
            float y = static_cast<float>(std::get<double>(vec[2]));
            float z = static_cast<float>(std::get<double>(vec[3]));
            bool ar = static_cast<bool>(std::get<long>(vec[4]));
            g->addAtom(std::make_shared<Atom>(name, x, y, z), 0, ar);
        }

        // edges
        for (const auto& vec : edges_vec) {
            int x = static_cast<int>(std::get<long>(vec[0]));
            int y = static_cast<int>(std::get<long>(vec[1]));
            std::string type = std::get<std::string>(vec[2]);
            g->addEdge(x, y, type);
            g->addEdge(y, x, type);
        }

        // polys
        for (const auto& vec : polys_vec) {
            int who = static_cast<int>(std::get<long>(vec[0]));
            float x = static_cast<float>(std::get<double>(vec[1]));
            float y = static_cast<float>(std::get<double>(vec[2]));
            float z = static_cast<float>(std::get<double>(vec[3]));
            g->addPoly(x, y, z, who);
        }

        if constexpr (std::is_same_v<T, Graph>) {
            // ring
            for (const auto& vec : ring_edges_vec) {
                int x = static_cast<int>(std::get<long>(vec[0]));
                int y = static_cast<int>(std::get<long>(vec[1]));
                g->addRingEdge(x, y);
            }
            g->calMainChain();
        }
    }

    std::shared_ptr<Graph> buildGraphFromPSmiles(const std::string& psmiles);
    std::shared_ptr<Graph> buildGraphFromMol2(const std::string& path);

    template <class T>
    bool buildCppTypeFromPSmiles(std::shared_ptr<T>& g, const std::string& psmiles) {
        if (psmiles.empty()) {
            throw exception::InvalidParameterException("The psmiles string is empty");
        }
        std::string smi;
        for (const auto& c : psmiles) {
            if (c == '*') {
                smi += "13C";
            }
            else smi += c;
        }

        std::cout << "psmiles = " << psmiles << std::endl;
        std::cout << "smi = " << smi << std::endl;

        // Py
        Py_Initialize();
        if (!Py_IsInitialized()) {
            std::cerr << "Python initialization failed" << std::endl;
            return false;
        }

        PyRun_SimpleString("import sys");
        PyRun_SimpleString("sys.path.append('./../script')");

        PyObject* module = PyImport_ImportModule("rdkit_helper");

        if (!module) {
            PyErr_Print();
            std::cerr << "Failed to load module 'rdkit_helper.py'" << std::endl;
            Py_Finalize();
            return false;
        }

        PyObject* func = PyObject_GetAttrString(module, "read_smi");
        if (!func || !PyCallable_Check(func)) {
            PyErr_Print();
            std::cerr << "Failed to load function 'read_smi'" << std::endl;
            Py_DECREF(module);
            Py_Finalize();
            return false;
        }

        PyObject* args = PyTuple_New(1);
        PyTuple_SET_ITEM(args, 0, Py_BuildValue("s", smi.c_str()));

        PyObject* result = PyObject_CallObject(func, args);
        if (!result) {
            PyErr_Print();
            std::cerr << "Function call failed" << std::endl;
            Py_DECREF(func);
            Py_DECREF(module);
            Py_Finalize();
            return false;
        }

        // 检查返回值是否为元组
        if (!PyTuple_Check(result)) {
            std::cerr << "Function did not return a tuple" << std::endl;
            Py_DECREF(result);
            Py_DECREF(func);
            Py_DECREF(module);
            Py_Finalize();
            return false;
        }

        // 提取返回的 4 个列表
        PyObject* atoms_cache = PyTuple_GetItem(result, 0);
        PyObject* edges_cache = PyTuple_GetItem(result, 1);
        PyObject* polys_cache = PyTuple_GetItem(result, 2);
        PyObject* ring_edges_cache = PyTuple_GetItem(result, 3);

        // 转换为 C++ 数据结构
        std::vector<std::vector<std::variant<long, double, std::string>>> atoms_vec = chemio::extract_list(atoms_cache);
        std::vector<std::vector<std::variant<long, double, std::string>>> edges_vec = chemio::extract_list(edges_cache);
        std::vector<std::vector<std::variant<long, double, std::string>>> polys_vec = chemio::extract_list(polys_cache);
        std::vector<std::vector<std::variant<long, double, std::string>>> ring_edges_vec = chemio::extract_list(ring_edges_cache);

        // 清理资源
        Py_DECREF(args);
        Py_DECREF(result);
        Py_DECREF(func);
        Py_DECREF(module);

        PyInfoConvertToCppType(g, atoms_vec, edges_vec, polys_vec, ring_edges_vec);
        return true;
    }

    template <class T>
    bool buildGraphFromMol2(std::shared_ptr<T>& g, const std::string& path) {
        if (!std::filesystem::exists(path)) {
            throw exception::IllegalStringException("buildGraphFromMol2: Illegal Mol2 Path");
        }

        // Py
        Py_Initialize();
        if (!Py_IsInitialized()) {
            std::cerr << "Python initialization failed" << std::endl;
            return false;
        }

        PyRun_SimpleString("import sys");
        PyRun_SimpleString("sys.path.append('./../script')");

        PyObject* module = PyImport_ImportModule("rdkit_helper");

        if (!module) {
            PyErr_Print();
            std::cerr << "Failed to load module 'rdkit_helper'" << std::endl;
            Py_Finalize();
            return false;
        }

        PyObject* func = PyObject_GetAttrString(module, "read_mol2");
        if (!func || !PyCallable_Check(func)) {
            PyErr_Print();
            std::cerr << "Failed to load function 'read_mol2'" << std::endl;
            Py_DECREF(module);
            Py_Finalize();
            return false;
        }

        PyObject* args = PyTuple_New(1);
        PyTuple_SET_ITEM(args, 0, Py_BuildValue("s", path.c_str()));

        PyObject* result = PyObject_CallObject(func, args);
        if (!result) {
            PyErr_Print();
            std::cerr << "Function call failed" << std::endl;
            Py_DECREF(func);
            Py_DECREF(module);
            Py_Finalize();
            return false;
        }

        // 检查返回值是否为元组
        if (!PyTuple_Check(result)) {
            std::cerr << "Function did not return a tuple" << std::endl;
            Py_DECREF(result);
            Py_DECREF(func);
            Py_DECREF(module);
            Py_Finalize();
            return false;
        }

        // 提取返回的 4 个列表
        PyObject* atoms_cache = PyTuple_GetItem(result, 0);
        PyObject* edges_cache = PyTuple_GetItem(result, 1);
        PyObject* polys_cache = PyTuple_GetItem(result, 2);
        PyObject* ring_edges_cache = PyTuple_GetItem(result, 3);

        // 转换为 C++ 数据结构
        std::vector<std::vector<std::variant<long, double, std::string>>> atoms_vec = chemio::extract_list(atoms_cache);
        std::vector<std::vector<std::variant<long, double, std::string>>> edges_vec = chemio::extract_list(edges_cache);
        std::vector<std::vector<std::variant<long, double, std::string>>> polys_vec = chemio::extract_list(polys_cache);
        std::vector<std::vector<std::variant<long, double, std::string>>> ring_edges_vec = chemio::extract_list(ring_edges_cache);

        // 清理资源
        Py_DECREF(args);
        Py_DECREF(result);
        Py_DECREF(func);
        Py_DECREF(module);

        PyInfoConvertToCppType(g, atoms_vec, edges_vec, polys_vec, ring_edges_vec);
        return true;
    }

    std::string getAtomType(const std::string& name, int bond_num, bool ar = false);
    void writeMol2File(const std::string& file_name, const std::string& adj_file_name, const std::shared_ptr<Graph>& g, const std::string& file_info);
    void writeMol2File(const std::string& file_name, const std::string& adj_file_name, const std::shared_ptr<CrosslinkingSystem>& cls, const std::string& file_info);
    void writeLoss2File(const std::string& loss_file_name, const std::unique_ptr<std::vector<std::pair<double, double>>>& loss_vec_ptr);
}

#endif //ATOM_SEARCH_CPP_CHEMIO_H
