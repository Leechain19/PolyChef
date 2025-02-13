//
// Created by AnthonyZhang on 2025/2/13.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_LOADLIB_H
#define ATOM_SEARCH_CPP_LOADLIB_H

#include <dlfcn.h>
#include <functional>
#include "atom.h"

typedef bool (*customFunction)(const Position &);

// RAII

class CustomFunctionLoader {
public:
    CustomFunctionLoader() = default;
    explicit CustomFunctionLoader(const std::string& lib_path, const std::string& funcname = "customFunction");
    CustomFunctionLoader(const CustomFunctionLoader& other) = delete;
    virtual ~CustomFunctionLoader();

    CustomFunctionLoader& operator=(const CustomFunctionLoader& other) = delete;

    void loadInfiSpace();
    void loadSharedLib(const std::string& lib_path, const std::string& funcname = "customFunction");
    bool checkPoint(const Position& point);

private:
    std::function<bool(const Position &)> f_;
    void* handle_ = nullptr;
};

#endif //ATOM_SEARCH_CPP_LOADLIB_H
