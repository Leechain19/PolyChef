//
// Created by AnthonyZhang on 2025/2/13.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_LOADLIB_H
#define ATOM_SEARCH_CPP_LOADLIB_H

#include <dlfcn.h>
#include <functional>
#include <array>
#include "atom.h"

typedef bool (*customFunction)(const Position &);
typedef int (*readinFunction)();

// RAII

class CustomFunctionLoader {
public:
    CustomFunctionLoader() = default;
    explicit CustomFunctionLoader(const std::string& lib_path, const std::string& funcname = "customFunction");
    CustomFunctionLoader(const CustomFunctionLoader& other) = delete;
    virtual ~CustomFunctionLoader();
    CustomFunctionLoader& operator=(const CustomFunctionLoader& other) = delete;
    bool checkPoint(const Position& point);
    [[nodiscard]] int x_min() const ;
    [[nodiscard]] int x_max() const ;
    [[nodiscard]] std::pair<int, int> x_min_max() const ;

    [[nodiscard]] int y_min() const ;
    [[nodiscard]] int y_max() const ;
    [[nodiscard]] std::pair<int, int> y_min_max() const ;

    [[nodiscard]] int z_min() const ;
    [[nodiscard]] int z_max() const ;
    [[nodiscard]] std::pair<int, int> z_min_max() const ;


private:
    void loadSharedLib(const std::string& lib_path, const std::string& funcname = "customFunction");
    std::function<bool(const Position &)> f_;
    void* handle_ = nullptr;
    int x_min_{}, x_max_{}, y_min_{}, y_max_{}, z_min_{}, z_max_{};
};

#endif //ATOM_SEARCH_CPP_LOADLIB_H
