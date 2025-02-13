//
// Created by AnthonyZhang on 2025/2/13.
//

#include "loadlib.h"

CustomFunctionLoader::CustomFunctionLoader(const std::string& lib_path, const std::string& funcname){
    if (lib_path.empty()) {
        loadInfiSpace();
    }
    else {
        loadSharedLib(lib_path, funcname);
    }
}

CustomFunctionLoader::~CustomFunctionLoader() {
    if (handle_) {
        dlclose(handle_);
    }
}

void CustomFunctionLoader::loadSharedLib(const std::string &lib_path, const std::string& funcname) {
    void* handle = dlopen(lib_path.c_str(), RTLD_LAZY);
    if (!handle) {
        std::cerr << "Error loading library: " << dlerror() << std::endl;
        throw std::runtime_error("Error loading library");
    }

    void* func_ptr = dlsym(handle, funcname.c_str());
    if (!func_ptr) {
        std::cerr << "Error loading function customFunction: " << dlerror() << std::endl;
        dlclose(handle);
        throw std::runtime_error("Error loading function customFunction");
    }
    f_ = reinterpret_cast<customFunction>(func_ptr);
}

void CustomFunctionLoader::loadInfiSpace() {
    f_ = [](const Position &point) {
        return true;
    };
}

bool CustomFunctionLoader::checkPoint(const Position &point) {
    return f_(point);
}
