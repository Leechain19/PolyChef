//
// Created by AnthonyZhang on 2025/2/13.
//

#include "loadlib.h"
#include "exception.h"

CustomFunctionLoader::CustomFunctionLoader(const std::string& lib_path, const std::string& funcname){
    if (lib_path.empty()) {
        throw exception::IllegalStringException("Error: empty libaray path");
    }
    loadSharedLib(lib_path, funcname);
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

    std::function<int()> read_f;

    auto readin = [&read_f, &handle](const char* funcname) -> int {
        void* func_ptr = dlsym(handle, funcname);
        if (!func_ptr) {
            std::cerr << "Error loading function customFunction: " << dlerror() << std::endl;
            dlclose(handle);
            throw std::runtime_error("Error loading function customFunction");
        }
        read_f = reinterpret_cast<readinFunction>(func_ptr);
        return read_f();
    };

    x_min_ = readin("xmin");
    x_max_ = readin("xmax");
    y_min_ = readin("ymin");
    y_max_ = readin("ymax");
    z_min_ = readin("zmin");
    z_max_ = readin("zmax");
}

bool CustomFunctionLoader::checkPoint(const Position &point) {
    return f_(point);
}

int CustomFunctionLoader::x_min() const {
    return x_min_;
}

int CustomFunctionLoader::x_max() const {
    return x_max_;
}

std::pair<int, int> CustomFunctionLoader::x_min_max() const {
    return {x_min_, x_max_};
}

int CustomFunctionLoader::y_min() const {
    return y_min_;
}

int CustomFunctionLoader::y_max() const {
    return y_max_;
}

std::pair<int, int> CustomFunctionLoader::y_min_max() const {
    return {y_min_, y_max_};
}

int CustomFunctionLoader::z_min() const {
    return z_min_;
}

int CustomFunctionLoader::z_max() const {
    return z_max_;
}

std::pair<int, int> CustomFunctionLoader::z_min_max() const {
    return {z_min_, z_max_};
}
