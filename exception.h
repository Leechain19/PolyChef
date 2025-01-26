//
// Created by AnthonyZhang on 2025/1/11.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_EXCEPTION_H
#define ATOM_SEARCH_CPP_EXCEPTION_H

#include <stdexcept>

namespace exception {
    class IllegalStringException : public std::exception {
    public:
        explicit IllegalStringException(std::string message) : msg_(std::move(message)) {}
        [[nodiscard]] const char* what() const noexcept override {
            return msg_.c_str();
        }
    private:
        std::string msg_;
    };

    class InvalidParameterException : public std::exception {
    public:
        explicit InvalidParameterException(std::string message) : msg_(std::move(message)) {}
        [[nodiscard]] const char* what() const noexcept override {
            return msg_.c_str();
        }
    private:
        std::string msg_;
    };
}


#endif //ATOM_SEARCH_CPP_EXCEPTION_H
