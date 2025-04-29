//
// Created by AnthonyZhang on 2025/1/11.
//
#pragma once
#ifndef ATOM_SEARCH_CPP_EXCEPTION_H
#define ATOM_SEARCH_CPP_EXCEPTION_H

#include <stdexcept>
#include <string>

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

    class MissingConfigError : public std::runtime_error {
        public:
            explicit MissingConfigError(const std::string& field_name)
            : std::runtime_error("Required configuration parameter '" + field_name + "' is missing in JSON file") {}
    };
}


#endif //ATOM_SEARCH_CPP_EXCEPTION_H
