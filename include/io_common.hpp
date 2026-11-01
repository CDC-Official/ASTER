//
//  io_common.hpp
//  Test_dylib
//
//  Created by 陈酊存 on 2023/8/30.
//

#ifndef CDC_IO_COMMON_HPP
#define CDC_IO_COMMON_HPP

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <cctype>
#include "utils_overwatch.hpp"

int line_count(std::ifstream& infile);

int line_count(std::string file_path);

void WriteString(std::ostream& out, std::string str, const int width = 10, bool ENDL = false);

void writeDouble(std::ostream& out, double number, const int width = 22, const int precision = 16, bool ENDL = false);

void writeInteger(std::ostream& out, int number, const int width = 6, bool ENDL = false);

void split(const std::string& str, const std::string& delimiter, std::vector<std::string>& tokens);

void split(const std::string& str, std::vector<std::string>& tokens);

std::string strip(const std::string& str);
#endif /* CDC_IO_COMMON_HPP */
