// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/utils.hpp"
#include "common/exception.hpp"
#include "common/maths.hpp"

namespace numav {

void trim_right_whitespace(std::string_view& sv) {
    constexpr std::string_view WHITE_SPACE = " \t\n\r\f\v";
    const size_t end = sv.find_last_not_of(WHITE_SPACE);
    sv = (end == std::string_view::npos) ? "" : sv.substr(0, end+1);
}

_FuncRealToCmplx convert_table_to_real_to_cmplx_func(
    const char* const impedance_text_file
) {
    // open file
    std::ifstream file(impedance_text_file);
    std::string line;
    if (!file.is_open()) {
        error("Could not open file: {}", impedance_text_file);
    }

    // first pass: count lines
    size_t line_count = 0;
    while (std::getline(file, line)) {
        ++line_count;
    }
    file.clear();
    file.seekg(0, std::ios::beg);
    std::vector<double> real_vec;
    real_vec.reserve(line_count);
    std::vector<_cmplx_t> cmplx_vec;
    cmplx_vec.reserve(line_count);

    // second pass: read each line
    while (std::getline(file, line))
    {
        // read frequency
        size_t first_comma_pos = line.find(',');
        if (first_comma_pos == std::string::npos) {
            continue; // malformed line is skipped
        }
        std::string col1_str = line.substr(0, first_comma_pos);
        std::istringstream col1_input_string(col1_str);
        double col1;
        col1_input_string >> col1;
        real_vec.push_back(col1);
        
        // read real part of complex vector
        size_t second_comma_pos = line.find(',', first_comma_pos+1);
        if (second_comma_pos == std::string::npos) {
            continue; // malformed line is skipped
        }
        std::string col2_str =
            line.substr(first_comma_pos + 1, second_comma_pos);
        std::istringstream col2_input_string(col2_str);
        double col2;
        col2_input_string >> col2;
        
        // read imaginary part of complex vector
        std::string col3_str = line.substr(second_comma_pos + 1);
        std::istringstream col3_input_string(col3_str);
        double col3;
        col3_input_string >> col3;
        
        // write complex vector
        cmplx_vec.push_back(_cmplx_t(col2, col3));
    }
    
    // create the _FuncRealToCmplx funciton
    auto func_real_to_cmplx = [real_vec, cmplx_vec](double real) {
        return interpolate(real_vec, cmplx_vec, real);
    };

    return func_real_to_cmplx;
}

} // namespace numav
