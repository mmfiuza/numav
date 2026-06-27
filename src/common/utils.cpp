// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/utils.hpp"
#include "common/exception.hpp"
#include "common/maths.hpp"

namespace numav {

void trim_right_whitespace(std::string_view& sv) {
    constexpr std::string_view WHITE_SPACE = " \t\n\r\f\v";
    const uint64_t end = sv.find_last_not_of(WHITE_SPACE);
    sv = (end == std::string_view::npos) ? "" : sv.substr(0UL, end + 1UL);
}

FuncFloatToCmplx const2func(
    const Cmplx constant
) {
    return [constant](Float f) -> Cmplx {
        (void) f;
        return constant;
    };
}

FuncFloatToCmplx table2func(
    const char* const pqv_text_file
) {
    // open file
    std::ifstream file(pqv_text_file);
    std::string line;
    if (!file.is_open()) {
        error("Could not open file: {}", pqv_text_file);
    }

    // first pass: count lines
    uint64_t line_count = 0UL;
    while (std::getline(file, line)) {
        ++line_count;
    }
    file.clear();
    file.seekg(0UL, std::ios::beg);
    std::vector<Float> real_vec;
    real_vec.reserve(line_count);
    std::vector<Cmplx> cmplx_vec;
    cmplx_vec.reserve(line_count);

    // second pass: read each line
    while (std::getline(file, line))
    {
        // read frequency
        uint64_t first_comma_pos = line.find(',');
        if (first_comma_pos == std::string::npos) {
            continue; // malformed line is skipped
        }
        std::string col1_str = line.substr(0UL, first_comma_pos);
        std::istringstream col1_input_string(col1_str);
        Float col1;
        col1_input_string >> col1;
        real_vec.push_back(col1);
        
        // read real part of complex vector
        uint64_t second_comma_pos = line.find(',', first_comma_pos + 1UL);
        if (second_comma_pos == std::string::npos) {
            continue; // malformed line is skipped
        }
        std::string col2_str =
            line.substr(first_comma_pos + 1UL, second_comma_pos);
        std::istringstream col2_input_string(col2_str);
        Float col2;
        col2_input_string >> col2;
        
        // read imaginary part of complex vector
        std::string col3_str = line.substr(second_comma_pos + 1UL);
        std::istringstream col3_input_string(col3_str);
        Float col3;
        col3_input_string >> col3;
        
        // write complex vector
        cmplx_vec.push_back(Cmplx(col2, col3));
    }
    
    // create the FuncFloatToCmplx funciton
    auto func_real_to_cmplx = [real_vec, cmplx_vec](Float real) {
        return interpolate(real_vec, cmplx_vec, real);
    };

    return func_real_to_cmplx;
}

} // namespace numav
