// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/utils.hpp"

void trim_right_whitespace(std::string_view& sv) {
    constexpr std::string_view WHITE_SPACE = " \t\n\r\f\v";
    const size_t end = sv.find_last_not_of(WHITE_SPACE);
    sv = (end == std::string_view::npos) ? "" : sv.substr(0, end+1);
}

