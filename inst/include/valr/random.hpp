// random.hpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef VALR_RANDOM_HPP
#define VALR_RANDOM_HPP

#include <random>

namespace valr {

// Random number generator types
using Engine = std::mt19937;
using UniformIntDist = std::uniform_int_distribution<int>;
using PiecewiseConstDist = std::piecewise_constant_distribution<>;

}  // namespace valr

#endif  // VALR_RANDOM_HPP
