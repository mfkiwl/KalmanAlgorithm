/*  __          _      __  __          _   _
| |/ /    /\   | |    |  \/  |   /\   | \ | |
| ' /    /  \  | |    | \  / |  /  \  |  \| |
|  <    / /\ \ | |    | |\/| | / /\ \ | . ` |
| . \  / ____ \| |____| |  | |/ ____ \| |\  |
|_|\_\/_/    \_\______|_|  |_/_/    \_\_| \_|

Kalman Filter
Version 0.4.0
https://github.com/FrancoisCarouge/Kalman

SPDX-License-Identifier: Unlicense

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <https://unlicense.org> */

#ifndef FCAROUGE_PHYSICAL_LINALG_HPP
#define FCAROUGE_PHYSICAL_LINALG_HPP

#include "fcarouge/linalg.hpp"
#include "fcarouge/unit.hpp"
#include "fcarouge/utility.hpp"

namespace fcarouge {

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
class physical_matrix {
public:
  inline constexpr physical_matrix() = default;

  inline constexpr physical_matrix(const physical_matrix &other) = default;

  inline constexpr physical_matrix &
  operator=(const physical_matrix &other) = default;

  inline constexpr physical_matrix(physical_matrix &&other) = default;

  inline constexpr physical_matrix &
  operator=(physical_matrix &&other) = default;

  explicit inline constexpr physical_matrix(const Matrix &value)
      : data{value} {}

  inline constexpr physical_matrix(const auto &...elements)
      // +requires scalar type compatibility
    requires(size<RowIndexes> != 1 && size<ColumnIndexes> == 1)
  {
    size_t<ColumnIndexes> i{0};
    // explicit conversion?
    ([&] { data[i++] = elements; }(), ...);
  }

  [[nodiscard]] static constexpr inline auto identity() {
    return identity_v<Matrix>;
  }

  [[nodiscard]] static constexpr inline auto zero() { return zero_v<Matrix>; }

  // private:
  Matrix data;
};

// Is int the right type for promotion rules in all cases?
using one_row = pack<int>;
using one_column = pack<int>;

template <typename Matrix, typename... RowIndexes>
using physical_column_vector =
    physical_matrix<Matrix, pack<RowIndexes...>, one_column>;

template <typename Matrix, typename... ColumnIndexes>
using physical_row_vector =
    physical_matrix<Matrix, one_row, pack<ColumnIndexes...>>;

template <typename Matrix1, typename Matrix2, typename Row, typename Size,
          typename Column>
[[nodiscard]] inline constexpr auto
operator*(const physical_matrix<Matrix1, Row, Size> &lhs,
          const physical_matrix<Matrix2, Size, Column> &rhs) {

  auto result{lhs.data * rhs.data};
  using matrix = decltype(result);
  return physical_matrix<matrix, Row, Column>{result};
}

template <typename Matrix, typename Row, typename Column>
[[nodiscard]] inline constexpr auto
transpose([[maybe_unused]] const physical_matrix<Matrix, Row, Column> &lhs) {
  internal::transposer t;
  auto result{t(lhs.data)};
  using matrix = decltype(result);
  return physical_matrix<matrix, Column, Row>{result};
}
} // namespace fcarouge

#endif // FCAROUGE_PHYSICAL_LINALG_HPP
