// dataframe.hpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef VALR_DATAFRAME_HPP
#define VALR_DATAFRAME_HPP

#include <cpp11.hpp>
#include <cpp11/R.hpp>

#include <string>
#include <vector>

namespace valr {

// Convert cpp11::integers to std::vector<int> with optional index adjustment
// If adjust_index is true, converts from R's 1-based to C++'s 0-based indexing
inline std::vector<int> to_std_vector(const cpp11::integers& x, bool adjust_index = false) {
  size_t n = x.size();
  std::vector<int> result(n);
  if (adjust_index) {
    for (size_t i = 0; i < n; ++i) {
      result[i] = x[i] - 1;
    }
  } else {
    for (size_t i = 0; i < n; ++i) {
      result[i] = x[i];
    }
  }
  return result;
}

// Convert std::vector<int> to cpp11::writable::integers
inline cpp11::writable::integers to_cpp11_integers(const std::vector<int>& x) {
  cpp11::writable::integers result(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    result[i] = x[i];
  }
  return result;
}

// Convert std::vector<double> to cpp11::writable::doubles
inline cpp11::writable::doubles to_cpp11_doubles(const std::vector<double>& x) {
  cpp11::writable::doubles result(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    result[i] = x[i];
  }
  return result;
}

// Helper to set factor levels on an integer vector
inline void set_factor_levels(SEXP x, SEXP levels) {
  if (TYPEOF(x) != INTSXP) {
    cpp11::stop("Internal error: Only integers can be made into factors");
  }

  SEXP factor_class = PROTECT(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(factor_class, 0, Rf_mkChar("factor"));

  Rf_setAttrib(x, R_LevelsSymbol, levels);
  Rf_setAttrib(x, R_ClassSymbol, factor_class);

  UNPROTECT(1);
}

// Subset a dataframe by row indices (0-based)
inline cpp11::writable::data_frame subset_dataframe(const cpp11::data_frame& df,
                                                    const std::vector<int>& row_indices) {
  int ncol = df.size();
  int nrow_out = row_indices.size();

  cpp11::writable::list output(ncol);

  // Copy column names (PROTECT to satisfy rchk)
  SEXP names = PROTECT(Rf_getAttrib(df, R_NamesSymbol));
  Rf_setAttrib(output, R_NamesSymbol, names);
  UNPROTECT(1);

  for (int j = 0; j < ncol; ++j) {
    SEXP col = VECTOR_ELT(df, j);
    SEXP new_col = PROTECT(Rf_allocVector(TYPEOF(col), nrow_out));

    for (int i = 0; i < nrow_out; ++i) {
      int idx = row_indices[i];

      switch (TYPEOF(new_col)) {
        case REALSXP:
          REAL(new_col)[i] = (idx == NA_INTEGER) ? NA_REAL : REAL(col)[idx];
          break;
        case INTSXP:
        case LGLSXP:
          INTEGER(new_col)[i] = (idx == NA_INTEGER) ? NA_INTEGER : INTEGER(col)[idx];
          break;
        case STRSXP:
          SET_STRING_ELT(new_col, i, (idx == NA_INTEGER) ? NA_STRING : STRING_ELT(col, idx));
          break;
        case VECSXP:
          SET_VECTOR_ELT(new_col, i, (idx == NA_INTEGER) ? R_NilValue : VECTOR_ELT(col, idx));
          break;
        default:
          UNPROTECT(1);
          cpp11::stop("Incompatible column type detected");
      }
    }

    // Preserve factor levels
    if (Rf_inherits(col, "factor")) {
      SEXP levels = PROTECT(Rf_getAttrib(col, R_LevelsSymbol));
      set_factor_levels(new_col, levels);
      UNPROTECT(1);
    }

    SET_VECTOR_ELT(output, j, new_col);
    UNPROTECT(1);
  }

  // Copy attributes
  Rf_copyMostAttrib(df, output);

  // Set row names
  SEXP row_names = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(row_names)[0] = NA_INTEGER;
  INTEGER(row_names)[1] = -nrow_out;
  Rf_setAttrib(output, R_RowNamesSymbol, row_names);
  UNPROTECT(1);

  // Set class
  Rf_setAttrib(output, R_ClassSymbol, Rf_mkString("data.frame"));

  return cpp11::writable::data_frame(output);
}

// Check if a dataframe is grouped
inline bool is_grouped(const cpp11::data_frame& df) {
  return Rf_inherits(df, "grouped_df");
}

// Handler for grouped dataframes from dplyr
class GroupedDataFrame {
 private:
  cpp11::data_frame data_;
  cpp11::data_frame groups_;

  static cpp11::data_frame check_is_grouped(const cpp11::data_frame& x) {
    if (!is_grouped(x)) {
      cpp11::stop("Input must be a grouped data frame");
    }
    return x;
  }

 public:
  explicit GroupedDataFrame(const cpp11::data_frame& x)
      : data_(check_is_grouped(x)),
        groups_(cpp11::as_cpp<cpp11::data_frame>(data_.attr("groups"))) {}

  const cpp11::data_frame& data() const { return data_; }

  int nrows() const { return static_cast<int>(data_.nrow()); }

  int ngroups() const { return static_cast<int>(groups_.nrow()); }

  const cpp11::data_frame& group_data() const { return groups_; }

  // Get the .rows column (list of integer vectors with row indices)
  cpp11::list indices() const {
    int last_col = groups_.size() - 1;
    return cpp11::as_cpp<cpp11::list>(groups_[last_col]);
  }

  // Strip grouping from a dataframe
  template <typename Data>
  static void strip_groups(Data& x) {
    x.attr("groups") = R_NilValue;
  }
};

// Builder for constructing output dataframes efficiently
class DataFrameBuilder {
 public:
  std::vector<std::string> names;
  cpp11::writable::list data;

  DataFrameBuilder() : data(static_cast<R_xlen_t>(0)) {}

  size_t size() const { return data.size(); }

  // Add a dataframe with optional suffix, optionally dropping chrom column
  void add_df(const cpp11::data_frame& df, const std::string& suffix = "", bool drop_chrom = true) {
    int ncol = df.size();
    cpp11::strings col_names(df.attr("names"));

    for (int i = 0; i < ncol; ++i) {
      std::string name(col_names[i]);

      if (name == "chrom" && drop_chrom) {
        continue;
      }

      if (!suffix.empty() && name != "chrom") {
        name += suffix;
      }

      names.push_back(name);
      data.push_back(df[i]);
    }
  }

  // Add a single column
  template <typename T>
  void add_column(const std::string& name, const T& col) {
    names.push_back(name);
    data.push_back(cpp11::as_sexp(col));
  }

  // Add a vector column
  void add_column(const std::string& name, const std::vector<int>& col) {
    names.push_back(name);
    cpp11::writable::integers vec(col.size());
    for (size_t i = 0; i < col.size(); ++i) {
      vec[i] = col[i];
    }
    data.push_back(vec);
  }

  void add_column(const std::string& name, const std::vector<double>& col) {
    names.push_back(name);
    cpp11::writable::doubles vec(col.size());
    for (size_t i = 0; i < col.size(); ++i) {
      vec[i] = col[i];
    }
    data.push_back(vec);
  }

  // Build the final dataframe
  // nrow: expected number of rows (used for validation and row.names)
  cpp11::writable::data_frame build(int nrow) {
    // Validate column lengths match expected row count
    for (R_xlen_t i = 0; i < data.size(); ++i) {
      SEXP col = data[i];
      R_xlen_t col_len = Rf_xlength(col);
      if (col_len != nrow) {
        cpp11::stop("Column length mismatch: column %d has %d rows, expected %d",
                    static_cast<int>(i), static_cast<int>(col_len), nrow);
      }
    }

    // Set names
    cpp11::writable::strings name_vec(names.size());
    for (size_t i = 0; i < names.size(); ++i) {
      name_vec[i] = names[i];
    }
    data.attr("names") = name_vec;

    // Set row names
    SEXP row_names = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(row_names)[0] = NA_INTEGER;
    INTEGER(row_names)[1] = -nrow;
    data.attr("row.names") = cpp11::sexp(row_names);
    UNPROTECT(1);

    // Create the data frame first, then set tibble class
    // (cpp11::writable::data_frame may reset class to just "data.frame")
    cpp11::writable::data_frame result(data);
    result.attr("class") = cpp11::writable::strings({"tbl_df", "tbl", "data.frame"});

    return result;
  }
};

// Helper to compare rows between two dataframes
inline bool compare_rows(const cpp11::data_frame& df_x, const cpp11::data_frame& df_y, int idx_x,
                         int idx_y) {
  int ncols = df_x.size();
  cpp11::strings names_x(df_x.attr("names"));
  cpp11::strings names_y(df_y.attr("names"));

  for (int i = 0; i < ncols; ++i) {
    std::string col_name(names_x[i]);
    if (col_name == ".rows")
      continue;

    // Find matching column in y
    int y_col = -1;
    for (int j = 0; j < df_y.size(); ++j) {
      if (std::string(names_y[j]) == col_name) {
        y_col = j;
        break;
      }
    }
    if (y_col < 0) {
      cpp11::stop("Column not found in y dataframe");
    }

    // Compare values based on column type
    SEXP col_x_sexp = df_x[i];
    SEXP col_y_sexp = df_y[y_col];

    switch (TYPEOF(col_x_sexp)) {
      case STRSXP: {
        cpp11::strings col_x(col_x_sexp);
        cpp11::strings col_y(col_y_sexp);
        if (std::string(col_x[idx_x]) != std::string(col_y[idx_y])) {
          return false;
        }
        break;
      }
      case REALSXP: {
        cpp11::doubles col_x(col_x_sexp);
        cpp11::doubles col_y(col_y_sexp);
        if (col_x[idx_x] != col_y[idx_y]) {
          return false;
        }
        break;
      }
      case INTSXP: {
        cpp11::integers col_x(col_x_sexp);
        cpp11::integers col_y(col_y_sexp);
        if (col_x[idx_x] != col_y[idx_y]) {
          return false;
        }
        break;
      }
      case LGLSXP: {
        cpp11::logicals col_x(col_x_sexp);
        cpp11::logicals col_y(col_y_sexp);
        if (static_cast<int>(col_x[idx_x]) != static_cast<int>(col_y[idx_y])) {
          return false;
        }
        break;
      }
      default:
        cpp11::stop("Unsupported column type in group comparison");
    }
  }
  return true;
}

}  // namespace valr

#endif  // VALR_DATAFRAME_HPP
