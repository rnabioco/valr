#ifndef dplyr_tools_GroupedDataFrame_H
#define dplyr_tools_GroupedDataFrame_H

#include <dplyr/registration.h>
#include <tools/SlicingIndex.h>
#include <tools/VectorView.h>

#include <tools/SymbolVector.h>
#include <tools/SymbolMap.h>
#include <tools/bad.h>

// namespace dplyr {
//
// class GroupedDataFrame;
//
// class GroupedDataFrameIndexIterator {
// public:
//   GroupedDataFrameIndexIterator(const GroupedDataFrame& gdf_);
//
//   GroupedDataFrameIndexIterator& operator++();
//
//   GroupedSlicingIndex operator*() const;
//
//   int i;
//   const GroupedDataFrame& gdf;
//   ListView indices;
// };
//
// class GroupedDataFrame {
// public:
//   typedef GroupedDataFrameIndexIterator group_iterator;
//   typedef GroupedSlicingIndex slicing_index;
//
//   GroupedDataFrame(DataFrame x);
//   GroupedDataFrame(DataFrame x, const GroupedDataFrame& model);
//
//   group_iterator group_begin() const {
//     return GroupedDataFrameIndexIterator(*this);
//   }
//
//   SymbolString symbol(int i) const {
//     return symbols.get_name(i);
//   }
//
//   DataFrame& data() {
//     return data_;
//   }
//   const DataFrame& data() const {
//     return data_;
//   }
//
//   inline int size() const {
//     return data_.size();
//   }
//
//   inline int ngroups() const {
//     return groups.nrow();
//   }
//
//   inline int nvars() const {
//     return nvars_ ;
//   }
//
//   inline int nrows() const {
//     return data_.nrows();
//   }
//
//   inline SEXP label(int i) const {
//     return groups[i];
//   }
//
//   inline bool has_group(const SymbolString& g) const {
//     return symbols.has(g);
//   }
//
//   inline SEXP indices() const {
//     return groups[groups.size() - 1] ;
//   }
//
//   inline SymbolVector get_vars() const {
//     return symbols.get_names();
//   }
//
//   static SymbolVector group_vars(SEXP x);
//
//   inline const DataFrame& group_data() const {
//     return groups;
//   }
//
//   template <typename Data>
//   static void strip_groups(Data& x) {
//     x.attr("groups") = R_NilValue;
//   }
//
//   template <typename Data>
//   static void set_groups(Data& x, SEXP groups) {
//     x.attr("groups") = groups;
//   }
//
//   template <typename Data1, typename Data2>
//   static void copy_groups(Data1& x, const Data2& y) {
//     x.attr("groups") = y.attr("groups");
//   }
//
//   static inline CharacterVector classes() {
//     return Rcpp::CharacterVector::create("grouped_df", "tbl_df", "tbl", "data.frame");
//   }
//
// private:
//
//   DataFrame data_;
//   SymbolMap symbols;
//   DataFrame groups;
//   int nvars_;
//
// };
//
// inline GroupedDataFrameIndexIterator::GroupedDataFrameIndexIterator(const GroupedDataFrame& gdf_) :
//   i(0), gdf(gdf_), indices(gdf.indices()) {}
//
// inline GroupedDataFrameIndexIterator& GroupedDataFrameIndexIterator::operator++() {
//   i++;
//   return *this;
// }
//
// inline GroupedSlicingIndex GroupedDataFrameIndexIterator::operator*() const {
//   return GroupedSlicingIndex(indices[i], i);
// }
//
// }
//



// wrapper to call dplyr DataFrameSubsetVisitors
// class GroupedSlicingIndex {
// public:
//   GroupedSlicingIndex(): data(), group_index(-1) {
//     R_PreserveObject(data);
//   }
//
//   ~GroupedSlicingIndex() {
//     if (group_index == -1) {
//       R_ReleaseObject(data);
//     }
//   }
//
//   GroupedSlicingIndex(SEXP data_, int group_) : data(data_), group_index(group_) {}
//
//   int size() const {
//     return data.size();
//   }
//
//   int operator[](int i) const {
//     return data[i] - 1;
//   }
//
//   int group() const {
//     return group_index;
//   }
//
//   inline operator SEXP() const {
//     return data;
//   }
//
// private:
//   // in general we don't need to protect data because
//   // it is already protected by the .rows column of the grouped_df
//   //
//   // but we do when using the default constructor, hence the
//   // R_PreserveObject / R_ReleaseObject above
//   Rcpp::IntegerVectorView data;
//
//   int group_index;
// };

// redefine groupeddataframe to simplify need for dplyr cpp functions
namespace dplyr {

class GroupedDataFrame;

class GroupedDataFrameIndexIterator {
public:
  GroupedDataFrameIndexIterator(const GroupedDataFrame& gdf_);

  GroupedDataFrameIndexIterator& operator++();

  GroupedSlicingIndex operator*() const;

  int i;
  const GroupedDataFrame& gdf;
  ListView indices;
};

class GroupedDataFrame {
public:
  typedef GroupedDataFrameIndexIterator group_iterator;
  typedef GroupedSlicingIndex slicing_index;

  GroupedDataFrame(DataFrame x) ;
  group_iterator group_begin() const {
    return GroupedDataFrameIndexIterator(*this);
  }

  DataFrame& data() {
    return data_;
  }
  const DataFrame& data() const {
    return data_;
  }

  inline int size() const {
    return data_.size();
  }

  inline int ngroups() const {
    return groups.nrows();
  }

  inline int nrows() const {
    return data_.nrows();
  }

  inline const DataFrame& group_data() const {
    return groups ;
  }

  inline SEXP indices() const {
    return groups[groups.size() - 1] ;
  }

  template <typename Data>
  static void strip_groups(Data& x) {
      x.attr("groups") = R_NilValue;
  }

private:
  DataFrame data_;
  DataFrame groups;
};


inline GroupedDataFrameIndexIterator::GroupedDataFrameIndexIterator(const GroupedDataFrame& gdf_) :
  i(0), gdf(gdf_), indices(gdf.indices()) {}

inline GroupedDataFrameIndexIterator& GroupedDataFrameIndexIterator::operator++() {
  i++;
  return *this;
}

inline GroupedSlicingIndex GroupedDataFrameIndexIterator::operator*() const {
  return GroupedSlicingIndex(indices[i], i);
}

}

namespace Rcpp {
using namespace dplyr;

template <>
inline bool is<GroupedDataFrame>(SEXP x) {
  return Rf_inherits(x, "grouped_df");
}

}

#endif
