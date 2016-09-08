#ifndef __INTERVAL_TREE_H
#define __INTERVAL_TREE_H

// Modified from IntervalTree.h from EKG
// https://github.com/ekg/intervaltree

#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

template <class T, typename K = int>
class Interval {
public:
    K start;
    K stop;
    T value;
    Interval(K s, K e, const T& v)
        : start(s)
        , stop(e)
        , value(v)
    { }
};

template <class T, typename K>
K intervalStart(const Interval<T,K>& i) {
    return i.start;
}

template <class T, typename K>
K intervalStop(const Interval<T,K>& i) {
    return i.stop;
}

template <class T, typename K>
  std::ostream& operator<<(std::ostream& out, Interval<T,K>& i) {
    out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
    return out;
}

template <class T, typename K = int>
class IntervalStartSorter {
public:
    bool operator() (const Interval<T,K>& a, const Interval<T,K>& b) {
        return a.start < b.start;
    }
};

template <class T, typename K = int>
class IntervalTree {

public:
    typedef Interval<T,K> interval;
    typedef std::vector<interval> intervalVector;
    typedef IntervalTree<T,K> intervalTree;

    intervalVector intervals;
    std::unique_ptr<intervalTree> left;
    std::unique_ptr<intervalTree> right;
    K center;

    IntervalTree<T,K>(void)
        : left(nullptr)
        , right(nullptr)
        , center(0) 
    { }

private:
    std::unique_ptr<intervalTree> copyTree(const intervalTree& orig){
        return std::unique_ptr<intervalTree>(new intervalTree(orig));
    }
  
public:

    IntervalTree<T,K>(const intervalTree& other)
    :   intervals(other.intervals),
        left(other.left ? copyTree(*other.left) : nullptr),
        right(other.right ? copyTree(*other.right) : nullptr),
        center(other.center)
    {
    }

public:

    IntervalTree<T,K>& operator=(const intervalTree& other) {
        center = other.center ;
        intervals = other.intervals;
        left = other.left ? copyTree(*other.left) : nullptr;
        right = other.right ? copyTree(*other.right) : nullptr;
        return *this;
    }

    // Note: changes the order of ivals
    IntervalTree<T,K>(
            intervalVector& ivals,
            std::size_t depth = 16,
            std::size_t minbucket = 64,
            K leftextent = 0,
            K rightextent = 0,
            std::size_t maxbucket = 512
            )
        : left(nullptr)
        , right(nullptr)
    {

        --depth;
        IntervalStartSorter<T,K> intervalStartSorter;
        if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
            std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
            intervals = ivals;
            center = ivals.at(ivals.size() / 2).start;
        } else {
            if (leftextent == 0 && rightextent == 0) {
                // sort intervals by start
              std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
            }

            K leftp = 0;
            K rightp = 0;
            K centerp = 0;

            if (leftextent || rightextent) {
                leftp = leftextent;
                rightp = rightextent;
            } else {
                leftp = ivals.front().start;
                std::vector<K> stops;
                stops.resize(ivals.size());
                transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T,K>);
                rightp = *max_element(stops.begin(), stops.end());
            }

            //centerp = ( leftp + rightp ) / 2;
            centerp = ivals.at(ivals.size() / 2).start;
            center = centerp;

            intervalVector lefts;
            intervalVector rights;

            for (typename intervalVector::const_iterator i = ivals.begin(); i != ivals.end(); ++i) {
                const interval& interval = *i;
                if (interval.stop < center) {
                    lefts.push_back(interval);
                } else if (interval.start > center) {
                    rights.push_back(interval);
                } else {
                    intervals.push_back(interval);
                }
            }

            if (!lefts.empty()) {
                left = std::unique_ptr<intervalTree>(new intervalTree(lefts, depth, minbucket, leftp, centerp));
            }
            if (!rights.empty()) {
                right = std::unique_ptr<intervalTree>(new intervalTree(rights, depth, minbucket, centerp, rightp));
            }
        }
    }

    intervalVector findOverlapping(K start, K stop) const {
    	intervalVector ov;
    	this->findOverlapping(start, stop, ov);
    	return ov;
    }

    void findOverlapping(K start, K stop, intervalVector& overlapping) const {
        if (!intervals.empty() && ! (stop < intervals.front().start)) {
            for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
                const interval& interval = *i;
                if (interval.stop >= start && interval.start <= stop) {
                    overlapping.push_back(interval);
                }
            }
        }

        if (left && start <= center) {
            left->findOverlapping(start, stop, overlapping);
        }

        if (right && stop >= center) {
            right->findOverlapping(start, stop, overlapping);
        }

    }

    intervalVector findContained(K start, K stop) const {
    	intervalVector contained;
    	this->findContained(start, stop, contained);
    	return contained;
    }

    void findContained(K start, K stop, intervalVector& contained) const {
        if (!intervals.empty() && ! (stop < intervals.front().start)) {
            for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
                const interval& interval = *i;
                if (interval.start >= start && interval.stop <= stop) {
                    contained.push_back(interval);
                }
            }
        }

        if (left && start <= center) {
            left->findContained(start, stop, contained);
        }

        if (right && stop >= center) {
            right->findContained(start, stop, contained);
        }

    }

  intervalVector findClosest(K start, K stop) const {
    intervalVector closest ;
    std::pair<int, intervalVector> min_dist_l, min_dist_r;
    // initiatialize left and right distances to minimize for closest
  //  min_dist_l = start
  //  min_dist_r = intervals.front().start)
    this->findClosest(start, stop, closest, min_dist_l, min_dist_r) ;
    return closest;
  }
  
  void findClosest(K start, K stop, intervalVector&  closest, 
                   std::pair<int, intervalVector>& min_dist_l, 
                   std::pair<int, intervalVector>& min_dist_r) const {
  //  Rcpp::Rcout << min_dist_l.first << " = min_dist_l" << std::endl ; 
  //    Rcpp::Rcout << stop << " = current_test_ivl_stop" << std::endl ;
  //    Rcpp::Rcout << intervals.front().start << " = current_ivl_tree_node_min" << std::endl ;
  //    Rcpp::Rcout << intervals.back().stop << " = current_ivl_tree_node_max" << std::endl ;
  //    Rcpp::Rcout << center << " = current_ivl_tree_node_center" << std::endl ; 
  //    Rcpp::Rcout << this << " = current_ivl_tree_node_address" << std::endl ;
  //    Rcpp::Rcout << (*this).center << " = current_ivl_tree_node_address_rereferenced_center" << std::endl ;
    // initialize distances from left and right not sure if necesary to have node specific max...
    int min_node_dist_l = min_dist_l.first ;
    int min_node_dist_r = min_dist_r.first ;
    min_dist_l.second.clear() ;
    
    //  Rcpp::Rcout << min_dist_l.first << " min dist" << std::endl;
    //  Rcpp::Rcout << start << " = ivl_test_start " << stop << " ivl_test_Stop " << center << " ivl center "  << std::endl ; 
    if (!intervals.empty() && ! (stop < intervals.front().start)) {
      for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
        const interval& interval = *i;
        //   Rcpp::Rcout << interval.start << " = ivl_tree_start " << interval.stop << " ivl_tree_Stop " << std::endl ; 
        if (interval.stop >= start && interval.start <= stop) {
          closest.push_back(interval);
          //     Rcpp::Rcout << "intersect" << std::endl ; 
        } else if (stop < interval.start) {
          // finddistance on left
          int ivl_dist_l = interval.start - stop ;
          //        Rcpp::Rcout << ivl_dist_l << " = ivl_dist_l smaller" << std::endl ; 
          // if smaller than best min dist found update dictionary with intervals
          if (ivl_dist_l < min_node_dist_l) {
            //      Rcpp::Rcout << ivl_dist_l << " = ivl_dist_l smaller" << std::endl ; 
            //        Rcpp::Rcout << interval.start << " " << interval.stop << " = interval smaller that " << start << " " << stop << std::endl ;
            //    Rcpp::Rcout << min_dist_l.first << " min dist loop before" << std::endl; 
            min_node_dist_l = ivl_dist_l;
            min_dist_l.first = min_node_dist_l ; 
            //    Rcpp::Rcout << min_dist_l.first << " min dist loop after" << std::endl; 
            min_dist_l.second.clear() ; 
            min_dist_l.second.push_back(interval) ;
            //    Rcpp::Rcout << min_dist_l.second[0] << " first interval loop after" << std::endl; 
          } else if (ivl_dist_l == min_node_dist_l) {
            min_dist_l.second.push_back(interval) ;
          }
        } else if (start > interval.stop) {
          // find distance on right
          int ivl_dist_r = start - interval.stop ;
          // if smaller than best min dist found update dictionary with intervals
          if (ivl_dist_r < min_node_dist_r) {
            //              Rcpp::Rcout << ivl_dist_r << " = ivl_dist_r smaller" << std::endl ; 
            //              Rcpp::Rcout << interval.start << " " << interval.stop << " = interval smaller that " << start << " " << stop << std::endl ;
            //Rcpp::Rcout << "right" << std::endl ;
            //Rcpp::Rcout << min_dist_r.first << " right min dist loop before " << std::endl; 
            min_node_dist_r = ivl_dist_r;
            min_dist_r.first = min_node_dist_r ; 
            //    Rcpp::Rcout << min_dist_r.first << " right min dist loop after " << std::endl; 
            min_dist_r.second.clear() ; 
            min_dist_r.second.push_back(interval) ;
            //   Rcpp::Rcout << min_dist_r.second[0] << " right first interval loop after" << std::endl;
          } else if (ivl_dist_r == min_node_dist_r) {
            min_dist_r.second.push_back(interval) ;
          }
        }
      }
        
        
    }
  
    
    
    if (left && start <= center) {
      // Rcpp::Rcout << "push to left tree" << std::endl ; 
      // Rcpp::Rcout  << (*left).center << " = left center" << std::endl ;
      //      Rcpp::Rcout  << start << " = start" << std::endl ;
      //        int left_val ;
      //      left_val = (*left).intervals.front().start ;
      //      Rcpp::Rcout << left_val << " left_tree_front_start" << std::endl ;
      //        Rcpp::Rcout << (*left).intervals.back().stop << " left_tree_end_stop" << std::endl ;
      left->findClosest(start, stop,  closest, min_dist_l, min_dist_r);
    }
    
    if (right && stop >= center) {
      // Rcpp::Rcout  << (right).center << " = right center" << std::endl ;
      //    Rcpp::Rcout  << stop << " = right stop" << std::endl ;
      //    int left_val ;
      //    left_val = (*right).intervals.front().start ;
      //    Rcpp::Rcout << left_val << " right_tree_front_start" << std::endl ;
      //        Rcpp::Rcout << (*right).intervals.back().stop << " right_tree_end_stop" << std::endl ;
      // Rcpp::Rcout << "push to right tree" << std::endl ; 
      right->findClosest(start, stop,  closest, min_dist_l, min_dist_r);
    }  
    
    //
    // Finally report all of the non-overlapping closest interva
   // Rcpp::Rcout << min_dist_l.second.size() << "left interval sizes" << std::endl ; 
   //  Rcpp::Rcout << min_dist_r.second.size() << "right interval sizes" << std::endl ; 
   // Rcpp::Rcout << min_dist_l.first << " left max " << min_dist_r.first << " right max " << std::endl ;
    if (min_dist_l.first < min_dist_r.first ){
      closest.insert(closest.end(), min_dist_l.second.begin(), min_dist_l.second.end())  ;
    } else if (min_dist_l.first == min_dist_r.first) {
      min_dist_l.second.insert(min_dist_l.second.end(), min_dist_r.second.begin(), min_dist_r.second.end()) ;
      closest.insert(closest.end(), min_dist_l.second.begin(), min_dist_l.second.end())  ;
    } else {
      closest.insert(closest.end(), min_dist_r.second.begin(), min_dist_r.second.end())  ;
    }
    //Rcpp::Rcout << closest.size()  << " final closeset size" << std::endl ; 
  }
    ~IntervalTree(void) = default;

};

#endif
