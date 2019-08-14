// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_UNION_FIND_HPP_INCLUDED
#define OVK_CORE_UNION_FIND_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/PointerIterator.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Set.hpp>

#include <utility>

namespace ovk {
namespace core {

class union_find {

public:

  // Elements aren't mutable, so iterator is const
  using iterator = pointer_iterator<union_find, const int *>;
  using const_iterator = iterator;

  union_find() = default;

  void Reserve(int Count) {
    Elements_.Reserve(Count);
    ParentIndices_.Reserve(Count);
    Ranks_.Reserve(Count);
  }

  void Insert(int Element) {
    auto Iter = Elements_.LowerBound(Element);
    if (Iter == Elements_.End() || *Iter > Element) {
      Iter = Elements_.Insert(Iter, Element);
      int iElement = int(Iter - Elements_.Begin());
      for (int iOtherElement = 0; iOtherElement < ParentIndices_.Count(); ++iOtherElement) {
        if (ParentIndices_(iOtherElement) >= iElement) {
          ++ParentIndices_(iOtherElement);
        }
      }
      ParentIndices_.Insert(iElement, iElement);
      Ranks_.Insert(iElement, 0);
    }
  }

  int Find(int Element) {
    auto Iter = Elements_.Find(Element);
    if (Iter != Elements_.End()) {
      int iElement = Iter - Elements_.Begin();
      int iRoot = Find_(iElement);
      return Elements_[iRoot];
    } else {
      return -1;
    }
  }

  int Union(int LeftElement, int RightElement) {
    auto LeftIter = Elements_.Find(LeftElement);
    auto RightIter = Elements_.Find(RightElement);
    int iLeft = -1;
    int iRight = -1;
    if (LeftIter != Elements_.End()) iLeft = LeftIter - Elements_.Begin();
    if (RightIter != Elements_.End()) iRight = RightIter - Elements_.Begin();
    if (iLeft >= 0 && iRight >= 0) {
      int iLeftRoot = Find_(iLeft);
      int iRightRoot = Find_(iRight);
      if (iLeftRoot != iRightRoot) {
        if (Ranks_(iLeftRoot) < Ranks_(iRightRoot)) {
          ParentIndices_(iLeftRoot) = iRightRoot;
          return Elements_[iRightRoot];
        } else if (Ranks_(iLeftRoot) > Ranks_(iRightRoot)) {
          ParentIndices_(iRightRoot) = iLeftRoot;
          return Elements_[iLeftRoot];
        } else {
          ParentIndices_(iRightRoot) = iLeftRoot;
          ++Ranks_(iLeftRoot);
          return Elements_[iLeftRoot];
        }
      } else {
        return Elements_[iLeftRoot];
      }
    } else if (iLeft >= 0) {
      int iLeftRoot = Find_(iLeft);
      return Elements_[iLeftRoot];
    } else if (iRight >= 0) {
      int iRightRoot = Find_(iRight);
      return Elements_[iRightRoot];
    } else {
      return -1;
    }
  }

  // Shuffle things around so that root of each set is the smallest contained element
  void Relabel() {
    map<int,int> OldRootIndexToNewRootIndex;
    for (int iElement = 0; iElement < Elements_.Count(); ++iElement) {
      int iOldRoot = Find_(iElement);
      int &iNewRoot = OldRootIndexToNewRootIndex.Fetch(iOldRoot, iOldRoot);
      if (Elements_[iElement] < Elements_[iOldRoot]) {
        iNewRoot = iElement;
      }
    }
    for (int iElement = 0; iElement < Elements_.Count(); ++iElement) {
      ParentIndices_(iElement) = OldRootIndexToNewRootIndex(ParentIndices_(iElement));
      Ranks_(iElement) = ParentIndices_(iElement) == iElement ? 1 : 0;
    }
  }

  int Count() const { return int(Elements_.Count()); }

  iterator Begin() const { return iterator(Elements_.Data()); }
  iterator End() const { return iterator(Elements_.Data()+Elements_.Count()); }

  // Google Test doesn't use free begin/end functions and instead expects container to have
  // lowercase begin/end methods
  iterator begin() const { return Begin(); }
  iterator end() const { return End(); }

private:

  set<int> Elements_;
  array<int> ParentIndices_;
  array<int> Ranks_;

  int Find_(int iElement) {
    if (ParentIndices_(iElement) != iElement) {
      ParentIndices_(iElement) = Find_(ParentIndices_(iElement));
    }
    return ParentIndices_(iElement);
  }

};

inline union_find::iterator begin(const union_find &UnionFind) {
  return UnionFind.Begin();
}

inline union_find::iterator end(const union_find &UnionFind) {
  return UnionFind.End();
}

}}

#endif
