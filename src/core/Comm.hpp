// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COMM_HPP_INCLUDED
#define OVK_CORE_COMM_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {

class comm {

public:

  comm() = default;
  explicit comm(MPI_Comm Comm, bool DuplicateComm=true);

  void Reset();

  explicit operator bool() const { return static_cast<bool>(Info_); }

  MPI_Comm Get() const { return Info_->Comm; }
  operator MPI_Comm() const { return Info_->Comm; }

  int Size() const { return Info_->Size; }
  int Rank() const { return Info_->Rank; }

private:

  struct comm_info {
    MPI_Comm Comm;
    int Size;
    int Rank;
  };

  // Use a shared pointer so we don't have to implement reference counting ourselves
  std::shared_ptr<comm_info> Info_;

};

inline bool operator==(const comm &Left, const comm &Right);
inline bool operator!=(const comm &Left, const comm &Right);

inline comm DuplicateComm(const comm &Comm);

inline comm CreateSubsetComm(const comm &Comm, bool InSubset);

inline comm CreateCartComm(const comm &Comm, int NumDims, const tuple<int> &Dims, const tuple<bool>
  &Periodic, bool AllowReorder=true);
inline bool IsCartComm(const comm &Comm);
inline int GetCartCommDimension(const comm &Comm);
inline tuple<int> GetCartCommDims(const comm &Comm);
inline tuple<bool> GetCartCommPeriodic(const comm &Comm);
inline tuple<int> GetCartCommCoords(const comm &Comm);

}}

#include <ovk/core/Comm.inl>

#endif
