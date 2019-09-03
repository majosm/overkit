// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COMM_HPP_INCLUDED
#define OVK_CORE_COMM_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {

class comm_view {

public:

  comm_view();
  template <typename T, OVK_FUNCDECL_REQUIRES(!std::is_same<T, comm_view>::value &&
    std::is_convertible<T, MPI_Comm>::value)> comm_view(const T &Comm);

  void Reset() { *this = comm_view(); }

  explicit operator bool() const { return Comm_ != MPI_COMM_NULL; }

  MPI_Comm Get() const { return Comm_; }
  operator MPI_Comm() const { return Comm_; }

  int Size() const;
  int Rank() const;

private:

  MPI_Comm Comm_;

  friend class core::test_helper<comm_view>;

};

class comm {

public:

  comm() = default;
  explicit comm(MPI_Comm Comm);

  comm(const comm &Other) = delete;
  comm(comm &&Other) noexcept;

  comm &operator=(comm Other) noexcept;

  void Reset();

  explicit operator bool() const { return static_cast<bool>(View_); }

  MPI_Comm Get() const { return View_.Get(); }
  operator MPI_Comm() const { return static_cast<MPI_Comm>(View_); }

  int Size() const { return View_.Size(); }
  int Rank() const { return View_.Rank(); }

  MPI_Comm Release();

private:

  struct resource {
    MPI_Comm Comm_;
    resource(MPI_Comm Comm):
      Comm_(Comm)
    {}
    ~resource() noexcept {
      if (Comm_ != MPI_COMM_NULL) {
        MPI_Comm_free(&Comm_);
      }
    }
  };

  std::unique_ptr<resource> Resource_;
  comm_view View_;

  friend class core::test_helper<comm>;

};

inline bool operator==(const comm_view &Left, const comm_view &Right);
inline bool operator!=(const comm_view &Left, const comm_view &Right);

comm DuplicateComm(comm_view Comm);

comm CreateSubsetComm(comm_view Comm, bool InSubset);

comm CreateCartComm(comm_view Comm, int NumDims, const tuple<int> &Dims, const tuple<bool>
  &Periodic, bool AllowReorder=true);
bool IsCartComm(comm_view Comm);
int GetCartCommDimension(comm_view Comm);
tuple<int> GetCartCommDims(comm_view Comm);
tuple<bool> GetCartCommPeriodic(comm_view Comm);
tuple<int> GetCartCommCoords(comm_view Comm);

}

#include <ovk/core/Comm.inl>

#endif
