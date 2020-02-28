// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename F, OVK_FUNCDEF_REQUIRES(IsCallableWith<F &&>())> auto Serialize(comm_view Comm,
  F &&Func) -> decltype(std::forward<F>(Func)()) {

  for (int OtherRank = 0; OtherRank < Comm.Rank(); ++OtherRank) {
    MPI_Barrier(Comm);
  }

  auto DoRemainingBarriers = OnScopeExit([&Comm] {
    for (int OtherRank = Comm.Rank(); OtherRank < Comm.Size(); ++OtherRank) {
      MPI_Barrier(Comm);
    }
  });

  return std::forward<F>(Func)();

}

}}
