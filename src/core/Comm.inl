// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

inline comm::comm(MPI_Comm Comm, bool DuplicateComm) {

  if (Comm != MPI_COMM_NULL) {
    auto FreeCommAndDelete = [](comm_info *Info) {
      MPI_Comm_free(&Info->Comm);
      delete Info;
    };
    Info_ = std::shared_ptr<comm_info>(new comm_info(), FreeCommAndDelete);
    if (DuplicateComm) {
      MPI_Comm_dup(Comm, &Info_->Comm);
    } else {
      Info_->Comm = Comm;
    }
    MPI_Comm_size(Info_->Comm, &Info_->Size);
    MPI_Comm_rank(Info_->Comm, &Info_->Rank);
  }

}

inline void comm::Reset() {

  Info_.reset();

}

inline comm DuplicateComm(const comm &Comm) {

  return comm(Comm.Get());

}

inline comm CreateSubsetComm(const comm &Comm, bool InSubset) {

  MPI_Comm CommRaw = Comm.Get();

  MPI_Comm SubsetCommRaw;
  MPI_Comm_split(CommRaw, InSubset ? 0 : MPI_UNDEFINED, Comm.Rank(), &SubsetCommRaw);

  return comm(SubsetCommRaw, false);

}

}}
