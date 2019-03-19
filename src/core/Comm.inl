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

inline bool operator==(const comm &Left, const comm &Right) {

  return Left.Get() == Right.Get();

}

inline bool operator!=(const comm &Left, const comm &Right) {

  return !(Left == Right);

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

inline comm CreateCartComm(const comm &Comm, int NumDims, const tuple<int> &Dims, const tuple<bool>
  &Periodic, bool AllowReorder) {

  tuple<int> PeriodicInt = tuple<int>(Periodic);
  MPI_Comm CartCommRaw;
  MPI_Cart_create(Comm, NumDims, Dims.Data(), PeriodicInt.Data(), AllowReorder, &CartCommRaw);

  return comm(CartCommRaw, false);

}

inline bool IsCartComm(const comm &Comm) {

  int TopologyType;
  MPI_Topo_test(Comm, &TopologyType);

  return TopologyType == MPI_CART;

}

inline int GetCartCommDimension(const comm &Comm) {

  int NumDims;
  MPI_Cartdim_get(Comm, &NumDims);

  return NumDims;

}

inline tuple<int> GetCartCommDims(const comm &Comm) {

  int NumDims = GetCartCommDimension(Comm);

  tuple<int> CartDims = MakeUniformTuple<int>(1);
  tuple<int> PeriodicInt;
  tuple<int> CartCoords;
  MPI_Cart_get(Comm, NumDims, CartDims.Data(), PeriodicInt.Data(), CartCoords.Data());

  return CartDims;

}

inline tuple<bool> GetCartCommPeriodic(const comm &Comm) {

  int NumDims = GetCartCommDimension(Comm);

  tuple<int> CartDims;
  tuple<int> PeriodicInt = MakeUniformTuple<int>(0);
  tuple<int> CartCoords;
  MPI_Cart_get(Comm, NumDims, CartDims.Data(), PeriodicInt.Data(), CartCoords.Data());

  return tuple<bool>(PeriodicInt);

}

inline tuple<int> GetCartCommCoords(const comm &Comm) {

  int NumDims = GetCartCommDimension(Comm);

  tuple<int> CartDims;
  tuple<int> PeriodicInt;
  tuple<int> CartCoords = MakeUniformTuple<int>(0);
  MPI_Cart_get(Comm, NumDims, CartDims.Data(), PeriodicInt.Data(), CartCoords.Data());

  return CartCoords;

}

}}
