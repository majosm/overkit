// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline comm_view::comm_view():
  Comm_(MPI_COMM_NULL)
{}

template <typename T, OVK_FUNCDEF_REQUIRES(!std::is_same<T,comm_view>::value &&
  std::is_convertible<T, MPI_Comm>::value)> comm_view::comm_view(const T &Comm):
  Comm_(Comm)
{}

inline int comm_view::Size() const {

  int Size;
  MPI_Comm_size(Comm_, &Size);

  return Size;

}

inline int comm_view::Rank() const {

  int Rank;
  MPI_Comm_rank(Comm_, &Rank);

  return Rank;

}

inline comm::comm(MPI_Comm Comm) {

  if (Comm != MPI_COMM_NULL) {
    Resource_.reset(new resource(Comm));
    View_ = Comm;
  }

}

inline comm::comm(comm &&Other) noexcept:
  Resource_(std::move(Other.Resource_)),
  View_(Other.View_)
{
  Other.View_.Reset();
}

inline comm &comm::operator=(comm Other) noexcept {

  using std::swap;

  swap(Resource_, Other.Resource_);
  swap(View_, Other.View_);

  return *this;

}

inline void comm::Reset() {

  Resource_.reset();
  View_.Reset();

}

inline MPI_Comm comm::Release() {

  MPI_Comm Comm = Resource_->Comm_;
  Resource_->Comm_ = MPI_COMM_NULL;

  Resource_.reset();
  View_.Reset();

  return Comm;

}

inline bool operator==(const comm_view &Left, const comm_view &Right) {
  return Left.Get() == Right.Get();
}
inline bool operator!=(const comm_view &Left, const comm_view &Right) {
  return !(Left == Right);
}

inline comm DuplicateComm(comm_view Comm) {

  MPI_Comm DuplicatedCommRaw;
  MPI_Comm_dup(Comm.Get(), &DuplicatedCommRaw);

  return comm(DuplicatedCommRaw);

}

inline comm CreateSubsetComm(comm_view Comm, bool InSubset) {

  MPI_Comm CommRaw = Comm.Get();

  MPI_Comm SubsetCommRaw;
  MPI_Comm_split(CommRaw, InSubset ? 0 : MPI_UNDEFINED, Comm.Rank(), &SubsetCommRaw);

  return comm(SubsetCommRaw);

}

inline comm CreateCartComm(comm_view Comm, int NumDims, const tuple<int> &Dims_, const tuple<bool>
  &Periodic, bool AllowReorder) {

  tuple<int> Dims = Dims_;
  MPI_Dims_create(Comm.Size(), NumDims, Dims.Data());

  tuple<int> PeriodicInt = tuple<int>(Periodic);

  MPI_Comm CartCommRaw;
  MPI_Cart_create(Comm, NumDims, Dims.Data(), PeriodicInt.Data(), AllowReorder, &CartCommRaw);

  return comm(CartCommRaw);

}

inline bool IsCartComm(comm_view Comm) {

  int TopologyType;
  MPI_Topo_test(Comm, &TopologyType);

  return TopologyType == MPI_CART;

}

inline int GetCartCommDimension(comm_view Comm) {

  int NumDims;
  MPI_Cartdim_get(Comm, &NumDims);

  return NumDims;

}

inline tuple<int> GetCartCommDims(comm_view Comm) {

  int NumDims = GetCartCommDimension(Comm);

  tuple<int> CartDims = MakeUniformTuple<int>(NumDims, 1, 1);
  tuple<int> PeriodicInt;
  tuple<int> CartCoords;
  MPI_Cart_get(Comm, NumDims, CartDims.Data(), PeriodicInt.Data(), CartCoords.Data());

  return CartDims;

}

inline tuple<bool> GetCartCommPeriodic(comm_view Comm) {

  int NumDims = GetCartCommDimension(Comm);

  tuple<int> CartDims;
  tuple<int> PeriodicInt = MakeUniformTuple<int>(NumDims, 0);
  tuple<int> CartCoords;
  MPI_Cart_get(Comm, NumDims, CartDims.Data(), PeriodicInt.Data(), CartCoords.Data());

  return tuple<bool>(PeriodicInt);

}

inline tuple<int> GetCartCommCoords(comm_view Comm) {

  int NumDims = GetCartCommDimension(Comm);

  tuple<int> CartDims;
  tuple<int> PeriodicInt;
  tuple<int> CartCoords = MakeUniformTuple<int>(NumDims, 0);
  MPI_Cart_get(Comm, NumDims, CartDims.Data(), PeriodicInt.Data(), CartCoords.Data());

  return CartCoords;

}

}
