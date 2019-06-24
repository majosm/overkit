// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

// Intel 17 didn't like this for some reason
// template <typename ArrayType, OVK_FUNCDEF_REQUIRES(IsArray<ArrayType>() && ArrayHasFootprint<
//   ArrayType, MAX_DIMS, array_layout::GRID>())> request halo::Exchange(ArrayType &Array) const {
template <typename ArrayType, OVK_FUNCDEF_REQUIRES(IsArray<ArrayType>() && ArrayRank<ArrayType>()
  == MAX_DIMS && ArrayLayout<ArrayType>() == array_layout::GRID)> request halo::Exchange(ArrayType
  &Array) const {

  using value_type = array_value_type<ArrayType>;

  OVK_DEBUG_ASSERT(IsSupportedDataType<value_type>(), "Unsupported data type.");

  core::profiler &Profiler = Context_->core_Profiler();

  Profiler.StartSync(TOTAL_TIME, Comm_);
  Profiler.Start(EXCHANGE_TIME);

  data_type DataType = GetDataType<value_type>();

  array<halo_exchanger> &HaloExchangersForType = HaloExchangers_.Get(int(DataType));

  int iHaloExchanger = 0;
  while (iHaloExchanger < HaloExchangersForType.Count() && HaloExchangersForType(iHaloExchanger).
    Active()) {
    ++iHaloExchanger;
  }
  if (iHaloExchanger == HaloExchangersForType.Count()) {
    Profiler.Stop(EXCHANGE_TIME);
    Profiler.Start(SETUP_TIME);
    HaloExchangersForType.Append(halo_internal::halo_exchanger_for_type<value_type>(*Context_,
      Comm_, HaloMap_));
    Profiler.Stop(SETUP_TIME);
    Profiler.Start(EXCHANGE_TIME);
  }
  halo_exchanger &HaloExchanger = HaloExchangersForType(iHaloExchanger);

  auto EndProfiles = OnScopeExit([&] {
    Profiler.Stop(EXCHANGE_TIME);
    Profiler.Stop(TOTAL_TIME);
  });

  return HaloExchanger.Exchange(ArrayData(Array));

}

namespace halo_internal {

template <typename T> halo_exchanger_for_type<T>::halo_exchanger_for_type(context &Context,
  comm_view Comm, const halo_map &HaloMap):
  FloatingRefGenerator_(*this),
  Context_(Context.GetFloatingRef()),
  Comm_(Comm),
  HaloMap_(HaloMap.GetFloatingRef())
{

  int NumNeighbors = HaloMap.NeighborRanks().Count();

  SendBuffers_.Resize({NumNeighbors});
  RecvBuffers_.Resize({NumNeighbors});
  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    SendBuffers_(iNeighbor).Resize({HaloMap.NeighborSendIndices(iNeighbor).Count()});
    RecvBuffers_(iNeighbor).Resize({HaloMap.NeighborRecvIndices(iNeighbor).Count()});
  }

  MPIRequests_.Reserve(2*NumNeighbors);

}

template <typename T> request halo_exchanger_for_type<T>::Exchange(value_type *ArrayData) {

  const halo_map &HaloMap = *HaloMap_;
  const array<int> &NeighborRanks = HaloMap.NeighborRanks();
  const array<long long> &LocalToLocalSourceIndices = HaloMap.LocalToLocalSourceIndices();
  const array<long long> &LocalToLocalDestIndices = HaloMap.LocalToLocalDestIndices();

  core::profiler &Profiler = Context_->core_Profiler();

  int NumNeighbors = HaloMap.NeighborRanks().Count();

  MPI_Datatype DataType = GetMPIDataType<mpi_value_type>();

  MPIRequests_.Clear();

  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    MPI_Request &Request = MPIRequests_.Append();
    Profiler.Start(MPI_TIME);
    MPI_Irecv(RecvBuffers_(iNeighbor).Data(), RecvBuffers_(iNeighbor).Count(), DataType,
      NeighborRanks(iNeighbor), 0, Comm_, &Request);
    Profiler.Stop(MPI_TIME);
  }

  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    Profiler.Start(PACK_TIME);
    const array<long long> &SendIndices = HaloMap.NeighborSendIndices(iNeighbor);
    for (long long iSendPoint = 0; iSendPoint < SendIndices.Count(); ++iSendPoint) {
      long long iPoint = SendIndices(iSendPoint);
      SendBuffers_(iNeighbor)(iSendPoint) = mpi_value_type(ArrayData[iPoint]);
    }
    Profiler.Stop(PACK_TIME);
    MPI_Request &Request = MPIRequests_.Append();
    Profiler.Start(MPI_TIME);
    MPI_Isend(SendBuffers_(iNeighbor).Data(), SendBuffers_(iNeighbor).Count(), DataType,
      NeighborRanks(iNeighbor), 0, Comm_, &Request);
    Profiler.Stop(MPI_TIME);
  }

  Profiler.Start(PACK_TIME);
  Profiler.Start(UNPACK_TIME);

  long long NumLocalToLocal = LocalToLocalSourceIndices.Count();
  for (long long iLocalToLocal = 0; iLocalToLocal < NumLocalToLocal; ++iLocalToLocal) {
    long long iSource = LocalToLocalSourceIndices(iLocalToLocal);
    long long iDest = LocalToLocalDestIndices(iLocalToLocal);
    ArrayData[iDest] = ArrayData[iSource];
  }

  Profiler.Stop(PACK_TIME);
  Profiler.Stop(UNPACK_TIME);

  Active_ = true;

  return exchange_request(*this, ArrayData);

}

template <typename T> halo_exchanger_for_type<T>::exchange_request::exchange_request(
  halo_exchanger_for_type &HaloExchanger, value_type *ArrayData):
  HaloExchanger_(HaloExchanger.FloatingRefGenerator_.Generate()),
  ArrayData_(ArrayData)
{}

template <typename T> void halo_exchanger_for_type<T>::exchange_request::OnMPIRequestComplete(int
  iMPIRequest) {

  halo_exchanger_for_type &HaloExchanger = *HaloExchanger_;
  const halo_map &HaloMap = *HaloExchanger.HaloMap_;

  profiler &Profiler = HaloExchanger.Context_->core_Profiler();

  if (iMPIRequest < HaloExchanger.RecvBuffers_.Count()) {

    Profiler.Start(UNPACK_TIME);

    int iNeighbor = iMPIRequest;
    const array<long long> &RecvIndices = HaloMap.NeighborRecvIndices(iNeighbor);
    for (long long iRecvPoint = 0; iRecvPoint < RecvIndices.Count(); ++iRecvPoint) {
      long long iPoint = RecvIndices(iRecvPoint);
      ArrayData_[iPoint] = value_type(HaloExchanger.RecvBuffers_(iNeighbor)(iRecvPoint));
    }

    Profiler.Stop(UNPACK_TIME);

  }

}

template <typename T> void halo_exchanger_for_type<T>::exchange_request::OnComplete() {

  HaloExchanger_->Active_ = false;

}

}

}}
