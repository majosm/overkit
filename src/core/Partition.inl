// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace partition_internal {

template <typename FieldType, OVK_FUNCDEF_REQUIRES(core::IsField<FieldType>())> request
  halo::Exchange(FieldType &Field) const {

  using value_type = core::array_value_type<FieldType>;

  field_view<value_type> View(Field);

  return Exchange(View);

}

template <typename T, OVK_FUNCDEF_REQUIRES(!std::is_const<T>::value)> request
  halo::Exchange(field_view<T> View) const {

  OVK_DEBUG_ASSERT(core::IsSupportedDataType<T>(), "Unsupported data type.");

  core::profiler &Profiler = Context_->core_Profiler();

  Profiler.StartSync(TOTAL_TIME, Comm_);
  Profiler.Start(EXCHANGE_TIME);

  data_type DataType = core::GetDataType<T>();

  array<halo_exchanger> &HaloExchangersForType = HaloExchangers_.Fetch(int(DataType));

  int iHaloExchanger = 0;
  while (iHaloExchanger < HaloExchangersForType.Count() && HaloExchangersForType(iHaloExchanger).
    Active()) {
    ++iHaloExchanger;
  }
  if (iHaloExchanger == HaloExchangersForType.Count()) {
    Profiler.Stop(EXCHANGE_TIME);
    Profiler.Start(SETUP_TIME);
    HaloExchangersForType.Append(halo_internal::halo_exchanger_for_type<T>(*Context_, Comm_,
      HaloMap_));
    Profiler.Stop(SETUP_TIME);
    Profiler.Start(EXCHANGE_TIME);
  }
  halo_exchanger &HaloExchanger = HaloExchangersForType(iHaloExchanger);

  auto EndProfiles = core::OnScopeExit([&] {
    Profiler.Stop(EXCHANGE_TIME);
    Profiler.Stop(TOTAL_TIME);
  });

  return HaloExchanger.Exchange(View.Data());

}

namespace halo_internal {

template <typename T> halo_exchanger_for_type<T>::halo_exchanger_for_type(context &Context,
  comm_view Comm, const halo_map &HaloMap):
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

template <typename T> request halo_exchanger_for_type<T>::Exchange(value_type *FieldData) {

  const halo_map &HaloMap = *HaloMap_;
  const array<int> &NeighborRanks = HaloMap.NeighborRanks();
  const array<long long> &LocalToLocalSourceIndices = HaloMap.LocalToLocalSourceIndices();
  const array<long long> &LocalToLocalDestIndices = HaloMap.LocalToLocalDestIndices();

  core::profiler &Profiler = Context_->core_Profiler();

  int NumNeighbors = HaloMap.NeighborRanks().Count();

  MPI_Datatype DataType = core::GetMPIDataType<mpi_value_type>();

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
      SendBuffers_(iNeighbor)(iSendPoint) = mpi_value_type(FieldData[iPoint]);
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
    FieldData[iDest] = FieldData[iSource];
  }

  Profiler.Stop(PACK_TIME);
  Profiler.Stop(UNPACK_TIME);

  Active_ = true;

  return exchange_request(*this, FieldData);

}

template <typename T> halo_exchanger_for_type<T>::exchange_request::exchange_request(
  halo_exchanger_for_type &HaloExchanger, value_type *FieldData):
  HaloExchanger_(HaloExchanger.FloatingRefGenerator_.Generate(HaloExchanger)),
  FieldData_(FieldData)
{}

template <typename T> void halo_exchanger_for_type<T>::exchange_request::OnMPIRequestComplete(int
  iMPIRequest) {

  halo_exchanger_for_type &HaloExchanger = *HaloExchanger_;
  const halo_map &HaloMap = *HaloExchanger.HaloMap_;

  core::profiler &Profiler = HaloExchanger.Context_->core_Profiler();

  if (iMPIRequest < HaloExchanger.RecvBuffers_.Count()) {

    Profiler.Start(UNPACK_TIME);

    int iNeighbor = iMPIRequest;
    const array<long long> &RecvIndices = HaloMap.NeighborRecvIndices(iNeighbor);
    for (long long iRecvPoint = 0; iRecvPoint < RecvIndices.Count(); ++iRecvPoint) {
      long long iPoint = RecvIndices(iRecvPoint);
      FieldData_[iPoint] = value_type(HaloExchanger.RecvBuffers_(iNeighbor)(iRecvPoint));
    }

    Profiler.Stop(UNPACK_TIME);

  }

}

template <typename T> void halo_exchanger_for_type<T>::exchange_request::OnComplete() {

  HaloExchanger_->Active_ = false;

}

}

}}
