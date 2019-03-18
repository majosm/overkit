// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

// Intel 17 didn't like this for some reason
// template <typename ArrayType, OVK_FUNCDEF_REQUIRES(IsArray<ArrayType>() && ArrayHasFootprint<
//   ArrayType, MAX_DIMS, array_layout::GRID>())> request halo::Exchange(ArrayType &Array) const {
template <typename ArrayType, OVK_FUNCDEF_REQUIRES(IsArray<ArrayType>() && ArrayRank<ArrayType>()
  == MAX_DIMS && ArrayLayout<ArrayType>() == array_layout::GRID)> request halo::Exchange(ArrayType &Array) const {

  StartProfile(*Profiler_, TotalTime_);

  array<exchanger> &ExchangersForType = Exchangers_[GetTypeID<array_value_type<ArrayType>>()];

  int iExchanger = 0;
  while (iExchanger < ExchangersForType.Count() && ExchangersForType(iExchanger).Active()) {
    ++iExchanger;
  }
  if (iExchanger == ExchangersForType.Count()) {
    StartProfile(*Profiler_, SetupTime_);
    ExchangersForType.Append(halo_internal::exchanger_for_type<array_value_type<ArrayType>>(
      Comm_, NeighborRanks_, NeighborSendIndices_, NeighborRecvIndices_, LocalToLocalSourceIndices_,
      LocalToLocalDestIndices_, *Profiler_));
    EndProfile(*Profiler_, SetupTime_);
  }
  exchanger &Exchanger = ExchangersForType(iExchanger);

  StartProfile(*Profiler_, ExchangeTime_);

  auto EndProfiles = OnScopeExit([&] {
    EndProfile(*Profiler_, ExchangeTime_);
    EndProfile(*Profiler_, TotalTime_);
  });

  return Exchanger.Exchange(ArrayData(Array));

}

namespace halo_internal {

template <typename T> exchanger_for_type<T>::exchanger_for_type(comm_view Comm, const array<int>
  &NeighborRanks, const array<array<long long>> &NeighborSendIndices, const array<array<long long>>
  &NeighborRecvIndices, const array<long long> &LocalToLocalSourceIndices, const array<long long>
  &LocalToLocalDestIndices, profiler &Profiler):
  Comm_(Comm),
  NeighborRanks_(NeighborRanks),
  NeighborSendIndices_(NeighborSendIndices),
  NeighborRecvIndices_(NeighborRecvIndices),
  LocalToLocalSourceIndices_(LocalToLocalSourceIndices),
  LocalToLocalDestIndices_(LocalToLocalDestIndices),
  Active_(new bool(false)),
  Profiler_(&Profiler),
  PackTime_(GetProfilerTimerID(Profiler, "Halo::Exchange::Pack")),
  UnpackTime_(GetProfilerTimerID(Profiler, "Halo::Exchange::Unpack")),
  MPITime_(GetProfilerTimerID(Profiler, "Halo::Exchange::MPI"))
{

  int NumNeighbors = NeighborRanks.Count();

  SendBuffers_.Resize({NumNeighbors});
  RecvBuffers_.Resize({NumNeighbors});
  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    SendBuffers_(iNeighbor).Resize({NeighborSendIndices_(iNeighbor).Count()});
    RecvBuffers_(iNeighbor).Resize({NeighborRecvIndices_(iNeighbor).Count()});
  }

  MPIRequests_.Reserve(2*NumNeighbors);

}

template <typename T> request exchanger_for_type<T>::Exchange(value_type *ArrayData) {

  int NumNeighbors = NeighborRanks_.Count();

  MPI_Datatype DataType = GetMPIDataType<mpi_value_type>();

  MPIRequests_.Clear();

  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    MPI_Request &Request = MPIRequests_.Append();
    StartProfile(*Profiler_, MPITime_);
    MPI_Irecv(RecvBuffers_(iNeighbor).Data(), RecvBuffers_(iNeighbor).Count(), DataType,
      NeighborRanks_(iNeighbor), 0, Comm_, &Request);
    EndProfile(*Profiler_, MPITime_);
  }

  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    StartProfile(*Profiler_, PackTime_);
    const array<long long> &SendIndices = NeighborSendIndices_(iNeighbor);
    for (long long iSendPoint = 0; iSendPoint < SendIndices.Count(); ++iSendPoint) {
      long long iPoint = SendIndices(iSendPoint);
      SendBuffers_(iNeighbor)(iSendPoint) = mpi_value_type(ArrayData[iPoint]);
    }
    EndProfile(*Profiler_, PackTime_);
    MPI_Request &Request = MPIRequests_.Append();
    StartProfile(*Profiler_, MPITime_);
    MPI_Isend(SendBuffers_(iNeighbor).Data(), SendBuffers_(iNeighbor).Count(), DataType,
      NeighborRanks_(iNeighbor), 0, Comm_, &Request);
    EndProfile(*Profiler_, MPITime_);
  }

  StartProfile(*Profiler_, PackTime_);
  StartProfile(*Profiler_, UnpackTime_);

  long long NumLocalToLocal = LocalToLocalSourceIndices_.Count();
  for (long long iLocalToLocal = 0; iLocalToLocal < NumLocalToLocal; ++iLocalToLocal) {
    long long iSource = LocalToLocalSourceIndices_(iLocalToLocal);
    long long iDest = LocalToLocalDestIndices_(iLocalToLocal);
    ArrayData[iDest] = ArrayData[iSource];
  }

  EndProfile(*Profiler_, PackTime_);
  EndProfile(*Profiler_, UnpackTime_);

  *Active_ = true;

  return exchanger_request<value_type>(ArrayData, NeighborRecvIndices_, RecvBuffers_, MPIRequests_,
    *Active_, *Profiler_);

}

template <typename T> exchanger_request<T>::exchanger_request(value_type *ArrayData, array_view<
  const array<long long>> NeighborRecvIndices, array_view<array<mpi_value_type>> RecvBuffers,
  array_view<MPI_Request> MPIRequests, bool &Active, profiler &Profiler):
  ArrayData_(ArrayData),
  NeighborRecvIndices_(NeighborRecvIndices),
  RecvBuffers_(RecvBuffers),
  MPIRequests_(MPIRequests),
  Active_(&Active),
  Profiler_(&Profiler),
  UnpackTime_(GetProfilerTimerID(Profiler, "Halo::Exchange::Unpack")),
  MPITime_(GetProfilerTimerID(Profiler, "Halo::Exchange::MPI"))
{}

template <typename T> void exchanger_request<T>::Finish(int iMPIRequest) {

  if (iMPIRequest < RecvBuffers_.Count()) {

    StartProfile(*Profiler_, UnpackTime_);

    int iNeighbor = iMPIRequest;
    const array<long long> &RecvIndices = NeighborRecvIndices_(iNeighbor);
    for (long long iRecvPoint = 0; iRecvPoint < RecvIndices.Count(); ++iRecvPoint) {
      long long iPoint = RecvIndices(iRecvPoint);
      ArrayData_[iPoint] = value_type(RecvBuffers_(iNeighbor)(iRecvPoint));
    }

    EndProfile(*Profiler_, UnpackTime_);

  }

}

template <typename T> void exchanger_request<T>::Wait() {

  while (true) {
    int iMPIRequest;
    StartProfile(*Profiler_, MPITime_);
    MPI_Waitany(MPIRequests_.Count(), MPIRequests_.Data(), &iMPIRequest, MPI_STATUSES_IGNORE);
    EndProfile(*Profiler_, MPITime_);
    if (iMPIRequest == MPI_UNDEFINED) {
      break;
    }
    Finish(iMPIRequest);
  }

  *Active_ = false;

}

}

}}
