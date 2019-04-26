// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Exchange.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Request.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

namespace ovk {

namespace core {

void CreateExchange(exchange &Exchange, const connectivity &Connectivity, logger &Logger,
  error_handler &ErrorHandler, profiler &Profiler) {

  Exchange.Connectivity_ = &Connectivity;

  Exchange.Logger_ = &Logger;
  Exchange.ErrorHandler_ = &ErrorHandler;

  Exchange.Profiler_ = &Profiler;

  GetConnectivityDimension(Connectivity, Exchange.NumDims_);

  Exchange.Comm_ = core::GetConnectivityComm(Connectivity);

  MPI_Barrier(Exchange.Comm_);

  Exchange.Logger_->LogStatus(Exchange.Comm_.Rank() == 0, 0, "Created exchange %s.",
    Connectivity.Name_);

}

void DestroyExchange(exchange &Exchange) {

  MPI_Barrier(Exchange.Comm_);

  const connectivity &Connectivity = *Exchange.Connectivity_;

  Exchange.Logger_->LogStatus(Exchange.Comm_.Rank() == 0, 0, "Destroyed exchange %s.",
    Connectivity.Name_);

  Exchange.Comm_.Reset();

}

void CreateExchangeInfo(exchange_info &Info, const exchange *Exchange, comm_view Comm) {

  bool IsLocal = Exchange != nullptr;
  bool IsRoot = false;
  if (IsLocal) {
    IsRoot = Exchange->Comm_.Rank() == 0;
  }

  const connectivity *Connectivity = nullptr;
  if (IsLocal) {
    Connectivity = Exchange->Connectivity_;
  }

  int RootRank;
  if (IsRoot) RootRank = Comm.Rank();
  core::BroadcastAnySource(&RootRank, 1, MPI_INT, IsRoot, Comm);

  if (IsRoot) {
    Info.DonorGridID_ = Connectivity->DonorGridID_;
    Info.ReceiverGridID_ = Connectivity->ReceiverGridID_;
    Info.NumDims_ = Connectivity->NumDims_;
  }
  MPI_Bcast(&Info.DonorGridID_, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info.ReceiverGridID_, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info.NumDims_, 1, MPI_INT, RootRank, Comm);

  int NameLength;
  if (IsRoot) NameLength = Connectivity->Name_.length();
  MPI_Bcast(&NameLength, 1, MPI_INT, RootRank, Comm);
  array<char> NameChars({NameLength});
  if (IsRoot) NameChars.Fill(Connectivity->Name_.begin());
  MPI_Bcast(NameChars.Data(), NameLength, MPI_CHAR, RootRank, Comm);
  Info.Name_.assign(NameChars.LinearBegin(), NameChars.LinearEnd());

  Info.RootRank_ = RootRank;

  Info.IsLocal_ = IsLocal;

}

void DestroyExchangeInfo(exchange_info &Info) {

  Info.Name_.clear();

}

}

bool RankHasExchangeDonorSide(const exchange &Exchange) {

  return RankHasConnectivityDonorSide(*Exchange.Connectivity_);

}

bool RankHasExchangeReceiverSide(const exchange &Exchange) {

  return RankHasConnectivityReceiverSide(*Exchange.Connectivity_);

}

namespace core {

void UpdateExchange(exchange &Exchange) {

  MPI_Barrier(Exchange.Comm_);

  const connectivity &Connectivity = *Exchange.Connectivity_;

  Exchange.Logger_->LogStatus(Exchange.Comm_.Rank() == 0, 0, "Updating exchange %s...",
    Connectivity.Name_);

  const connectivity::edits *Edits;
  core::GetConnectivityEdits(Connectivity, Edits);

  MPI_Barrier(Exchange.Comm_);

  Exchange.Logger_->LogStatus(Exchange.Comm_.Rank() == 0, 0, "Done updating exchange %s.",
    Connectivity.Name_);

}

}

void GetExchangeInfoDonorGridID(const exchange_info &Info, int &DonorGridID) {

  DonorGridID = Info.DonorGridID_;

}

void GetExchangeInfoReceiverGridID(const exchange_info &Info, int &ReceiverGridID) {

  ReceiverGridID = Info.ReceiverGridID_;

}

void GetExchangeInfoName(const exchange_info &Info, std::string &Name) {

  Name = Info.Name_;

}

void GetExchangeInfoDimension(const exchange_info &Info, int &NumDims) {

  NumDims = Info.NumDims_;

}

void GetExchangeInfoRootRank(const exchange_info &Info, int &RootRank) {

  RootRank = Info.RootRank_;

}

void GetExchangeInfoIsLocal(const exchange_info &Info, bool &IsLocal) {

  IsLocal = Info.IsLocal_;

}

}
