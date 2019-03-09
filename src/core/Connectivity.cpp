// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Connectivity.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/ConnectivityD.hpp"
#include "ovk/core/ConnectivityR.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/ErrorHandler.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <memory>
#include <string>

namespace ovk {

namespace {

bool EditingDonorSide(const connectivity &Connectivity);
bool EditingReceiverSide(const connectivity &Connectivity);

void CreateDonorInfo(connectivity::donor_info &DonorInfo, const grid *DonorGrid,
  int DestinationGridID, const core::comm &Comm);
void DestroyDonorInfo(connectivity::donor_info &DonorInfo);

void CreateReceiverInfo(connectivity::receiver_info &ReceiverInfo, const grid *ReceiverGrid,
  int SourceGridID, const core::comm &Comm);
void DestroyReceiverInfo(connectivity::receiver_info &ReceiverInfo);

void EditDonorSideGlobal(connectivity &Connectivity, connectivity_d *&Donors, bool IsLocal);
void ReleaseDonorSideGlobal(connectivity &Connectivity, connectivity_d *&Donors, bool IsLocal);
void EditReceiverSideGlobal(connectivity &Connectivity, connectivity_r *&Receivers, bool IsLocal);
void ReleaseReceiverSideGlobal(connectivity &Connectivity, connectivity_r *&Receivers, bool IsLocal);

void DefaultEdits(connectivity::edits &Edits);

}

namespace core {

void CreateConnectivity(connectivity &Connectivity, int NumDims, core::comm Comm, const grid
  *DonorGrid, const grid *ReceiverGrid, core::logger &Logger, core::error_handler &ErrorHandler) {

  Connectivity.Comm_ = std::move(Comm);

  MPI_Barrier(Connectivity.Comm_);

  Connectivity.Logger_ = &Logger;
  Connectivity.ErrorHandler_ = &ErrorHandler;

  // Don't know source/destination grid IDs globally yet, so use -1 as placeholder
  CreateDonorInfo(Connectivity.DonorInfo_, DonorGrid, -1, Connectivity.Comm_);
  CreateReceiverInfo(Connectivity.ReceiverInfo_, ReceiverGrid, -1, Connectivity.Comm_);

  GetGridInfoID(Connectivity.DonorInfo_, Connectivity.DonorGridID_);
  GetGridInfoID(Connectivity.ReceiverInfo_, Connectivity.ReceiverGridID_);

  // Now fill in the source/destination grid IDs
  Connectivity.DonorInfo_.DestinationGridID_ = Connectivity.ReceiverGridID_;
  Connectivity.ReceiverInfo_.SourceGridID_ = Connectivity.DonorGridID_;

  std::string DonorGridName, ReceiverGridName;
  GetGridInfoName(Connectivity.DonorInfo_, DonorGridName);
  GetGridInfoName(Connectivity.ReceiverInfo_, ReceiverGridName);

  Connectivity.Name_ = StringPrint("(%s,%s)", DonorGridName, ReceiverGridName);

  Connectivity.NumDims_ = NumDims;

  DefaultEdits(Connectivity.Edits_);

  if (DonorGrid) {
    Connectivity.Donors_.reset(new connectivity_d());
    core::CreateConnectivityDonorSide(*Connectivity.Donors_, *DonorGrid,
      Connectivity.ReceiverGridID_, Logger, ErrorHandler);
  }

  if (ReceiverGrid) {
    Connectivity.Receivers_.reset(new connectivity_r());
    core::CreateConnectivityReceiverSide(*Connectivity.Receivers_, *ReceiverGrid,
      Connectivity.DonorGridID_, Logger, ErrorHandler);
  }

  MPI_Barrier(Connectivity.Comm_);

  core::LogStatus(*Connectivity.Logger_, Connectivity.Comm_.Rank() == 0, 0, "Created connectivity "
    "%s.", Connectivity.Name_);

}

void DestroyConnectivity(connectivity &Connectivity) {

  MPI_Barrier(Connectivity.Comm_);

  if (Connectivity.Donors_) {
    core::DestroyConnectivityDonorSide(*Connectivity.Donors_);
    Connectivity.Donors_.reset();
  }

  if (Connectivity.Receivers_) {
    core::DestroyConnectivityReceiverSide(*Connectivity.Receivers_);
    Connectivity.Receivers_.reset();
  }

  DestroyDonorInfo(Connectivity.DonorInfo_);
  DestroyReceiverInfo(Connectivity.ReceiverInfo_);

  MPI_Barrier(Connectivity.Comm_);

  LogStatus(*Connectivity.Logger_, Connectivity.Comm_.Rank() == 0, 0, "Destroyed connectivity %s.",
    Connectivity.Name_);

  Connectivity.Comm_.Reset();

}

void CreateConnectivityInfo(connectivity_info &Info, const connectivity *Connectivity, const comm
  &Comm) {

  bool IsLocal = Connectivity != nullptr;
  bool IsRoot = false;
  if (IsLocal) {
    IsRoot = Connectivity->Comm_.Rank() == 0;
  }

  int RootRank;
  if (IsRoot) RootRank = Comm.Rank();
  BroadcastAnySource(&RootRank, 1, MPI_INT, IsRoot, Comm);

  if (IsRoot) {
    Info.DonorGridID_ = Connectivity->DonorGridID_;
    Info.ReceiverGridID_ = Connectivity->ReceiverGridID_;
  }
  MPI_Bcast(&Info.DonorGridID_, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info.ReceiverGridID_, 1, MPI_INT, RootRank, Comm);

  int NameLength;
  if (IsRoot) NameLength = Connectivity->Name_.length();
  MPI_Bcast(&NameLength, 1, MPI_INT, RootRank, Comm);
  array<char> NameChars({NameLength});
  if (IsRoot) NameChars.Fill(Connectivity->Name_.begin());
  MPI_Bcast(NameChars.Data(), NameLength, MPI_CHAR, RootRank, Comm);
  Info.Name_.assign(NameChars.LinearBegin(), NameChars.LinearEnd());

  if (IsRoot) Info.NumDims_ = Connectivity->NumDims_;
  MPI_Bcast(&Info.NumDims_, 1, MPI_INT, RootRank, Comm);

  Info.RootRank_ = RootRank;

  Info.IsLocal_ = IsLocal;

}

void DestroyConnectivityInfo(connectivity_info &Info) {

  Info.Name_.clear();

}

}

namespace {

void CreateDonorInfo(connectivity::donor_info &DonorInfo, const grid *DonorGrid,
  int DestinationGridID, const core::comm &Comm) {

  core::CreateGridInfo(DonorInfo, DonorGrid, Comm);

  DonorInfo.DestinationGridID_ = DestinationGridID;
  DonorInfo.EditRefCount_ = 0;

}

void DestroyDonorInfo(connectivity::donor_info &DonorInfo) {

  core::DestroyGridInfo(DonorInfo);

}

void CreateReceiverInfo(connectivity::receiver_info &ReceiverInfo, const grid *ReceiverGrid,
  int SourceGridID, const core::comm &Comm) {

  core::CreateGridInfo(ReceiverInfo, ReceiverGrid, Comm);

  ReceiverInfo.SourceGridID_ = SourceGridID;
  ReceiverInfo.EditRefCount_ = 0;

}

void DestroyReceiverInfo(connectivity::receiver_info &ReceiverInfo) {

  core::DestroyGridInfo(ReceiverInfo);

}

}

void GetConnectivityDonorGridID(const connectivity &Connectivity, int &DonorGridID) {

  DonorGridID = Connectivity.DonorGridID_;

}

void GetConnectivityReceiverGridID(const connectivity &Connectivity, int &ReceiverGridID) {

  ReceiverGridID = Connectivity.ReceiverGridID_;

}

void GetConnectivityName(const connectivity &Connectivity, std::string &Name) {

  Name = Connectivity.Name_;

}

void GetConnectivityDimension(const connectivity &Connectivity, int &NumDims) {

  NumDims = Connectivity.NumDims_;

}

void GetConnectivityComm(const connectivity &Connectivity, MPI_Comm &Comm) {

  Comm = Connectivity.Comm_.Get();

}

void GetConnectivityCommSize(const connectivity &Connectivity, int &CommSize) {

  CommSize = Connectivity.Comm_.Size();

}

void GetConnectivityCommRank(const connectivity &Connectivity, int &CommRank) {

  CommRank = Connectivity.Comm_.Rank();

}

namespace core {

const comm &GetConnectivityComm(const connectivity &Connectivity) {

  return Connectivity.Comm_;

}

}

void GetConnectivityDonorGridInfo(const connectivity &Connectivity, const grid_info *&DonorGridInfo)
  {

  DonorGridInfo = &Connectivity.DonorInfo_;

}

void GetConnectivityReceiverGridInfo(const connectivity &Connectivity, const grid_info
  *&ReceiverGridInfo) {

  ReceiverGridInfo = &Connectivity.ReceiverInfo_;

}

bool RankHasConnectivityDonorSide(const connectivity &Connectivity) {

  return Connectivity.Donors_ != nullptr;

}

void GetConnectivityDonorSide(const connectivity &Connectivity, const connectivity_d *&Donors) {

  OVK_DEBUG_ASSERT(Connectivity.Donors_, "Connectivity %s does not have donor-side data on "
    "rank @rank@.", Connectivity.Name_);

  Donors = Connectivity.Donors_.get();

}

void EditConnectivityDonorSideLocal(connectivity &Connectivity, connectivity_d *&Donors) {

  EditDonorSideGlobal(Connectivity, Donors, true);

}

void EditConnectivityDonorSideRemote(connectivity &Connectivity) {

  connectivity_d *IgnoredDonors = nullptr;
  EditDonorSideGlobal(Connectivity, IgnoredDonors, false);
  
}

void ReleaseConnectivityDonorSideLocal(connectivity &Connectivity, connectivity_d *&Donors) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");

  ReleaseDonorSideGlobal(Connectivity, Donors, true);

}

void ReleaseConnectivityDonorSideRemote(connectivity &Connectivity) {

  connectivity_d *IgnoredDonors = nullptr;
  ReleaseDonorSideGlobal(Connectivity, IgnoredDonors, false);
  
}

namespace {

void EditDonorSideGlobal(connectivity &Connectivity, connectivity_d *&Donors, bool IsLocal) {

  if (OVK_DEBUG && IsLocal) {
    OVK_DEBUG_ASSERT(Connectivity.Donors_, "Connectivity %s does not have donor-side data on "
      "rank @rank@.", Connectivity.Name_);
  }

  bool StartEdit = Connectivity.DonorInfo_.EditRefCount_ == 0;
  ++Connectivity.DonorInfo_.EditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Connectivity.Comm_);
  }

  if (IsLocal) {
    Donors = Connectivity.Donors_.get();
  }

}

void ReleaseDonorSideGlobal(connectivity &Connectivity, connectivity_d *&Donors, bool IsLocal) {

  OVK_DEBUG_ASSERT(EditingDonorSide(Connectivity), "Unable to release connectivity %s donor-side "
    "data; not currently being edited.", Connectivity.Name_);

  if (OVK_DEBUG && IsLocal) {
    OVK_DEBUG_ASSERT(Connectivity.Donors_, "Connectivity %s does not have donor-side data on "
      "rank @rank@.", Connectivity.Name_);
    OVK_DEBUG_ASSERT(Donors == Connectivity.Donors_.get(), "Invalid donors pointer.");
  }

  --Connectivity.DonorInfo_.EditRefCount_;
  bool EndEdit = Connectivity.DonorInfo_.EditRefCount_ == 0;

  if (IsLocal) {
    Donors = nullptr;
  }

  if (EndEdit) {

    MPI_Barrier(Connectivity.Comm_);

    const connectivity_d::edits *DonorEdits;
    if (IsLocal) {
      core::GetConnectivityDonorSideEdits(*Connectivity.Donors_, DonorEdits);
    }

    connectivity::edits &Edits = Connectivity.Edits_;

    int RootRank;
    GetGridInfoRootRank(Connectivity.DonorInfo_, RootRank);

    bool IsGridRoot = Connectivity.Comm_.Rank() == RootRank;

    int EditedNumDonors = 0;
    if (IsGridRoot) EditedNumDonors = DonorEdits->Count_;
    MPI_Bcast(&EditedNumDonors, 1, MPI_INT, RootRank, Connectivity.Comm_);
    Edits.NumDonors_ = Edits.NumDonors_ || EditedNumDonors;

    int EditedExtents = 0;
    if (IsGridRoot) EditedExtents = DonorEdits->Extents_;
    MPI_Bcast(&EditedExtents, 1, MPI_INT, RootRank, Connectivity.Comm_);
    Edits.DonorExtents_ = Edits.DonorExtents_ || EditedExtents;

    int EditedCoords = 0;
    if (IsGridRoot) EditedCoords = DonorEdits->Coords_;
    MPI_Bcast(&EditedCoords, 1, MPI_INT, RootRank, Connectivity.Comm_);
    Edits.DonorCoords_ = Edits.DonorCoords_ || EditedCoords;

    int EditedInterpCoefs = 0;
    if (IsGridRoot) EditedInterpCoefs = DonorEdits->InterpCoefs_;
    MPI_Bcast(&EditedInterpCoefs, 1, MPI_INT, RootRank, Connectivity.Comm_);
    Edits.DonorInterpCoefs_ = Edits.DonorInterpCoefs_ || EditedInterpCoefs;

    int EditedDestinations = 0;
    if (IsGridRoot) EditedDestinations = DonorEdits->Destinations_;
    MPI_Bcast(&EditedDestinations, 1, MPI_INT, RootRank, Connectivity.Comm_);
    Edits.DonorDestinations_ = Edits.DonorDestinations_ || EditedDestinations;

  }

}

}

bool RankHasConnectivityReceiverSide(const connectivity &Connectivity) {

  return Connectivity.Receivers_ != nullptr;

}

void GetConnectivityReceiverSide(const connectivity &Connectivity, const connectivity_r *&Receivers)
  {

  OVK_DEBUG_ASSERT(Connectivity.Receivers_, "Connectivity %s does not have receiver-side data on "
    "rank @rank@.", Connectivity.Name_);

  Receivers = Connectivity.Receivers_.get();

}

void EditConnectivityReceiverSideLocal(connectivity &Connectivity, connectivity_r *&Receivers) {

  EditReceiverSideGlobal(Connectivity, Receivers, true);

}

void EditConnectivityReceiverSideRemote(connectivity &Connectivity) {

  connectivity_r *IgnoredReceivers = nullptr;
  EditReceiverSideGlobal(Connectivity, IgnoredReceivers, false);
  
}

void ReleaseConnectivityReceiverSideLocal(connectivity &Connectivity, connectivity_r *&Receivers) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");

  ReleaseReceiverSideGlobal(Connectivity, Receivers, true);

}

void ReleaseConnectivityReceiverSideRemote(connectivity &Connectivity) {

  connectivity_r *IgnoredReceivers = nullptr;
  ReleaseReceiverSideGlobal(Connectivity, IgnoredReceivers, false);
  
}

namespace {

void EditReceiverSideGlobal(connectivity &Connectivity, connectivity_r *&Receivers, bool IsLocal) {

  if (OVK_DEBUG && IsLocal) {
    OVK_DEBUG_ASSERT(Connectivity.Receivers_, "Connectivity %s does not have receiver-side data "
      "on rank @rank@.", Connectivity.Name_);
  }

  bool StartEdit = Connectivity.ReceiverInfo_.EditRefCount_ == 0;
  ++Connectivity.ReceiverInfo_.EditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Connectivity.Comm_);
  }

  if (IsLocal) {
    Receivers = Connectivity.Receivers_.get();
  }

}

void ReleaseReceiverSideGlobal(connectivity &Connectivity, connectivity_r *&Receivers,
  bool IsLocal) {

  OVK_DEBUG_ASSERT(EditingReceiverSide(Connectivity), "Unable to release connectivity %s "
    "receiver-side data; not currently being edited.", Connectivity.Name_);

  if (IsLocal) {
    OVK_DEBUG_ASSERT(Connectivity.Receivers_, "Connectivity %s does not have receiver-side data "
      "on rank @rank@.", Connectivity.Name_);
    OVK_DEBUG_ASSERT(Receivers == Connectivity.Receivers_.get(), "Invalid receivers pointer.");
  }

  --Connectivity.ReceiverInfo_.EditRefCount_;
  bool EndEdit = Connectivity.ReceiverInfo_.EditRefCount_ == 0;

  if (IsLocal) {
    Receivers = nullptr;
  }

  if (EndEdit) {

    MPI_Barrier(Connectivity.Comm_);

    const connectivity_r::edits *ReceiverEdits;
    if (IsLocal) {
      core::GetConnectivityReceiverSideEdits(*Connectivity.Receivers_, ReceiverEdits);
    }

    connectivity::edits &Edits = Connectivity.Edits_;

    int RootRank;
    GetGridInfoRootRank(Connectivity.ReceiverInfo_, RootRank);

    bool IsGridRoot = Connectivity.Comm_.Rank() == RootRank;

    int EditedNumReceivers = 0;
    if (IsGridRoot) EditedNumReceivers = ReceiverEdits->Count_;
    MPI_Bcast(&EditedNumReceivers, 1, MPI_INT, RootRank, Connectivity.Comm_);
    Edits.NumReceivers_ = Edits.NumReceivers_ || EditedNumReceivers;

    int EditedPoints = 0;
    if (IsGridRoot) EditedPoints = ReceiverEdits->Points_;
    MPI_Bcast(&EditedPoints, 1, MPI_INT, RootRank, Connectivity.Comm_);
    Edits.ReceiverPoints_ = Edits.ReceiverPoints_ || EditedPoints;

    int EditedSources = 0;
    if (IsGridRoot) EditedSources = ReceiverEdits->Sources_;
    MPI_Bcast(&EditedSources, 1, MPI_INT, RootRank, Connectivity.Comm_);
    Edits.ReceiverSources_ = Edits.ReceiverSources_ || EditedSources;

  }

}

bool EditingDonorSide(const connectivity &Connectivity) {

  return Connectivity.DonorInfo_.EditRefCount_ > 0;

}

bool EditingReceiverSide(const connectivity &Connectivity) {

  return Connectivity.ReceiverInfo_.EditRefCount_ > 0;

}

}

namespace core {

void GetConnectivityEdits(const connectivity &Connectivity, const connectivity::edits *&Edits) {

  Edits = &Connectivity.Edits_;

}

void ResetConnectivityEdits(connectivity &Connectivity) {

  DefaultEdits(Connectivity.Edits_);

  if (Connectivity.Donors_) {
    ResetConnectivityDonorSideEdits(*Connectivity.Donors_);
  }

  if (Connectivity.Receivers_) {
    ResetConnectivityReceiverSideEdits(*Connectivity.Receivers_);
  }

}

}

namespace {

void DefaultEdits(connectivity::edits &Edits) {

  Edits.NumDonors_ = false;
  Edits.DonorExtents_ = false;
  Edits.DonorCoords_ = false;
  Edits.DonorInterpCoefs_ = false;
  Edits.DonorDestinations_ = false;
  Edits.NumReceivers_ = false;
  Edits.ReceiverPoints_ = false;
  Edits.ReceiverSources_ = false;

}

}

void GetConnectivityInfoDonorGridID(const connectivity_info &Info, int &DonorGridID) {

  DonorGridID = Info.DonorGridID_;

}

void GetConnectivityInfoReceiverGridID(const connectivity_info &Info, int &ReceiverGridID) {

  ReceiverGridID = Info.ReceiverGridID_;

}

void GetConnectivityInfoName(const connectivity_info &Info, std::string &Name) {

  Name = Info.Name_;

}

void GetConnectivityInfoDimension(const connectivity_info &Info, int &NumDims) {

  NumDims = Info.NumDims_;

}

void GetConnectivityInfoRootRank(const connectivity_info &Info, int &RootRank) {

  RootRank = Info.RootRank_;

}

void GetConnectivityInfoIsLocal(const connectivity_info &Info, bool &IsLocal) {

  IsLocal = Info.IsLocal_;

}

}
