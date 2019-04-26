// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONNECTIVITY_HPP_INCLUDED
#define OVK_CORE_CONNECTIVITY_HPP_INCLUDED

#include <ovk/core/Comm.hpp>
#include <ovk/core/ConnectivityD.hpp>
#include <ovk/core/ConnectivityR.hpp>
#include <ovk/core/ErrorHandler.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Logger.hpp>

#include <mpi.h>

#include <memory>
#include <string>

namespace ovk {

struct connectivity {

  struct edits {
    bool NumDonors_;
    bool DonorExtents_;
    bool DonorCoords_;
    bool DonorInterpCoefs_;
    bool DonorDestinations_;
    bool NumReceivers_;
    bool ReceiverPoints_;
    bool ReceiverSources_;
  };

  struct donor_info : grid_info {
    int DestinationGridID_;
    int EditRefCount_;
  };

  struct receiver_info : grid_info {
    int SourceGridID_;
    int EditRefCount_;
  };

  mutable core::logger *Logger_;
  mutable core::error_handler *ErrorHandler_;
  int DonorGridID_;
  int ReceiverGridID_;
  std::string Name_;
  int NumDims_;
  core::comm Comm_;
  edits Edits_;
  donor_info DonorInfo_;
  receiver_info ReceiverInfo_;
  std::unique_ptr<connectivity_d> Donors_;
  std::unique_ptr<connectivity_r> Receivers_;
  // optional<connectivity_d> Donors_;
  // optional<connectivity_r> Receivers_;

};

struct connectivity_info {
  int DonorGridID_;
  int ReceiverGridID_;
  std::string Name_;
  int NumDims_;
  int RootRank_;
  int IsLocal_;
};

namespace core {

void CreateConnectivity(connectivity &Connectivity, int NumDims, core::comm Comm, const grid
  *DonorGrid, const grid *ReceiverGrid, core::logger &Logger, core::error_handler &ErrorHandler);
void DestroyConnectivity(connectivity &Connectivity);

void CreateConnectivityInfo(connectivity_info &Info, const connectivity *Connectivity, comm_view
  Comm);
void DestroyConnectivityInfo(connectivity_info &Info);

}

void GetConnectivityDonorGridID(const connectivity &Connectivity, int &DonorGridID);
void GetConnectivityReceiverGridID(const connectivity &Connectivity, int &ReceiverGridID);
void GetConnectivityName(const connectivity &Connectivity, std::string &Name);
void GetConnectivityDimension(const connectivity &Connectivity, int &NumDims);
void GetConnectivityComm(const connectivity &Connectivity, MPI_Comm &Comm);
void GetConnectivityCommSize(const connectivity &Connectivity, int &CommSize);
void GetConnectivityCommRank(const connectivity &Connectivity, int &CommRank);

namespace core {
const comm &GetConnectivityComm(const connectivity &Connectivity);
}

void GetConnectivityDonorGridInfo(const connectivity &Connectivity, const grid_info
  *&DonorGridInfo);
void GetConnectivityReceiverGridInfo(const connectivity &Connectivity, const grid_info
  *&ReceiverGridInfo);

bool RankHasConnectivityDonorSide(const connectivity &Connectivity);
void GetConnectivityDonorSide(const connectivity &Connectivity, const connectivity_d *&Donors);
void EditConnectivityDonorSideLocal(connectivity &Connectivity, connectivity_d *&Donors);
void EditConnectivityDonorSideRemote(connectivity &Connectivity);
void ReleaseConnectivityDonorSideLocal(connectivity &Connectivity, connectivity_d *&Donors);
void ReleaseConnectivityDonorSideRemote(connectivity &Connectivity);

bool RankHasConnectivityReceiverSide(const connectivity &Connectivity);
void GetConnectivityReceiverSide(const connectivity &Connectivity, const connectivity_r
  *&Receivers);
void EditConnectivityReceiverSideLocal(connectivity &Connectivity, connectivity_r *&Receivers);
void EditConnectivityReceiverSideRemote(connectivity &Connectivity);
void ReleaseConnectivityReceiverSideLocal(connectivity &Connectivity, connectivity_r *&Receivers);
void ReleaseConnectivityReceiverSideRemote(connectivity &Connectivity);

namespace core {
void GetConnectivityEdits(const connectivity &Connectivity, const connectivity::edits *&Edits);
void ResetConnectivityEdits(connectivity &Connectivity);
}

void GetConnectivityInfoDonorGridID(const connectivity_info &Info, int &DonorGridID);
void GetConnectivityInfoReceiverGridID(const connectivity_info &Info, int &ReceiverGridID);
void GetConnectivityInfoName(const connectivity_info &Info, std::string &Name);
void GetConnectivityInfoDimension(const connectivity_info &Info, int &NumDims);
void GetConnectivityInfoRootRank(const connectivity_info &Info, int &RootRank);
void GetConnectivityInfoIsLocal(const connectivity_info &Info, bool &IsLocal);

}

#endif
