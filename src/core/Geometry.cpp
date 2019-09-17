// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Geometry.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Event.hpp"
#include "ovk/core/Field.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Partition.hpp"
#include "ovk/core/Request.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

namespace geometry_internal {

geometry_base::geometry_base(std::shared_ptr<context> &&Context, const grid &Grid):
  Context_(std::move(Context)),
  Grid_(&Grid),
  Comm_(Grid.Comm())
{
  MPI_Barrier(Comm_);
}

geometry_base::~geometry_base() noexcept {

  if (Context_) {
    MPI_Barrier(Comm_);
    core::logger &Logger = Context_->core_Logger();
    Logger.LogDebug(Comm_.Rank() == 0, 0, "Destroyed geometry %s.", Grid_->Name());
  }

}

}

geometry::geometry(std::shared_ptr<context> &&Context, const grid &Grid, params &&Params):
  geometry_base(std::move(Context), Grid),
  NumDims_(Grid.Dimension()),
  Type_(Params.Type_),
  PeriodicLength_(Params.PeriodicLength_),
  Coords_({MAX_DIMS})
{

  if (OVK_DEBUG) {
    for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
      OVK_DEBUG_ASSERT(PeriodicLength_(iDim) == 0., "Periodic length has incorrect dimension.");
    }
  }

  // In 1D, geometry type is always uniform or rectilinear
  if (NumDims_ == 1) {
    switch (Type_) {
    case geometry_type::UNIFORM:
    case geometry_type::ORIENTED_UNIFORM:
      Type_ = geometry_type::UNIFORM;
      break;
    case geometry_type::RECTILINEAR:
    case geometry_type::ORIENTED_RECTILINEAR:
    case geometry_type::CURVILINEAR:
      Type_ = geometry_type::RECTILINEAR;
      break;
    default:
      OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
      break;
    }
  }

  const range &ExtendedRange = Grid.ExtendedRange();

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Coords_(iDim).Assign(Grid.SharedPartition());
  }

  for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
    for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
      for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
        Coords_(0)(i,j,k) = double(i);
        Coords_(1)(i,j,k) = double(j);
        Coords_(2)(i,j,k) = double(k);
      }
    }
  }

  Volumes_.Assign(Grid.SharedPartition(), 1.);
  CellVolumes_.Assign(Grid.SharedCellPartition(), 1.);

  MPI_Barrier(Comm_);

  core::logger &Logger = Context_->core_Logger();
  Logger.LogDebug(Comm_.Rank() == 0, 0, "Created geometry %s.", Grid.Name());

}

geometry::~geometry() noexcept {

  if (Context_) {
    // Barrier before cleaning up
    MPI_Barrier(Comm_);
  }

}

geometry geometry::internal_Create(std::shared_ptr<context> &&Context, const grid &Grid, params
  &&Params) {

  return {std::move(Context), Grid, std::move(Params)};

}

namespace core {

geometry CreateGeometry(std::shared_ptr<context> Context, const grid &Grid, geometry::params Params)
  {

  OVK_DEBUG_ASSERT(Context, "Invalid context.");

  return geometry::internal_Create(std::move(Context), Grid, std::move(Params));

}

}

void geometry::SetPeriodicLength(const tuple<double> &PeriodicLength) {

  if (OVK_DEBUG) {
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      OVK_DEBUG_ASSERT(PeriodicLength(iDim) >= 0., "Periodic length must be nonnegative.");
    }
    for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
      OVK_DEBUG_ASSERT(PeriodicLength(iDim) == 0., "Periodic length has incorrect dimension.");
    }
  }

  PeriodicLength_ = PeriodicLength;

  CoordsEvent_.Trigger();

}

bool geometry::EditingCoords() const {

  return CoordsEditor_.Active();

}

edit_handle<array<distributed_field<double>>> geometry::EditCoords() {

  if (!CoordsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<geometry> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      geometry &Geometry = *FloatingRef;
      Geometry.OnCoordsEndEdit_();
      MPI_Barrier(Geometry.Comm_);
    };
    CoordsEditor_.Activate(std::move(DeactivateFunc));
  }

  return CoordsEditor_.Edit(Coords_);

}

void geometry::RestoreCoords() {

  OVK_DEBUG_ASSERT(CoordsEditor_.Active(), "Unable to restore coords; not currently being edited.");

  CoordsEditor_.Restore();

}

void geometry::OnCoordsEndEdit_() {

  core::logger &Logger = Context_->core_Logger();

  const grid &Grid = *Grid_;
  int NumDims = Grid.Dimension();
  const partition &Partition = Grid.Partition();
  const cart &Cart = Grid.Cart();
  const range &GlobalRange = Grid.GlobalRange();
  const range &LocalRange = Grid.LocalRange();
  const range &ExtendedRange = Grid.ExtendedRange();
  const range &CellLocalRange = Grid.CellLocalRange();
  const range &CellExtendedRange = Grid.CellExtendedRange();

  request Requests[MAX_DIMS];

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Requests[iDim] = Partition.Exchange(Coords_(iDim));
  }

  WaitAll(Requests);

  for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
    for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
      for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        if (!GlobalRange.Contains(Point)) {
          tuple<int> Period = Cart.GetPeriod(Point);
          for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
            Coords_(iDim)(Point) += double(Period(iDim)) * PeriodicLength_(iDim);
          }
        }
      }
    }
  }

  MPI_Barrier(Comm_);

  CoordsEvent_.Trigger();

  Logger.LogDebug(Comm_.Rank() == 0, 0, "Updating auxiliary data for geometry %s...", Grid.Name());

  for (int k = CellLocalRange.Begin(2); k < CellLocalRange.End(2); ++k) {
    for (int j = CellLocalRange.Begin(1); j < CellLocalRange.End(1); ++j) {
      for (int i = CellLocalRange.Begin(0); i < CellLocalRange.End(0); ++i) {
        tuple<int> Cell = {i,j,k};
        CellVolumes_(Cell) = core::CellVolume(NumDims, Coords_, Type_, Cell);
      }
    }
  }

  CellVolumes_.Exchange();

  // Compute the volume at each point by averaging the volumes of neighboring cells
  for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
    for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        range NeighborCellRange;
        for (int iDim = 0; iDim < NumDims; ++iDim) {
          NeighborCellRange.Begin(iDim) = Point(iDim)-1;
          NeighborCellRange.End(iDim) = Point(iDim)+1;
        }
        for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
          NeighborCellRange.Begin(iDim) = 0;
          NeighborCellRange.End(iDim) = 1;
        }
        NeighborCellRange = IntersectRanges(CellExtendedRange, NeighborCellRange);
        Volumes_(Point) = 0.;
        for (int o = NeighborCellRange.Begin(2); o < NeighborCellRange.End(2); ++o) {
          for (int n = NeighborCellRange.Begin(1); n < NeighborCellRange.End(1); ++n) {
            for (int m = NeighborCellRange.Begin(0); m < NeighborCellRange.End(0); ++m) {
              tuple<int> NeighborCell = {m,n,o};
              Volumes_(Point) += CellVolumes_(NeighborCell);
            }
          }
        }
        Volumes_(Point) /= double(NeighborCellRange.Count());
      }
    }
  }

  Volumes_.Exchange();

  MPI_Barrier(Comm_);

  Logger.LogDebug(Comm_.Rank() == 0, 0, "Done updating auxiliary data for geometry %s.",
    Grid.Name());

  VolumesEvent_.Trigger();
  CellVolumesEvent_.Trigger();

}

geometry::params &geometry::params::SetType(geometry_type Type) {

  OVK_DEBUG_ASSERT(ValidGeometryType(Type), "Invalid type.");

  Type_ = Type;

  return *this;

}

geometry::params &geometry::params::SetPeriodicLength(const tuple<double> &PeriodicLength) {

  if (OVK_DEBUG) {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      OVK_DEBUG_ASSERT(PeriodicLength(iDim) >= 0., "Periodic length must be nonnegative.");
    }
  }

  PeriodicLength_ = PeriodicLength;

  return *this;

}

}
