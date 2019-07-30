// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Geometry.hpp"

#include "ovk/core/Array.hpp"
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
#include "ovk/core/Request.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

namespace geometry_internal {

geometry_base::geometry_base(std::shared_ptr<context> &&Context, const grid &Grid):
  Context_(std::move(Context)),
  Grid_(&Grid),
  Comm_(Grid_->Comm())
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
  NumDims_(Grid_->Dimension()),
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
    }
  }

  const range &ExtendedRange = Grid_->ExtendedRange();

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Coords_(iDim).Resize(ExtendedRange);
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

  MPI_Barrier(Comm_);

  core::logger &Logger = Context_->core_Logger();
  Logger.LogDebug(Comm_.Rank() == 0, 0, "Created geometry %s.", Grid_->Name());

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

edit_handle<array<field<double>>> geometry::EditCoords() {

  if (!CoordsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<geometry> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      geometry &Geometry = *FloatingRef;
      Geometry.OnCoordsEndEdit_();
      MPI_Barrier(Geometry.Comm_);
      Geometry.CoordsEvent_.Trigger();
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

  const grid &Grid = *Grid_;
  const range &GlobalRange = Grid.GlobalRange();
  const range &ExtendedRange = Grid.ExtendedRange();
  const core::halo &Halo = Grid.core_Partition().Halo();

  request Requests[MAX_DIMS];

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Requests[iDim] = Halo.Exchange(Coords_(iDim));
  }

  WaitAll(Requests);

  for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
    for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
      for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          if (Point(iDim) > GlobalRange.End(iDim)) {
            Coords_(iDim)(Point) += PeriodicLength_(iDim);
          } else if (Point(iDim) < GlobalRange.Begin(iDim)) {
            Coords_(iDim)(Point) -= PeriodicLength_(iDim);
          }
        }
      }
    }
  }

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
