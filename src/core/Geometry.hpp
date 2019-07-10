// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GEOMETRY_HPP_INCLUDED
#define OVK_CORE_GEOMETRY_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Editor.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Geometry.h>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

enum class geometry_type {
  UNIFORM = OVK_GEOMETRY_TYPE_UNIFORM,
  ORIENTED_UNIFORM = OVK_GEOMETRY_TYPE_ORIENTED_UNIFORM,
  RECTILINEAR = OVK_GEOMETRY_TYPE_RECTILINEAR,
  ORIENTED_RECTILINEAR = OVK_GEOMETRY_TYPE_ORIENTED_RECTILINEAR,
  CURVILINEAR = OVK_GEOMETRY_TYPE_CURVILINEAR
};

inline bool ValidGeometryType(geometry_type GeometryType) {
  return ovkValidGeometryType(ovk_geometry_type(GeometryType));
}

namespace geometry_internal {

// For doing stuff before creation and after destruction
class geometry_base {

protected:

  geometry_base(std::shared_ptr<context> &&Context, const grid &Grid);

  geometry_base(const geometry_base &Other) = delete;
  geometry_base(geometry_base &&Other) noexcept = default;

  geometry_base &operator=(const geometry_base &Other) = delete;
  geometry_base &operator=(geometry_base &&Other) noexcept = default;

  ~geometry_base() noexcept;

  std::shared_ptr<context> Context_;

  const grid *Grid_;

  comm_view Comm_;

};

}

class geometry : private geometry_internal::geometry_base {

public:

  class params {
  public:
    params() = default;
    geometry_type GeometryType() const { return GeometryType_; }
    params &SetGeometryType(geometry_type GeometryType);
    const tuple<double> &PeriodicLength() const { return PeriodicLength_; }
    params &SetPeriodicLength(const tuple<double> &PeriodicLength);
  private:
    geometry_type GeometryType_ = geometry_type::CURVILINEAR;
    tuple<double> PeriodicLength_ = {0., 0., 0.};
    friend class geometry;
  };

  geometry(const geometry &Other) = delete;
  geometry(geometry &&Other) noexcept = default;

  geometry &operator=(const geometry &Other) = delete;
  geometry &operator=(geometry &&Other) = default;

  ~geometry() noexcept;

  floating_ref<const geometry> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<geometry> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

  const grid &Grid() const { return *Grid_; }

  int Dimension() const { return NumDims_; }

  comm_view Comm() const { return Comm_; }

  geometry_type GeometryType() const { return GeometryType_; }

  const tuple<double> &PeriodicLength() const { return PeriodicLength_; }
  double PeriodicLength(int iDim) const { return PeriodicLength_[iDim]; }
  void SetPeriodicLength(const tuple<double> &PeriodicLength);

  const array<field<double>> &Coords() const { return Coords_; }
  bool EditingCoords() const;
  edit_handle<array<field<double>>> EditCoords();
  void RestoreCoords();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddCoordsEventListener(F Listener) const {
    return CoordsEvent_.AddListener(std::move(Listener));
  }

  static geometry internal_Create(std::shared_ptr<context> &&Context, const grid &Grid, params
    &&Params);

private:

  floating_ref_generator FloatingRefGenerator_;

  int NumDims_;

  geometry_type GeometryType_;

  tuple<double> PeriodicLength_;

  array<field<double>> Coords_;
  editor CoordsEditor_;
  mutable event<void()> CoordsEvent_;

  geometry(std::shared_ptr<context> &&Context, const grid &Grid, params &&Params);

};

namespace core {
geometry CreateGeometry(std::shared_ptr<context> Context, const grid &Grid, geometry::params
  Params={});
}

}

#endif
