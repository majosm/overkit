// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GEOMETRY_COMPONENT_HPP_INCLUDED
#define OVK_CORE_GEOMETRY_COMPONENT_HPP_INCLUDED

#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/DomainBase.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Geometry.hpp>
#include <ovk/core/GeometryComponent.h>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/StringWrapper.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>

namespace ovk {

enum class geometry_event_flags : typename std::underlying_type<ovk_geometry_event_flags>::type {
  NONE = OVK_GEOMETRY_EVENT_FLAGS_NONE,
  CREATE = OVK_GEOMETRY_EVENT_FLAGS_CREATE,
  DESTROY = OVK_GEOMETRY_EVENT_FLAGS_DESTROY,
  EDIT_COORDS = OVK_GEOMETRY_EVENT_FLAGS_EDIT_COORDS,
  ALL = OVK_GEOMETRY_EVENT_FLAGS_ALL
};

constexpr inline geometry_event_flags operator|(geometry_event_flags Left, geometry_event_flags
  Right) {
  return geometry_event_flags(int(Left) | int(Right));
}
constexpr inline geometry_event_flags operator&(geometry_event_flags Left, geometry_event_flags
  Right) {
  return geometry_event_flags(int(Left) & int(Right));
}
constexpr inline geometry_event_flags operator^(geometry_event_flags Left, geometry_event_flags
  Right) {
  return geometry_event_flags(int(Left) ^ int(Right));
}
constexpr inline geometry_event_flags operator~(geometry_event_flags EventFlags) {
  return geometry_event_flags(~int(EventFlags));
}
inline geometry_event_flags operator|=(geometry_event_flags &Left, geometry_event_flags Right) {
  return Left = Left | Right;
}
inline geometry_event_flags operator&=(geometry_event_flags &Left, geometry_event_flags Right) {
  return Left = Left & Right;
}
inline geometry_event_flags operator^=(geometry_event_flags &Left, geometry_event_flags Right) {
  return Left = Left ^ Right;
}

namespace geometry_component_internal {

// For doing stuff before creation and after destruction
class geometry_component_base {

public:

  geometry_component_base(const core::domain_base &Domain, std::string &&Name);

  geometry_component_base(const geometry_component_base &Other) = delete;
  geometry_component_base(geometry_component_base &&Other) noexcept = default;

  geometry_component_base &operator=(const geometry_component_base &Other) = delete;
  geometry_component_base &operator=(geometry_component_base &&Other) noexcept = default;

  ~geometry_component_base() noexcept;

  floating_ref<const context> Context_;
  floating_ref<const core::domain_base> Domain_;

  core::string_wrapper Name_;

};

}

class geometry_component : private geometry_component_internal::geometry_component_base {

public:

  class params {
  public:
    params() {}
    const std::string &Name() const { return Name_; }
    params &SetName(std::string Name);
  private:
    core::string_wrapper Name_ = "GeometryComponent";
    friend class geometry_component;
  };

  geometry_component(const core::domain_base &Domain, params Params={});

  geometry_component(const geometry_component &Other) = delete;
  geometry_component(geometry_component &&Other) noexcept = default;

  geometry_component &operator=(const geometry_component &Other) = delete;
  geometry_component &operator=(geometry_component &&Other) noexcept = default;

  floating_ref<const geometry_component> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<geometry_component> GetFloatingRef() {
    return FloatingRefGenerator_.Generate(*this);
  }

  const std::string &Name() const { return *Name_; }

  int GeometryCount() const;

  const set<int> &GeometryIDs() const;

  bool GeometryExists(int GridID) const;

  void CreateGeometry(int GridID, optional<geometry::params> MaybeParams={});
  void CreateGeometries(array_view<const int> GridIDs, array<optional<geometry::params>>
    MaybeParams={});

  void DestroyGeometry(int GridID);
  void DestroyGeometries(array_view<const int> GridIDs);

  void ClearGeometries();

  int LocalGeometryCount() const;

  const set<int> &LocalGeometryIDs() const;

  const geometry &Geometry(int GridID) const;
  bool EditingGeometry(int GridID) const;
  edit_handle<geometry> EditGeometry(int GridID);
  void RestoreGeometry(int GridID);

  void StartEdit();
  void EndEdit();

  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F, int, geometry_event_flags,
    bool>())> event_listener_handle AddGeometryEventListener(F Listener) const {
    return GeometryEvent_.AddListener(std::move(Listener));
  }

private:

  struct geometry_record {};

  struct local {
    geometry Geometry;
    geometry_event_flags EventFlags;
    floating_ref_generator FloatingRefGenerator;
    event_listener_handle CoordsEventListener;
    editor Editor;
    explicit local(geometry Geometry);
  };

  floating_ref_generator FloatingRefGenerator_;

  map<int,grid_event_flags> GridEventFlags_;
  event_listener_handle GridEventListener_;

  map<int,geometry_record> GeometryRecords_;
  map_noncontig<int,local> Locals_;

  mutable event<void(int, geometry_event_flags, bool)> GeometryEvent_;

  void OnGridEvent_();
  void DestroyGeometriesForDyingGrids_();

  void SyncEdits_();

};

}

#endif
