// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DOMAIN_BASE_HPP_INCLUDED
#define OVK_CORE_DOMAIN_BASE_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Domain.h>
#include <ovk/core/Event.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Logger.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/StringWrapper.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

namespace ovk {

enum class grid_event_flags : typename std::underlying_type<ovk_grid_event_flags>::type {
  NONE = OVK_GRID_EVENT_FLAGS_NONE,
  CREATE = OVK_GRID_EVENT_FLAGS_CREATE,
  DESTROY = OVK_GRID_EVENT_FLAGS_DESTROY,
  ALL = OVK_GRID_EVENT_FLAGS_ALL
};

constexpr inline grid_event_flags operator|(grid_event_flags Left, grid_event_flags Right) {
  return grid_event_flags(int(Left) | int(Right));
}
constexpr inline grid_event_flags operator&(grid_event_flags Left, grid_event_flags Right) {
  return grid_event_flags(int(Left) & int(Right));
}
constexpr inline grid_event_flags operator^(grid_event_flags Left, grid_event_flags Right) {
  return grid_event_flags(int(Left) ^ int(Right));
}
constexpr inline grid_event_flags operator~(grid_event_flags EventFlags) {
  return grid_event_flags(~int(EventFlags));
}
inline grid_event_flags operator|=(grid_event_flags &Left, grid_event_flags Right) {
  return Left = Left | Right;
}
inline grid_event_flags operator&=(grid_event_flags &Left, grid_event_flags Right) {
  return Left = Left & Right;
}
inline grid_event_flags operator^=(grid_event_flags &Left, grid_event_flags Right) {
  return Left = Left ^ Right;
}

namespace core {

namespace domain_internal {

// For doing stuff before creation and after destruction
class domain_base_1 {

protected:

  domain_base_1(std::shared_ptr<context> &&Context, std::string &&Name, MPI_Comm Comm);

  domain_base_1(const domain_base_1 &Other) = delete;
  domain_base_1(domain_base_1 &&Other) noexcept = default;

  domain_base_1 &operator=(const domain_base_1 &Other) = delete;
  domain_base_1 &operator=(domain_base_1 &&Other) noexcept = default;

  ~domain_base_1() noexcept;

  std::shared_ptr<context> Context_;

  core::logger::status_level_and_indent_handle Level1_;
  core::logger::status_level_and_indent_handle Suppress_;

  core::string_wrapper Name_;

  comm Comm_;

};

}

class domain_base : protected domain_internal::domain_base_1 {

public:

  domain_base(std::shared_ptr<context> &&Context, std::string &&Name, int NumDims, MPI_Comm Comm);

  domain_base(const domain_base &Other) = delete;
  domain_base(domain_base &&Other) noexcept = default;

  domain_base &operator=(const domain_base &Other) = delete;
  domain_base &operator=(domain_base &&Other) noexcept = default;

  floating_ref<const domain_base> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<domain_base> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

  const std::string &Name() const { return *Name_; }
  int Dimension() const { return NumDims_; }

  const comm &Comm() const { return Comm_; }

  int GridCount() const { return GridRecords_.Count(); }

  const set<int> &GridIDs() const { return GridRecords_.Keys(); }

  bool GridExists(int GridID) const;

  grid::params MakeGridParams() const;

  void CreateGrid(int GridID, optional<grid::params> MaybeParams);
  void CreateGrids(array_view<const int> GridIDs, array<optional<grid::params>> MaybeParams);

  void DestroyGrid(int GridID);
  void DestroyGrids(array_view<const int> GridIDs);

  const grid_info &GridInfo(int GridID) const;
  int LocalGridCount() const { return LocalGrids_.Count(); }
  const set<int> &LocalGridIDs() const { return LocalGrids_.Keys(); }
  bool GridIsLocal(int GridID) const;
  const grid &Grid(int GridID) const;

  int GridIDToIndex(int GridID) const;
  int GridIDToLocalIndex(int GridID) const;
  int GridIndexToID(int GridIndex) const;
  int GridLocalIndexToID(int LocalGridIndex) const;
  const grid_info &GridInfoByIndex(int GridIndex) const;
  const grid &GridByIndex(int LocalGridIndex) const;

  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F, int, grid_event_flags,
    bool>())> event_listener_handle AddGridEventListener(F Listener) const {
    return GridEvent_.AddListener(std::move(Listener));
  }

protected:

  struct grid_record {
    grid_info GridInfo;
    explicit grid_record(grid_info &&GridInfo_):
      GridInfo(std::move(GridInfo_))
    {}
  };

  floating_ref_generator FloatingRefGenerator_;

  int NumDims_;

  map_noncontig<int,grid_record> GridRecords_;
  map_noncontig<int,grid> LocalGrids_;

  mutable event<void(int, grid_event_flags, bool)> GridEvent_;

};

}

}

#endif
