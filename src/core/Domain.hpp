// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DOMAIN_HPP_INCLUDED
#define OVK_CORE_DOMAIN_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Component.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/DomainBase.hpp>
#include <ovk/core/Editor.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/IDMap.hpp>
#include <ovk/core/IDSet.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/Requires.hpp>

#include <mpi.h>

#include <algorithm>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

namespace ovk {

enum class component_event_flags : int {
  NONE = OVK_COMPONENT_EVENT_FLAGS_NONE,
  CREATE = OVK_COMPONENT_EVENT_FLAGS_CREATE,
  DESTROY = OVK_COMPONENT_EVENT_FLAGS_DESTROY,
  EDIT = OVK_COMPONENT_EVENT_FLAGS_EDIT,
  ALL = OVK_COMPONENT_EVENT_FLAGS_ALL
};

constexpr inline component_event_flags operator|(component_event_flags Left, component_event_flags
  Right) {
  return component_event_flags(int(Left) | int(Right));
}
constexpr inline component_event_flags operator&(component_event_flags Left, component_event_flags
  Right) {
  return component_event_flags(int(Left) & int(Right));
}
constexpr inline component_event_flags operator^(component_event_flags Left, component_event_flags
  Right) {
  return component_event_flags(int(Left) ^ int(Right));
}
constexpr inline component_event_flags operator~(component_event_flags EventFlags) {
  return component_event_flags(~int(EventFlags));
}
inline component_event_flags operator|=(component_event_flags &Left, component_event_flags Right) {
  return Left = Left | Right;
}
inline component_event_flags operator&=(component_event_flags &Left, component_event_flags Right) {
  return Left = Left & Right;
}
inline component_event_flags operator^=(component_event_flags &Left, component_event_flags Right) {
  return Left = Left ^ Right;
}

class domain : public core::domain_base {

public:

  class params {
  public:
    params() = default;
    const std::string &Name() const { return *Name_; }
    params &SetName(std::string Name);
    int Dimension() const { return NumDims_; }
    params &SetDimension(int NumDims);
    MPI_Comm Comm() const { return Comm_; }
    params &SetComm(MPI_Comm Comm);
  private:
    core::string_wrapper Name_ = "Domain";
    int NumDims_ = 2;
    MPI_Comm Comm_ = MPI_COMM_NULL;
    friend class domain;
  };

  domain(const domain &Other) = delete;
  domain(domain &&Other) noexcept = default;

  domain &operator=(const domain &Other) = delete;
  domain &operator=(domain &&Other) noexcept = default;

  ~domain() noexcept;

  floating_ref<const domain> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate<const domain>();
  }
  floating_ref<domain> GetFloatingRef() { return FloatingRefGenerator_.Generate<domain>(); }

  const id_set<1> &ComponentIDs() const { return Components_.Keys(); }

  bool ComponentExists(int ComponentID) const;

  template <typename T, typename... Args, OVK_FUNCDECL_REQUIRES(std::is_constructible<T, const
    core::domain_base &, Args &&...>::value)> void CreateComponent(int ComponentID, Args &&...
    Arguments);
  void DestroyComponent(int ComponentID);

  template <typename T> const T &Component(int ComponentID) const;
  bool EditingComponent(int ComponentID) const;
  template <typename T> edit_handle<T> EditComponent(int ComponentID);
  void RestoreComponent(int ComponentID);

  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F, int, component_event_flags>()
    )> event_listener_handle AddComponentEventListener(F Listener) const {
    return ComponentEvent_.AddListener(std::move(Listener));
  }

  static domain internal_Create(std::shared_ptr<context> &&Context, params &&Params);

private:

  struct component_data {
    core::component Component;
    editor Editor;
    explicit component_data(core::component &&Component_):
      Component(std::move(Component_))
    {}
  };

  id_map<1,component_data,false> Components_;

  mutable event<void(int, component_event_flags)> ComponentEvent_;

  domain(std::shared_ptr<context> &&Context, params &&Params);

  void DestroyExchangesForGrid_(int GridID);

};

domain CreateDomain(std::shared_ptr<context> Context, domain::params Params);

}

#include <ovk/core/Domain.inl>

#endif
