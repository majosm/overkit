// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

template <typename T, typename... Args, OVK_FUNCDEF_REQUIRES(std::is_constructible<T, const
  core::domain_base &, Args &&...>::value)> void domain::CreateComponent(int ComponentID, Args &&...
  Arguments) {

  MPI_Barrier(Comm_);

  OVK_DEBUG_ASSERT(ComponentID >= 0, "Invalid component ID.");
  OVK_DEBUG_ASSERT(!ComponentExists(ComponentID), "Component %i already exists.", ComponentID);

  T Component(static_cast<const core::domain_base &>(*this), std::forward<Args>(Arguments)...);
  Components_.Insert(ComponentID, std::move(Component));

  MPI_Barrier(Comm_);

  ComponentEvent_.Trigger(ComponentID, component_event_flags::CREATE);

  MPI_Barrier(Comm_);

}

template <typename T> const T &domain::Component(int ComponentID) const {

  OVK_DEBUG_ASSERT(ComponentID >= 0, "Invalid component ID.");
  OVK_DEBUG_ASSERT(ComponentExists(ComponentID), "Component %i does not exist.", ComponentID);

  const component_data &Data = Components_(ComponentID);
  const core::component &Component = Data.Component;

  return Component.Get<T>();

}

template <typename T> edit_handle<T> domain::EditComponent(int ComponentID) {

  OVK_DEBUG_ASSERT(ComponentID >= 0, "Invalid component ID.");
  OVK_DEBUG_ASSERT(ComponentExists(ComponentID), "Component %i does not exist.", ComponentID);

  component_data &Data = Components_(ComponentID);
  core::component &Component = Data.Component;
  editor &Editor = Data.Editor;

  if (!Editor.Active()) {
    MPI_Barrier(Comm_);
    Component.StartEdit();
    floating_ref<domain> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef, ComponentID] {
      domain &Domain = *FloatingRef;
      component_data &Data = Domain.Components_(ComponentID);
      core::component &Component = Data.Component;
      Component.EndEdit();
      MPI_Barrier(Domain.Comm_);
      Domain.ComponentEvent_.Trigger(ComponentID, component_event_flags::EDIT);
      MPI_Barrier(Domain.Comm_);
    };
    Editor.Activate(std::move(DeactivateFunc));
  }

  return Editor.Edit(Component.Get<T>());

}

}
