// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ASSEMBLER_HPP_INCLUDED
#define OVK_CORE_ASSEMBLER_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Assembler.h>
#include <ovk/core/Comm.hpp>
#include <ovk/core/ConnectivityComponent.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Domain.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/GeometryComponent.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/IDMap.hpp>
#include <ovk/core/IDSet.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/OverlapComponent.hpp>
#include <ovk/core/StateComponent.hpp>
#include <ovk/core/StringWrapper.hpp>

#include <mpi.h>

#include <memory>
#include <string>

namespace ovk {

enum class occludes {
  NONE = OVK_OCCLUDES_NONE,
  ALL = OVK_OCCLUDES_ALL,
  COARSE = OVK_OCCLUDES_COARSE
};

inline bool ValidOccludes(occludes Occludes) {
  return ovkValidOccludes(ovk_occludes(Occludes));
}

enum class connection_type {
  NONE = OVK_CONNECTION_NONE,
  NEAREST = OVK_CONNECTION_NEAREST,
  LINEAR = OVK_CONNECTION_LINEAR,
  CUBIC = OVK_CONNECTION_CUBIC
};

inline bool ValidConnectionType(connection_type ConnectionType) {
  return ovkValidConnectionType(ovk_connection_type(ConnectionType));
}

class assembler {

public:

  class params {
  public:
    params() = default;
    const std::string &Name() const { return Name_; }
    params &SetName(std::string Name);
  private:
    core::string_wrapper Name_ = "Assembler";
    friend class assembler;
  };

  class bindings {
  public:
    bindings() = default;
    int GeometryComponentID() const { return GeometryComponentID_; }
    bindings &SetGeometryComponentID(int GeometryComponentID);
    int StateComponentID() const { return StateComponentID_; }
    bindings &SetStateComponentID(int StateComponentID);
    int OverlapComponentID() const { return OverlapComponentID_; }
    bindings &SetOverlapComponentID(int OverlapComponentID);
    int ConnectivityComponentID() const { return ConnectivityComponentID_; }
    bindings &SetConnectivityComponentID(int ConnectivityComponentID);
  private:
    int GeometryComponentID_ = -1;
    int StateComponentID_ = -1;
    int OverlapComponentID_ = -1;
    int ConnectivityComponentID_ = -1;
    friend class assembler;
  };

  class options {
  public:
    bool Overlappable(int MGridID, int NGridID) const;
    options &SetOverlappable(int MGridID, int NGridID, bool Overlappable);
    options &ResetOverlappable(int MGridID, int NGridID);
    double OverlapTolerance(int MGridID, int NGridID) const;
    options &SetOverlapTolerance(int MGridID, int NGridID, double OverlapTolerance);
    options &ResetOverlapTolerance(int MGridID, int NGridID);
    double OverlapAccelDepthAdjust(int MGridID) const;
    options &SetOverlapAccelDepthAdjust(int MGridID, double OverlapAccelDepthAdjust);
    options &ResetOverlapAccelDepthAdjust(int MGridID);
    double OverlapAccelResolutionAdjust(int MGridID) const;
    options &SetOverlapAccelResolutionAdjust(int MGridID, double OverlapAccelResolutionAdjust);
    options &ResetOverlapAccelResolutionAdjust(int MGridID);
    bool InferBoundaries(int GridID) const;
    options &SetInferBoundaries(int GridID, bool InferBoundaries);
    options &ResetInferBoundaries(int GridID);
    bool CutBoundaryHoles(int MGridID, int NGridID) const;
    options &SetCutBoundaryHoles(int MGridID, int NGridID, bool CutBoundaryHoles);
    options &ResetCutBoundaryHoles(int MGridID, int NGridID);
    occludes Occludes(int MGridID, int NGridID) const;
    options &SetOccludes(int MGridID, int NGridID, occludes Occludes);
    options &ResetOccludes(int MGridID, int NGridID);
    int EdgePadding(int MGridID, int NGridID) const;
    options &SetEdgePadding(int MGridID, int NGridID, int EdgePadding);
    options &ResetEdgePadding(int MGridID, int NGridID);
    int EdgeSmoothing(int NGridID) const;
    options &SetEdgeSmoothing(int NGridID, int EdgeSmoothing);
    options &ResetEdgeSmoothing(int NGridID);
    connection_type ConnectionType(int MGridID, int NGridID) const;
    options &SetConnectionType(int MGridID, int NGridID, connection_type ConnectionType);
    options &ResetConnectionType(int MGridID, int NGridID);
    int FringeSize(int NGridID) const;
    options &SetFringeSize(int NGridID, int FringeSize);
    options &ResetFringeSize(int NGridID);
    bool MinimizeOverlap(int MGridID, int NGridID) const;
    options &SetMinimizeOverlap(int MGridID, int NGridID, bool MinimizeOverlap);
    options &ResetMinimizeOverlap(int MGridID, int NGridID);
  private:
    id_set<1> GridIDs_;
    id_map<2,bool> Overlappable_;
    id_map<2,double> OverlapTolerance_;
    id_map<1,double> OverlapAccelDepthAdjust_;
    id_map<1,double> OverlapAccelResolutionAdjust_;
    id_map<1,bool> InferBoundaries_;
    id_map<2,bool> CutBoundaryHoles_;
    id_map<2,occludes> Occludes_;
    id_map<2,int> EdgePadding_;
    id_map<1,int> EdgeSmoothing_;
    id_map<2,connection_type> ConnectionType_;
    id_map<1,int> FringeSize_;
    id_map<2,bool> MinimizeOverlap_;
    options() = default;
    void AddGrids(const id_set<1> &GridIDs);
    void RemoveGrids(const id_set<1> &GridIDs);
    template <typename T> T GetOption_(const id_map<1,T> &Option, int GridID, T DefaultValue) const;
    template <typename T> T GetOption_(const id_map<2,T> &Option, int MGridID, int NGridID, T
      DefaultValue) const;
    template <typename T> void SetOption_(id_map<1,T> &Option, int GridID, T Value, T DefaultValue);
    template <typename T> void SetOption_(id_map<2,T> &Option, int MGridID, int NGridID, T Value, T
      DefaulValue);
    void PrintOptions_();
    friend class assembler;
  };

  assembler(const assembler &Other) = delete;
  assembler(assembler &&Other) noexcept = default;

  assembler &operator=(const assembler &Other) = delete;
  assembler &operator=(assembler &&Other) noexcept = default;

  ~assembler() noexcept;

  const context &Context() const;
  context &Context();
  const std::shared_ptr<context> &SharedContext() const;

  const std::string &Name() const { return *Name_; }

  bool Bound() const;
  void Bind(domain &Domain, bindings Bindings);
  void Unbind();

  const domain &Domain() const;
  domain &Domain();

  const options &Options() const { return Options_; }
  bool EditingOptions() const;
  edit_handle<options> EditOptions();
  void RestoreOptions();

  void Assemble();

  static assembler internal_Create(std::shared_ptr<context> &&Context, params &&Params);

private:

  struct update_manifest {
    id_set<1> AddGridsToOptions;
    id_set<1> RemoveGridsFromOptions;
    id_set<1> RemoveAssemblyManifestEntries;
  };

  struct assembly_manifest {
    id_set<2> DetectOverlap;
    id_set<1> InferBoundaries;
    id_set<2> CutBoundaryHoles;
    id_set<2> ComputeOcclusion;
    id_set<2> ApplyPadding;
    id_set<1> ApplySmoothing;
    id_set<2> MinimizeOverlap;
    id_set<2> GenerateConnectivity;
  };

  floating_ref_generator FloatingRefGenerator_;

  std::shared_ptr<context> Context_;

  core::string_wrapper Name_;

  floating_ref<domain> Domain_;
  event_listener_handle GridEventListener_;
  event_listener_handle ComponentEventListener_;

  int GeometryComponentID_ = -1;
  event_listener_handle GeometryEventListener_;

  int StateComponentID_ = -1;
  event_listener_handle StateEventListener_;

  int OverlapComponentID_ = -1;
  event_listener_handle OverlapEventListener_;

  int ConnectivityComponentID_ = -1;
  event_listener_handle ConnectivityEventListener_;

  options Options_;
  options CachedOptions_;
  editor OptionsEditor_;

  update_manifest UpdateManifest_;

  assembly_manifest AssemblyManifest_;

  struct local_grid_aux_data {
    field<bool> ActiveMask;
    field<bool> CellActiveMask;
    field<bool> DomainBoundaryMask;
    field<bool> InternalBoundaryMask;
  };

  struct assembly_data {
    id_map<1,local_grid_aux_data> LocalGridAuxData;
    assembly_data(int NumDims, comm_view Comm);
  };

  optional<assembly_data> AssemblyData_;

  assembler(std::shared_ptr<context> &&Context, params &&Params);

  void OnGridEvent_(int GridID, grid_event_flags Flags, bool LastInSequence);
  void OnComponentEvent_(int ComponentID, component_event_flags Flags);
  void OnGeometryEvent_(int GridID, geometry_event_flags Flags, bool LastInSequence);
  void OnStateEvent_(int GridID, state_event_flags Flags, bool LastInSequence);
  void OnOverlapEvent_(int MGridID, int NGridID, overlap_event_flags Flags, bool LastInSequence);
  void OnConnectivityEvent_(int MGrid, int NGridID, connectivity_event_flags Flags, bool
    LastInSequence);

  void OnOptionsStartEdit_();
  void OnOptionsEndEdit_();

  void Update_();

  void AddGridsToOptions_();
  void RemoveGridsFromOptions_();
  void RemoveAssemblyManifestEntries_();

  void InitializeAssembly_();

};

assembler CreateAssembler(std::shared_ptr<context> Context, assembler::params Params={});

}

#endif
