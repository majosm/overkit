// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ASSEMBLER_HPP_INCLUDED
#define OVK_CORE_ASSEMBLER_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Assembler.h>
#include <ovk/core/CollectMap.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/ConnectivityComponent.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/DisperseMap.hpp>
#include <ovk/core/DistributedField.hpp>
#include <ovk/core/DistributedRegionHash.hpp>
#include <ovk/core/Domain.hpp>
#include <ovk/core/ElemMap.hpp>
#include <ovk/core/ElemSet.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/GeometryComponent.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/OverlapComponent.hpp>
#include <ovk/core/Partition.hpp>
#include <ovk/core/RecvMap.hpp>
#include <ovk/core/SendMap.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/StateComponent.hpp>
#include <ovk/core/StringWrapper.hpp>

#include <mpi.h>

#include <memory>
#include <string>
#include <type_traits>

namespace ovk {

enum class occludes : typename std::underlying_type<ovk_occludes>::type {
  NONE = OVK_OCCLUDES_NONE,
  ALL = OVK_OCCLUDES_ALL,
  COARSE = OVK_OCCLUDES_COARSE
};

inline bool ValidOccludes(occludes Occludes) {
  return ovkValidOccludes(ovk_occludes(Occludes));
}

enum class connection_type : typename std::underlying_type<ovk_connection_type>::type {
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
    bool Overlappable(const elem<int,2> &GridIDPair) const;
    options &SetOverlappable(const elem<int,2> &GridIDPair, bool Overlappable);
    options &ResetOverlappable(const elem<int,2> &GridIDPair);
    double OverlapTolerance(const elem<int,2> &GridIDPair) const;
    options &SetOverlapTolerance(const elem<int,2> &GridIDPair, double OverlapTolerance);
    options &ResetOverlapTolerance(const elem<int,2> &GridIDPair);
    double OverlapAccelDepthAdjust(int MGridID) const;
    options &SetOverlapAccelDepthAdjust(int MGridID, double OverlapAccelDepthAdjust);
    options &ResetOverlapAccelDepthAdjust(int MGridID);
    double OverlapAccelResolutionAdjust(int MGridID) const;
    options &SetOverlapAccelResolutionAdjust(int MGridID, double OverlapAccelResolutionAdjust);
    options &ResetOverlapAccelResolutionAdjust(int MGridID);
    bool InferBoundaries(int GridID) const;
    options &SetInferBoundaries(int GridID, bool InferBoundaries);
    options &ResetInferBoundaries(int GridID);
    bool CutBoundaryHoles(const elem<int,2> &GridIDPair) const;
    options &SetCutBoundaryHoles(const elem<int,2> &GridIDPair, bool CutBoundaryHoles);
    options &ResetCutBoundaryHoles(const elem<int,2> &GridIDPair);
    occludes Occludes(const elem<int,2> &GridIDPair) const;
    options &SetOccludes(const elem<int,2> &GridIDPair, occludes Occludes);
    options &ResetOccludes(const elem<int,2> &GridIDPair);
    int EdgePadding(const elem<int,2> &GridIDPair) const;
    options &SetEdgePadding(const elem<int,2> &GridIDPair, int EdgePadding);
    options &ResetEdgePadding(const elem<int,2> &GridIDPair);
    int EdgeSmoothing(int NGridID) const;
    options &SetEdgeSmoothing(int NGridID, int EdgeSmoothing);
    options &ResetEdgeSmoothing(int NGridID);
    connection_type ConnectionType(const elem<int,2> &GridIDPair) const;
    options &SetConnectionType(const elem<int,2> &GridIDPair, connection_type ConnectionType);
    options &ResetConnectionType(const elem<int,2> &GridIDPair);
    int FringeSize(int NGridID) const;
    options &SetFringeSize(int NGridID, int FringeSize);
    options &ResetFringeSize(int NGridID);
    bool MinimizeOverlap(const elem<int,2> &GridIDPair) const;
    options &SetMinimizeOverlap(const elem<int,2> &GridIDPair, bool MinimizeOverlap);
    options &ResetMinimizeOverlap(const elem<int,2> &GridIDPair);
    bool DisjointConnections(const elem<int,2> &GridIDPair) const;
    options &SetDisjointConnections(const elem<int,2> &GridIDPair, bool DisjointConnections);
    options &ResetDisjointConnections(const elem<int,2> &GridIDPair);
  private:
    set<int> GridIDs_;
    elem_map<int,2,bool> Overlappable_;
    elem_map<int,2,double> OverlapTolerance_;
    map<int,double> OverlapAccelDepthAdjust_;
    map<int,double> OverlapAccelResolutionAdjust_;
    map<int,bool> InferBoundaries_;
    elem_map<int,2,bool> CutBoundaryHoles_;
    elem_map<int,2,occludes> Occludes_;
    elem_map<int,2,int> EdgePadding_;
    map<int,int> EdgeSmoothing_;
    elem_map<int,2,connection_type> ConnectionType_;
    map<int,int> FringeSize_;
    elem_map<int,2,bool> MinimizeOverlap_;
    elem_map<int,2,bool> DisjointConnections_;
    options() = default;
    void AddGrids(const set<int> &GridIDs);
    void RemoveGrids(const set<int> &GridIDs);
    template <typename T> T GetOption_(const map<int,T> &Option, int GridID, T DefaultValue) const;
    template <typename T> T GetOption_(const elem_map<int,2,T> &Option, const elem<int,2>
      &GridIDPair, T DefaultValue) const;
    template <typename T> void SetOption_(map<int,T> &Option, int GridID, T Value, T DefaultValue);
    template <typename T> void SetOption_(elem_map<int,2,T> &Option, const elem<int,2> &GridIDPair,
      T Value, T DefaulValue);
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
    set<int> AddGridsToOptions;
    set<int> RemoveGridsFromOptions;
    set<int> RemoveAssemblyManifestEntries;
  };

  struct assembly_manifest {
    elem_set<int,2> DetectOverlap;
    set<int> InferBoundaries;
    elem_set<int,2> CutBoundaryHoles;
    elem_set<int,2> ComputeOcclusion;
    elem_set<int,2> ApplyPadding;
    set<int> ApplySmoothing;
    elem_set<int,2> MinimizeOverlap;
    elem_set<int,2> GenerateConnectivity;
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
    core::partition_pool PartitionPool;
    distributed_field<bool> ActiveMask;
    distributed_field<bool> CellActiveMask;
    distributed_field<bool> DomainBoundaryMask;
    distributed_field<bool> InternalBoundaryMask;
    explicit local_grid_aux_data(core::partition_pool PartitionPool_):
      PartitionPool(std::move(PartitionPool_))
    {}
  };

  using bounding_box_hash = core::distributed_region_hash<double>;
  using bounding_box_hash_region_data = core::distributed_region_data<double>;
  using bounding_box_hash_retrieved_bins = core::distributed_region_hash_retrieved_bins<double>;

  struct local_overlap_m_aux_data {
    core::collect_map CollectMap;
    core::send_map SendMap;
  };

  struct local_overlap_n_aux_data {
    core::recv_map RecvMap;
    core::disperse_map DisperseMap;
    distributed_field<bool> OverlapMask;
    array<double> Volumes;
  };

  struct assembly_data {
    map<int,local_grid_aux_data> LocalGridAuxData;
    bounding_box_hash BoundingBoxHash;
    elem_map<int,2,local_overlap_m_aux_data> LocalOverlapMAuxData;
    elem_map<int,2,local_overlap_n_aux_data> LocalOverlapNAuxData;
    elem_map<int,2,distributed_field<bool>> ProjectedBoundaryMasks;
    map<int,distributed_field<bool>> OuterFringeMasks;
    elem_map<int,2,distributed_field<bool>> PairwiseOcclusionMasks;
    map<int,distributed_field<bool>> OcclusionMasks;
    map<int,distributed_field<bool>> OverlapMinimizationMasks;
    map<int,distributed_field<bool>> InnerFringeMasks;
    assembly_data(int NumDims, comm_view Comm);
  };

  optional<assembly_data> AssemblyData_;

  assembler(std::shared_ptr<context> &&Context, params &&Params);

  void OnGridEvent_(int GridID, grid_event_flags Flags, bool LastInSequence);
  void OnComponentEvent_(int ComponentID, component_event_flags Flags);
  void OnGeometryEvent_(int GridID, geometry_event_flags Flags, bool LastInSequence);
  void OnStateEvent_(int GridID, state_event_flags Flags, bool LastInSequence);
  void OnOverlapEvent_(const elem<int,2> &OverlapID, overlap_event_flags Flags, bool
    LastInSequence);
  void OnConnectivityEvent_(const elem<int,2> &ConnectivityID, connectivity_event_flags Flags, bool
    LastInSequence);

  void OnOptionsStartEdit_();
  void OnOptionsEndEdit_();

  void Update_();

  void AddGridsToOptions_();
  void RemoveGridsFromOptions_();
  void RemoveAssemblyManifestEntries_();

  void InitializeAssembly_();
  void ValidateOptions_();
  void DetectOverlap_();
  void InferBoundaries_();
  void CutBoundaryHoles_();
  void LocateOuterFringe_();
  void DetectOccluded_();
  void MinimizeOverlap_();
  void GenerateConnectivityData_();

  static constexpr int OVERLAP_TIME = core::profiler::ASSEMBLER_OVERLAP_TIME;
  static constexpr int OVERLAP_BB_TIME = core::profiler::ASSEMBLER_OVERLAP_BB_TIME;
  static constexpr int OVERLAP_BB_SUBDIVIDE_TIME = core::profiler::ASSEMBLER_OVERLAP_BB_SUBDIVIDE_TIME;
  static constexpr int OVERLAP_BB_HASH_CREATE_TIME = core::profiler::ASSEMBLER_OVERLAP_BB_HASH_CREATE_TIME;
  static constexpr int OVERLAP_BB_HASH_MAP_TIME = core::profiler::ASSEMBLER_OVERLAP_BB_HASH_MAP_TIME;
  static constexpr int OVERLAP_BB_HASH_RETRIEVE_TIME = core::profiler::ASSEMBLER_OVERLAP_BB_HASH_RETRIEVE_TIME;
  static constexpr int OVERLAP_CONNECT_TIME = core::profiler::ASSEMBLER_OVERLAP_CONNECT_TIME;
  static constexpr int OVERLAP_TRANSFER_TIME = core::profiler::ASSEMBLER_OVERLAP_TRANSFER_TIME;
  static constexpr int OVERLAP_SEARCH_TIME = core::profiler::ASSEMBLER_OVERLAP_SEARCH_TIME;
  static constexpr int OVERLAP_SEARCH_BUILD_ACCEL_TIME = core::profiler::ASSEMBLER_OVERLAP_SEARCH_BUILD_ACCEL_TIME;
  static constexpr int OVERLAP_SEARCH_QUERY_ACCEL_TIME = core::profiler::ASSEMBLER_OVERLAP_SEARCH_QUERY_ACCEL_TIME;
  static constexpr int OVERLAP_FILL_TIME = core::profiler::ASSEMBLER_OVERLAP_FILL_TIME;
  static constexpr int OVERLAP_CREATE_EXCHANGE_TIME = core::profiler::ASSEMBLER_OVERLAP_CREATE_EXCHANGE_TIME;
  static constexpr int OVERLAP_CREATE_AUX_TIME = core::profiler::ASSEMBLER_OVERLAP_CREATE_AUX_TIME;
  static constexpr int INFER_BOUNDARIES_TIME = core::profiler::ASSEMBLER_INFER_BOUNDARIES_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_PROJECT_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_PROJECT_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_DETECT_EXTERIOR_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_DETECT_EXTERIOR_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_UPDATE_AUX_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_UPDATE_AUX_TIME;
  static constexpr int LOCATE_OUTER_FRINGE_TIME = core::profiler::ASSEMBLER_LOCATE_OUTER_FRINGE_TIME;
  static constexpr int OCCLUSION_TIME = core::profiler::ASSEMBLER_OCCLUSION_TIME;
  static constexpr int OCCLUSION_PAIRWISE_TIME = core::profiler::ASSEMBLER_OCCLUSION_PAIRWISE_TIME;
  static constexpr int OCCLUSION_PAD_SMOOTH_TIME = core::profiler::ASSEMBLER_OCCLUSION_PAD_SMOOTH_TIME;
  static constexpr int OCCLUSION_ACCUMULATE_TIME = core::profiler::ASSEMBLER_OCCLUSION_ACCUMULATE_TIME;
  static constexpr int MINIMIZE_OVERLAP_TIME = core::profiler::ASSEMBLER_MINIMIZE_OVERLAP_TIME;
  static constexpr int CONNECTIVITY_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_TIME;
  static constexpr int CONNECTIVITY_LOCATE_RECEIVERS_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_LOCATE_RECEIVERS_TIME;
  static constexpr int CONNECTIVITY_CHOOSE_DONORS_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_CHOOSE_DONORS_TIME;
  static constexpr int CONNECTIVITY_FILL_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_FILL_TIME;

};

assembler CreateAssembler(std::shared_ptr<context> Context, assembler::params Params={});

}

#endif
