// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ASSEMBLER_HPP_INCLUDED
#define OVK_CORE_ASSEMBLER_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Assembler.h>
#include <ovk/core/Box.hpp>
#include <ovk/core/CollectMap.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/ConnectivityComponent.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/DataTypeOps.hpp>
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
#include <ovk/core/HashableRegionTraits.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/MPISerializableTraits.hpp>
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

namespace assembler_internal {
struct fragment {
  int GridID = -1;
  int ID = -1;
  range CellRange;
//   static_array<double,MAX_DIMS*MAX_DIMS,2> Orientation;
  tuple<tuple<double>> Orientation;
  box Box;
  fragment() = default;
};
}

namespace core {
template <> struct hashable_region_traits<assembler_internal::fragment> {
  using coord_type = double;
  static box ComputeExtents(int NumDims, const assembler_internal::fragment &Region) {
    box Extents;
    switch (NumDims) {
    case 1:
      Extents = ComputeExtents1D(Region);
      break;
    case 2:
      Extents = ComputeExtents2D(Region);
      break;
    default:
      Extents = ComputeExtents3D(Region);
      break;
    }
    return Extents;
  }
  static box ComputeExtents1D(const assembler_internal::fragment &Region) {
    return Region.Box;
  }
  static box ComputeExtents2D(const assembler_internal::fragment &Region) {
    auto OrientInv = [&Region](const elem<double,2> &Vec) -> elem<double,2> {
      const tuple<tuple<double>> &O = Region.Orientation;
      return {
        O(0)(0)*Vec(0) + O(0)(1)*Vec(1),
        O(1)(0)*Vec(0) + O(1)(1)*Vec(1)
      };
    };
    elem<double,2> RegionCorners[4] = {
      OrientInv({Region.Box.Begin(0),Region.Box.Begin(1)}),
      OrientInv({  Region.Box.End(0),Region.Box.Begin(1)}),
      OrientInv({Region.Box.Begin(0),  Region.Box.End(1)}),
      OrientInv({  Region.Box.End(0),  Region.Box.End(1)})
    };
    box Extents = MakeEmptyBox(2);
    for (int l = 0; l < 4; ++l) {
      Extents = ExtendBox(Extents, {RegionCorners[l](0), RegionCorners[l](1), 0.});
    }
    return Extents;
  }
  static box ComputeExtents3D(const assembler_internal::fragment &Region) {
    auto OrientInv = [&Region](const elem<double,3> &Vec) -> elem<double,3> {
      const tuple<tuple<double>> &O = Region.Orientation;
      return {
        O(0)(0)*Vec(0) + O(0)(1)*Vec(1) + O(0)(2)*Vec(2),
        O(1)(0)*Vec(0) + O(1)(1)*Vec(1) + O(1)(2)*Vec(2),
        O(2)(0)*Vec(0) + O(2)(1)*Vec(1) + O(2)(2)*Vec(2)
      };
    };
    elem<double,3> RegionCorners[8] = {
      OrientInv({Region.Box.Begin(0),Region.Box.Begin(1), Region.Box.Begin(2)}),
      OrientInv({  Region.Box.End(0),Region.Box.Begin(1), Region.Box.Begin(2)}),
      OrientInv({Region.Box.Begin(0),  Region.Box.End(1), Region.Box.Begin(2)}),
      OrientInv({  Region.Box.End(0),  Region.Box.End(1), Region.Box.Begin(2)}),
      OrientInv({Region.Box.Begin(0),Region.Box.Begin(1),   Region.Box.End(2)}),
      OrientInv({  Region.Box.End(0),Region.Box.Begin(1),   Region.Box.End(2)}),
      OrientInv({Region.Box.Begin(0),  Region.Box.End(1),   Region.Box.End(2)}),
      OrientInv({  Region.Box.End(0),  Region.Box.End(1),   Region.Box.End(2)})
    };
    box Extents = MakeEmptyBox(3);
    for (int l = 0; l < 8; ++l) {
      Extents = ExtendBox(Extents, RegionCorners[l]);
    }
    return Extents;
  }
  template <typename IndexerType> static set<typename IndexerType::index_type> MapToBins(int
    NumDims, const range &BinRange, const IndexerType &BinIndexer, const tuple<double> &LowerCorner,
    const tuple<double> &BinSize, const assembler_internal::fragment &Region) {
    set<typename IndexerType::index_type> Bins;
    switch (NumDims) {
    case 1:
      Bins = MapToBins1D(BinRange, BinIndexer, LowerCorner, BinSize, Region);
      break;
    case 2:
      Bins = MapToBins2D(BinRange, BinIndexer, LowerCorner, BinSize, Region);
      break;
    default:
      Bins = MapToBins3D(BinRange, BinIndexer, LowerCorner, BinSize, Region);
      break;
    }
    return Bins;
  }
  template <typename IndexerType> static set<typename IndexerType::index_type> MapToBins1D(const
    range &BinRange, const IndexerType &BinIndexer, const tuple<double> &LowerCorner, const
    tuple<double> &BinSize, const assembler_internal::fragment &Region) {
    return hashable_region_traits<box>::MapToBins(1, BinRange, BinIndexer, LowerCorner, BinSize,
      Region.Box);
  }
  template <typename IndexerType> static set<typename IndexerType::index_type> MapToBins2D(const
    range &BinRange, const IndexerType &BinIndexer, const tuple<double> &LowerCorner, const
    tuple<double> &BinSize, const assembler_internal::fragment &Region) {
    using index_type = typename IndexerType::index_type;
    auto OrientInv = [&Region](const elem<double,2> &Vec) -> elem<double,2> {
      const tuple<tuple<double>> &O = Region.Orientation;
      return {
        O(0)(0)*Vec(0) + O(0)(1)*Vec(1),
        O(1)(0)*Vec(0) + O(1)(1)*Vec(1)
      };
    };
    elem<double,2> RegionCorners[4] = {
      OrientInv({Region.Box.Begin(0),Region.Box.Begin(1)}),
      OrientInv({  Region.Box.End(0),Region.Box.Begin(1)}),
      OrientInv({Region.Box.Begin(0),  Region.Box.End(1)}),
      OrientInv({  Region.Box.End(0),  Region.Box.End(1)})
    };
    box RegionExtents = MakeEmptyBox(2);
    for (int l = 0; l < 4; ++l) {
      RegionExtents = ExtendBox(RegionExtents, {RegionCorners[l](0), RegionCorners[l](1), 0.});
    }
    set<index_type> CandidateBins = hashable_region_traits<box>::MapToBins(2, BinRange, BinIndexer,
      LowerCorner, BinSize, RegionExtents);
    set<index_type> Bins;
    Bins.Reserve(CandidateBins.Count());
    elem<double,2> SeparatingLineAxes[4] = {
      {1., 0.},
      {0., 1.},
      OrientInv({1.,0.}),
      OrientInv({0.,1.})
    };
    auto OverlapsAlongAxis = [](array_view<const elem<double,2>> Corners1, array_view<const
      elem<double,2>> Corners2, const elem<double,2> &Axis) -> bool {
      auto IntervalAlongAxis = [](const array_view<const elem<double,2>> &Corners, const
        elem<double,2> &Axis) -> interval<double> {
        interval<double> Interval;
        Interval.Begin(0) = std::numeric_limits<double>::max();
        Interval.End(0) = std::numeric_limits<double>::min();
        for (auto &Corner : Corners) {
          double Dot = Corner(0)*Axis(0) + Corner(1)*Axis(1);
          Interval.Begin(0) = Min(Interval.Begin(0), Dot);
          Interval.End(0) = Max(Interval.End(0), Dot);
        }
        return Interval;
      };
      interval<double> Interval1 = IntervalAlongAxis(Corners1, Axis);
      interval<double> Interval2 = IntervalAlongAxis(Corners2, Axis);
      return Interval2.End(0) >= Interval1.Begin(0) && Interval1.End(0) >= Interval2.Begin(0);
    };
    for (index_type iBin : CandidateBins) {
      tuple<int> BinLoc = BinIndexer.ToTuple(iBin);
      box BinExtents = MakeEmptyBox(2);
      for (int iDim = 0; iDim < 2; ++iDim) {
        BinExtents.Begin(iDim) = LowerCorner(iDim) + double(BinLoc(iDim))*BinSize(iDim);
        BinExtents.End(iDim) = LowerCorner(iDim) + double(BinLoc(iDim)+1)*BinSize(iDim);
      }
      elem<double,2> BinCorners[4] = {
        {BinExtents.Begin(0),BinExtents.Begin(1)},
        {  BinExtents.End(0),BinExtents.Begin(1)},
        {BinExtents.Begin(0),  BinExtents.End(1)},
        {  BinExtents.End(0),  BinExtents.End(1)}
      };
      bool Overlaps = true;
      for (auto &Axis : SeparatingLineAxes) {
        Overlaps = Overlaps && OverlapsAlongAxis(RegionCorners, BinCorners, Axis);
        if (!Overlaps) break;
      }
      if (Overlaps) {
        Bins.Insert(iBin);
      }
    }
    return Bins;
  }
  template <typename IndexerType> static set<typename IndexerType::index_type> MapToBins3D(const
    range &BinRange, const IndexerType &BinIndexer, const tuple<double> &LowerCorner, const
    tuple<double> &BinSize, const assembler_internal::fragment &Region) {
    using index_type = typename IndexerType::index_type;
    auto OrientInv = [&Region](const elem<double,3> &Vec) -> elem<double,3> {
      const tuple<tuple<double>> &O = Region.Orientation;
      return {
        O(0)(0)*Vec(0) + O(0)(1)*Vec(1) + O(0)(2)*Vec(2),
        O(1)(0)*Vec(0) + O(1)(1)*Vec(1) + O(1)(2)*Vec(2),
        O(2)(0)*Vec(0) + O(2)(1)*Vec(1) + O(2)(2)*Vec(2)
      };
    };
    elem<double,3> RegionCorners[8] = {
      OrientInv({Region.Box.Begin(0),Region.Box.Begin(1), Region.Box.Begin(2)}),
      OrientInv({  Region.Box.End(0),Region.Box.Begin(1), Region.Box.Begin(2)}),
      OrientInv({Region.Box.Begin(0),  Region.Box.End(1), Region.Box.Begin(2)}),
      OrientInv({  Region.Box.End(0),  Region.Box.End(1), Region.Box.Begin(2)}),
      OrientInv({Region.Box.Begin(0),Region.Box.Begin(1),   Region.Box.End(2)}),
      OrientInv({  Region.Box.End(0),Region.Box.Begin(1),   Region.Box.End(2)}),
      OrientInv({Region.Box.Begin(0),  Region.Box.End(1),   Region.Box.End(2)}),
      OrientInv({  Region.Box.End(0),  Region.Box.End(1),   Region.Box.End(2)})
    };
    box RegionExtents = MakeEmptyBox(3);
    for (int l = 0; l < 8; ++l) {
      RegionExtents = ExtendBox(RegionExtents, RegionCorners[l]);
    }
    set<index_type> CandidateBins = hashable_region_traits<box>::MapToBins(3, BinRange, BinIndexer,
      LowerCorner, BinSize, RegionExtents);
    set<index_type> Bins;
    Bins.Reserve(CandidateBins.Count());
    elem<double,3> RegionAxisI = OrientInv({1.,0.,0.});
    elem<double,3> RegionAxisJ = OrientInv({0.,1.,0.});
    elem<double,3> RegionAxisK = OrientInv({0.,0.,1.});
    elem<double,3> SeparatingPlaneAxes[15] = {
      {1., 0., 0.},
      {0., 1., 0.},
      {0., 0., 1.},
      RegionAxisI,
      RegionAxisJ,
      RegionAxisK,
      {0., -RegionAxisI(2), RegionAxisI(1)},
      {0., -RegionAxisJ(2), RegionAxisJ(1)},
      {0., -RegionAxisK(2), RegionAxisK(1)},
      {RegionAxisI(2), 0., -RegionAxisI(0)},
      {RegionAxisJ(2), 0., -RegionAxisJ(0)},
      {RegionAxisK(2), 0., -RegionAxisK(0)},
      {-RegionAxisI(1), RegionAxisI(0), 0.},
      {-RegionAxisJ(1), RegionAxisJ(0), 0.},
      {-RegionAxisJ(1), RegionAxisK(0), 0.}
    };
    auto OverlapsAlongAxis = [](array_view<const elem<double,3>> Corners1, array_view<const
      elem<double,3>> Corners2, const elem<double,3> &Axis) -> bool {
      auto IntervalAlongAxis = [](const array_view<const elem<double,3>> &Corners, const
        elem<double,3> &Axis) -> interval<double> {
        interval<double> Interval;
        Interval.Begin(0) = std::numeric_limits<double>::max();
        Interval.End(0) = std::numeric_limits<double>::min();
        for (auto &Corner : Corners) {
          double Dot = Corner(0)*Axis(0) + Corner(1)*Axis(1) + Corner(2)*Axis(2);
          Interval.Begin(0) = Min(Interval.Begin(0), Dot);
          Interval.End(0) = Max(Interval.End(0), Dot);
        }
        return Interval;
      };
      interval<double> Interval1 = IntervalAlongAxis(Corners1, Axis);
      interval<double> Interval2 = IntervalAlongAxis(Corners2, Axis);
      return Interval2.End(0) >= Interval1.Begin(0) && Interval1.End(0) >= Interval2.Begin(0);
    };
    for (index_type iBin : CandidateBins) {
      tuple<int> BinLoc = BinIndexer.ToTuple(iBin);
      box BinExtents = MakeEmptyBox(3);
      for (int iDim = 0; iDim < 3; ++iDim) {
        BinExtents.Begin(iDim) = LowerCorner(iDim) + double(BinLoc(iDim))*BinSize(iDim);
        BinExtents.End(iDim) = LowerCorner(iDim) + double(BinLoc(iDim)+1)*BinSize(iDim);
      }
      elem<double,3> BinCorners[8] = {
        {BinExtents.Begin(0),BinExtents.Begin(1), BinExtents.Begin(2)},
        {  BinExtents.End(0),BinExtents.Begin(1), BinExtents.Begin(2)},
        {BinExtents.Begin(0),  BinExtents.End(1), BinExtents.Begin(2)},
        {  BinExtents.End(0),  BinExtents.End(1), BinExtents.Begin(2)},
        {BinExtents.Begin(0),BinExtents.Begin(1),   BinExtents.End(2)},
        {  BinExtents.End(0),BinExtents.Begin(1),   BinExtents.End(2)},
        {BinExtents.Begin(0),  BinExtents.End(1),   BinExtents.End(2)},
        {  BinExtents.End(0),  BinExtents.End(1),   BinExtents.End(2)}
      };
      bool Overlaps = true;
      for (auto &Axis : SeparatingPlaneAxes) {
        Overlaps = Overlaps && OverlapsAlongAxis(RegionCorners, BinCorners, Axis);
        if (!Overlaps) break;
      }
      if (Overlaps) {
        Bins.Insert(iBin);
      }
    }
    return Bins;
  }
};
template <> struct mpi_serializable_traits<assembler_internal::fragment> {
  using packed_type = assembler_internal::fragment;
  static handle<MPI_Datatype> CreateMPIType() {
    auto RangeMPIType = mpi_serializable_traits<range>::CreateMPIType();
    auto OrientationMPIType = mpi_serializable_traits<tuple<tuple<double>>>::CreateMPIType();
    auto BoxMPIType = mpi_serializable_traits<box>::CreateMPIType();
    assembler_internal::fragment Dummy;
    return CreateMPIStructType(sizeof(Dummy), {
      {GetByteOffset(Dummy,Dummy.GridID),2,MPI_INT},
      {GetByteOffset(Dummy,Dummy.CellRange),1,RangeMPIType.Get()},
      {GetByteOffset(Dummy,Dummy.Orientation),1,OrientationMPIType.Get()},
      {GetByteOffset(Dummy,Dummy.Box),1,BoxMPIType.Get()}
    });
  }
  static packed_type Pack(const assembler_internal::fragment &Fragment) { return Fragment; }
  static assembler_internal::fragment Unpack(const packed_type &PackedFragment) {
    return PackedFragment;
  }
};
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

  using fragment = assembler_internal::fragment;

  using fragment_hash = core::distributed_region_hash<fragment>;
  using fragment_hash_region_data = core::distributed_region_data<fragment>;
  using fragment_hash_retrieved_bins = core::distributed_region_hash_retrieved_bins<fragment>;

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
    fragment_hash FragmentHash;
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
  static constexpr int OVERLAP_FRAGMENT_TIME = core::profiler::ASSEMBLER_OVERLAP_FRAGMENT_TIME;
  static constexpr int OVERLAP_HASH_TIME = core::profiler::ASSEMBLER_OVERLAP_HASH_TIME;
  static constexpr int OVERLAP_HASH_CREATE_TIME = core::profiler::ASSEMBLER_OVERLAP_HASH_CREATE_TIME;
  static constexpr int OVERLAP_HASH_MAP_TO_BINS_TIME = core::profiler::ASSEMBLER_OVERLAP_HASH_MAP_TO_BINS_TIME;
  static constexpr int OVERLAP_HASH_RETRIEVE_BINS_TIME = core::profiler::ASSEMBLER_OVERLAP_HASH_RETRIEVE_BINS_TIME;
  static constexpr int OVERLAP_CONNECT_TIME = core::profiler::ASSEMBLER_OVERLAP_CONNECT_TIME;
  static constexpr int OVERLAP_SEARCH_TIME = core::profiler::ASSEMBLER_OVERLAP_SEARCH_TIME;
  static constexpr int OVERLAP_SEARCH_BUILD_ACCEL_TIME = core::profiler::ASSEMBLER_OVERLAP_SEARCH_BUILD_ACCEL_TIME;
  static constexpr int OVERLAP_SEARCH_QUERY_ACCEL_TIME = core::profiler::ASSEMBLER_OVERLAP_SEARCH_QUERY_ACCEL_TIME;
  static constexpr int OVERLAP_SYNC_TIME = core::profiler::ASSEMBLER_OVERLAP_SYNC_TIME;
  static constexpr int OVERLAP_CREATE_TIME = core::profiler::ASSEMBLER_OVERLAP_CREATE_TIME;
  static constexpr int OVERLAP_FILL_TIME = core::profiler::ASSEMBLER_OVERLAP_FILL_TIME;
  static constexpr int OVERLAP_CREATE_EXCHANGE_TIME = core::profiler::ASSEMBLER_OVERLAP_CREATE_EXCHANGE_TIME;
  static constexpr int OVERLAP_CREATE_AUX_TIME = core::profiler::ASSEMBLER_OVERLAP_CREATE_AUX_TIME;
  static constexpr int INFER_BOUNDARIES_TIME = core::profiler::ASSEMBLER_INFER_BOUNDARIES_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_PROJECT_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_PROJECT_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_PROJECT_CREATE_EXCHANGE_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_PROJECT_CREATE_EXCHANGE_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_PROJECT_EXCHANGE_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_PROJECT_EXCHANGE_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_PROJECT_GEN_COVER_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_PROJECT_GEN_COVER_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_PROJECT_GEN_BOUNDARY_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_PROJECT_GEN_BOUNDARY_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_DETECT_EXTERIOR_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_DETECT_EXTERIOR_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_DETECT_EXTERIOR_SEED_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_DETECT_EXTERIOR_SEED_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_DETECT_EXTERIOR_FLOOD_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_DETECT_EXTERIOR_FLOOD_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_UPDATE_AUX_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_UPDATE_AUX_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_UPDATE_AUX_GRID_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_UPDATE_AUX_GRID_TIME;
  static constexpr int CUT_BOUNDARY_HOLES_UPDATE_AUX_OVERLAP_TIME = core::profiler::ASSEMBLER_CUT_BOUNDARY_HOLES_UPDATE_AUX_OVERLAP_TIME;
  static constexpr int LOCATE_OUTER_FRINGE_TIME = core::profiler::ASSEMBLER_LOCATE_OUTER_FRINGE_TIME;
  static constexpr int OCCLUSION_TIME = core::profiler::ASSEMBLER_OCCLUSION_TIME;
  static constexpr int OCCLUSION_PAIRWISE_TIME = core::profiler::ASSEMBLER_OCCLUSION_PAIRWISE_TIME;
  static constexpr int OCCLUSION_PAD_SMOOTH_TIME = core::profiler::ASSEMBLER_OCCLUSION_PAD_SMOOTH_TIME;
  static constexpr int OCCLUSION_ACCUMULATE_TIME = core::profiler::ASSEMBLER_OCCLUSION_ACCUMULATE_TIME;
  static constexpr int MINIMIZE_OVERLAP_TIME = core::profiler::ASSEMBLER_MINIMIZE_OVERLAP_TIME;
  static constexpr int CONNECTIVITY_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_TIME;
  static constexpr int CONNECTIVITY_LOCATE_RECEIVERS_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_LOCATE_RECEIVERS_TIME;
  static constexpr int CONNECTIVITY_DONOR_EDGE_DISTANCE_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_DONOR_EDGE_DISTANCE_TIME;
  static constexpr int CONNECTIVITY_DONOR_EDGE_DISTANCE_COMPUTE_DISTANCES_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_DONOR_EDGE_DISTANCE_COMPUTE_DISTANCES_TIME;
  static constexpr int CONNECTIVITY_DONOR_EDGE_DISTANCE_CREATE_EXCHANGE_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_DONOR_EDGE_DISTANCE_CREATE_EXCHANGE_TIME;
  static constexpr int CONNECTIVITY_DONOR_EDGE_DISTANCE_EXCHANGE_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_DONOR_EDGE_DISTANCE_EXCHANGE_TIME;
  static constexpr int CONNECTIVITY_CHOOSE_DONORS_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_CHOOSE_DONORS_TIME;
  static constexpr int CONNECTIVITY_SYNC_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_SYNC_TIME;
  static constexpr int CONNECTIVITY_SYNC_CREATE_EXCHANGE_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_SYNC_CREATE_EXCHANGE_TIME;
  static constexpr int CONNECTIVITY_SYNC_EXCHANGE_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_SYNC_EXCHANGE_TIME;
  static constexpr int CONNECTIVITY_SYNC_FINALIZE_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_SYNC_FINALIZE_TIME;
  static constexpr int CONNECTIVITY_CREATE_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_CREATE_TIME;
  static constexpr int CONNECTIVITY_FILL_TIME = core::profiler::ASSEMBLER_CONNECTIVITY_FILL_TIME;

};

assembler CreateAssembler(std::shared_ptr<context> Context, assembler::params Params={});

}

#endif
