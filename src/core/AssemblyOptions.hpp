// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ASSEMBLY_OPTIONS_HPP_INCLUDED
#define OVK_CORE_ASSEMBLY_OPTIONS_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>

#include <map>

namespace ovk {

struct assembly_options {
  int NumDims_;
  int NumGrids_;
  std::map<int, std::map<int, bool>> Overlappable_;
  std::map<int, std::map<int, double>> OverlapTolerance_;
  std::map<int, double> OverlapAccelDepthAdjust_;
  std::map<int, double> OverlapAccelResolutionAdjust_;
  std::map<int, bool> InferBoundaries_;
  std::map<int, std::map<int, bool>> CutBoundaryHoles_;
  std::map<int, std::map<int, occludes>> Occludes_;
  std::map<int, std::map<int, int>> EdgePadding_;
  std::map<int, int> EdgeSmoothing_;
  std::map<int, std::map<int, connection_type>> ConnectionType_;
  std::map<int, int> FringeSize_;
  std::map<int, std::map<int, bool>> MinimizeOverlap_;
};

void CreateAssemblyOptions(assembly_options &Options, int NumDims, int NumGrids, int *GridIDs);
void DestroyAssemblyOptions(assembly_options &Options);

void GetAssemblyOptionsDimension(const assembly_options &Options, int &NumDims);
void GetAssemblyOptionsGridCount(const assembly_options &Options, int &NumGrids);

void GetAssemblyOptionOverlappable(const assembly_options &Options, int OverlappingGridID,
  int OverlappedGridID, bool &Overlappable);
void SetAssemblyOptionOverlappable(assembly_options &Options, int OverlappingGridID,
  int OverlappedGridID, bool Overlappable);

void GetAssemblyOptionOverlapTolerance(const assembly_options &Options, int OverlappingGridID,
  int OverlappedGridID, double &OverlapTolerance);
void SetAssemblyOptionOverlapTolerance(assembly_options &Options, int OverlappingGridID,
  int OverlappedGridID, double OverlapTolerance);

void GetAssemblyOptionOverlapAccelDepthAdjust(const assembly_options &Options,
  int OverlappingGridID, double &OverlapAccelDepthAdjust);
void SetAssemblyOptionOverlapAccelDepthAdjust(assembly_options &Options,
  int OverlappingGridID, double OverlapAccelDepthAdjust);

void GetAssemblyOptionOverlapAccelResolutionAdjust(const assembly_options &Options,
  int OverlappingGridID, double &OverlapAccelResolutionAdjust);
void SetAssemblyOptionOverlapAccelResolutionAdjust(assembly_options &Options,
  int OverlappingGridID, double OverlapAccelResolutionAdjust);

void GetAssemblyOptionInferBoundaries(const assembly_options &Options, int GridID,
  bool &InferBoundaries);
void SetAssemblyOptionInferBoundaries(assembly_options &Options, int GridID,
  bool InferBoundaries);

void GetAssemblyOptionCutBoundaryHoles(const assembly_options &Options, int CuttingGridID,
  int CutGridID, bool &CutBoundaryHoles);
void SetAssemblyOptionCutBoundaryHoles(assembly_options &Options, int CuttingGridID,
  int CutGridID, bool CutBoundaryHoles);

void GetAssemblyOptionOccludes(const assembly_options &Options, int OccludingGridID,
  int OccludedGridID, occludes &Occludes);
void SetAssemblyOptionOccludes(assembly_options &Options, int OccludingGridID,
  int OccludedGridID, occludes Occludes);

void GetAssemblyOptionEdgePadding(const assembly_options &Options, int OccludingGridID,
  int OccludedGridID, int &EdgePadding);
void SetAssemblyOptionEdgePadding(assembly_options &Options, int OccludingGridID,
  int OccludedGridID, int EdgePadding);

void GetAssemblyOptionEdgeSmoothing(const assembly_options &Options, int OccludedGridID,
  int &EdgeSmoothing);
void SetAssemblyOptionEdgeSmoothing(assembly_options &Options, int OccludedGridID,
  int EdgeSmoothing);

void GetAssemblyOptionConnectionType(const assembly_options &Options, int DonorGridID,
  int ReceiverGridID, connection_type &ConnectionType);
void SetAssemblyOptionConnectionType(assembly_options &Options, int DonorGridID,
  int ReceiverGridID, connection_type ConnectionType);

void GetAssemblyOptionFringeSize(const assembly_options &Options, int GridID, int &FringeSize);
void SetAssemblyOptionFringeSize(assembly_options &Options, int GridID, int FringeSize);

void GetAssemblyOptionMinimizeOverlap(const assembly_options &Options, int DonorGridID,
  int ReceiverGridID, bool &MinimizeOverlap);
void SetAssemblyOptionMinimizeOverlap(assembly_options &Options, int DonorGridID,
  int ReceiverGridID, bool MinimizeOverlap);

}

#endif
