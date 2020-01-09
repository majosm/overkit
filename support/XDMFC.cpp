// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "support/XDMFC.h"

#include "support/XDMF.hpp"

#include <ovk/core/Array.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Global.hpp>

#include <algorithm>
#include <cstring>
#include <string>
#include <utility>

#ifdef __cplusplus
extern "C" {
#endif

void support_CreateXDMFGridMeta(support_xdmf_grid_meta **Grid, const char *Name, const int *Size) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto GridCPPPtr = new support::xdmf_grid_meta(Name, Size);
  *Grid = reinterpret_cast<support_xdmf_grid_meta *>(GridCPPPtr);

}

void support_DestroyXDMFGridMeta(support_xdmf_grid_meta **Grid) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(*Grid, "Invalid grid pointer.");

  auto GridCPPPtr = reinterpret_cast<support::xdmf_grid_meta *>(*Grid);

  delete GridCPPPtr;

  *Grid = nullptr;

}

void support_GetXDMFGridMetaName(const support_xdmf_grid_meta *Grid, char *Name) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &GridCPP = *reinterpret_cast<const support::xdmf_grid_meta *>(Grid);
  std::strcpy(Name, GridCPP.Name().c_str());

}

void support_GetXDMFGridMetaSize(const support_xdmf_grid_meta *Grid, int *Size) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Size, "Invalid size pointer.");

  auto &GridCPP = *reinterpret_cast<const support::xdmf_grid_meta *>(Grid);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    Size[iDim] = GridCPP.Size()(iDim);
  }

}

void support_CreateXDMFAttributeMeta(support_xdmf_attribute_meta **Attribute, const char *Name,
  support_xdmf_attribute_type Type) {

  OVK_DEBUG_ASSERT(Attribute, "Invalid attribute pointer.");

  auto AttributeCPPPtr = new support::xdmf_attribute_meta(Name, support::xdmf_attribute_type(Type));
  *Attribute = reinterpret_cast<support_xdmf_attribute_meta *>(AttributeCPPPtr);

}

void support_DestroyXDMFAttributeMeta(support_xdmf_attribute_meta **Attribute) {

  OVK_DEBUG_ASSERT(Attribute, "Invalid attribute pointer.");
  OVK_DEBUG_ASSERT(*Attribute, "Invalid attribute pointer.");

  auto AttributeCPPPtr = reinterpret_cast<support::xdmf_attribute_meta *>(*Attribute);

  delete AttributeCPPPtr;

  *Attribute = nullptr;

}

void support_GetXDMFAttributeMetaName(const support_xdmf_attribute_meta *Attribute, char *Name) {

  OVK_DEBUG_ASSERT(Attribute, "Invalid attribute pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &AttributeCPP = *reinterpret_cast<const support::xdmf_attribute_meta *>(Attribute);
  std::strcpy(Name, AttributeCPP.Name().c_str());

}

void support_GetXDMFAttributeMetaType(const support_xdmf_attribute_meta *Attribute,
  support_xdmf_attribute_type *Type) {

  OVK_DEBUG_ASSERT(Attribute, "Invalid attribute pointer.");
  OVK_DEBUG_ASSERT(Type, "Invalid type pointer.");

  auto &AttributeCPP = *reinterpret_cast<const support::xdmf_attribute_meta *>(Attribute);
  *Type = support_xdmf_attribute_type(AttributeCPP.Type());

}

void support_CreateXDMF(support_xdmf **XDMF, const char *Path, int NumDims, MPI_Comm Comm,
  int NumGrids, support_xdmf_grid_meta **Grids, int NumAttributes, support_xdmf_attribute_meta
  **Attributes, support_xdmf_error *Error) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");
  OVK_DEBUG_ASSERT(Grids, "Invalid grids pointer.");
  OVK_DEBUG_ASSERT(Attributes, "Invalid attributes pointer.");
  OVK_DEBUG_ASSERT(Error, "Invalid error pointer.");

  ovk::array<support::xdmf_grid_meta> GridsCPP;
  GridsCPP.Reserve(NumGrids);

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    auto GridCPPPtr = reinterpret_cast<support::xdmf_grid_meta *>(Grids[iGrid]);
    GridsCPP.Append(std::move(*GridCPPPtr));
    delete GridCPPPtr;
    Grids[iGrid] = nullptr;
  }

  ovk::array<support::xdmf_attribute_meta> AttributesCPP;
  AttributesCPP.Reserve(NumAttributes);

  for (int iAttribute = 0; iAttribute < NumAttributes; ++iAttribute) {
    auto AttributeCPPPtr = reinterpret_cast<support::xdmf_attribute_meta *>(Attributes[iAttribute]);
    AttributesCPP.Append(std::move(*AttributeCPPPtr));
    delete AttributeCPPPtr;
    Attributes[iAttribute] = nullptr;
  }

  support::captured_xdmf_error ErrorCPP;
  auto MaybeXDMFCPP = support::CreateXDMF(Path, NumDims, Comm, std::move(GridsCPP),
    std::move(AttributesCPP), ErrorCPP);

  if (!ErrorCPP) {
    auto XDMFCPPPtr = new support::xdmf(MaybeXDMFCPP.Release());
    *XDMF = reinterpret_cast<support_xdmf *>(XDMFCPPPtr);
  } else {
    *XDMF = nullptr;
  }

  *Error = support_xdmf_error(ErrorCPP.Code());

}

void support_OpenXDMF(support_xdmf **XDMF, const char *Path, MPI_Comm Comm, support_xdmf_error
  *Error) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");
  OVK_DEBUG_ASSERT(Path, "Invalid path pointer.");
  OVK_DEBUG_ASSERT(Error, "Invalid error pointer.");

  support::captured_xdmf_error ErrorCPP;
  auto MaybeXDMFCPP = support::OpenXDMF(Path, Comm, ErrorCPP);

  if (!ErrorCPP) {
    auto XDMFCPPPtr = new support::xdmf(MaybeXDMFCPP.Release());
    *XDMF = reinterpret_cast<support_xdmf *>(XDMFCPPPtr);
  } else {
    *XDMF = nullptr;
  }

  *Error = support_xdmf_error(ErrorCPP.Code());

}

void support_CloseXDMF(support_xdmf **XDMF) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");
  OVK_DEBUG_ASSERT(*XDMF, "Invalid xdmf pointer.");

  auto XDMFCPPPtr = reinterpret_cast<support::xdmf *>(*XDMF);

  delete XDMFCPPPtr;

  *XDMF = nullptr;

}

void support_WriteXDMFGeometry(support_xdmf *XDMF, const char *GridName, int Dimension, const double
  *Data, const int *DataBegin, const int *DataEnd) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");

  auto &XDMFCPP = *reinterpret_cast<support::xdmf *>(XDMF);
  XDMFCPP.WriteGeometry(GridName, Dimension, {Data, {DataBegin, DataEnd}});

}

void support_WriteXDMFGeometryRange(support_xdmf *XDMF, const char *GridName, int Dimension, const
  double *Data, const int *DataBegin, const int *DataEnd, const int *WriteRangeBegin, const int
  *WriteRangeEnd) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");

  auto &XDMFCPP = *reinterpret_cast<support::xdmf *>(XDMF);
  XDMFCPP.WriteGeometry(GridName, Dimension, {Data, {DataBegin, DataEnd}}, {WriteRangeBegin,
    WriteRangeEnd});

}

void support_WriteXDMFAttribute(support_xdmf *XDMF, const char *GridName, const char *AttributeName,
  const void *Data, const int *DataBegin, const int *DataEnd) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");

  auto &XDMFCPP = *reinterpret_cast<support::xdmf *>(XDMF);

  const ovk::array<support::xdmf_attribute_meta> &Attributes = XDMFCPP.Attributes();

  auto Iter = std::find_if(Attributes.Begin(), Attributes.End(), [AttributeName](const
    support::xdmf_attribute_meta &Attribute) -> bool {
    return Attribute.Name() == AttributeName;
  });

  OVK_DEBUG_ASSERT(Iter != Attributes.End(), "Invalid attribute '%s'.\n", AttributeName);

  const support::xdmf_attribute_meta &Attribute = *Iter;

  switch (Attribute.Type()) {
  case support::xdmf_attribute_type::INT:
    XDMFCPP.WriteAttribute(GridName, AttributeName, {static_cast<const int *>(Data),
      {DataBegin, DataEnd}});
    break;
  case support::xdmf_attribute_type::LONG_LONG:
    XDMFCPP.WriteAttribute(GridName, AttributeName, {static_cast<const long long *>(Data),
      {DataBegin, DataEnd}});
    break;
  case support::xdmf_attribute_type::DOUBLE:
    XDMFCPP.WriteAttribute(GridName, AttributeName, {static_cast<const double *>(Data),
      {DataBegin, DataEnd}});
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
  }

}

void support_WriteXDMFAttributeRange(support_xdmf *XDMF, const char *GridName, const char
  *AttributeName, const void *Data, const int *DataBegin, const int *DataEnd, const int
  *WriteRangeBegin, const int *WriteRangeEnd) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");

  auto &XDMFCPP = *reinterpret_cast<support::xdmf *>(XDMF);

  const ovk::array<support::xdmf_attribute_meta> &Attributes = XDMFCPP.Attributes();

  auto Iter = std::find_if(Attributes.Begin(), Attributes.End(), [AttributeName](const
    support::xdmf_attribute_meta &Attribute) -> bool {
    return Attribute.Name() == AttributeName;
  });

  OVK_DEBUG_ASSERT(Iter != Attributes.End(), "Invalid attribute '%s'\n", AttributeName);

  const support::xdmf_attribute_meta &Attribute = *Iter;

  switch (Attribute.Type()) {
  case support::xdmf_attribute_type::INT:
    XDMFCPP.WriteAttribute(GridName, AttributeName, {static_cast<const int *>(Data),
      {DataBegin, DataEnd}}, {WriteRangeBegin, WriteRangeEnd});
    break;
  case support::xdmf_attribute_type::LONG_LONG:
    XDMFCPP.WriteAttribute(GridName, AttributeName, {static_cast<const long long *>(Data),
      {DataBegin, DataEnd}}, {WriteRangeBegin, WriteRangeEnd});
    break;
  case support::xdmf_attribute_type::DOUBLE:
    XDMFCPP.WriteAttribute(GridName, AttributeName, {static_cast<const double *>(Data),
      {DataBegin, DataEnd}}, {WriteRangeBegin, WriteRangeEnd});
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
  }

}

void support_GetXDMFPath(const support_xdmf *XDMF, char *Path) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");
  OVK_DEBUG_ASSERT(Path, "Invalid path pointer.");

  auto &XDMFCPP = *reinterpret_cast<const support::xdmf *>(XDMF);
  std::strcpy(Path, XDMFCPP.Path().c_str());

}

int support_GetXDMFDimension(const support_xdmf *XDMF) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");

  auto &XDMFCPP = *reinterpret_cast<const support::xdmf *>(XDMF);
  return XDMFCPP.Dimension();

}

void support_GetXDMFComm(const support_xdmf *XDMF, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &XDMFCPP = *reinterpret_cast<const support::xdmf *>(XDMF);
  *Comm = XDMFCPP.Comm();

}

int support_GetXDMFGridCount(const support_xdmf *XDMF) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");

  auto &XDMFCPP = *reinterpret_cast<const support::xdmf *>(XDMF);
  return XDMFCPP.Grids().Count();


}

void support_GetXDMFGrid(const support_xdmf *XDMF, int Index, const support_xdmf_grid_meta **Grid) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &XDMFCPP = *reinterpret_cast<const support::xdmf *>(XDMF);
  *Grid = reinterpret_cast<const support_xdmf_grid_meta *>(&XDMFCPP.Grids()(Index));

}

int support_GetXDMFAttributeCount(const support_xdmf *XDMF) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");

  auto &XDMFCPP = *reinterpret_cast<const support::xdmf *>(XDMF);
  return XDMFCPP.Attributes().Count();

}

void support_GetXDMFAttribute(const support_xdmf *XDMF, int Index, const support_xdmf_attribute_meta
  **Attribute) {

  OVK_DEBUG_ASSERT(XDMF, "Invalid xdmf pointer.");
  OVK_DEBUG_ASSERT(Attribute, "Invalid attribute pointer.");

  auto &XDMFCPP = *reinterpret_cast<const support::xdmf *>(XDMF);
  *Attribute = reinterpret_cast<const support_xdmf_attribute_meta *>(&XDMFCPP.Attributes()(Index));

}

#ifdef __cplusplus
}
#endif
