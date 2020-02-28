// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_SUPPORT_XDMF_C_H_LOADED
#define OVK_SUPPORT_XDMF_C_H_LOADED

#include <support/XDMF.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef OVK_HAVE_XDMF

struct support_xdmf_grid_meta;
typedef struct support_xdmf_grid_meta support_xdmf_grid_meta;

struct support_xdmf_attribute_meta;
typedef struct support_xdmf_attribute_meta support_xdmf_attribute_meta;

struct support_xdmf;
typedef struct support_xdmf support_xdmf;

void support_CreateXDMFGridMeta(support_xdmf_grid_meta **Grid, const char *Name, const int *Size);
void support_DestroyXDMFGridMeta(support_xdmf_grid_meta **Grid);

void support_GetXDMFGridMetaName(const support_xdmf_grid_meta *Grid, char *Name);
void support_GetXDMFGridMetaSize(const support_xdmf_grid_meta *Grid, int *Size);

void support_CreateXDMFAttributeMeta(support_xdmf_attribute_meta **Attribute, const char *Name,
  support_xdmf_attribute_type Type);
void support_DestroyXDMFAttributeMeta(support_xdmf_attribute_meta **Attribute);

void support_GetXDMFAttributeMetaName(const support_xdmf_attribute_meta *Attribute, char *Name);
void support_GetXDMFAttributeMetaType(const support_xdmf_attribute_meta *Attribute,
  support_xdmf_attribute_type *Type);

void support_CreateXDMF(support_xdmf **XDMF, const char *Path, int NumDims, MPI_Comm Comm,
  int NumGrids, support_xdmf_grid_meta **Grids, int NumAttributes, support_xdmf_attribute_meta
  **Attributes, support_xdmf_error *Error);
void support_OpenXDMF(support_xdmf **XDMF, const char *Path, MPI_Comm Comm, support_xdmf_error
  *Error);
void support_CloseXDMF(support_xdmf **XDMF);

void support_WriteXDMFGeometry(support_xdmf *XDMF, const char *GridName, int Dimension, const double
  *Data, const int *DataBegin, const int *DataEnd);
void support_WriteXDMFGeometryRange(support_xdmf *XDMF, const char *GridName, int Dimension, const
  double *Data, const int *DataBegin, const int *DataEnd, const int *WriteRangeBegin, const int
  *WriteRangeEnd);

void support_WriteXDMFAttribute(support_xdmf *XDMF, const char *GridName, const char *AttributeName,
  const void *Data, const int *DataBegin, const int *DataEnd);
void support_WriteXDMFAttributeRange(support_xdmf *XDMF, const char *GridName, const char
  *AttributeName, const void *Data, const int *DataBegin, const int *DataEnd, const int
  *WriteRangeBegin, const int *WriteRangeEnd);

void support_GetXDMFPath(const support_xdmf *XDMF, char *Path);
int support_GetXDMFDimension(const support_xdmf *XDMF);
void support_GetXDMFComm(const support_xdmf *XDMF, MPI_Comm *Comm);
int support_GetXDMFGridCount(const support_xdmf *XDMF);
void support_GetXDMFGrid(const support_xdmf *XDMF, int Index, const support_xdmf_grid_meta **Grid);
int support_GetXDMFAttributeCount(const support_xdmf *XDMF);
void support_GetXDMFAttribute(const support_xdmf *XDMF, int Index, const support_xdmf_attribute_meta
  **Attribute);

#endif

#ifdef __cplusplus
}
#endif

#endif
