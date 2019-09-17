// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_OVERLAP_COMPONENT_H_INCLUDED
#define OVK_CORE_C_OVERLAP_COMPONENT_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/Global.h>
#include <ovk/core-c/OverlapM.h>
#include <ovk/core-c/OverlapN.h>
#include <ovk/core/OverlapComponent.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_overlap_component;
typedef struct ovk_overlap_component ovk_overlap_component;

struct ovk_overlap_component_params;
typedef struct ovk_overlap_component_params ovk_overlap_component_params;

int ovkOverlapCount(const ovk_overlap_component *OverlapComponent);

bool ovkOverlapExists(const ovk_overlap_component *OverlapComponent, int MGridID, int NGridID);

void ovkCreateOverlap(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID);
void ovkCreateOverlaps(ovk_overlap_component *OverlapComponent, int Count, const int *MGridIDs,
  const int *NGridIDs);

void ovkDestroyOverlap(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID);
void ovkDestroyOverlaps(ovk_overlap_component *OverlapComponent, int Count, const int *MGridIDs,
  const int *NGridIDs);

int ovkLocalOverlapMCount(const ovk_overlap_component *OverlapComponent);
int ovkLocalOverlapMCountForGrid(const ovk_overlap_component *OverlapComponent, int MGridID);

void ovkGetOverlapM(const ovk_overlap_component *OverlapComponent, int MGridID, int NGridID, const
  ovk_overlap_m **OverlapM);
bool ovkEditingOverlapM(const ovk_overlap_component *OverlapComponent, int MGridID, int NGridID);
void ovkEditOverlapM(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID,
  ovk_overlap_m **OverlapM);
void ovkRestoreOverlapM(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID,
  ovk_overlap_m **OverlapM);

int ovkLocalOverlapNCount(const ovk_overlap_component *OverlapComponent);
int ovkLocalOverlapNCountForGrid(const ovk_overlap_component *OverlapComponent, int NGridID);

void ovkGetOverlapN(const ovk_overlap_component *OverlapComponent, int MGridID, int NGridID, const
  ovk_overlap_n **OverlapN);
bool ovkEditingOverlapN(const ovk_overlap_component *OverlapComponent, int MGridID, int NGridID);
void ovkEditOverlapN(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID,
  ovk_overlap_n **OverlapN);
void ovkRestoreOverlapN(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID,
  ovk_overlap_n **OverlapN);

void ovkCreateOverlapComponentParams(ovk_overlap_component_params **Params);
void ovkDestroyOverlapComponentParams(ovk_overlap_component_params **Params);
void ovkGetOverlapComponentParamName(const ovk_overlap_component_params *Params, char *Name);
void ovkSetOverlapComponentParamName(ovk_overlap_component_params *Params, const char *Name);

#ifdef __cplusplus
}
#endif

#endif
