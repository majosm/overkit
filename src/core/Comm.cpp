// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Comm.hpp"

#include <ovk/core/Global.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {

comm DuplicateComm(comm_view Comm) {

  MPI_Comm DuplicatedCommRaw;
  MPI_Comm_dup(Comm.Get(), &DuplicatedCommRaw);

  return comm(DuplicatedCommRaw);

}

comm CreateSubsetComm(comm_view Comm, bool InSubset) {

  MPI_Comm CommRaw = Comm.Get();

  MPI_Comm SubsetCommRaw;
  MPI_Comm_split(CommRaw, InSubset ? 0 : MPI_UNDEFINED, Comm.Rank(), &SubsetCommRaw);

  return comm(SubsetCommRaw);

}

comm CreateCartComm(comm_view Comm, int NumDims, const tuple<int> &Dims_, const tuple<bool>
  &Periodic, bool AllowReorder) {

  tuple<int> Dims = Dims_;
  MPI_Dims_create(Comm.Size(), NumDims, Dims.Data());

  tuple<int> PeriodicInt = tuple<int>(Periodic);

  MPI_Comm CartCommRaw;
  MPI_Cart_create(Comm, NumDims, Dims.Data(), PeriodicInt.Data(), AllowReorder, &CartCommRaw);

  return comm(CartCommRaw);

}

bool IsCartComm(comm_view Comm) {

  int TopologyType;
  MPI_Topo_test(Comm, &TopologyType);

  return TopologyType == MPI_CART;

}

int GetCartCommDimension(comm_view Comm) {

  int NumDims;
  MPI_Cartdim_get(Comm, &NumDims);

  return NumDims;

}

tuple<int> GetCartCommDims(comm_view Comm) {

  int NumDims = GetCartCommDimension(Comm);

  tuple<int> CartDims = MakeUniformTuple<int>(NumDims, 1, 1);
  tuple<int> PeriodicInt;
  tuple<int> CartCoords;
  MPI_Cart_get(Comm, NumDims, CartDims.Data(), PeriodicInt.Data(), CartCoords.Data());

  return CartDims;

}

tuple<bool> GetCartCommPeriodic(comm_view Comm) {

  int NumDims = GetCartCommDimension(Comm);

  tuple<int> CartDims;
  tuple<int> PeriodicInt = MakeUniformTuple<int>(NumDims, 0);
  tuple<int> CartCoords;
  MPI_Cart_get(Comm, NumDims, CartDims.Data(), PeriodicInt.Data(), CartCoords.Data());

  return tuple<bool>(PeriodicInt);

}

tuple<int> GetCartCommCoords(comm_view Comm) {

  int NumDims = GetCartCommDimension(Comm);

  tuple<int> CartDims;
  tuple<int> PeriodicInt;
  tuple<int> CartCoords = MakeUniformTuple<int>(NumDims, 0);
  MPI_Cart_get(Comm, NumDims, CartDims.Data(), PeriodicInt.Data(), CartCoords.Data());

  return CartCoords;

}

}
