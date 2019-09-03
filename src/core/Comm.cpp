// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Comm.hpp"

#include <ovk/core/Array.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <cmath>
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

namespace core {

tuple<int> CreateCartesianDecompDims(int Size, int NumDims, const tuple<int> &InputDims) {

  tuple<int> Dims = InputDims;

  MPI_Dims_create(Size, NumDims, Dims.Data());

  // MPI_Dims_create doesn't always give the desired decompositions (e.g., OpenMPI returns (N*N,1)
  // instead of (N,N) when N is a prime, and MVAPICH returns (2*M,N) instead of (M,2*N) when M is
  // a prime and N is small). Below is a naive attempt to improve the results.

  constexpr int SmallPrimes[] = {
      2,   3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,  61,  67,
     71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163,
    167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
    271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383,
    389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
    503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619,
    631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751,
    757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881,
    883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997
  };

  if (OVK_DEBUG) {
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      int ReducedDim = Dims(iDim);
      for (int Prime : SmallPrimes) {
        while (ReducedDim % Prime == 0) {
          ReducedDim /= Prime;
        }
        if (ReducedDim == 1) break;
      }
      OVK_DEBUG_ASSERT(ReducedDim == 1, "Cartesian communicator dimension has prime factors larger "
        "than those defined; add more to list.");
    }
  }

  auto DimsVariance = [NumDims](const tuple<int> &Dims) -> int {
    int Mean = 0;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      Mean += Dims(iDim);
    }
    Mean = std::round(double(Mean)/double(NumDims));
    int Variance = 0;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      int Deviation = Dims(iDim) - Mean;
      Variance += Deviation*Deviation;
    }
    return Variance;
  };

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    if (InputDims(iDim) != 0) continue;
    array<int> PrimeFactors;
    PrimeFactors.Reserve(32);
    int ReducedDim = Dims(iDim);
    for (int Prime : SmallPrimes) {
      if (Prime > ReducedDim) break;
      while (ReducedDim % Prime == 0) {
        PrimeFactors.Append(Prime);
        ReducedDim /= Prime;
      }
    }
    for (int Factor : PrimeFactors) {
      tuple<int> BestDims = Dims;
      int BestVariance = DimsVariance(Dims);
      for (int iCyclicDim = iDim+1; iCyclicDim < iDim+NumDims; ++iCyclicDim) {
        int iOtherDim = iCyclicDim % NumDims;
        if (InputDims(iOtherDim) != 0) continue;
        tuple<int> NewDims = Dims;
        NewDims(iDim) /= Factor;
        NewDims(iOtherDim) *= Factor;
        int NewVariance = DimsVariance(NewDims);
        if (NewVariance < BestVariance) {
          BestDims = NewDims;
          BestVariance = NewVariance;
        }
      }
      Dims = BestDims;
    }
  }

  if (InputDims(0) == 0 && InputDims(1) == 0 && Dims(0) < Dims(1)) {
    std::swap(Dims(0), Dims(1));
  }
  if (InputDims(0) == 0 && InputDims(2) == 0 && Dims(0) < Dims(2)) {
    std::swap(Dims(0), Dims(2));
  }
  if (InputDims(1) == 0 && InputDims(2) == 0 && Dims(1) < Dims(2)) {
    std::swap(Dims(1), Dims(2));
  }

  return Dims;

}

}

comm CreateCartComm(comm_view Comm, int NumDims, const tuple<int> &Dims_, const tuple<bool>
  &Periodic, bool AllowReorder) {

  tuple<int> Dims = core::CreateCartesianDecompDims(Comm.Size(), NumDims, Dims_);

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
