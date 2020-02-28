// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DATA_TYPE_OPS_HPP_INCLUDED
#define OVK_CORE_DATA_TYPE_OPS_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Handle.hpp>

#include <mpi.h>

#include <cstdint>
#include <initializer_list>

namespace ovk {
namespace core {

inline handle<MPI_Datatype> CreateMPIContiguousType(int Count, MPI_Datatype ElementMPIType) {

  MPI_Datatype MPITypeRaw;
  MPI_Type_contiguous(Count, ElementMPIType, &MPITypeRaw);

  return {MPITypeRaw, &MPI_Type_free};

}

struct mpi_struct_block {
  std::ptrdiff_t Offset;
  int Count;
  MPI_Datatype Type;
};

inline handle<MPI_Datatype> CreateMPIStructType(std::size_t Size, std::initializer_list<
  mpi_struct_block> Blocks) {

  int NumBlocks = int(Blocks.size());

  array<MPI_Aint> BlockOffsets({NumBlocks});
  array<int> BlockCounts({NumBlocks});
  array<MPI_Datatype> BlockTypes({NumBlocks});

  int iBlock = 0;
  for (auto &Block : Blocks) {
    BlockOffsets(iBlock) = MPI_Aint(Block.Offset);
    BlockCounts(iBlock) = Block.Count;
    BlockTypes(iBlock) = Block.Type;
    ++iBlock;
  }

  MPI_Datatype BlocksMPITypeRaw;
  MPI_Type_create_struct(NumBlocks, BlockCounts.Data(), BlockOffsets.Data(), BlockTypes.Data(),
    &BlocksMPITypeRaw);
  handle<MPI_Datatype> BlocksMPIType(BlocksMPITypeRaw, &MPI_Type_free);

  MPI_Datatype MPITypeRaw;
  MPI_Type_create_resized(BlocksMPIType, 0, Size, &MPITypeRaw);

  return {MPITypeRaw, &MPI_Type_free};

}

template <typename T1, typename T2> long long GetByteOffset(const T1 &Left, const T2 &Right) {
  return reinterpret_cast<const byte *>(&Right) - reinterpret_cast<const byte *>(&Left);
}

}}

#endif
