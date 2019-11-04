// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
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

struct mpi_struct_block {
  std::ptrdiff_t Offset;
  int Count;
  MPI_Datatype Type;
};

inline handle<MPI_Datatype> CreateStructMPIDatatype(std::size_t Size, std::initializer_list<
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

  MPI_Datatype ElementsTypeRaw;
  MPI_Type_create_struct(NumBlocks, BlockCounts.Data(), BlockOffsets.Data(), BlockTypes.Data(),
    &ElementsTypeRaw);
  handle<MPI_Datatype> ElementsType(ElementsTypeRaw, &MPI_Type_free);

  MPI_Datatype StructTypeRaw;
  MPI_Type_create_resized(ElementsType, 0, Size, &StructTypeRaw);
  MPI_Type_commit(&StructTypeRaw);

  return {StructTypeRaw, &MPI_Type_free};

}

}}

#endif
