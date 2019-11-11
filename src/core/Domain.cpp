// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Domain.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Component.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Optional.hpp"
#include "ovk/core/Set.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <memory>
#include <string>
#include <utility>

namespace ovk {

namespace core {

namespace domain_internal {

domain_base_1::domain_base_1(std::shared_ptr<context> &&Context, std::string &&Name, MPI_Comm Comm):
  Context_(std::move(Context)),
  Name_(std::move(Name)),
  Comm_(DuplicateComm(Comm))
{
  core::logger &Logger = Context_->core_Logger();
  Logger.LogStatus(Comm_.Rank() == 0, "Creating domain %s...", *Name_);
  Level1_ = Logger.IncreaseStatusLevelAndIndent();
}

domain_base_1::~domain_base_1() noexcept {

  if (Context_) {
    MPI_Barrier(Comm_);
    core::logger &Logger = Context_->core_Logger();
    Level1_.Reset();
    Logger.LogStatus(Comm_.Rank() == 0, "Done destroying domain %s.", *Name_);
  }

}

}

domain_base::domain_base(std::shared_ptr<context> &&Context, std::string &&Name, int NumDims,
  MPI_Comm Comm):
  domain_base_1(std::move(Context), std::move(Name), Comm),
  NumDims_(NumDims)
{}

bool domain_base::GridExists(int GridID) const {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");

  return GridRecords_.Contains(GridID);

}

grid::params domain_base::MakeGridParams() const {

  return grid::params().SetDimension(NumDims_);

}

void domain_base::CreateGrid(int GridID, optional<grid::params> MaybeParams) {

  MPI_Barrier(Comm_);

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(!GridExists(GridID), "Grid %i already exists.", GridID);

  core::logger &Logger = Context_->core_Logger();

  bool IsLocal = MaybeParams.Present();

  if (OVK_DEBUG) {
    int IsLocalInt = IsLocal ? 1 : 0;
    int AtLeastOneLocal;
    MPI_Allreduce(&IsLocalInt, &AtLeastOneLocal, 1, MPI_INT, MPI_LOR, Comm_);
    OVK_DEBUG_ASSERT(AtLeastOneLocal, "Grid must be local to at least one rank.");
  }

  if (Logger.LoggingStatus()) {
    bool IsRoot = false; 
    if (IsLocal) {
      int Rank;
      MPI_Comm_rank(MaybeParams.Get().Comm(), &Rank);
      IsRoot = Rank == 0;
    }
    std::string GridName;
    if (IsRoot) GridName = MaybeParams.Get().Name();
    core::BroadcastStringAnySource(GridName, IsRoot, Comm_);
    Logger.LogStatus(Comm_.Rank() == 0, "Creating grid %s.%s...", *Name_, GridName);
  }
  auto Level1 = Logger.IncreaseStatusLevelAndIndent();

  grid *MaybeGrid = nullptr;

  if (IsLocal) {
    grid Grid = core::CreateGrid(Context_, MaybeParams.Release());
    MaybeGrid = &LocalGrids_.Insert(GridID, std::move(Grid));
  }

  grid_info GridInfo = core::CreateGridInfo(MaybeGrid, Comm_);
  GridRecords_.Insert(GridID, std::move(GridInfo));

  MPI_Barrier(Comm_);

  Level1.Reset();
  grid_record &Record = GridRecords_(GridID);
  Logger.LogStatus(Comm_.Rank() == 0, "Done creating grid %s.%s.", *Name_, Record.GridInfo.Name());

  GridEvent_.Trigger(GridID, grid_event_flags::CREATE, true);

  MPI_Barrier(Comm_);

}

void domain_base::CreateGrids(array_view<const int> GridIDs, array<optional<grid::params>>
  MaybeParams) {

  MPI_Barrier(Comm_);

  int NumCreates = GridIDs.Count();

  OVK_DEBUG_ASSERT(MaybeParams.Count() == NumCreates, "Incorrect params array size.");

  core::logger &Logger = Context_->core_Logger();

  array<bool> IsLocal({NumCreates}, false);
  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    IsLocal(iCreate) = MaybeParams(iCreate).Present();
  }

  if (OVK_DEBUG) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      int GridID = GridIDs(iCreate);
      OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
      OVK_DEBUG_ASSERT(!GridExists(GridID), "Grid %i already exists.", GridID);
      int IsLocalInt = IsLocal(iCreate) ? 1 : 0;
      int AtLeastOneLocal;
      MPI_Allreduce(&IsLocalInt, &AtLeastOneLocal, 1, MPI_INT, MPI_LOR, Comm_);
      OVK_DEBUG_ASSERT(AtLeastOneLocal, "Grid must be local to at least one rank.");
    }
  }

  if (Logger.LoggingStatus()) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      bool IsRoot = false; 
      if (IsLocal(iCreate)) {
        int Rank;
        MPI_Comm_rank(MaybeParams(iCreate).Get().Comm(), &Rank);
        IsRoot = Rank == 0;
      }
      std::string GridName;
      if (IsRoot) GridName = MaybeParams(iCreate).Get().Name();
      core::BroadcastStringAnySource(GridName, IsRoot, Comm_);
      Logger.LogStatus(Comm_.Rank() == 0, "Creating grid %s.%s...", *Name_, GridName);
    }
  }
  auto Level1 = Logger.IncreaseStatusLevelAndIndent();

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int GridID = GridIDs(iCreate);
    if (IsLocal(iCreate)) {
      grid Grid = core::CreateGrid(Context_, MaybeParams(iCreate).Release());
      LocalGrids_.Insert(GridID, std::move(Grid));
    }
  }

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int GridID = GridIDs(iCreate);
    grid *MaybeGrid = nullptr;
    if (IsLocal(iCreate)) {
      MaybeGrid = &LocalGrids_(GridID);
    }
    grid_info GridInfo = core::CreateGridInfo(MaybeGrid, Comm_);
    GridRecords_.Insert(GridID, std::move(GridInfo));
  }

  MPI_Barrier(Comm_);

  Level1.Reset();
  if (Logger.LoggingStatus()) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      int GridID = GridIDs(iCreate);
      grid_record &Record = GridRecords_(GridID);
      Logger.LogStatus(Comm_.Rank() == 0, "Done creating grid %s.%s.", *Name_,
        Record.GridInfo.Name());
    }
  }

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int GridID = GridIDs(iCreate);
    GridEvent_.Trigger(GridID, grid_event_flags::CREATE, iCreate == NumCreates-1);
  }

  MPI_Barrier(Comm_);

}

void domain_base::DestroyGrid(int GridID) {

  MPI_Barrier(Comm_);

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridExists(GridID), "Grid %i does not exist.", GridID);

  GridEvent_.Trigger(GridID, grid_event_flags::DESTROY, true);

  MPI_Barrier(Comm_);

  core::logger &Logger = Context_->core_Logger();

  std::string GridName;
  if (Logger.LoggingStatus()) {
    grid_record &Record = GridRecords_(GridID);
    GridName = Record.GridInfo.Name();
    Logger.LogStatus(Comm_.Rank() == 0, "Destroying grid %s.%s...", *Name_, GridName);
  }
  auto Level1 = Logger.IncreaseStatusLevelAndIndent();

  LocalGrids_.Erase(GridID);

  GridRecords_.Erase(GridID);

  MPI_Barrier(Comm_);

  Level1.Reset();
  if (Logger.LoggingStatus()) {
    Logger.LogStatus(Comm_.Rank() == 0, "Done destroying grid %s.%s.", *Name_, GridName);
  }

}

void domain_base::DestroyGrids(array_view<const int> GridIDs) {

  MPI_Barrier(Comm_);

  int NumDestroys = GridIDs.Count();

  if (OVK_DEBUG) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int GridID = GridIDs(iDestroy);
      OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
      OVK_DEBUG_ASSERT(GridExists(GridID), "Grid %i does not exist.", GridID);
    }
  }

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int GridID = GridIDs(iDestroy);
    GridEvent_.Trigger(GridID, grid_event_flags::DESTROY, iDestroy == NumDestroys-1);
  }

  MPI_Barrier(Comm_);

  core::logger &Logger = Context_->core_Logger();

  array<std::string> GridNames;
  if (Logger.LoggingStatus()) {
    GridNames.Resize({NumDestroys});
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int GridID = GridIDs(iDestroy);
      grid_record &Record = GridRecords_(GridID);
      GridNames(iDestroy) = Record.GridInfo.Name();
      Logger.LogStatus(Comm_.Rank() == 0, "Destroying grid %s.%s...", *Name_, GridNames(iDestroy));
    }
  }
  auto Level1 = Logger.IncreaseStatusLevelAndIndent();

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int GridID = GridIDs(iDestroy);
    LocalGrids_.Erase(GridID);
  }

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int GridID = GridIDs(iDestroy);
    GridRecords_.Erase(GridID);
  }

  MPI_Barrier(Comm_);

  Level1.Reset();
  if (Logger.LoggingStatus()) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      Logger.LogStatus(Comm_.Rank() == 0, "Done destroying grid %s.%s.", *Name_,
        GridNames(iDestroy));
    }
  }

}

const grid_info &domain_base::GridInfo(int GridID) const {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridExists(GridID), "Grid %i does not exist.", GridID);

  return GridRecords_(GridID).GridInfo;

}

bool domain_base::GridIsLocal(int GridID) const {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridExists(GridID), "Grid %i does not exist.", GridID);

  return LocalGrids_.Contains(GridID);

}

const grid &domain_base::Grid(int GridID) const {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridExists(GridID), "Grid %i does not exist.", GridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = GridRecords_(GridID).GridInfo;
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "Grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  return LocalGrids_(GridID);

}

}

domain::domain(std::shared_ptr<context> &&Context, params &&Params):
  domain_base(std::move(Context), std::move(*Params.Name_), Params.NumDims_, Params.Comm_)
{

  core::logger &Logger = Context_->core_Logger();

  MPI_Barrier(Comm_);

  if (Comm_.Rank() == 0 && Logger.LoggingStatus()) {
    std::string ProcessesString = core::FormatNumber(Comm_.Size(), "processes", "process");
    Logger.LogStatus(true, "Created %1iD domain %s on %s.", NumDims_, *Name_, ProcessesString);
  }

  Level1_.Reset();
  Logger.LogStatus(Comm_.Rank() == 0, "Done creating domain %s.", *Name_);

}

domain::~domain() noexcept {

  if (Context_) {
    MPI_Barrier(Comm_);
    core::logger &Logger = Context_->core_Logger();
    Logger.LogStatus(Comm_.Rank() == 0, "Destroying domain %s...", *Name_);
    Level1_ = Logger.IncreaseStatusLevelAndIndent();
  }

}

domain domain::internal_Create(std::shared_ptr<context> &&Context, params &&Params) {

  return {std::move(Context), std::move(Params)};

}

domain CreateDomain(std::shared_ptr<context> Context, domain::params Params) {

  OVK_DEBUG_ASSERT(Context, "Invalid context.");

  return domain::internal_Create(std::move(Context), std::move(Params));

}

bool domain::ComponentExists(int ComponentID) const {

  OVK_DEBUG_ASSERT(ComponentID >= 0, "Invalid component ID.");

  return Components_.Contains(ComponentID);

}

void domain::DestroyComponent(int ComponentID) {

  MPI_Barrier(Comm_);

  OVK_DEBUG_ASSERT(ComponentID >= 0, "Invalid component ID.");
  OVK_DEBUG_ASSERT(ComponentExists(ComponentID), "Component %i does not exist.", ComponentID);

  ComponentEvent_.Trigger(ComponentID, component_event_flags::DESTROY);

  MPI_Barrier(Comm_);

  Components_.Erase(ComponentID);

  MPI_Barrier(Comm_);

}

bool domain::EditingComponent(int ComponentID) const {

  OVK_DEBUG_ASSERT(ComponentID >= 0, "Invalid component ID.");
  OVK_DEBUG_ASSERT(ComponentExists(ComponentID), "Component %i does not exist.", ComponentID);

  const component_data &Data = Components_(ComponentID);
  const editor &Editor = Data.Editor;

  return Editor.Active();

}

void domain::RestoreComponent(int ComponentID) {

  OVK_DEBUG_ASSERT(ComponentID >= 0, "Invalid component ID.");
  OVK_DEBUG_ASSERT(ComponentExists(ComponentID), "Component %i does not exist.", ComponentID);

  component_data &Data = Components_(ComponentID);
  editor &Editor = Data.Editor;

  OVK_DEBUG_ASSERT(Editor.Active(), "Unable to restore component %i; not currently being edited.",
    ComponentID);

  Editor.Restore();

}

domain::params &domain::params::SetName(std::string Name) {

  Name_ = std::move(Name);

  return *this;

}

domain::params &domain::params::SetDimension(int NumDims) {

  NumDims_ = NumDims;

  return *this;

}

domain::params &domain::params::SetComm(MPI_Comm Comm) {

  Comm_ = Comm;

  return *this;

}

}
