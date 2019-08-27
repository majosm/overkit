// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_EXCHANGER_HPP_INCLUDED
#define OVK_CORE_EXCHANGER_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Collect.hpp>
#include <ovk/core/CollectMap.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/ConnectivityComponent.hpp>
#include <ovk/core/ConnectivityM.hpp>
#include <ovk/core/ConnectivityN.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Disperse.hpp>
#include <ovk/core/DisperseMap.hpp>
#include <ovk/core/Domain.hpp>
#include <ovk/core/ElemMap.hpp>
#include <ovk/core/ElemSet.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/Exchanger.h>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Recv.hpp>
#include <ovk/core/RecvMap.hpp>
#include <ovk/core/Send.hpp>
#include <ovk/core/SendMap.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/StringWrapper.hpp>

#include <mpi.h>

#include <memory>
#include <string>

namespace ovk {

enum class collect_op {
  NONE = OVK_COLLECT_NONE,
  ANY = OVK_COLLECT_ANY,
  NOT_ALL = OVK_COLLECT_NOT_ALL,
  ALL = OVK_COLLECT_ALL,
  INTERPOLATE = OVK_COLLECT_INTERPOLATE
};

inline bool ValidCollectOp(collect_op CollectOp) {
  return ovkValidCollectOp(ovk_collect_op(CollectOp));
}

enum class disperse_op {
  OVERWRITE = OVK_DISPERSE_OVERWRITE,
  APPEND = OVK_DISPERSE_APPEND
};

inline bool ValidDisperseOp(disperse_op DisperseOp) {
  return ovkValidDisperseOp(ovk_disperse_op(DisperseOp));
}

class exchanger {

public:

  class params {
  public:
    params() = default;
    const std::string &Name() const { return Name_; }
    params &SetName(std::string Name);
  private:
    core::string_wrapper Name_ = "Exchanger";
    friend class exchanger;
  };

  class bindings {
  public:
    bindings() = default;
    int ConnectivityComponentID() const { return ConnectivityComponentID_; }
    bindings &SetConnectivityComponentID(int ConnectivityComponentID);
  private:
    int ConnectivityComponentID_ = -1;
    friend class exchanger;
  };

  exchanger(const exchanger &Other) = delete;
  exchanger(exchanger &&Other) noexcept = default;

  exchanger &operator=(const exchanger &Other) = delete;
  exchanger &operator=(exchanger &&Other) noexcept = default;

  ~exchanger() noexcept;

  const context &Context() const;
  context &Context();
  const std::shared_ptr<context> &SharedContext() const;

  const std::string &Name() const { return *Name_; }

  bool Bound() const;
  void Bind(const domain &Domain, bindings Bindings);
  void Unbind();

  const domain &Domain() const;

  const set<int> &CollectIDs(int MGridID, int NGridID) const;
  bool CollectExists(int MGridID, int NGridID, int CollectID) const;
  void CreateCollect(int MGridID, int NGridID, int CollectID, collect_op CollectOp,
    data_type ValueType, int Count, const range &GridValuesRange, array_layout GridValuesLayout);
  void DestroyCollect(int MGridID, int NGridID, int CollectID);
  // "GridValues" actual type is const T * const *
  // "DonorValues" actual type is T **
  void Collect(int MGridID, int NGridID, int CollectID, const void *GridValues, void *DonorValues);

  const set<int> &SendIDs(int MGridID, int NGridID) const;
  bool SendExists(int MGridID, int NGridID, int SendID) const;
  void CreateSend(int MGridID, int NGridID, int SendID, data_type ValueType, int Count, int Tag);
  void DestroySend(int MGridID, int NGridID, int SendID);
  // "DonorValues" actual type is const T * const *
  request Send(int MGridID, int NGridID, int SendID, const void *DonorValues);

  const set<int> &ReceiveIDs(int MGridID, int NGridID) const;
  bool ReceiveExists(int MGridID, int NGridID, int RecvID) const;
  void CreateReceive(int MGridID, int NGridID, int RecvID, data_type ValueType, int Count, int Tag);
  void DestroyReceive(int MGridID, int NGridID, int RecvID);
  // "ReceiverValues" actual type is T **
  request Receive(int MGridID, int NGridID, int RecvID, void *ReceiverValues);

  const set<int> &DisperseIDs(int MGridID, int NGridID) const;
  bool DisperseExists(int MGridID, int NGridID, int DisperseID) const;
  void CreateDisperse(int MGridID, int NGridID, int DisperseID, disperse_op DisperseOp,
    data_type ValueType, int Count, const range &GridValuesRange, array_layout GridValuesLayout);
  void DestroyDisperse(int MGridID, int NGridID, int DisperseID);
  // "ReceiverValues" actual type is const T * const *
  // "GridValues" actual type is T **
  void Disperse(int MGridID, int NGridID, int DisperseID, const void *ReceiverValues, void
    *GridValues);

  static exchanger internal_Create(std::shared_ptr<context> &&Context, params &&Params);

private:

  struct local_m {
    const connectivity_m *Connectivity;
    array<int> DestinationRanks;
    core::collect_map CollectMap;
    core::send_map SendMap;
    map<int,core::collect> Collects;
    map<int,core::send> Sends;
  };

  struct local_n {
    const connectivity_n *Connectivity;
    array<int> SourceRanks;
    core::recv_map RecvMap;
    core::disperse_map DisperseMap;
    map<int,core::recv> Recvs;
    map<int,core::disperse> Disperses;
  };

  struct update_manifest {
    elem_set<int,2> CreateLocal;
    elem_set<int,2> DestroyLocal;
    elem_set<int,2> UpdateSourceDestRanks;
    elem_set<int,2> ResetExchanges;
  };

  floating_ref_generator FloatingRefGenerator_;

  std::shared_ptr<context> Context_;

  core::string_wrapper Name_;

  floating_ref<const domain> Domain_;
  event_listener_handle ComponentEventListener_;

  int ConnectivityComponentID_ = -1;
  floating_ref<const connectivity_component> ConnectivityComponent_;
  event_listener_handle ConnectivityEventListener_;

  elem_map_noncontig<int,2,local_m> LocalMs_;
  elem_map_noncontig<int,2,local_n> LocalNs_;

  update_manifest UpdateManifest_;

  exchanger(std::shared_ptr<context> &&Context, params &&Params);

  void OnComponentEvent_(int ComponentID, component_event_flags Flags);
  void OnConnectivityEvent_(int MGridID, int NGridID, connectivity_event_flags Flags, bool
    LastInSequence);

  void Update_();

  void CreateLocals_();
  void DestroyLocals_();
  void UpdateSourceDestRanks_();
  void ResetExchanges_();

  static constexpr int COLLECT_TIME = core::profiler::EXCHANGER_COLLECT_TIME;
  static constexpr int SEND_RECV_TIME = core::profiler::EXCHANGER_SEND_RECV_TIME;
  static constexpr int DISPERSE_TIME = core::profiler::EXCHANGER_DISPERSE_TIME;

};

exchanger CreateExchanger(std::shared_ptr<context> Context, exchanger::params Params={});

}

#endif
