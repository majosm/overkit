// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
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
#include <ovk/core/Constants.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Disperse.hpp>
#include <ovk/core/Domain.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/IDMap.hpp>
#include <ovk/core/IDSet.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/Recv.hpp>
#include <ovk/core/RecvMap.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/Send.hpp>
#include <ovk/core/SendMap.hpp>
#include <ovk/core/StringWrapper.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <memory>
#include <string>

namespace ovk {

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
    bindings &SetConnectivityComponentID(int ConnectivityComponentID) {
      ConnectivityComponentID_ = ConnectivityComponentID;
      return *this;
    }
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
  void Bind(std::shared_ptr<domain> Domain, bindings Bindings);
  void Unbind();

  const domain &Domain() const;
  domain &Domain();
  const std::shared_ptr<domain> &SharedDomain();

  const id_set<1> &CollectIDs(int MGridID, int NGridID) const;
  bool CollectExists(int MGridID, int NGridID, int CollectID) const;
  void CreateCollect(int MGridID, int NGridID, int CollectID, collect_op CollectOp,
    data_type ValueType, int Count, const range &GridValuesRange, array_layout GridValuesLayout);
  void DestroyCollect(int MGridID, int NGridID, int CollectID);
  void Collect(int MGridID, int NGridID, int CollectID, const void * const *GridValues,
    void **DonorValues);

  const id_set<1> &SendIDs(int MGridID, int NGridID) const;
  bool SendExists(int MGridID, int NGridID, int SendID) const;
  void CreateSend(int MGridID, int NGridID, int SendID, data_type ValueType, int Count, int Tag);
  void DestroySend(int MGridID, int NGridID, int SendID);
  request Send(int MGridID, int NGridID, int SendID, const void * const *DonorValues);

  const id_set<1> &ReceiveIDs(int MGridID, int NGridID) const;
  bool ReceiveExists(int MGridID, int NGridID, int RecvID) const;
  void CreateReceive(int MGridID, int NGridID, int RecvID, data_type ValueType, int Count, int Tag);
  void DestroyReceive(int MGridID, int NGridID, int RecvID);
  request Receive(int MGridID, int NGridID, int RecvID, void **ReceiverValues);

  void Wait(request &Request);
  void WaitAll(array_view<request> Requests);
  void WaitAny(array_view<request> Requests, int &Index);
  // Needed for C API
  void core_WaitAll(array_view<request *> Requests);
  void core_WaitAny(array_view<request *> Requests, int &Index);

  const id_set<1> &DisperseIDs(int MGridID, int NGridID) const;
  bool DisperseExists(int MGridID, int NGridID, int DisperseID) const;
  void CreateDisperse(int MGridID, int NGridID, int DisperseID, disperse_op DisperseOp,
    data_type ValueType, int Count, const range &GridValuesRange, array_layout GridValuesLayout);
  void DestroyDisperse(int MGridID, int NGridID, int DisperseID);
  void Disperse(int MGridID, int NGridID, int DisperseID, const void * const *ReceiverValues,
    void **GridValues);

  static exchanger internal_Create(std::shared_ptr<context> &&Context, params &&Params);

private:

  struct local_m {
    const connectivity_m *Connectivity;
    array<int> DestinationRanks;
    core::collect_map CollectMap;
    core::send_map SendMap;
    id_map<1,core::collect> Collects;
    id_map<1,core::send> Sends;
  };

  struct local_n {
    const connectivity_n *Connectivity;
    array<int> SourceRanks;
    core::recv_map RecvMap;
    id_map<1,core::recv> Recvs;
    id_map<1,core::disperse> Disperses;
  };

  floating_ref_generator<exchanger> FloatingRefGenerator_;

  std::shared_ptr<context> Context_;

  core::string_wrapper Name_;

  std::shared_ptr<domain> Domain_;
  floating_ref<const connectivity_component> ConnectivityComponent_;

  event_listener_handle ComponentEventListener_;

  id_map<2,connectivity_event_flags> ConnectivityEventFlags_;
  event_listener_handle ConnectivityEventListener_;

  id_map<2,local_m,false> LocalMs_;
  id_map<2,local_n,false> LocalNs_;

  exchanger(std::shared_ptr<context> &&Context, params &&Params);

  void OnConnectivityEvent_();

  void DestroyExchangesForDyingConnectivities_();
  void CreateExchangesForNewConnectivities_();
  void UpdateExchangesForModifiedConnectivities_();

  void UpdateSourceDestRanks_();
  void PurgeExchanges_();

  static constexpr int COLLECT_TIME = core::profiler::EXCHANGER_COLLECT_TIME;
  static constexpr int SEND_RECV_TIME = core::profiler::EXCHANGER_SEND_RECV_TIME;
  static constexpr int DISPERSE_TIME = core::profiler::EXCHANGER_DISPERSE_TIME;

};

exchanger CreateExchanger(std::shared_ptr<context> Context, exchanger::params Params={});

}

#endif
