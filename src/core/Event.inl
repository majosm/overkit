// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

template <typename... Args> event<void(Args...)>::event():
  FloatingRefGenerator_(*this)
{}

template <typename... Args> template <typename F, OVK_FUNCDEF_REQUIRES(core::IsCallableWith<F,
  Args...>())> event_listener_handle event<void(Args...)>::AddListener(F Listener) {

  int ID = Listeners_.NextAvailableKey();

  Listeners_.Insert(ID, std::move(Listener));

  return {*this, ID};

}

template <typename... Args> void event<void(Args...)>::Trigger(Args... Arguments) {

  for (auto &Entry : Listeners_) {
    listener &Listener = Entry.Value();
    Listener(Arguments...);
  }

}

template <typename FuncSignature> event_listener_handle::event_listener_handle(event<FuncSignature>
  &Event, int ID):
  ID_(ID)
{

  floating_ref<event<FuncSignature>> EventRef = Event.FloatingRefGenerator_.Generate();

  RemoveFunc_.reset(new std::function<void(int)>([EventRef](int ID) {
    EventRef->Listeners_.Erase(ID);
  }));

}

inline event_listener_handle::event_listener_handle(event_listener_handle &&Other) noexcept:
  ID_(Other.ID_),
  RemoveFunc_(std::move(Other.RemoveFunc_))
{
  Other.ID_ = -1;
}

inline event_listener_handle &event_listener_handle::operator=(event_listener_handle Other)
  noexcept {

  using std::swap;

  swap(ID_, Other.ID_);
  swap(RemoveFunc_, Other.RemoveFunc_);

  return *this;

}

inline event_listener_handle::~event_listener_handle() noexcept {

  Reset();

}

inline void event_listener_handle::Reset() {

  if (ID_ >= 0) {
    (*RemoveFunc_)(ID_);
    RemoveFunc_.reset();
    ID_ = -1;
  }

}

}
