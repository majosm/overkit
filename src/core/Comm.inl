// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline comm_view::comm_view():
  Comm_(MPI_COMM_NULL)
{}

template <typename T, OVK_FUNCDEF_REQUIRES(!std::is_same<T,comm_view>::value &&
  std::is_convertible<T, MPI_Comm>::value)> comm_view::comm_view(const T &Comm):
  Comm_(Comm)
{}

inline int comm_view::Size() const {

  int Size;
  MPI_Comm_size(Comm_, &Size);

  return Size;

}

inline int comm_view::Rank() const {

  int Rank;
  MPI_Comm_rank(Comm_, &Rank);

  return Rank;

}

inline comm::comm(MPI_Comm Comm) {

  if (Comm != MPI_COMM_NULL) {
    Resource_.reset(new resource(Comm));
    View_ = Comm;
  }

}

inline comm::comm(comm &&Other) noexcept:
  Resource_(std::move(Other.Resource_)),
  View_(Other.View_)
{
  Other.View_.Reset();
}

inline comm &comm::operator=(comm Other) noexcept {

  using std::swap;

  swap(Resource_, Other.Resource_);
  swap(View_, Other.View_);

  return *this;

}

inline void comm::Reset() {

  Resource_.reset();
  View_.Reset();

}

inline MPI_Comm comm::Release() {

  MPI_Comm Comm = Resource_->Comm_;
  Resource_->Comm_ = MPI_COMM_NULL;

  Resource_.reset();
  View_.Reset();

  return Comm;

}

inline bool operator==(const comm_view &Left, const comm_view &Right) {
  return Left.Get() == Right.Get();
}
inline bool operator!=(const comm_view &Left, const comm_view &Right) {
  return !(Left == Right);
}

}
