// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline cart::cart(int NumDims):
  NumDims_(NumDims),
  Range_(MakeEmptyRange(NumDims)),
  Periodic_(false,false,false),
  PeriodicStorage_(periodic_storage::UNIQUE)
{}

inline cart::cart(int NumDims, const range &Range, const tuple<bool> &Periodic, periodic_storage
  PeriodicStorage):
  NumDims_(NumDims),
  Range_(Range),
  Periodic_(Periodic),
  PeriodicStorage_(PeriodicStorage)
{}

inline tuple<int> cart::GetPeriodSize() const {

  tuple<int> PeriodSize;

  switch (PeriodicStorage_) {
  case periodic_storage::UNIQUE:
    PeriodSize = Range_.Size();
    break;
  case periodic_storage::DUPLICATED:
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      PeriodSize(iDim) = Range_.Size(iDim) - (Periodic_(iDim) ? 1 : 0);
    }
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    PeriodSize = Range_.Size();
    break;
  }

  return PeriodSize;

}

inline tuple<int> cart::GetPeriod(const tuple<int> &Tuple) const {

  tuple<int> PeriodSize = GetPeriodSize();

  tuple<int> Period;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    if (Periodic_(iDim)) {
      int Offset = Tuple(iDim) - Range_.Begin(iDim);
      Period(iDim) = Offset/PeriodSize(iDim) - int(Offset % PeriodSize(iDim) < 0);
    } else {
      Period(iDim) = 0;
    }
  }

  return Period;

}

inline tuple<int> cart::PeriodicAdjust(const tuple<int> &Tuple) const {

  tuple<int> PeriodSize = GetPeriodSize();

  tuple<int> AdjustedTuple = Tuple;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    if (Periodic_(iDim)) {
      int Mod = (Tuple(iDim) - Range_.Begin(iDim)) % PeriodSize(iDim);
      AdjustedTuple(iDim) = Range_.Begin(iDim) + Mod + PeriodSize(iDim) * (Mod < 0);
    }
  }

  return AdjustedTuple;

}

inline optional<tuple<int>> cart::MapToRange(const range &Range, const tuple<int> &Tuple) const {

  tuple<int> PeriodSize = GetPeriodSize();

  tuple<int> BeginPeriod = GetPeriod(Range.Begin());
  tuple<int> TuplePeriod = GetPeriod(Tuple);

  tuple<int> MappedTuple;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    MappedTuple(iDim) = Tuple(iDim) + PeriodSize(iDim) * (BeginPeriod(iDim) - TuplePeriod(iDim));
    if (MappedTuple(iDim) < Range.Begin(iDim)) MappedTuple(iDim) += PeriodSize(iDim);
  }

  optional<tuple<int>> MaybeMappedTuple;

  if (Range.Contains(MappedTuple)) {
    MaybeMappedTuple = MappedTuple;
  }

  return MaybeMappedTuple;

}

inline optional<range> cart::MapToRange(const range &Range, const range &OtherRange) const {

  optional<range> MaybeMappedOtherRange;

  if (RangesOverlap(Range, OtherRange)) {
    MaybeMappedOtherRange = OtherRange;
  }

  // Prefer to map lower corner of OtherRange into Range if possible
  if (!MaybeMappedOtherRange) {
    auto MaybeMappedOtherLower = MapToRange(Range, OtherRange.Begin());
    if (MaybeMappedOtherLower) {
      range MappedOtherRange;
      MappedOtherRange.Begin() = *MaybeMappedOtherLower;
      MappedOtherRange.End() = MappedOtherRange.Begin() + OtherRange.Size();
      MaybeMappedOtherRange = MappedOtherRange;
    }
  }

  if (!MaybeMappedOtherRange) {
    tuple<int> PeriodSize = GetPeriodSize();
    range AdjustedRange;
    AdjustedRange.Begin() = *MapToRange(Range_, Range.Begin());
    AdjustedRange.End() = AdjustedRange.Begin() + Range.Size();
    range AdjustedOtherRange;
    AdjustedOtherRange.Begin() = *MapToRange(Range_, OtherRange.Begin());
    AdjustedOtherRange.End() = AdjustedOtherRange.Begin() + OtherRange.Size();
    range PeriodSearchRange;
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      PeriodSearchRange.Begin(iDim) = -1;
      PeriodSearchRange.End(iDim) = 1;
    }
    for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
      PeriodSearchRange.Begin(iDim) = 0;
      PeriodSearchRange.End(iDim) = 1;
    }
    for (int k = PeriodSearchRange.Begin(2); k < PeriodSearchRange.End(2); ++k) {
      for (int j = PeriodSearchRange.Begin(1); j < PeriodSearchRange.End(1); ++j) {
        for (int i = PeriodSearchRange.Begin(0); i < PeriodSearchRange.End(0); ++i) {
          tuple<int> Period = {i,j,k};
          range ShiftedOtherRange;
          for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
            ShiftedOtherRange.Begin(iDim) = AdjustedOtherRange.Begin(iDim) + PeriodSize(iDim) *
              Period(iDim);
            ShiftedOtherRange.End(iDim) = AdjustedOtherRange.End(iDim) + PeriodSize(iDim) *
              Period(iDim);
          }
          if (RangesOverlap(AdjustedRange, ShiftedOtherRange)) {
            MaybeMappedOtherRange = ShiftedOtherRange;
          }
        }
      }
    }
  }

  return MaybeMappedOtherRange;

}

inline bool operator==(const cart &Left, const cart &Right) {

  return
    Left.NumDims_ == Right.NumDims_ &&
    Left.Range_ == Right.Range_ &&
    Left.Periodic_ == Right.Periodic_ &&
    Left.PeriodicStorage_ == Right.PeriodicStorage_;

}

inline bool operator!=(const cart &Left, const cart &Right) {

  return !(Left == Right);

}

inline cart MakeEmptyCart(int NumDims) {

  return {NumDims, MakeEmptyRange(NumDims), MakeUniformTuple<int>(NumDims, false),
    periodic_storage::UNIQUE};

}

}
