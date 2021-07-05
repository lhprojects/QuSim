#include "SplittingUtils.h"

const Real tpow1t = pow(2, 1 / 3.0);

template<> const Real SplitingConstants<Real>::C1 = 1 / (2 * (2 - tpow1t));
template<> const Real SplitingConstants<Real>::C2 = (1 - tpow1t) / (2 * (2 - tpow1t));
template<> const Real SplitingConstants<Real>::D1 = 1 / (2 - tpow1t);
template<> const Real SplitingConstants<Real>::D2 = -tpow1t / (2 - tpow1t);

template<> const Real32 SplitingConstants<Real32>::C1 = SplitingConstants<Real>::C1;
template<> const Real32 SplitingConstants<Real32>::C2 = SplitingConstants<Real>::C2;
template<> const Real32 SplitingConstants<Real32>::D1 = SplitingConstants<Real>::D1;
template<> const Real32 SplitingConstants<Real32>::D2 = SplitingConstants<Real>::D2;

