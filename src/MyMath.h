#ifndef MYMATH_H_
#define MYMATH_H_


#include "biodynamo.h"

namespace bdm {

inline real_t mod(Real3 a) {
  real_t mod = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  return mod;
}

inline real_t DotProduct(Real3 a, Real3 b) {
  real_t product = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return product;
}

}; // namespace bdm

#endif  // MOVE_H_

