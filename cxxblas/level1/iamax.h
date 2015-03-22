#ifndef ULMBLAS_CXXBLAS_LEVEL1_IAMAX_H
#define ULMBLAS_CXXBLAS_LEVEL1_IAMAX_H 1

namespace cxxblas {

template <typename IndexType, typename VX>
    IndexType
    iamax(IndexType      n,
          const VX       *x,
          IndexType      incX);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/iamax.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_IAMAX_H
