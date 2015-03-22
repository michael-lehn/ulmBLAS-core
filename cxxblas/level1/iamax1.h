#ifndef ULMBLAS_CXXBLAS_LEVEL1_IAMAX1_H
#define ULMBLAS_CXXBLAS_LEVEL1_IAMAX1_H 1

namespace cxxblas {

template <typename IndexType, typename VX>
    IndexType
    iamax1(IndexType      n,
           const VX       *x,
           IndexType      incX);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/iamax1.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_IAMAX1_H
