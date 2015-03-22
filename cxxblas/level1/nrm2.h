#ifndef ULMBLAS_CXXBLAS_LEVEL1_NRM2_H
#define ULMBLAS_CXXBLAS_LEVEL1_NRM2_H 1

namespace cxxblas {

template <typename IndexType, typename VX, typename Result>
    void
    nrm2(IndexType  n,
         const VX   *x,
         IndexType  incX,
         Result     &result);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/nrm2.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_NRM2_H
