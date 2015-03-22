#ifndef ULMBLAS_CXXBLAS_LEVEL1_CONJ_H
#define ULMBLAS_CXXBLAS_LEVEL1_CONJ_H 1

namespace cxxblas {

template <typename IndexType, typename VX>
    void
    conj(IndexType      n,
         VX             *x,
         IndexType      incX);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/conj.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_CONJ_H 1
