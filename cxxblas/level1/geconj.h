#ifndef ULMBLAS_CXXBLAS_LEVEL1_GECONJ_H
#define ULMBLAS_CXXBLAS_LEVEL1_GECONJ_H 1

namespace cxxblas {

template <typename IndexType, typename MX, typename MY>
    void
    geconj(IndexType      m,
           IndexType      n,
           MX             *X,
           IndexType      incRowX,
           IndexType      incColX);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/geconj.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_GECONJ_H
