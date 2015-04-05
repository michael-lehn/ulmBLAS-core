#ifndef ULMBLAS_CXXBLAS_LEVEL1_GECOTR_H
#define ULMBLAS_CXXBLAS_LEVEL1_GECOTR_H 1

namespace cxxblas {

template <typename IndexType, typename MX>
    void
    gecotr(IndexType      n,
           bool           transX,
           bool           conjX,
           MX             *X,
           IndexType      incRowX,
           IndexType      incColX);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/gecotr.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_GECOTR_H
