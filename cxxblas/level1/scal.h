#ifndef ULMBLAS_CXXBLAS_LEVEL1_SCAL_H
#define ULMBLAS_CXXBLAS_LEVEL1_SCAL_H 1

#ifdef USE_EXTERNAL_BLAS
#include <ulmblas/external/blis/scal.h>
#endif // USE_EXTERNAL_BLAS

namespace cxxblas {

template <typename IndexType, typename Alpha, typename VX>
    void
    scal(IndexType      n,
         const Alpha    &alpha,
         VX             *x,
         IndexType      incX);

template <typename IndexType, typename Alpha, typename MA>
    void
    gescal(IndexType    m,
           IndexType    n,
           const Alpha  &alpha,
           MA           *A,
           IndexType    incRowA,
           IndexType    incColA);

template <typename IndexType, typename Alpha, typename MA>
    void
    trscal(IndexType    m,
           IndexType    n,
           const Alpha  &alpha,
           bool         lowerA,
           bool         unitDiagA,
           MA           *A,
           IndexType    incRowA,
           IndexType    incColA);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/scal.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_SCAL_H
