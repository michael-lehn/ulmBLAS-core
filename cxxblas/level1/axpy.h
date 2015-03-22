#ifndef ULMBLAS_CXXBLAS_LEVEL1_AXPY_H
#define ULMBLAS_CXXBLAS_LEVEL1_AXPY_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TY>
    void
    axpy(IndexType    n,
         const Alpha  &alpha,
         const TX     *x,
         IndexType    incX,
         TY           *y,
         IndexType    incY);

template <typename IndexType, typename Alpha, typename TX, typename TY>
    void
    acxpy(IndexType    n,
          const Alpha  &alpha,
          const TX     *x,
          IndexType    incX,
          TY           *y,
          IndexType    incY);

template <typename IndexType, typename Alpha, typename MX, typename MY>
    void
    geaxpy(IndexType      m,
           IndexType      n,
           const Alpha    &alpha,
           bool           transX,
           bool           conjX,
           const MX       *X,
           IndexType      incRowX,
           IndexType      incColX,
           MY             *Y,
           IndexType      incRowY,
           IndexType      incColY);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/axpy.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_AXPY_H
