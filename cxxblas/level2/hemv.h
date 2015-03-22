#ifndef ULMBLAS_CXXBLAS_LEVEL2_HEMV_H
#define ULMBLAS_CXXBLAS_LEVEL2_HEMV_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
    void
    hemv(IndexType    n,
         const Alpha  &alpha,
         bool         lowerA,
         const TA     *A,
         IndexType    incRowA,
         IndexType    incColA,
         const TX     *x,
         IndexType    incX,
         const Beta   &beta,
         TY           *y,
         IndexType    incY);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/hemv.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_HEMV_H
