#ifndef ULMBLAS_CXXBLAS_LEVEL2_SYR2_H
#define ULMBLAS_CXXBLAS_LEVEL2_SYR2_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TY,
          typename TA>
    void
    syr2(IndexType    n,
         const Alpha  &alpha,
         const TX     *x,
         IndexType    incX,
         const TA     *y,
         IndexType    incY,
         bool         lowerA,
         TA           *A,
         IndexType    incRowA,
         IndexType    incColA);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/syr2.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_SYR2_H
