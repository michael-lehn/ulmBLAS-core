#ifndef ULMBLAS_CXXBLAS_LEVEL2_SYMV_H
#define ULMBLAS_CXXBLAS_LEVEL2_SYMV_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
    void
    symv(IndexType    n,
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

#include <ulmblas/cxxblas/level2/symv.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_SYMV_H
