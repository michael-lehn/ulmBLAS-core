#ifndef ULMBLAS_CXXBLAS_LEVEL1_ASUM_H
#define ULMBLAS_CXXBLAS_LEVEL1_ASUM_H 1

namespace cxxblas {

template <typename IndexType, typename TX, typename Result>
    void
    asum(IndexType    n,
         const TX     *x,
         IndexType    incX,
         Result       &result);

template <typename IndexType, typename TX, typename Result>
    void
    asum1(IndexType    n,
          const TX     *x,
          IndexType    incX,
          Result       &result);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/asum.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_ASUM_H
