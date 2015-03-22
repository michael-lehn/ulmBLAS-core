#ifndef ULMBLAS_CXXBLAS_LEVEL1_ASUM_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_ASUM_TCC 1

#include <ulmblas/cxxblas/level1/asum.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename TX, typename Result>
void
asum(IndexType    n,
     const TX     *x,
     IndexType    incX,
     Result       &result)
{
    if (incX<0) {
        x -= incX*(n-1);
    }
    ulmBLAS::asum(n, x, incX, result);
}

template <typename IndexType, typename TX, typename Result>
void
asum1(IndexType    n,
      const TX     *x,
      IndexType    incX,
      Result       &result)
{
    if (incX<0) {
        x -= incX*(n-1);
    }
    ulmBLAS::asum1(n, x, incX, result);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_ASUM_TCC
