#ifndef ULMBLAS_CXXBLAS_LEVEL1_ASUM_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_ASUM_TCC 1

#include <ulmblas/cxxblas/level1/asum.h>
#include <ulmblas/cxxblas/auxiliary/macros.h>
#include <ulmblas/impl/level1/asum.h>
#include <ulmblas/impl/level1extensions/asum1.h>

namespace cxxblas {

template <typename IndexType, typename TX, typename Result>
void
asum(IndexType    n,
     const TX     *x,
     IndexType    incX,
     Result       &result)
{
    CXXBLAS_DEBUG_OUT("asum");

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
    CXXBLAS_DEBUG_OUT("asum1");

    if (incX<0) {
        x -= incX*(n-1);
    }
    ulmBLAS::asum1(n, x, incX, result);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_ASUM_TCC
