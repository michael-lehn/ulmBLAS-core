#ifndef ULMBLAS_CXXBLAS_LEVEL1_RAXPY_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_RAXPY_TCC 1

#include <ulmblas/cxxblas/level1/axpy.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TY>
void
raxpy(IndexType    n,
      const Alpha  &alpha,
      const TX     *x,
      IndexType    incX,
      const TY     *y,
      IndexType    incY)
{
    ulmBLAS::raxpy(n, alpha, x, incX, y, incY);
}

template <typename IndexType, typename Alpha, typename TX, typename TY>
void
racxpy(IndexType    n,
       const Alpha  &alpha,
       const TX     *x,
       IndexType    incX,
       const TY     *y,
       IndexType    incY)
{
    ulmBLAS::racxpy(n, alpha, x, incX, y, incY);
}

template <typename IndexType, typename Alpha, typename MX, typename MY>
void
geraxpy(IndexType      m,
        IndexType      n,
        const Alpha    &alpha,
        bool           transX,
        bool           conjX,
        const MX       *X,
        IndexType      incRowX,
        IndexType      incColX,
        MY             *Y,
        IndexType      incRowY,
        IndexType      incColY)
{
    if (!transX) {
        if (!conjX) {
            ulmBLAS::geraxpy(m, n,
                             alpha,
                             X, incRowX, incColX,
                             Y, incRowY, incColY);
        } else {
            ulmBLAS::geracxpy(m, n,
                              alpha,
                              X, incRowX, incColX,
                              Y, incRowY, incColY);
        }
    } else {
        if (!conjX) {
            ulmBLAS::geraxpy(m, n,
                             alpha,
                             X, incColX, incRowX,
                             Y, incRowY, incColY);
        } else {
            ulmBLAS::geracxpy(m, n,
                              alpha,
                              X, incColX, incRowX,
                              Y, incRowY, incColY);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_RAXPY_TCC
