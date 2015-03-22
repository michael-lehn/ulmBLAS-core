#ifndef ULMBLAS_CXXBLAS_LEVEL1_AXPY_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_AXPY_TCC 1

#include <ulmblas/cxxblas/level1/axpy.h>
#include <ulmblas/impl/level1/axpy.h>
#include <ulmblas/impl/level1extensions/geaxpy.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TY>
void
axpy(IndexType    n,
     const Alpha  &alpha,
     const TX     *x,
     IndexType    incX,
     TY           *y,
     IndexType    incY)
{
    if (incX<0) {
        x -= incX*(n-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }
    ulmBLAS::axpy(n, alpha, x, incX, y, incY);
}

template <typename IndexType, typename Alpha, typename TX, typename TY>
void
acxpy(IndexType    n,
      const Alpha  &alpha,
      const TX     *x,
      IndexType    incX,
      TY           *y,
      IndexType    incY)
{
    if (incX<0) {
        x -= incX*(n-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }
    ulmBLAS::acxpy(n, alpha, x, incX, y, incY);
}

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
       IndexType      incColY)
{
    if (!transX) {
        if (!conjX) {
            ulmBLAS::geaxpy(m, n,
                            alpha,
                            X, incRowX, incColX,
                            Y, incRowY, incColY);
        } else {
            ulmBLAS::geacxpy(m, n,
                             alpha,
                             X, incRowX, incColX,
                             Y, incRowY, incColY);
        }
    } else {
        if (!conjX) {
            ulmBLAS::geaxpy(m, n,
                            alpha,
                            X, incColX, incRowX,
                            Y, incRowY, incColY);
        } else {
            ulmBLAS::geacxpy(m, n,
                             alpha,
                             X, incColX, incRowX,
                             Y, incRowY, incColY);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_AXPY_TCC
