#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GEAXPY_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GEAXPY_TCC 1

#include <ulmblas/impl/level1extensions/geaxpy.h>
#include <ulmblas/impl/level1/axpy.h>
#include <ulmblas/impl/auxiliary/conjugate.h>

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename MX, typename MY>
void
geaxpy(IndexType      m,
       IndexType      n,
       const Alpha    &alpha,
       const MX       *X,
       IndexType      incRowX,
       IndexType      incColX,
       MY             *Y,
       IndexType      incRowY,
       IndexType      incColY)
{
    const IndexType    UnitStride(1);

    if (m<=0 || n<=0 || alpha==Alpha(0)) {
        return;
    }

    if (incRowX==UnitStride && incRowY==UnitStride) {
//
//      X and Y are both column major
//
        for (IndexType j=0; j<n; ++j) {
            axpy(m, alpha,
                 &X[j*incColX], UnitStride,
                 &Y[j*incColY], UnitStride);
        }
    } else if (incColX==UnitStride && incColY==UnitStride) {
//
//      X and Y are both row major
//
        for (IndexType i=0; i<m; ++i) {
            axpy(n, alpha,
                 &X[i*incRowX], UnitStride,
                 &Y[i*incRowY], UnitStride);
        }
    } else {
//
//      General case
//
        for (IndexType j=0; j<n; ++j) {
            for (IndexType i=0; i<m; ++i) {
                Y[i*incRowY+j*incColY] += alpha*X[i*incRowX+j*incColX];
            }
        }
    }
}

template <typename IndexType, typename Alpha, typename MX, typename MY>
void
geacxpy(IndexType      m,
        IndexType      n,
        const Alpha    &alpha,
        const MX       *X,
        IndexType      incRowX,
        IndexType      incColX,
        MY             *Y,
        IndexType      incRowY,
        IndexType      incColY)
{
    const IndexType    UnitStride(1);

    if (m<=0 || n<=0 || alpha==Alpha(0)) {
        return;
    }

    if (incRowX==UnitStride && incRowY==UnitStride) {
//
//      X and Y are both column major
//
        for (IndexType j=0; j<n; ++j) {
            acxpy(m, alpha,
                  &X[j*incColX], UnitStride,
                  &Y[j*incColY], UnitStride);
        }
    } else if (incColX==UnitStride && incColY==UnitStride) {
//
//      X and Y are both row major
//
        for (IndexType i=0; i<m; ++i) {
            acxpy(n, alpha,
                  &X[i*incRowX], UnitStride,
                  &Y[i*incRowY], UnitStride);
        }
    } else {
//
//      General case
//
        for (IndexType j=0; j<n; ++j) {
            for (IndexType i=0; i<m; ++i) {
                Y[i*incRowY+j*incColY] += alpha
                                          *conjugate(X[i*incRowX+j*incColX]);
            }
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GEAXPY_TCC 1
