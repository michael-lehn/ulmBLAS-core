#ifndef ULMBLAS_CXXBLAS_LEVEL1_COPY_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_COPY_TCC 1

#include <ulmblas/cxxblas/level1/copy.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename TX, typename TY>
void
copy(IndexType    n,
     bool         conjX,
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
    ulmBLAS::copy(n, conjX, x, incX, y, incY);
}

template <typename IndexType, typename TX, typename TY>
void
copy(IndexType    n,
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
    ulmBLAS::copy(n, false, x, incX, y, incY);
}

template <typename IndexType, typename MX, typename MY>
void
gecopy(IndexType      m,
       IndexType      n,
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
        ulmBLAS::gecopy(m, n, conjX, X, incRowX, incColX, Y, incRowY, incColY);
    } else {
        ulmBLAS::gecopy(m, n, conjX, X, incColX, incRowX, Y, incRowY, incColY);
    }
}

template <typename IndexType, typename MX, typename MY>
void
gecopy(IndexType      m,
       IndexType      n,
       const MX       *X,
       IndexType      incRowX,
       IndexType      incColX,
       MY             *Y,
       IndexType      incRowY,
       IndexType      incColY)
{
    gecopy(m, n,
           false, false,
           X, incRowX, incColX,
           Y, incRowY, incColY);
}

template <typename IndexType, typename MA, typename MB>
void
trcopy(IndexType      m,
       IndexType      n,
       bool           lowerA,
       bool           transA,
       bool           conjA,
       bool           unitDiagA,
       const MA       *A,
       IndexType      incRowA_,
       IndexType      incColA_,
       MB             *B,
       IndexType      incRowB,
       IndexType      incColB)
{
    IndexType incRowA = (!transA) ? incRowA_ : incColA_;
    IndexType incColA = (!transA) ? incColA_ : incRowA_;

    if (!transA) {
        if (lowerA) {
            //
            // Copy lower part of A to lower part of B
            //
            ulmBLAS::trlcopy(m, n,
                             unitDiagA, conjA, A, incRowA, incColA,
                             B, incRowB, incColB);
        } else {
            //
            // Copy upper part of A to upper part of B
            //
            ulmBLAS::trlcopy(n, m,
                             unitDiagA, conjA, A, incColA, incRowA,
                             B, incColB, incRowB);
        }
    } else {
        if (lowerA) {
            //
            // Copy lower part of A to upper part of B
            //
            ulmBLAS::trlcopy(m, n,
                             unitDiagA, conjA, A, incRowA, incColA,
                             B, incColB, incRowB);
        } else {
            //
            // Copy upper part of A to lower part of B
            //
            ulmBLAS::trlcopy(n, m,
                             unitDiagA, conjA, A, incColA, incRowA,
                             B, incRowB, incColB);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_COPY_TCC
