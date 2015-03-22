#ifndef ULMBLAS_CXXBLAS_LEVEL1_COPY_H
#define ULMBLAS_CXXBLAS_LEVEL1_COPY_H 1

namespace cxxblas {

template <typename IndexType, typename TX, typename TY>
    void
    copy(IndexType    n,
         bool         conjX,
         const TX     *x,
         IndexType    incX,
         TY           *y,
         IndexType    incY);

template <typename IndexType, typename TX, typename TY>
    void
    copy(IndexType    n,
         const TX     *x,
         IndexType    incX,
         TY           *y,
         IndexType    incY);

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
           IndexType      incColY);

template <typename IndexType, typename MX, typename MY>
    void
    gecopy(IndexType      m,
           IndexType      n,
           const MX       *X,
           IndexType      incRowX,
           IndexType      incColX,
           MY             *Y,
           IndexType      incRowY,
           IndexType      incColY);

template <typename IndexType, typename MA, typename MB>
    void
    trcopy(IndexType      m,
           IndexType      n,
           bool           lowerA,
           bool           transA,
           bool           conjA,
           bool           unitDiagA,
           const MA       *A,
           IndexType      incRowA,
           IndexType      incColA,
           MB             *B,
           IndexType      incRowB,
           IndexType      incColB);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/copy.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_COPY_H
