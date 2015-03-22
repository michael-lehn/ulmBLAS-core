#ifndef ULMBLAS_CXXBLAS_LEVEL1_SCAL_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_SCAL_TCC 1

#include <ulmblas/cxxblas/level1/scal.h>
#include <ulmblas/impl/level1extensions/gescal.h>
#include <ulmblas/impl/level1/scal.h>
#include <ulmblas/impl/level1extensions/trlscal.h>
#include <ulmblas/impl/level1extensions/truscal.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename VX>
void
scal(IndexType      n,
     const Alpha    &alpha,
     VX             *x,
     IndexType      incX)
{
    ulmBLAS::scal(n, alpha, x, incX);
}

template <typename IndexType, typename Alpha, typename MA>
void
gescal(IndexType    m,
       IndexType    n,
       const Alpha  &alpha,
       MA           *A,
       IndexType    incRowA,
       IndexType    incColA)
{
    ulmBLAS::gescal(m, n, alpha, A, incRowA, incColA);
}

template <typename IndexType, typename Alpha, typename MA>
void
trscal(IndexType    m,
       IndexType    n,
       const Alpha  &alpha,
       bool         lowerA,
       bool         unitDiagA,
       MA           *A,
       IndexType    incRowA,
       IndexType    incColA)
{
    if (lowerA) {
        ulmBLAS::trlscal(m, n, unitDiagA, alpha, A, incRowA, incColA);
    } else {
        ulmBLAS::truscal(m, n, unitDiagA, alpha, A, incRowA, incColA);
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_SCAL_TCC
