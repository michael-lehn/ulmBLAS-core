#ifndef ULMBLAS_IMPL_LAPACK_GETF2_TCC
#define ULMBLAS_IMPL_LAPACK_GETF2_TCC 1

#include <algorithm>
#include <ulmblas/impl/lapack/getf2.h>
#include <ulmblas/impl/lapack/safemin.h>
#include <ulmblas/impl/level1/iamax.h>
#include <ulmblas/impl/level1/scal.h>
#include <ulmblas/impl/level1/swap.h>
#include <ulmblas/impl/level2/ger.h>

namespace ulmBLAS {

template <typename IndexType, typename T>
IndexType
getf2(IndexType    m,
      IndexType    n,
      T            *A,
      IndexType    incRowA,
      IndexType    incColA,
      IndexType    *piv,
      IndexType    incPiv)
{
    const T sMin = safeMin<T>();

    IndexType info = 0;

    if (m==0 || n==0) {
        return info;
    }

    for (IndexType j=0; j<std::min(m,n); ++j) {

        IndexType jp  = j+iamax(m-j, &A[j*incRowA+j*incColA], incRowA);
        piv[j*incPiv] = jp;

        if (A[jp*incRowA+j*incColA]!=T(0)) {
            if (jp!=j) {
                swap(n, &A[j*incRowA], incColA, &A[jp*incRowA], incColA);
            }
            if (j<m-1) {
                if (fabs(A[j*incRowA+j*incColA])>= sMin) {
                    scal(m-j-1, T(1)/A[j*incRowA+j*incColA],
                         &A[(j+1)*incRowA+j*incColA], incRowA);
                } else {
                    for (IndexType i=1; i<m-j; ++i) {
                        A[(j+i)*incRowA+j*incColA] /= A[j*incRowA+j*incColA];
                    }
                }
            }
        } else if (info==0) {
            info = j+1;
        }
        if (j<std::min(m,n)-1) {
            ger(m-j-1, n-j-1, -T(1),
                &A[(j+1)*incRowA +  j   *incColA], incRowA,
                &A[ j   *incRowA + (j+1)*incColA], incColA,
                &A[(j+1)*incRowA + (j+1)*incColA], incRowA, incColA);
        }
    }
    return info;
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LAPACK_GETF2_TCC
