#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESCAL_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESCAL_TCC 1

#include <ulmblas/impl/level1/scal.h>
#include <ulmblas/impl/level1extensions/gescal.h>

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename MA>
void
gescal(IndexType    m,
       IndexType    n,
       const Alpha  &alpha,
       MA           *A,
       IndexType    incRowA,
       IndexType    incColA)
{
    const IndexType    UnitStride(1);

    if (alpha==Alpha(1) || m<=0 || n<=0) {
        return;
    }

    if (incRowA==UnitStride) {
        for (IndexType j=0; j<n; ++j) {
            scal(m, alpha, &A[j*incColA], UnitStride);
        }
    } else if (incColA==UnitStride) {
        for (IndexType i=0; i<m; ++i) {
            scal(n, alpha, &A[i*incRowA], UnitStride);
        }
    } else {
        for (IndexType j=0; j<n; ++j) {
            for (IndexType i=0; i<m; ++i) {
                A[i*incRowA+j*incColA] *= alpha;
            }
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESCAL_TCC 1
