#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOTR_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOTR_TCC 1

#include <ulmblas/impl/level1extensions/gecotr.h>
#include <ulmblas/impl/level1/swap.h>
#include <ulmblas/impl/auxiliary/conjugate.h>

namespace ulmBLAS {

template <typename IndexType, typename MX>
void
gecotr(IndexType      n,
       bool           transX,
       bool           conjX,
       MX             *X,
       IndexType      incRowX,
       IndexType      incColX)
{
    if (!transX && !conjX) {
        return;
    }

    if (transX) {
        if (!conjX) {
            for (IndexType i=0; i<n; ++i) {
                swap(n-1-i,
                     &X[(i+1)*incRowX+i*incColX], incRowX,
                     &X[i*incRowX+(i+1)*incColX], incColX);
            }
        } else {
            for (IndexType i=0; i<n; ++i) {
                swapc(n-1-i,
                      &X[(i+1)*incRowX+i*incColX], incRowX,
                      &X[i*incRowX+(i+1)*incColX], incColX);
                X[i*(incRowX+incColX)] = conjugate(X[i*(incRowX+incColX)]);
            }
        }
    } else {
        for (IndexType j=0; j<n; ++j) {
            for (IndexType i=0; i<n; ++i) {
                X[i*incRowX+j*incColX] = conjugate(X[i*incRowX+j*incColX]);
            }
        }
    }
}


} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOTR_TCC 1
