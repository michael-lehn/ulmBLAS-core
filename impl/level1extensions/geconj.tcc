#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECONJ_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECONJ_TCC 1

#include <ulmblas/impl/level1extensions/geconj.h>

namespace ulmBLAS {

template <typename IndexType, typename MX>
void
geconj(IndexType      m,
       IndexType      n,
       MX             *X,
       IndexType      incRowX,
       IndexType      incColX)
{
    if (incRowX<incColX) {
        for (IndexType j=0; j<n; ++j) {
            conj(m, &X[j*incColX], incRowX);
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            conj(n, &X[i*incRowX], incColX);
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECONJ_TCC 1
