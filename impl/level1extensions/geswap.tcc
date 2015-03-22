#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESWAP_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESWAP_TCC 1

#include <ulmblas/impl/level1extensions/geswap.h>

namespace ulmBLAS {

template <typename IndexType, typename MX, typename MY>
void
geswap(IndexType      m,
       IndexType      n,
       const MX       *X,
       IndexType      incRowX,
       IndexType      incColX,
       MY             *Y,
       IndexType      incRowY,
       IndexType      incColY)
{
    if (incRowX<incColX) {
        for (IndexType j=0; j<n; ++j) {
            swap(m, &X[j*incColX], incRowX, &Y[j*incColY], incRowY);
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            swap(n, &X[i*incRowX], incColX, &Y[i*incRowY], incColY);
        }
    }
}

template <typename IndexType, typename MX, typename MY>
void
geswapc(IndexType      m,
        IndexType      n,
        const MX       *X,
        IndexType      incRowX,
        IndexType      incColX,
        MY             *Y,
        IndexType      incRowY,
        IndexType      incColY)
{
    if (incRowX<incColX) {
        for (IndexType j=0; j<n; ++j) {
            swapc(m, &X[j*incColX], incRowX, &Y[j*incColY], incRowY);
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            swapc(n, &X[i*incRowX], incColX, &Y[i*incRowY], incColY);
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESWAP_TCC 1
