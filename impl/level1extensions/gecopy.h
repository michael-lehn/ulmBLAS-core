#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOPY_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOPY_H 1

namespace ulmBLAS {

template <typename IndexType, typename MX, typename MY>
    void
    gecopy(IndexType      m,
           IndexType      n,
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


} // namespace ulmBLAS

#include <ulmblas/impl/level1extensions/gecopy.tcc>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOPY_H 1
