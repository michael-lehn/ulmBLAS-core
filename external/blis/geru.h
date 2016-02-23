#ifndef ULMBLAS_EXTERNAL_BLIS_GERU_H
#define ULMBLAS_EXTERNAL_BLIS_GERU_H 1

namespace cxxblas {

template <typename IndexType, typename T>
    void
    geru(IndexType    m_,
         IndexType    n_,
         const T      &alpha,
         const T      *x,
         IndexType    incX,
         const T      *y,
         IndexType    incY,
         T            *A,
         IndexType    incRowA,
         IndexType    incColA);

} // namespace cxxblas

#endif // ULMBLAS_EXTERNAL_BLIS_GERU_H

#include <ulmblas/external/blis/geru.tcc>
