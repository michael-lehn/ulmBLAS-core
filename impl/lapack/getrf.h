#ifndef ULMBLAS_IMPL_LAPACK_GETRF_H
#define ULMBLAS_IMPL_LAPACK_GETRF_H 1

namespace ulmBLAS {

template <typename IndexType, typename T>
    IndexType
    getrf(IndexType    m,
          IndexType    n,
          T            *A,
          IndexType    incRowA,
          IndexType    incColA,
          IndexType    *piv,
          IndexType    incPiv);

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LAPACK_GETRF_H

#include <ulmblas/impl/lapack/getrf.tcc>
