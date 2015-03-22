#ifndef ULMBLAS_IMPL_LEVEL3_UKERNEL_USYURK_H
#define ULMBLAS_IMPL_LEVEL3_UKERNEL_USYURK_H 1

namespace ulmBLAS {

template <typename IndexType, typename T, typename Beta, typename TC>
    void
    usyurk(IndexType    mr,
           IndexType    nr,
           IndexType    kc,
           const T      &alpha,
           const T      *A,
           const T      *B,
           const Beta   &beta,
           TC           *C,
           IndexType    incRowC,
           IndexType    incColC,
           const T      *nextA,
           const T      *nextB);

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_UKERNEL_USYURK_H

#include <ulmblas/impl/level3/ukernel/usyurk.tcc>
