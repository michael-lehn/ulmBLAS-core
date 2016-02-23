#ifndef ULMBLAS_EXTERNAL_BLIS_TRSM_H
#define ULMBLAS_EXTERNAL_BLIS_TRSM_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TB>
    void
    trsm(bool         leftA,
         IndexType    m,
         IndexType    n,
         const Alpha  &alpha,
         bool         lowerA,
         bool         transA,
         bool         conjA,
         bool         unitDiagA,
         const TA     *A,
         IndexType    incRowA,
         IndexType    incColA,
         TB           *B,
         IndexType    incRowB,
         IndexType    incColB);

} // namespace cxxblas

#endif // ULMBLAS_EXTERNAL_BLIS_TRSM_H

#include <ulmblas/external/blis/trsm.tcc>

