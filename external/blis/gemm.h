#ifndef ULMBLAS_EXTERNAL_BLIS_GEMM_H
#define ULMBLAS_EXTERNAL_BLIS_GEMM_H 1

namespace cxxblas {

template <typename IndexType, typename T>
    void
    gemm(IndexType    m_,
         IndexType    n_,
         IndexType    k_,
         const T      alpha,
         bool         transA,
         bool         conjA,
         const T      *A,
         IndexType    incRowA,
         IndexType    incColA,
         bool         transB,
         bool         conjB,
         const T      *B,
         IndexType    incRowB,
         IndexType    incColB,
         const T      &beta,
         T            *C,
         IndexType    incRowC,
         IndexType    incColC);

} // namespace cxxblas

#endif // ULMBLAS_EXTERNAL_BLIS_GEMM_H

#include <ulmblas/external/blis/gemm.tcc>
