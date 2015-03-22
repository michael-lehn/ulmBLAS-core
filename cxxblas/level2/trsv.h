#ifndef ULMBLAS_CXXBLAS_LEVEL2_TRSV_H
#define ULMBLAS_CXXBLAS_LEVEL2_TRSV_H 1

namespace cxxblas {

template <typename IndexType, typename TA, typename TX>
    void
    trsv(IndexType    n,
         bool         lowerA,
         bool         transA,
         bool         conjA,
         bool         unitDiagA,
         const TA     *A,
         IndexType    incRowA,
         IndexType    incColA,
         TX           *x,
         IndexType    incX);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/trsv.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_TRSV_H
