#ifndef ULMBLAS_CXXBLAS_LEVEL2_SYR_H
#define ULMBLAS_CXXBLAS_LEVEL2_SYR_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TA>
    void
    syr(IndexType    n,
        const Alpha  &alpha,
        const TX     *x,
        IndexType    incX,
        bool         lowerA,
        TA           *A,
        IndexType    incRowA,
        IndexType    incColA);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/syr.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_SYR_H
