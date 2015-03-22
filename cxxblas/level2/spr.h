#ifndef ULMBLAS_CXXBLAS_LEVEL2_SPR_H
#define ULMBLAS_CXXBLAS_LEVEL2_SPR_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TA>
    void
    spr(IndexType    n,
        const Alpha  &alpha,
        const TX     *x,
        IndexType    incX,
        bool         colMajorA,
        bool         lowerA,
        TA           *A);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/spr.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_SPR_H
