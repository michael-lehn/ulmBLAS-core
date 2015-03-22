#ifndef ULMBLAS_CXXBLAS_LEVEL2_HER_H
#define ULMBLAS_CXXBLAS_LEVEL2_HER_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TA>
    void
    her(IndexType    n,
        const Alpha  &alpha,
        const TX     *x,
        IndexType    incX,
        bool         lowerA,
        TA           *A,
        IndexType    incRowA,
        IndexType    incColA);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/her.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_HER_H
