#ifndef ULMBLAS_EXTERNAL_BLIS_GERU_TCC
#define ULMBLAS_EXTERNAL_BLIS_GERU_TCC 1

#include <ulmblas/external/blis/geru.h>
#include <blis/blis.h>
#include <complex>

namespace cxxblas {

// bli_sger
void
bli_ger(conj_t conjx, conj_t conjy,
        dim_t m, dim_t n,
        float alpha,
        const float *x, inc_t incx,
        const float *y, inc_t incy,
        float *A, inc_t rs_A, inc_t cs_A)
{
    bli_sger(conjx, conjy,
             m, n,
             (float *)&alpha,
             (float *)x, incx,
             (float *)y, incy,
             A, rs_A, cs_A);
}

// bli_dger
void
bli_ger(conj_t conjx, conj_t conjy,
        dim_t m, dim_t n,
        double alpha,
        const double *x, inc_t incx,
        const double *y, inc_t incy,
        double *A, inc_t rs_A, inc_t cs_A)
{
    bli_dger(conjx, conjy,
             m, n,
             (double *)&alpha,
             (double *)x, incx,
             (double *)y, incy,
             A, rs_A, cs_A);
}

// bli_cger
void
bli_ger(conj_t conjx, conj_t conjy,
        dim_t m, dim_t n,
        const std::complex<float> &alpha,
        const std::complex<float> *x, inc_t incx,
        const std::complex<float> *y, inc_t incy,
        std::complex<float> *A, inc_t rs_A, inc_t cs_A)
{
    bli_cger(conjx, conjy,
             m, n,
             (scomplex *)&alpha,
             (scomplex *)x, incx,
             (scomplex *)y, incy,
             (scomplex *)A, rs_A, cs_A);
}

// bli_zger
void
bli_ger(conj_t conjx, conj_t conjy,
        dim_t m, dim_t n,
        const std::complex<double> &alpha,
        const std::complex<double> *x, inc_t incx,
        const std::complex<double> *y, inc_t incy,
        std::complex<double> *A, inc_t rs_A, inc_t cs_A)
{
    bli_zger(conjx, conjy,
             m, n,
             (dcomplex *)&alpha,
             (dcomplex *)x, incx,
             (dcomplex *)y, incy,
             (dcomplex *)A, rs_A, cs_A);
}

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
     IndexType    incColA)
{
#   ifdef CXXBLAS_DEBUG
    std::cout << "BLIS GERU" << std::endl;
#   endif

    conj_t    conjx  = BLIS_NO_CONJUGATE;
    conj_t    conjy  = BLIS_NO_CONJUGATE;

    dim_t     m      = m_;
    dim_t     n      = n_;

    inc_t     incx   = incX;
    inc_t     incy   = incY;

    inc_t     rs_A   = incRowA;
    inc_t     cs_A   = incColA;

    err_t     init_result;

    bli_init_auto(&init_result);
    bli_ger(conjx, conjy,
            m, n,
            alpha,
            x, incx,
            y, incy,
            A, rs_A, cs_A);
    bli_finalize_auto(init_result);
}

} // namespace cxxblas

#endif // ULMBLAS_EXTERNAL_BLIS_GERU_TCC

