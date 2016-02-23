#ifndef ULMBLAS_EXTERNAL_BLIS_SCAL_TCC
#define ULMBLAS_EXTERNAL_BLIS_SCAL_TCC 1

#include <ulmblas/external/blis/scal.h>
#include <blis/blis.h>
#include <complex>

namespace cxxblas {

// bli_sscalv
void
bli_scalv(conj_t conj,
          dim_t  n,
          float  alpha,
          float  *x,
          inc_t  incx)
{
    bli_sscalv(conj, n, (float *)&alpha, x, incx);
}

// bli_dscalv
void
bli_scalv(conj_t conj,
          dim_t  n,
          double alpha,
          double *x,
          inc_t  incx)
{
    bli_dscalv(conj, n, (double *)&alpha, x, incx);
}

// bli_cscalv
void
bli_scalv(conj_t conj,
          dim_t  n,
          const std::complex<float>  &alpha,
          std::complex<float>  *x,
          inc_t  incx)
{
    bli_cscalv(conj, n, (scomplex *)&alpha, (scomplex *)x, incx);
}

// bli_zscalv
void
bli_scalv(conj_t conj,
          dim_t  n,
          const std::complex<double>  &alpha,
          std::complex<double>  *x,
          inc_t  incx)
{
    bli_zscalv(conj, n, (dcomplex *)&alpha, (dcomplex *)x, incx);
}


INSERT_GENTPROT_BASIC( scalv )

template <typename IndexType, typename T>
void
scal(IndexType      n_,
     const T        &alpha,
     T              *x,
     IndexType      incX)
{
#   ifdef CXXBLAS_DEBUG
    std::cout << "BLIS SCAL" << std::endl;
#   endif

    conj_t    conj = BLIS_NO_CONJUGATE;
    dim_t     n      = n_;
    inc_t     incx   = incX;

    err_t     init_result;

    bli_init_auto(&init_result);
    bli_scalv(conj, n, alpha, x, incx);
    bli_finalize_auto(init_result);

}

} // namespace cxxblas

#endif // ULMBLAS_EXTERNAL_BLIS_SCAL_TCC
