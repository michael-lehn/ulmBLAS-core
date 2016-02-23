#ifndef ULMBLAS_EXTERNAL_BLIS_TRSM_TCC
#define ULMBLAS_EXTERNAL_BLIS_TRSM_TCC 1

#include <ulmblas/external/blis/trsm.h>
#include <blis/blis.h>
#include <complex>

namespace cxxblas {

// bli_strsm
void
bli_trsm(side_t side, uplo_t uploa, trans_t transa, diag_t diaga,
         dim_t m, dim_t n,
         float alpha,
         const float *A, inc_t rs_A, inc_t cs_A,
         float *B, inc_t rs_B, inc_t cs_B)
{
    bli_strsm(side, uploa, transa, diaga,
              m, n,
              &alpha,
              (float *)A, rs_A, cs_A,
              (float *)B, rs_B, cs_B);
}

// bli_dtrsm
void
bli_trsm(side_t side, uplo_t uploa, trans_t transa, diag_t diaga,
         dim_t m, dim_t n,
         double alpha,
         const double *A, inc_t rs_A, inc_t cs_A,
         double *B, inc_t rs_B, inc_t cs_B)
{
    bli_dtrsm(side, uploa, transa, diaga,
              m, n,
              &alpha,
              (double *)A, rs_A, cs_A,
              (double *)B, rs_B, cs_B);
}

// bli_ctrsm
void
bli_trsm(side_t side, uplo_t uploa, trans_t transa, diag_t diaga,
         dim_t m, dim_t n,
         std::complex<float> alpha,
         const std::complex<float> *A, inc_t rs_A, inc_t cs_A,
         std::complex<float> *B, inc_t rs_B, inc_t cs_B)
{
    bli_ctrsm(side, uploa, transa, diaga,
              m, n,
              (scomplex *)&alpha,
              (scomplex *)A, rs_A, cs_A,
              (scomplex *)B, rs_B, cs_B);
}

// bli_ztrsm
void
bli_trsm(side_t side, uplo_t uploa, trans_t transa, diag_t diaga,
         dim_t m, dim_t n,
         std::complex<double> &alpha,
         const std::complex<double> *A, inc_t rs_A, inc_t cs_A,
         std::complex<double> *B, inc_t rs_B, inc_t cs_B)
{
    bli_ztrsm(side, uploa, transa, diaga,
              m, n,
              (dcomplex *)&alpha,
              (dcomplex *)A, rs_A, cs_A,
              (dcomplex *)B, rs_B, cs_B);
}

template <typename IndexType, typename T>
void
trsm(bool         leftA,
     IndexType    m_,
     IndexType    n_,
     const T      &alpha,
     bool         lowerA,
     bool         transA,
     bool         conjA,
     bool         unitDiagA,
     const T      *A,
     IndexType    incRowA,
     IndexType    incColA,
     T            *B,
     IndexType    incRowB,
     IndexType    incColB)
{
#   ifdef CXXBLAS_DEBUG
    std::cout << "BLIS TRSM" << std::endl;
#   endif

    side_t  side  = (leftA)     ? BLIS_LEFT
                                : BLIS_RIGHT;
    uplo_t  uplo  = (lowerA)    ? BLIS_LOWER
                                : BLIS_UPPER;
    trans_t trans = (!transA)   ? (!conjA) ? BLIS_NO_TRANSPOSE
                                           : BLIS_CONJ_NO_TRANSPOSE
                                : (!conjA) ? BLIS_TRANSPOSE
                                           : BLIS_CONJ_TRANSPOSE;
    diag_t  diag  = (unitDiagA) ? BLIS_UNIT_DIAG
                                : BLIS_NONUNIT_DIAG;
    dim_t   m     = m_;
    dim_t   n     = n_;
    inc_t   rs_A  = incRowA;
    inc_t   cs_A  = incColA;
    inc_t   rs_B  = incRowB;
    inc_t   cs_B  = incColB;

    err_t     init_result;

    bli_init_auto(&init_result);
    bli_trsm(side, uplo, trans, diag,
             m, n,
             alpha,
             A, rs_A, cs_A,
             B, rs_B, cs_B);
    bli_finalize_auto(init_result);
}

} // namespace cxxblas

#endif // ULMBLAS_EXTERNAL_BLIS_TRSM_TCC
