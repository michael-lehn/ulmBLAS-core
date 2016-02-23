#ifndef ULMBLAS_EXTERNAL_BLIS_GEMM_TCC
#define ULMBLAS_EXTERNAL_BLIS_GEMM_TCC 1

#include <ulmblas/external/blis/gemm.h>
#include <blis/blis.h>
#include <complex>

namespace cxxblas {

// bli_sgemm
void
bli_gemm(trans_t transA, trans_t transB,
         dim_t m, dim_t n, dim_t k,
         float alpha,
         const float *A, inc_t rs_A, inc_t cs_A,
         const float *B, inc_t rs_B, inc_t cs_B,
         float beta,
         float *C, inc_t rs_C, inc_t cs_C)
{
    bli_sgemm(transA, transB,
              m, n, k,
              &alpha,
              (float *)A, rs_A, cs_A,
              (float *)B, rs_B, cs_B,
              &beta,
              C, rs_C, cs_C);
}

// bli_dgemm
void
bli_gemm(trans_t transA, trans_t transB,
         dim_t m, dim_t n, dim_t k,
         double alpha,
         const double *A, inc_t rs_A, inc_t cs_A,
         const double *B, inc_t rs_B, inc_t cs_B,
         double beta,
         double *C, inc_t rs_C, inc_t cs_C)
{
    bli_dgemm(transA, transB,
              m, n, k,
              &alpha,
              (double *)A, rs_A, cs_A,
              (double *)B, rs_B, cs_B,
              &beta,
              C, rs_C, cs_C);
}

// bli_cgemm
void
bli_gemm(trans_t transA, trans_t transB,
         dim_t m, dim_t n, dim_t k,
         const std::complex<float> &alpha,
         const std::complex<float> *A, inc_t rs_A, inc_t cs_A,
         const std::complex<float> *B, inc_t rs_B, inc_t cs_B,
         const std::complex<float> &beta,
         std::complex<float> *C, inc_t rs_C, inc_t cs_C)
{
    bli_cgemm(transA, transB,
              m, n, k,
              (scomplex *)&alpha,
              (scomplex *)A, rs_A, cs_A,
              (scomplex *)B, rs_B, cs_B,
              (scomplex *)&beta,
              (scomplex *)C, rs_C, cs_C);
}

// bli_zgemm
void
bli_gemm(trans_t transA, trans_t transB,
         dim_t m, dim_t n, dim_t k,
         const std::complex<double> &alpha,
         const std::complex<double> *A, inc_t rs_A, inc_t cs_A,
         const std::complex<double> *B, inc_t rs_B, inc_t cs_B,
         const std::complex<double> &beta,
         std::complex<double> *C, inc_t rs_C, inc_t cs_C)
{
    bli_zgemm(transA, transB,
              m, n, k,
              (dcomplex *)&alpha,
              (dcomplex *)const_cast<std::complex<double> *>(A), rs_A, cs_A,
              (dcomplex *)const_cast<std::complex<double> *>(B), rs_B, cs_B,
              (dcomplex *)&beta,
              (dcomplex *)C, rs_C, cs_C);
}

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
     IndexType    incColC)
{
#   ifdef CXXBLAS_DEBUG
    std::cout << "BLIS GEMM" << std::endl;
#   endif

    trans_t   transa = (!transA)
                             ? (!conjA) ? BLIS_NO_TRANSPOSE
                                        : BLIS_CONJ_NO_TRANSPOSE
                             : (!conjA) ? BLIS_TRANSPOSE
                                        : BLIS_CONJ_TRANSPOSE;
    trans_t   transb = (!transB)
                             ? (!conjB) ? BLIS_NO_TRANSPOSE
                                        : BLIS_CONJ_NO_TRANSPOSE
                             : (!conjB) ? BLIS_TRANSPOSE
                                        : BLIS_CONJ_TRANSPOSE;
    dim_t     m      = m_;
    dim_t     n      = n_;
    dim_t     k      = k_;
    inc_t     rs_A   = incRowA;
    inc_t     cs_A   = incColA;
    inc_t     rs_B   = incRowB;
    inc_t     cs_B   = incColB;
    inc_t     rs_C   = incRowC;
    inc_t     cs_C   = incColC;

    err_t     init_result;

    bli_init_auto(&init_result);
    bli_gemm(transa, transb, m, n, k, alpha,
             A, rs_A, cs_A,
             B, rs_B, cs_B,
             beta,
             C, rs_C, cs_C);
    bli_finalize_auto(init_result);
}

} // namespace cxxblas

#endif // ULMBLAS_EXTERNAL_BLIS_GEMM_TCC

