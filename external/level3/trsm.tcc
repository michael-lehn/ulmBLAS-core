/*
 *   Copyright (c) 2009, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ULMBLAS_EXTERNAL_LEVEL3_TRSM_TCC
#define ULMBLAS_EXTERNAL_LEVEL3_TRSM_TCC 1

namespace cblas {

// strsm
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
trsm(bool colMajor,
     bool leftA, bool lowerA, bool transA, bool conjA, bool unitDiagA,
     IndexType m, IndexType n,
     float alpha,
     const float *A, IndexType ldA,
     float *B, IndexType ldB)
{
    CBLAS_DEBUG_OUT("[" BLAS_IMPL "]: cblas_strsm");

    cblas_strsm(Blas<>::order(colMajor),
                Blas<>::side(leftA),
                Blas<>::uplo(lowerA),
                Blas<>::trans(transA, conjA),
                Blas<>::diag(unitDiagA),
                m, n,
                alpha,
                A, ldA,
                B, ldB);
}

// dtrsm
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
trsm(bool colMajor,
     bool leftA, bool lowerA, bool transA, bool conjA, bool unitDiagA,
     IndexType m, IndexType n,
     double alpha,
     const double *A, IndexType ldA,
     double *B, IndexType ldB)
{
    CBLAS_DEBUG_OUT("[" BLAS_IMPL "]: cblas_dtrsm");

    cblas_dtrsm(Blas<>::order(colMajor),
                Blas<>::side(leftA),
                Blas<>::uplo(lowerA),
                Blas<>::trans(transA, conjA),
                Blas<>::diag(unitDiagA),
                m, n,
                alpha,
                A, ldA,
                B, ldB);
}

// ctrsm
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
trsm(bool colMajor,
     bool leftA, bool lowerA, bool transA, bool conjA, bool unitDiagA,
     IndexType m, IndexType n,
     const ComplexFloat &alpha,
     const ComplexFloat *A, IndexType ldA,
     ComplexFloat *B, IndexType ldB)
{
    CBLAS_DEBUG_OUT("[" BLAS_IMPL "]: cblas_ctrsm");

    cblas_ctrsm(Blas<>::order(colMajor),
                Blas<>::side(leftA),
                Blas<>::uplo(lowerA),
                Blas<>::trans(transA, conjA),
                Blas<>::diag(unitDiagA),
                m, n,
                reinterpret_cast<const float *>(&alpha),
                reinterpret_cast<const float *>(A), ldA,
                reinterpret_cast<float *>(B), ldB);
}

// ztrsm
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
trsm(bool colMajor,
     bool leftA, bool lowerA, bool transA, bool conjA, bool unitDiagA,
     IndexType m, IndexType n,
     const ComplexDouble &alpha,
     const ComplexDouble *A, IndexType ldA,
     ComplexDouble *B, IndexType ldB)
{
    CBLAS_DEBUG_OUT("[" BLAS_IMPL "]: cblas_ztrsm");

    cblas_ztrsm(Blas<>::order(colMajor),
                Blas<>::side(leftA),
                Blas<>::uplo(lowerA),
                Blas<>::trans(transA, conjA),
                Blas<>::diag(unitDiagA),
                m, n,
                reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double *>(A), ldA,
                reinterpret_cast<double *>(B), ldB);
}

} // namespace cblas

#endif // ULMBLAS_EXTERNAL_LEVEL3_TRSM_TCC 1
