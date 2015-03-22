/*
 * Copyright (C) 2014, The University of Texas at Austin
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  - Neither the name of The University of Texas at Austin nor the names
 *    of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*
 * Copyright (C) 2014-2015, Michael Lehn
 *
 * ulmBLAS adopted general ideas from BLIS.  Using micro kernels from BLIS
 * only requires minor modifications,
 *
 */

#ifndef ULMBLAS_IMPL_LEVEL2_HELMV_TCC
#define ULMBLAS_IMPL_LEVEL2_HELMV_TCC 1

#include <ulmblas/impl/auxiliary/real.h>
#include <ulmblas/impl/level1/scal.h>
#include <ulmblas/impl/level1/axpy.h>
#include <ulmblas/impl/level1extensions/dotaxpy.h>
#include <ulmblas/impl/level1extensions/dotxaxpyf.h>
#include <ulmblas/impl/level2/helmv.h>

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
void
helmv(IndexType    n,
      const Alpha  &alpha,
      bool         conjA,
      const TA     *A,
      IndexType    incRowA,
      IndexType    incColA,
      const TX     *x,
      IndexType    incX,
      const Beta   &beta,
      TY           *y,
      IndexType    incY)
{
    typedef decltype(Alpha(0)*TA(0)*TX(0)+Beta(0)*TY(0))  T;

    const IndexType    UnitStride(1);
    static const bool  homogeneousTypes = std::is_same<T,Alpha>::value
                                       && std::is_same<T,TA>::value
                                       && std::is_same<T,TX>::value
                                       && std::is_same<T,TY>::value
                                       && false;

    scal(n, beta, y, incY);

    if (homogeneousTypes && incRowA==UnitStride) {
        const IndexType bf = FuseFactor<T>::dotxaxpyf;
        const IndexType nb = (n/bf)*bf;

        T rho, rho_[bf];

        for (IndexType j=0; j<nb; j+=bf) {

            dotxaxpyf(n-j-bf, conjA, !conjA, false,
                      alpha,
                      &x[j*incX], incX,
                      &A[(j+bf)*incRowA+j*incColA], incRowA, incColA,
                      &x[(j+bf)*incX], incX,
                      &y[(j+bf)*incY], incY,
                      rho_, 1);

            for (IndexType l=0; l<bf; ++l) {
                dotaxpy(bf-1-l, conjA, !conjA, false,
                        alpha*x[(j+l)*incX],
                        &A[(j+l+1)*incRowA+(j+l)*incColA], incRowA,
                        &x[(j+l+1)*incX], incX,
                        &y[(j+l+1)*incY], incY,
                        rho);
                y[(j+l)*incY] += alpha*(rho+rho_[l]);
                y[(j+l)*incY] += alpha*real(A[(j+l)*incRowA+(j+l)*incColA])
                                      *x[(j+l)*incX];
            }

        }
        for (IndexType j=nb; j<n; ++j) {
            dotaxpy(n-1-j, conjA, !conjA, false,
                    alpha*x[j*incX],
                    &A[(j+1)*incRowA+j*incColA], incRowA,
                    &x[(j+1)*incX], incX,
                    &y[(j+1)*incY], incY,
                    rho);
            y[j*incY] += alpha*(real(A[j*incRowA+j*incColA])*x[j*incX]+rho);
        }
    } else if (homogeneousTypes && incColA==UnitStride) {
        const IndexType bf = FuseFactor<T>::dotxaxpyf;
        const IndexType nb = (n/bf)*bf;

        T rho, rho_[bf];

        for (IndexType i=0; i<nb; i+=bf) {

            dotxaxpyf(i, conjA, !conjA, false,
                      alpha, &x[i*incX], incX,
                      &A[i*incRowA], incColA, incRowA,
                      &x[0*incX], incX,
                      &y[0*incY], incY,
                      rho_, 1);

            for (IndexType l=0; l<bf; ++l) {
                dotaxpy(l, conjA, !conjA, false,
                        alpha*x[(i+l)*incX],
                        &A[(i+l)*incRowA+i*incColA], incColA,
                        &x[i*incX], incX,
                        &y[i*incY], incY,
                        rho);
                y[(i+l)*incY] += alpha*(rho+rho_[l]);
                y[(i+l)*incY] += alpha*real(A[(i+l)*incRowA+(i+l)*incColA])
                                      *x[(i+l)*incX];
            }
        }
        for (IndexType i=nb; i<n; ++i) {
            dotaxpy(i, conjA, !conjA, false,
                    alpha*x[i*incX],
                    &A[i*incRowA], incColA,
                    &x[0*incX], incX,
                    &y[0*incY], incY,
                    rho);
            y[i*incY] += alpha*(real(A[i*incRowA+i*incColA])*x[i*incX]+rho);
        }
    } else {
//
//      Otherwise we just use a variant with axpy as reference implementation.
//      TODO: packing, switching between dot/axpy variant depending on
//            abs(incRowA) and abs(incColA)
//
        if (!conjA) {
            for (IndexType i=0; i<n; ++i) {
                acxpy(i, alpha*x[i*incX],
                      &A[i*incRowA], incColA,
                      &y[0*incY], incY);
                y[i*incY] += alpha*real(A[i*incRowA+i*incColA])*x[i*incX];
                axpy(n-1-i, alpha*x[i*incX],
                     &A[(i+1)*incRowA+i*incColA], incRowA,
                     &y[(i+1)*incY], incY);
            }
        } else {
            for (IndexType i=0; i<n; ++i) {
                axpy(i, alpha*x[i*incX],
                     &A[i*incRowA], incColA,
                     &y[0*incY], incY);
                y[i*incY] += alpha*real(A[i*incRowA+i*incColA])*x[i*incX];
                acxpy(n-1-i, alpha*x[i*incX],
                      &A[(i+1)*incRowA+i*incColA], incRowA,
                      &y[(i+1)*incY], incY);
            }
         }
    }
}

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
void
helmv(IndexType    n,
      const Alpha  &alpha,
      const TA     *A,
      IndexType    incRowA,
      IndexType    incColA,
      const TX     *x,
      IndexType    incX,
      const Beta   &beta,
      TY           *y,
      IndexType    incY)
{
    helmv(n, alpha, false, A, incRowA, incColA, x, incX, beta, y, incY);
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL2_HELMV_TCC
