/*
 * Copyright (C) 2014, The University of Texas at Austin
 * Copyright (C) 2014-2015, Michael Lehn
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

#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_SSE_DOTXAXPYF_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_SSE_DOTXAXPYF_TCC 1

#include <immintrin.h>
#include <ulmblas/impl/auxiliary/isaligned.h>
#include <ulmblas/impl/config/fusefactor.h>
#include <ulmblas/impl/level1extensions/kernel/sse/dotxaxpyf.h>
#include <ulmblas/impl/level1extensions/kernel/ref/dotxaxpyf.h>

namespace ulmBLAS { namespace sse {

//
// ----------------
// Double Precision
// ----------------
//

template <typename IndexType>
typename std::enable_if<std::is_integral<IndexType>::value
                     && FuseFactor<double>::dotuxf==2,
void>::type
dotxaxpyf(IndexType      n,
          bool           conjX,
          bool           conjXt,
          bool           conjY,
          const double   &alpha,
          const double   *a,
          IndexType      incA,
          const double   *X,
          IndexType      incRowX,
          IndexType      incColX,
          const double   *y,
          IndexType      incY,
          double         *z,
          IndexType      incZ,
          double         *rho,
          IndexType      incRho)
{
    const IndexType bf = FuseFactor<double>::dotxaxpyf;

    for (IndexType l=0; l<bf; ++l) {
        rho[l*incRho] = 0;
    }

    if (n<=0) {
        return;
    }

    if (incRowX!=1 || incY!=1 || incZ!=1 || conjX || conjXt || conjY) {
        ref::dotxaxpyf(n, conjX, conjXt, conjY, alpha, a, incA,
                       X, incRowX, incColX, y, incY,
                       z, incZ, rho, incRho);
        return;
    }

    const double alpha0 = alpha*a[0*incA];
    const double alpha1 = alpha*a[1*incA];

    const double *x0    = &X[0*incColX];
    const double *x1    = &X[1*incColX];

    double &rho0        = rho[0*incRho];
    double &rho1        = rho[1*incRho];

    bool x0Aligned      = isAligned(x0, 16);
    bool x1Aligned      = isAligned(x1, 16);
    bool yAligned       = isAligned(y, 16);
    bool zAligned       = isAligned(z, 16);

    rho0 = rho1 = 0;

    if (!x0Aligned && !x1Aligned && !yAligned && !zAligned) {
        z[0] += alpha0*x0[0] + alpha1*x1[0];
        rho0 += x0[0]*y[0];
        rho1 += x1[0]*y[0];

        ++x0;
        ++x1;
        ++y;
        ++z;
        --n;
        x0Aligned = x1Aligned = yAligned = zAligned = true;
    }
    if (x0Aligned && x1Aligned && yAligned && zAligned) {
        IndexType nb = n / 4;
        IndexType nl = n % 4;

        __m128d rho00, rho11;
        __m128d alpha0_11, alpha1_11;
        __m128d z12, z34;
        __m128d y12, y34;
        __m128d tmp0, tmp1;
        __m128d x0_12, x0_34;
        __m128d x1_12, x1_34;

        alpha0_11 = _mm_loaddup_pd(&alpha0);
        alpha1_11 = _mm_loaddup_pd(&alpha1);

        rho00     = _mm_setzero_pd();
        rho11     = _mm_setzero_pd();

        double rho00_[2], rho11_[2];

        for (IndexType i=0; i<nb; ++i) {
            x0_12 = _mm_load_pd(x0);
            x0_34 = _mm_load_pd(x0+2);

            x1_12 = _mm_load_pd(x1);
            x1_34 = _mm_load_pd(x1+2);

            y12   = _mm_load_pd(y);
            y34   = _mm_load_pd(y+2);

            z12   = _mm_load_pd(z);
            z34   = _mm_load_pd(z+2);

            tmp0  = y12;
            tmp1  = y12;

            tmp0  = tmp0 * x0_12;
            tmp1  = tmp1 * x1_12;
            x0_12 = alpha0_11 * x0_12;
            x1_12 = alpha1_11 * x1_12;
            z12   = z12 + x0_12;
            z12   = z12 + x1_12;
            rho00 = rho00 + tmp0;
            rho11 = rho11 + tmp1;
            _mm_store_pd(z, z12);

            tmp0  = y34;
            tmp1  = y34;

            tmp0  = tmp0 * x0_34;
            tmp1  = tmp1 * x1_34;
            x0_34 = alpha0_11 * x0_34;
            x1_34 = alpha1_11 * x1_34;
            z34   = z34 + x0_34;
            z34   = z34 + x1_34;
            rho00 = rho00 + tmp0;
            rho11 = rho11 + tmp1;
            _mm_store_pd(z+2, z34);

            x0 += 4;
            x1 += 4;
            y  += 4;
            z  += 4;
        }
        _mm_store_pd(rho00_, rho00);
        _mm_store_pd(rho11_, rho11);

        rho0 += rho00_[0] + rho00_[1];
        rho1 += rho11_[0] + rho11_[1];

        for (IndexType i=0; i<nl; ++i) {
            z[i] += alpha0*x0[i] + alpha1*x1[i];
            rho0 += x0[i]*y[i];
            rho1 += x1[i]*y[i];
        }

    } else {
        ref::dotxaxpyf(n, conjX, conjXt, conjY, alpha, a, incA,
                       X, incRowX, incColX, y, incY,
                       z, incZ, rho, incRho);
    }
}

} } // namespace sse, ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_SSE_DOTXAXPYF_TCC 1
