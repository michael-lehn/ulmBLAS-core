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

#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_SSE_AXPY2V_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_SSE_AXPY2V_TCC 1

#include <immintrin.h>
#include <ulmblas/impl/auxiliary/isaligned.h>
#include <ulmblas/impl/level1extensions/kernel/ref/axpy2v.h>
#include <ulmblas/impl/level1extensions/kernel/sse/axpy2v.h>

namespace ulmBLAS { namespace sse {

//
// ----------------
// Double Precision
// ----------------
//

template <typename IndexType>
void
axpy2v(IndexType      n,
       const double   &alpha0,
       const double   &alpha1,
       const double   *x0,
       IndexType      incX0,
       const double   *x1,
       IndexType      incX1,
       double         *y,
       IndexType      incY)
{
    if (n<=0) {
        return;
    }

    if (incX0!=1 || incX1!=1 || incY!=1) {
        ref::axpy2v(n, alpha0, alpha1, x0, incX0, x1, incX1, y, incY);
        return;
    }

    bool x0Aligned = isAligned(x0, 16);
    bool x1Aligned = isAligned(x1, 16);
    bool yAligned  = isAligned(y, 16);

    if (!x0Aligned && !x1Aligned && !yAligned) {
        y[0] += alpha0*x0[0];
        y[0] += alpha1*x1[0];
        ++x0;
        ++x1;
        ++y;
        --n;
        x0Aligned = x1Aligned = yAligned = true;
    }
    if (x0Aligned && x1Aligned && yAligned) {
        IndexType nb = n / 6;
        IndexType nl = n % 6;

        __m128d y12, y34, y56;
        __m128d alpha0_11, x0_12, x0_34, x0_56;
        __m128d alpha1_11, x1_12, x1_34, x1_56;

        alpha0_11 = _mm_loaddup_pd(&alpha0);
        alpha1_11 = _mm_loaddup_pd(&alpha1);

        for (IndexType i=0; i<nb; ++i) {
            x0_12 = _mm_load_pd(x0);
            x1_12 = _mm_load_pd(x1);
            y12   = _mm_load_pd(y);

            x0_12 = x0_12 * alpha0_11;
            y12  += x0_12;
            x1_12 = x1_12 * alpha1_11;
            y12  += x1_12;
            _mm_store_pd(y, y12);

            x0_34 = _mm_load_pd(x0+2);
            x1_34 = _mm_load_pd(x1+2);
            y34   = _mm_load_pd(y +2);

            x0_34 = x0_34 * alpha0_11;
            y34  += x0_34;
            x1_34 = x1_34 * alpha1_11;
            y34  += x1_34;
            _mm_store_pd(y+2, y34);

            x0_56 = _mm_load_pd(x0+4);
            x1_56 = _mm_load_pd(x1+4);
            y56   = _mm_load_pd(y +4);

            x0_56 = x0_56 * alpha0_11;
            y56  += x0_56;
            x1_56 = x1_56 * alpha1_11;
            y56  += x1_56;
            _mm_store_pd(y+4, y56);

            x0 += 6;
            x1 += 6;
            y  += 6;
        }
        for (IndexType i=0; i<nl; ++i) {
            y[i] += alpha0*x0[i] + alpha1*x1[i];
        }
    } else {
        ref::axpy2v(n, alpha0, alpha1, x0, incX0, x1, incX1, y, incY);
    }
}

} } // namespace sse, ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_SSE_AXPY2V_TCC 1
