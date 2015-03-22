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

#ifndef ULMBLAS_IMPL_LEVEL1_KERNEL_SSE_DOT_TCC
#define ULMBLAS_IMPL_LEVEL1_KERNEL_SSE_DOT_TCC 1

#include <immintrin.h>

#include <ulmblas/impl/auxiliary/isaligned.h>
#include <ulmblas/impl/level1/kernel/ref/dot.h>
#include <ulmblas/impl/level1/kernel/sse/dot.h>

namespace ulmBLAS { namespace sse {

//
// ----------------
// Double Precision
// ----------------
//

template <typename IndexType>
void
dotu(IndexType      n,
     const double   *x,
     IndexType      incX,
     const double   *y,
     IndexType      incY,
     double         &result)
{
    if (n<=0) {
        result = 0;
        return;
    }

    if (incX!=1 || incY!=1) {
        ref::dotu(n, x, incX, y, incY, result);
        return;
    }

    bool xAligned = isAligned(x, 16);
    bool yAligned = isAligned(y, 16);

    double result_ = 0;

    if (!xAligned && !yAligned) {
        result_ = y[0]*x[0];
        ++x;
        ++y;
        --n;
        xAligned = yAligned = true;
    }
    if (xAligned && yAligned) {
        IndexType nb = n / 2;
        IndexType nl = n % 2;

        __m128d x12, y12, result12;
        double  result12_[2];

        result12 = _mm_setzero_pd();
        for (IndexType i=0; i<nb; ++i) {
            x12 = _mm_load_pd(x);
            y12 = _mm_load_pd(y);

            x12      = x12*y12;
            result12 = result12 + x12;

            x += 2;
            y += 2;
        }
        _mm_store_pd(result12_, result12);

        result = result_ + result12_[0] + result12_[1];

        for (IndexType i=0; i<nl; ++i) {
            result += x[i]*y[i];
        }

    } else {
        ref::dotu(n, x, incX, y, incY, result);
    }
}

} } // namespace ref, ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1_KERNEL_SSE_DOT_TCC 1
