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

#ifndef ULMBLAS_IMPL_LEVEL2_HPLR2_TCC
#define ULMBLAS_IMPL_LEVEL2_HPLR2_TCC 1

#include <ulmblas/impl/auxiliary/conjugate.h>
#include <ulmblas/impl/auxiliary/real.h>
#include <ulmblas/impl/level1extensions/axpy2v.h>
#include <ulmblas/impl/level2/hplr2.h>

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename TX, typename TY,
          typename TA>
void
hplr2(IndexType    n,
      bool         conj,
      const Alpha  &alpha,
      const TX     *x,
      IndexType    incX,
      const TY     *y,
      IndexType    incY,
      TA           *A)
{
    if (n==0 || alpha==Alpha(0)) {
        return;
    }

    if (!conj) {
        for (IndexType j=0; j<n; ++j) {
            axpy2v(n-j,
                   alpha*conjugate(y[j*incY]),
                   conjugate(alpha*x[j*incX]),
                   &x[j*incX], incX,
                   &y[j*incY], incY,
                   A, IndexType(1));
            A[0] = real(A[0]);
            A += n-j;
        }
    } else {
        for (IndexType j=0; j<n; ++j) {
            acxpy(n-j,
                  conjugate(alpha)*y[j*incY],
                  &x[j*incX], incX,
                  A, IndexType(1));
            acxpy(n-j,
                  alpha*x[j*incX],
                  &y[j*incY], incY,
                  A, IndexType(1));
            A[0] = real(A[0]);
            A += n-j;
        }
    }
}

template <typename IndexType, typename Alpha, typename TX, typename TY,
          typename TA>
void
hplr2(IndexType    n,
      const Alpha  &alpha,
      const TX     *x,
      IndexType    incX,
      const TY     *y,
      IndexType    incY,
      TA           *A)
{
    hplr2(n, false, alpha, x, incX, y, incY, A);
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL2_HPLR2_TCC
