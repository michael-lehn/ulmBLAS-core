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

#ifndef ULMBLAS_IMPL_LEVEL2_GBMV_TCC
#define ULMBLAS_IMPL_LEVEL2_GBMV_TCC 1

#include <ulmblas/impl/level1/axpy.h>
#include <ulmblas/impl/level1/scal.h>
#include <ulmblas/impl/level2/gbmv.h>

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
void
gbmv(IndexType    m,
     IndexType    n,
     IndexType    kl,
     IndexType    ku,
     const Alpha  &alpha,
     bool         conjA,
     const TA     *A,
     IndexType    ldA,
     const TX     *x,
     IndexType    incX,
     const Beta   &beta,
     TY           *y,
     IndexType    incY)
{
    if (m==0 || n==0 || (alpha==Alpha(0) && beta==Beta(1))) {
        return;
    }

    scal(m, beta, y, incY);

    if (alpha==Alpha(0)) {
        return;
    }

    if (!conjA) {
        for (IndexType j=0; j<n; ++j) {
            IndexType i0  = std::max(IndexType(0), ku-j);
            IndexType i1  = ku+1+kl - std::max(IndexType(0), (j+1+kl)-m);
            IndexType len = std::max(IndexType(0), i1-i0);

            IndexType iY  = std::max(IndexType(0), j-ku);

            axpy(len, alpha*x[j*incX],
                 &A[i0], IndexType(1),
                 &y[iY*incY], incY);
            A += ldA;
        }
    } else {
        for (IndexType j=0; j<n; ++j) {
            IndexType i0  = std::max(IndexType(0), ku-j);
            IndexType i1  = ku+1+kl - std::max(IndexType(0), (j+1+kl)-m);
            IndexType len = std::max(IndexType(0), i1-i0);

            IndexType iY  = std::max(IndexType(0), j-ku);

            acxpy(len, alpha*x[j*incX],
                  &A[i0], IndexType(1),
                  &y[iY*incY], incY);
            A += ldA;
        }
    }
}

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
void
gbmv(IndexType    m,
     IndexType    n,
     IndexType    kl,
     IndexType    ku,
     const Alpha  &alpha,
     const TA     *A,
     IndexType    ldA,
     const TX     *x,
     IndexType    incX,
     const Beta   &beta,
     TY           *y,
     IndexType    incY)
{
    gbmv(m, n, kl, ku, alpha, false,  A, ldA, x, incX, beta, y, incY);
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL2_GBMV_TCC
