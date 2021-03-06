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

#ifndef ULMBLAS_IMPL_LEVEL2_TBLSV_TCC
#define ULMBLAS_IMPL_LEVEL2_TBLSV_TCC 1

#include <ulmblas/impl/auxiliary/conjugate.h>
#include <ulmblas/impl/level1/axpy.h>
#include <ulmblas/impl/level2/tblsv.h>

namespace ulmBLAS {

template <typename IndexType, typename TA, typename TX>
void
tblsv(IndexType    n,
      IndexType    k,
      bool         unitDiag,
      bool         conjA,
      const TA     *A,
      IndexType    ldA,
      TX           *x,
      IndexType    incX)
{
    if (n==0) {
        return;
    }

    if (!conjA) {
        for (IndexType j=0; j<n; ++j) {
            IndexType i1  = k - std::max(IndexType(0), (j+1+k)-n);
            IndexType len = std::max(IndexType(0), i1+1);

            if (!unitDiag) {
                x[j*incX] /= A[0];
            }
            axpy(len-1, -x[j*incX], &A[1], IndexType(1), &x[(j+1)*incX], incX);
            A += ldA;
        }
    } else {
        for (IndexType j=0; j<n; ++j) {
            IndexType i1  = k - std::max(IndexType(0), (j+1+k)-n);
            IndexType len = std::max(IndexType(0), i1+1);

            if (!unitDiag) {
                x[j*incX] /= conjugate(A[0]);
            }
            acxpy(len-1, -x[j*incX], &A[1], IndexType(1), &x[(j+1)*incX], incX);
            A += ldA;
        }
    }
}

template <typename IndexType, typename TA, typename TX>
void
tblsv(IndexType    n,
      IndexType    k,
      bool         unitDiag,
      const TA     *A,
      IndexType    ldA,
      TX           *x,
      IndexType    incX)
{
    tblsv(n, k, unitDiag, false, A, ldA, x, incX);
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL2_TBLSV_TCC
