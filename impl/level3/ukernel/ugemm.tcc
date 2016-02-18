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

#ifndef ULMBLAS_IMPL_LEVEL3_UKERNEL_UGEMM_TCC
#define ULMBLAS_IMPL_LEVEL3_UKERNEL_UGEMM_TCC 1

#include <ulmblas/impl/auxiliary/printmatrix.h>
#include <ulmblas/impl/level1extensions/geaxpy.h>
#include <ulmblas/impl/level1extensions/gescal.h>
#include <ulmblas/impl/level3/ukernel/ugemm.h>
#include <algorithm>

namespace ulmBLAS {

//
//  Buffered variant.  Used for zero padded panels.
//
template <typename IndexType, typename T, typename Beta, typename TC>
void
ugemm(IndexType    mr,
      IndexType    nr,
      IndexType    kc,
      const T      &alpha,
      const T      *A,
      const T      *B,
      const Beta   &beta,
      TC           *C,
      IndexType    incRowC,
      IndexType    incColC,
      const T      *nextA,
      const T      *nextB)
{
    const IndexType MR = BlockSize<T>::MR;
    const IndexType NR = BlockSize<T>::NR;

    T   C_[MR*NR];

    // make sure there is no NaN in buffer
    std::fill_n(C_, MR*NR, T(0));

    ugemm(kc, alpha, A, B, T(0), C_, IndexType(1), MR, nextA, nextB);
    gescal(mr, nr, beta, C, incRowC, incColC);
    geaxpy(mr, nr, T(1), C_, IndexType(1), MR, C, incRowC, incColC);
}

//
//  Buffered variant.  Used if the result alpha*A*B needs to be upcasted for
//  computing C <- beta*C + (alpha*A*B)
//
template <typename IndexType, typename T, typename Beta, typename TC>
void
ugemm(IndexType   kc,
      const T     &alpha,
      const T     *A,
      const T     *B,
      const Beta  &beta,
      TC          *C,
      IndexType   incRowC,
      IndexType   incColC,
      const T     *nextA,
      const T     *nextB)
{
    const IndexType MR = BlockSize<T>::MR;
    const IndexType NR = BlockSize<T>::NR;

    ugemm(MR, NR, kc, alpha, A, B, beta, C, incRowC, incColC, nextA, nextB);
}

//
//  Unbuffered variant.
//
template <typename IndexType, typename T>
void
ugemm(IndexType   kc,
      const T     &alpha,
      const T     *A,
      const T     *B,
      const T     &beta,
      T           *C,
      IndexType   incRowC,
      IndexType   incColC,
      const T     *nextA,
      const T     *nextB)
{
    SELECT_UGEMM_KERNEL::ugemm(kc, alpha, A, B, beta, C, incRowC, incColC,
                               nextA, nextB);
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_UKERNEL_UGEMM_TCC
