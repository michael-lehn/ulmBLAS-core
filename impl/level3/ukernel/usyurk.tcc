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

#ifndef ULMBLAS_IMPL_LEVEL3_UKERNEL_USYURK_TCC
#define ULMBLAS_IMPL_LEVEL3_UKERNEL_USYURK_TCC 1

#include <ulmblas/impl/level1extensions/truaxpy.h>
#include <ulmblas/impl/level1extensions/truscal.h>
#include <ulmblas/impl/level3/ukernel/ugemm.h>
#include <ulmblas/impl/level3/ukernel/usyurk.h>

namespace ulmBLAS {

template <typename IndexType, typename T, typename Beta, typename TC>
void
usyurk(IndexType    mr,
       IndexType    nr,
       IndexType    kc,
       IndexType    ic,
       IndexType    jc,
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
    const IndexType MR = BlockSizeUGemm<T>::MR;
    const IndexType NR = BlockSizeUGemm<T>::NR;

    T   C_[MR*NR];

    ugemm(kc, alpha, A, B, T(0), C_, IndexType(1), MR, nextA, nextB);

    if (jc>ic) {
        gescal(jc-ic, nr, beta, C, incRowC, incColC);
        geaxpy(jc-ic, nr, Beta(1), C_, IndexType(1), MR, C, incRowC, incColC);
        truscal(mr-(jc-ic), nr, false, beta,
                &C[(jc-ic)*incRowC], incRowC, incColC);
        truaxpy(mr-(jc-ic), nr, false, Beta(1),
                &C_[jc-ic], IndexType(1), MR,
                &C[(jc-ic)*incRowC], incRowC, incColC);
    } else {
        truscal(mr, nr-(ic-jc), false, beta,
                &C[(ic-jc)*incColC], incRowC, incColC);
        truaxpy(mr, nr-(ic-jc), false, Beta(1),
                &C_[(ic-jc)*MR], IndexType(1), MR,
                &C[(ic-jc)*incColC], incRowC, incColC);
   }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_UKERNEL_USYURK_TCC
