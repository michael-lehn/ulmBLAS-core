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

#ifndef ULMBLAS_IMPL_LEVEL3_UKERNEL_UTRUSM_TCC
#define ULMBLAS_IMPL_LEVEL3_UKERNEL_UTRUSM_TCC 1

#include <ulmblas/impl/level1/scal.h>
#include <ulmblas/impl/level1extensions/gecopy.h>
#include <ulmblas/impl/level3/ukernel/ugemm.h>
#include <ulmblas/impl/level3/ukernel/utrusm.h>

//
//  Selected optimized micro kernel
//
#if defined(USE_SSE)
#   define  SELECT_UTRUSM_KERNEL    sse
#   include <ulmblas/impl/level3/ukernel/sse/utrusm.h>
#else
#   define  SELECT_UTRUSM_KERNEL    ref
#   include <ulmblas/impl/level3/ukernel/ref/utrusm.h>
#endif


namespace ulmBLAS {

//
//  Buffered variant.  Used for zero padded panels.
//
template <typename IndexType, typename T, typename TC>
void
utrusm(IndexType    mr,
       IndexType    nr,
       const T      *A,
       const T      *B,
       TC           *C,
       IndexType    incRowC,
       IndexType    incColC)
{
    const IndexType MR = BlockSize<T>::MR;
    const IndexType NR = BlockSize<T>::NR;

    T   A_[MR*MR];
    T   B_[MR*NR];
    T   C_[MR*NR];

    scal(MR*MR, T(0), A_, IndexType(1));
    scal(MR*NR, T(0), B_, IndexType(1));

    gecopy(mr, mr, false, A, IndexType(1), MR, A_, IndexType(1), MR);
    gecopy(mr, nr, false, B, NR, IndexType(1), B_, NR, IndexType(1));

    utrusm(A_, B_, C_, IndexType(1), MR);
    gecopy(mr, nr, false, C_, IndexType(1), MR, C, incRowC, incColC);
}

//
//  Buffered variant.  Used if the result A^(-1)*B needs to be upcasted for
//  computing C <- A^(-1)*B
//
template <typename T, typename TC, typename IndexType>
void
utrusm(const T     *A,
       const T     *B,
       TC          *C,
       IndexType   incRowC,
       IndexType   incColC)
{
    const IndexType MR = BlockSize<T>::MR;
    const IndexType NR = BlockSize<T>::NR;

    utrusm(MR, NR, A, B, C, incRowC, incColC);
}

//
//  Unbuffered variant.
//
template <typename IndexType, typename T>
void
utrusm(const T     *A,
       const T     *B,
       T           *C,
       IndexType   incRowC,
       IndexType   incColC)
{
    SELECT_UTRUSM_KERNEL::utrusm(A, B, C, incRowC, incColC);
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_UKERNEL_UTRUSM_TCC
