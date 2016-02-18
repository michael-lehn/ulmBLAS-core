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

#ifndef ULMBLAS_IMPL_LEVEL3_UKERNEL_UGEMM_H
#define ULMBLAS_IMPL_LEVEL3_UKERNEL_UGEMM_H 1

//
//  Selected optimized micro kernel
//
#if defined(USE_SSE)
#   define  SELECT_UGEMM_KERNEL     sse
#   include <ulmblas/impl/level3/ukernel/sse/ugemm.h>
#elif defined (USE_AVX)
#   define  SELECT_UGEMM_KERNEL     avx
#   include <ulmblas/impl/level3/ukernel/avx/ugemm.h>
#elif defined (USE_FMA)
#   define  SELECT_UGEMM_KERNEL     fma
#   include <ulmblas/impl/level3/ukernel/fma/ugemm.h>
#else
#   define  SELECT_UGEMM_KERNEL     ref
#   include <ulmblas/impl/level3/ukernel/ref/ugemm.h>
#endif

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
          const T      *nextB);

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
          const T     *nextB);

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
          const T     *nextB);

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_UKERNEL_UGEMM_H

#include <ulmblas/impl/level3/ukernel/ugemm.tcc>
