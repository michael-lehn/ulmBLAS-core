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

#ifndef ULMBLAS_IMPL_LEVEL3_UKERNEL_SSE_UGEMM_H
#define ULMBLAS_IMPL_LEVEL3_UKERNEL_SSE_UGEMM_H 1

#include <type_traits>
#include <ulmblas/impl/level3/ukernel/ref/ugemm.h>

namespace ulmBLAS { namespace sse {

using ref::ugemm;

template <typename IndexType>
typename std::enable_if<std::is_convertible<IndexType, std::int64_t>::value
                     && BlockSize<double>::MR==4
                     && BlockSize<double>::NR==4
                     && BlockSize<double>::align==16,
             void>::type
    ugemm(IndexType      kc_,
          const double   &alpha,
          const double   *A,
          const double   *B,
          const double   &beta,
          double         *C,
          IndexType      incRowC_,
          IndexType      incColC_,
          const double   *nextA,
          const double   *nextB);

} } // namespace sse, ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_UKERNEL_SSE_UGEMM_H

#include <ulmblas/impl/level3/ukernel/sse/ugemm.tcc>
