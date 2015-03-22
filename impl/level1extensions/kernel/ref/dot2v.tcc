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

#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_DOT2V_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_DOT2V_TCC 1

#include <ulmblas/impl/level1extensions/kernel/ref/dot2v.h>

namespace ulmBLAS { namespace ref {

template <typename IndexType, typename VX0, typename VX1, typename VY,
          typename Result>
void
dotu2v(IndexType      n,
       const VX0      *x0,
       IndexType      incX0,
       const VX1      *x1,
       IndexType      incX1,
       VY             *y,
       IndexType      incY,
       Result         *result,
       IndexType      resultInc)
{
    Result &result0 = result[0*resultInc];
    Result &result1 = result[1*resultInc];

    result0 = result1 = Result(0);

    for (IndexType i=0; i<n; ++i) {
        result0 += x0[i*incX0]*y[i*incY];
        result1 += x1[i*incX1]*y[i*incY];
    }
}

} } // namespace ref, ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_DOT2V_TCC 1
