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

#ifndef ULMBLAS_IMPL_CONFIG_FUSEFACTOR_H
#define ULMBLAS_IMPL_CONFIG_FUSEFACTOR_H 1

#include <complex>
#include <type_traits>

#include <ulmblas/impl/config/simd.h>

namespace ulmBLAS {

#if defined(USE_TESTPARAM)

#define DAXPYF_FUSEFACTOR   4
#define ZAXPYF_FUSEFACTOR   4
#define DDOTUXF_FUSEFACTOR  4
#define ZDOTUXF_FUSEFACTOR  4

template <typename T>
struct FuseFactor
{
    typedef std::complex<float>   fcomplex;
    typedef std::complex<double>  dcomplex;

    static const int axpyf = std::is_same<T,double>::value   ? DAXPYF_FUSEFACTOR
                           : std::is_same<T,dcomplex>::value ? ZAXPYF_FUSEFACTOR
                           : 1;

    static const int acxpyf = axpyf;

    static const int dotuxf = std::is_same<T,double>::value ? DDOTUXF_FUSEFACTOR
                           : std::is_same<T,dcomplex>::value? ZDOTUXF_FUSEFACTOR
                           : 1;

    static const int dotcxf = dotuxf;

    static const int dotxaxpyf = axpyf;
};


#elif defined(USE_SSE)

template <typename T>
struct FuseFactor
{
    typedef std::complex<float>   fcomplex;
    typedef std::complex<double>  dcomplex;

    static const int axpyf = std::is_same<T,double>::value   ? 2
                           : std::is_same<T,dcomplex>::value ? 4
                           : 1;

    static const int acxpyf = axpyf;

    static const int dotuxf = std::is_same<T,double>::value ? 4
                           : std::is_same<T,dcomplex>::value? 4
                           : 1;
    static const int dotcxf = dotuxf;

    static const int dotxaxpyf = axpyf;
};

#else

template <typename T>
struct FuseFactor
{
    typedef std::complex<float>   fcomplex;
    typedef std::complex<double>  dcomplex;

    static const int axpyf = std::is_same<T,double>::value   ? 4
                           : std::is_same<T,dcomplex>::value ? 4
                           : 1;

    static const int acxpyf = axpyf;

    static const int dotuxf = std::is_same<T,double>::value   ? 4
                            : std::is_same<T,dcomplex>::value ? 4
                            : 1;

    static const int dotxaxpyf = axpyf;
};


#endif



} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_CONFIG_BLOCKSIZE_H 1
