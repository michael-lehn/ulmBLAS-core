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

#ifndef ULMBLAS_IMPL_AUXILIARY_POW_TCC
#define ULMBLAS_IMPL_AUXILIARY_POW_TCC 1

#include <ulmblas/impl/auxiliary/pow.h>

#ifdef WITH_MPFR
#include <external/real.hpp>
#endif

namespace ulmBLAS {

template <typename T>
typename std::enable_if<std::is_same<T, int>::value,
         T>::type
pow(const T &base, const T &exponent)
{
    ASSERT( exponent>=0 );

    T result = T(1);
    for (int i=1; i<=exponent; ++i) {
        result *= base;
    }
    return result;
}

#ifdef WITH_MPFR
template <typename T>
typename std::enable_if<! std::is_same<T, int>::value
                     && IsNotComplex<T>::value
                     && !IsMpfrReal<T>::value,
         T>::type
pow(const T &base, int exponent)
{
//
//  TODO: Make this more general and call an external Fortran routine
//        that computes 'pow(base, exponent)' for comparison
//
#   ifdef CHECK_CXXLAPACK
    if (exponent==2) {
        return base*base;
    }
#   endif
    return std::pow(base, T(exponent));
}
#else
template <typename T>
typename std::enable_if<!std::is_same<T,int>::value
                     && !IsComplex<T>::value,
         T>::type
pow(const T &base, int exponent)
{
//
//  TODO: Make this more general and call an external Fortran routine
//        that computes 'pow(base, exponent)' for comparison
//
#   ifdef CHECK_CXXLAPACK
    if (exponent==2) {
        return base*base;
    }
#   endif
    return std::pow(base, T(exponent));
}
#endif

template <typename T>
std::complex<T>
pow(const std::complex<T> &base, int exponent)
{
//
//  TODO: Make this more general and call an external Fortran routine
//        that computes 'pow(base, exponent)' for comparison
//
#   ifdef CHECK_CXXLAPACK
    if (exponent==2) {
        return base*base;
    }
#   endif
    return std::pow(base, T(exponent));
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_AUXILIARY_POW_TCC
