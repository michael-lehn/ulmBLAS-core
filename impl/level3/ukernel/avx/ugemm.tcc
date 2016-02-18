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

#ifndef ULMBLAS_IMPL_LEVEL3_UKERNEL_AVX_UGEMM_TCC
#define ULMBLAS_IMPL_LEVEL3_UKERNEL_AVX_UGEMM_TCC 1

#include <type_traits>
#include <ulmblas/impl/level3/ukernel/ref/ugemm.h>

namespace ulmBLAS { namespace avx {

template <typename Index>
typename std::enable_if<std::is_convertible<Index, std::int64_t>::value
                     && BlockSize<double>::MR==4
                     && BlockSize<double>::NR==8
                     && BlockSize<double>::align==32,
         void>::type
ugemm(Index kc_, double alpha,
      const double *A, const double *B,
      double beta,
      double *C, Index incRowC_, Index incColC_,
      const double *, const double *)
{
    int64_t kc      = kc_;
    int64_t incRowC = incRowC_;
    int64_t incColC = incColC_;

    double *pAlpha  = &alpha;
    double *pBeta   = &beta;

//
//  Compute AB = A*B
//
    __asm__ volatile
    (
    "movq      %0,           %%rdi    \n\t"  // kc
    "movq      %1,           %%rsi    \n\t"  // A
    "movq      %2,           %%rdx    \n\t"  // B
    "movq      %5,           %%rcx    \n\t"  // C
    "movq      %6,           %%r8     \n\t"  // incRowC
    "movq      %7,           %%r9     \n\t"  // incColC

    "vmovapd           0 * 32(%%rdx),         %%ymm4\n\t"

    "vbroadcastsd       0 * 8(%%rsi),         %%ymm0\n\t"
    "vbroadcastsd       1 * 8(%%rsi),         %%ymm1\n\t"
    "vbroadcastsd       2 * 8(%%rsi),         %%ymm2\n\t"
    "vbroadcastsd       3 * 8(%%rsi),         %%ymm3\n\t"

    "vxorpd                  %%ymm8,          %%ymm8,          %%ymm8\n\t"
    "vxorpd                  %%ymm9,          %%ymm9,          %%ymm9\n\t"
    "vxorpd                  %%ymm10,         %%ymm10,         %%ymm10\n\t"
    "vxorpd                  %%ymm11,         %%ymm11,         %%ymm11\n\t"
    "vxorpd                  %%ymm12,         %%ymm12,         %%ymm12\n\t"
    "vxorpd                  %%ymm13,         %%ymm13,         %%ymm13\n\t"
    "vxorpd                  %%ymm14,         %%ymm14,         %%ymm14\n\t"
    "vxorpd                  %%ymm15,         %%ymm15,         %%ymm15\n\t"

    "jmp                     check%=\n\t"

    "loop%=:\n\t"

    "vmovapd           1 * 32(%%rdx),         %%ymm5\n\t"

    "vmulpd                  %%ymm0,          %%ymm4,          %%ymm6\n\t"
    "vaddpd                  %%ymm6,          %%ymm8,          %%ymm8\n\t"
    "vmulpd                  %%ymm1,          %%ymm4,          %%ymm7\n\t"
    "vaddpd                  %%ymm7,          %%ymm9,          %%ymm9\n\t"
    "vmulpd                  %%ymm2,          %%ymm4,          %%ymm6\n\t"
    "vaddpd                  %%ymm6,          %%ymm10,         %%ymm10\n\t"
    "vmulpd                  %%ymm3,          %%ymm4,          %%ymm7\n\t"
    "vaddpd                  %%ymm7,          %%ymm11,         %%ymm11\n\t"

    "vmovapd           2 * 32(%%rdx),         %%ymm4\n\t"

    "vmulpd                  %%ymm0,          %%ymm5,          %%ymm6\n\t"
    "vaddpd                  %%ymm6,          %%ymm12,         %%ymm12\n\t"
    "vbroadcastsd       4 * 8(%%rsi),         %%ymm0\n\t"
    "vmulpd                  %%ymm1,          %%ymm5,          %%ymm7\n\t"
    "vaddpd                  %%ymm7,          %%ymm13,         %%ymm13\n\t"
    "vbroadcastsd       5 * 8(%%rsi),         %%ymm1\n\t"
    "vmulpd                  %%ymm2,          %%ymm5,          %%ymm6\n\t"
    "vaddpd                  %%ymm6,          %%ymm14,         %%ymm14\n\t"
    "vbroadcastsd       6 * 8(%%rsi),         %%ymm2\n\t"
    "vmulpd                  %%ymm3,          %%ymm5,          %%ymm7\n\t"
    "vaddpd                  %%ymm7,          %%ymm15,         %%ymm15\n\t"
    "vbroadcastsd       7 * 8(%%rsi),         %%ymm3\n\t"

    "addq                    $32,            %%rsi\n\t"
    "addq                    $2*32,          %%rdx\n\t"
    "decq                    %%rdi\n\t"

    "check%=:\n\t"
    "testq                   %%rdi,           %%rdi\n\t"
    "jg                      loop%=\n\t"

    "movq      %3,           %%rdi                  \n\t"  // alpha
    "movq      %4,           %%rsi                  \n\t"  // beta
    "vbroadcastsd           (%%rdi),          %%ymm6\n\t"
    "vbroadcastsd           (%%rsi),          %%ymm7\n\t"


    "vmulpd                  %%ymm6,          %%ymm8,          %%ymm8\n\t"
    "vmulpd                  %%ymm6,          %%ymm9,          %%ymm9\n\t"
    "vmulpd                  %%ymm6,          %%ymm10,         %%ymm10\n\t"
    "vmulpd                  %%ymm6,          %%ymm11,         %%ymm11\n\t"
    "vmulpd                  %%ymm6,          %%ymm12,         %%ymm12\n\t"
    "vmulpd                  %%ymm6,          %%ymm13,         %%ymm13\n\t"
    "vmulpd                  %%ymm6,          %%ymm14,         %%ymm14\n\t"
    "vmulpd                  %%ymm6,          %%ymm15,         %%ymm15\n\t"

    "leaq                    (,%%r8,8),       %%r8\n\t"
    "leaq                    (,%%r9,8),       %%r9\n\t"

    "leaq                    (,%%r9,2),       %%r10\n\t"
    "leaq                    (%%r10,%%r9),    %%r11\n\t"
    "leaq                    (%%rcx,%%r10,2), %%rdx\n\t"

    "#\n\t"
    "#       Update C(0,:)\n\t"
    "#\n\t"
    "vmovlpd                 (%%rcx),         %%xmm0,          %%xmm0\n\t"
    "vmovhpd                 (%%rcx,%%r9),    %%xmm0,          %%xmm0\n\t"
    "vmovlpd                 (%%rcx,%%r10),   %%xmm1,          %%xmm1\n\t"
    "vmovhpd                 (%%rcx,%%r11),   %%xmm1,          %%xmm1\n\t"
    "vmovlpd                 (%%rdx),         %%xmm2,          %%xmm2\n\t"
    "vmovhpd                 (%%rdx,%%r9),    %%xmm2,          %%xmm2\n\t"
    "vmovlpd                 (%%rdx,%%r10),   %%xmm3,          %%xmm3\n\t"
    "vmovhpd                 (%%rdx,%%r11),   %%xmm3,          %%xmm3\n\t"

    "vmulpd                  %%xmm7,          %%xmm0,          %%xmm0\n\t"
    "vmulpd                  %%xmm7,          %%xmm1,          %%xmm1\n\t"
    "vmulpd                  %%xmm7,          %%xmm2,          %%xmm2\n\t"
    "vmulpd                  %%xmm7,          %%xmm3,          %%xmm3\n\t"

    "vextractf128            $1,              %%ymm8,          %%xmm4\n\t"
    "vextractf128            $1,              %%ymm12,         %%xmm5\n\t"

    "vaddpd                  %%xmm0,          %%xmm8,          %%xmm0\n\t"
    "vaddpd                  %%xmm1,          %%xmm4,          %%xmm1\n\t"
    "vaddpd                  %%xmm2,          %%xmm12,         %%xmm2\n\t"
    "vaddpd                  %%xmm3,          %%xmm5,          %%xmm3\n\t"

    "vmovlpd                 %%xmm0,          (%%rcx)\n\t"
    "vmovhpd                 %%xmm0,          (%%rcx,%%r9)\n\t"
    "vmovlpd                 %%xmm1,          (%%rcx,%%r10)\n\t"
    "vmovhpd                 %%xmm1,          (%%rcx,%%r11)\n\t"
    "vmovlpd                 %%xmm2,          (%%rdx)\n\t"
    "vmovhpd                 %%xmm2,          (%%rdx,%%r9)\n\t"
    "vmovlpd                 %%xmm3,          (%%rdx,%%r10)\n\t"
    "vmovhpd                 %%xmm3,          (%%rdx,%%r11)\n\t"

    "#\n\t"
    "#       Update C(1,:)\n\t"
    "#\n\t"
    "addq                    %%r8,            %%rcx\n\t"
    "addq                    %%r8,            %%rdx\n\t"

    "vmovlpd                 (%%rcx),         %%xmm0,          %%xmm0\n\t"
    "vmovhpd                 (%%rcx,%%r9),    %%xmm0,          %%xmm0\n\t"
    "vmovlpd                 (%%rcx,%%r10),   %%xmm1,          %%xmm1\n\t"
    "vmovhpd                 (%%rcx,%%r11),   %%xmm1,          %%xmm1\n\t"
    "vmovlpd                 (%%rdx),         %%xmm2,          %%xmm2\n\t"
    "vmovhpd                 (%%rdx,%%r9),    %%xmm2,          %%xmm2\n\t"
    "vmovlpd                 (%%rdx,%%r10),   %%xmm3,          %%xmm3\n\t"
    "vmovhpd                 (%%rdx,%%r11),   %%xmm3,          %%xmm3\n\t"

    "vmulpd                  %%xmm7,          %%xmm0,          %%xmm0\n\t"
    "vmulpd                  %%xmm7,          %%xmm1,          %%xmm1\n\t"
    "vmulpd                  %%xmm7,          %%xmm2,          %%xmm2\n\t"
    "vmulpd                  %%xmm7,          %%xmm3,          %%xmm3\n\t"

    "vextractf128            $1,              %%ymm9,          %%xmm4\n\t"
    "vextractf128            $1,              %%ymm13,         %%xmm5\n\t"

    "vaddpd                  %%xmm0,          %%xmm9,          %%xmm0\n\t"
    "vaddpd                  %%xmm1,          %%xmm4,          %%xmm1\n\t"
    "vaddpd                  %%xmm2,          %%xmm13,         %%xmm2\n\t"
    "vaddpd                  %%xmm3,          %%xmm5,          %%xmm3\n\t"

    "vmovlpd                 %%xmm0,          (%%rcx)\n\t"
    "vmovhpd                 %%xmm0,          (%%rcx,%%r9)\n\t"
    "vmovlpd                 %%xmm1,          (%%rcx,%%r10)\n\t"
    "vmovhpd                 %%xmm1,          (%%rcx,%%r11)\n\t"
    "vmovlpd                 %%xmm2,          (%%rdx)\n\t"
    "vmovhpd                 %%xmm2,          (%%rdx,%%r9)\n\t"
    "vmovlpd                 %%xmm3,          (%%rdx,%%r10)\n\t"
    "vmovhpd                 %%xmm3,          (%%rdx,%%r11)\n\t"

    "#\n\t"
    "#       Update C(2,:)\n\t"
    "#\n\t"
    "addq                    %%r8,            %%rcx\n\t"
    "addq                    %%r8,            %%rdx\n\t"

    "vmovlpd                 (%%rcx),         %%xmm0,          %%xmm0\n\t"
    "vmovhpd                 (%%rcx,%%r9),    %%xmm0,          %%xmm0\n\t"
    "vmovlpd                 (%%rcx,%%r10),   %%xmm1,          %%xmm1\n\t"
    "vmovhpd                 (%%rcx,%%r11),   %%xmm1,          %%xmm1\n\t"
    "vmovlpd                 (%%rdx),         %%xmm2,          %%xmm2\n\t"
    "vmovhpd                 (%%rdx,%%r9),    %%xmm2,          %%xmm2\n\t"
    "vmovlpd                 (%%rdx,%%r10),   %%xmm3,          %%xmm3\n\t"
    "vmovhpd                 (%%rdx,%%r11),   %%xmm3,          %%xmm3\n\t"

    "vmulpd                  %%xmm7,          %%xmm0,          %%xmm0\n\t"
    "vmulpd                  %%xmm7,          %%xmm1,          %%xmm1\n\t"
    "vmulpd                  %%xmm7,          %%xmm2,          %%xmm2\n\t"
    "vmulpd                  %%xmm7,          %%xmm3,          %%xmm3\n\t"

    "vextractf128            $1,              %%ymm10,         %%xmm4\n\t"
    "vextractf128            $1,              %%ymm14,         %%xmm5\n\t"

    "vaddpd                  %%xmm0,          %%xmm10,         %%xmm0\n\t"
    "vaddpd                  %%xmm1,          %%xmm4,          %%xmm1\n\t"
    "vaddpd                  %%xmm2,          %%xmm14,         %%xmm2\n\t"
    "vaddpd                  %%xmm3,          %%xmm5,          %%xmm3\n\t"

    "vmovlpd                 %%xmm0,          (%%rcx)\n\t"
    "vmovhpd                 %%xmm0,          (%%rcx,%%r9)\n\t"
    "vmovlpd                 %%xmm1,          (%%rcx,%%r10)\n\t"
    "vmovhpd                 %%xmm1,          (%%rcx,%%r11)\n\t"
    "vmovlpd                 %%xmm2,          (%%rdx)\n\t"
    "vmovhpd                 %%xmm2,          (%%rdx,%%r9)\n\t"
    "vmovlpd                 %%xmm3,          (%%rdx,%%r10)\n\t"
    "vmovhpd                 %%xmm3,          (%%rdx,%%r11)\n\t"

    "#\n\t"
    "#       Update C(3,:)\n\t"
    "#\n\t"
    "addq                    %%r8,            %%rcx\n\t"
    "addq                    %%r8,            %%rdx\n\t"

    "vmovlpd                 (%%rcx),         %%xmm0,          %%xmm0\n\t"
    "vmovhpd                 (%%rcx,%%r9),    %%xmm0,          %%xmm0\n\t"
    "vmovlpd                 (%%rcx,%%r10),   %%xmm1,          %%xmm1\n\t"
    "vmovhpd                 (%%rcx,%%r11),   %%xmm1,          %%xmm1\n\t"
    "vmovlpd                 (%%rdx),         %%xmm2,          %%xmm2\n\t"
    "vmovhpd                 (%%rdx,%%r9),    %%xmm2,          %%xmm2\n\t"
    "vmovlpd                 (%%rdx,%%r10),   %%xmm3,          %%xmm3\n\t"
    "vmovhpd                 (%%rdx,%%r11),   %%xmm3,          %%xmm3\n\t"

    "vmulpd                  %%xmm7,          %%xmm0,          %%xmm0\n\t"
    "vmulpd                  %%xmm7,          %%xmm1,          %%xmm1\n\t"
    "vmulpd                  %%xmm7,          %%xmm2,          %%xmm2\n\t"
    "vmulpd                  %%xmm7,          %%xmm3,          %%xmm3\n\t"

    "vextractf128            $1,              %%ymm11,         %%xmm4\n\t"
    "vextractf128            $1,              %%ymm15,         %%xmm5\n\t"

    "vaddpd                  %%xmm0,          %%xmm11,         %%xmm0\n\t"
    "vaddpd                  %%xmm1,          %%xmm4,          %%xmm1\n\t"
    "vaddpd                  %%xmm2,          %%xmm15,         %%xmm2\n\t"
    "vaddpd                  %%xmm3,          %%xmm5,          %%xmm3\n\t"

    "vmovlpd                 %%xmm0,          (%%rcx)\n\t"
    "vmovhpd                 %%xmm0,          (%%rcx,%%r9)\n\t"
    "vmovlpd                 %%xmm1,          (%%rcx,%%r10)\n\t"
    "vmovhpd                 %%xmm1,          (%%rcx,%%r11)\n\t"
    "vmovlpd                 %%xmm2,          (%%rdx)\n\t"
    "vmovhpd                 %%xmm2,          (%%rdx,%%r9)\n\t"
    "vmovlpd                 %%xmm3,          (%%rdx,%%r10)\n\t"
    "vmovhpd                 %%xmm3,          (%%rdx,%%r11)\n\t"

    : // output
    : // input
        "m" (kc),       // 0
        "m" (A),        // 1
        "m" (B),        // 2
        "m" (pAlpha),   // 3
        "m" (pBeta),    // 4
        "m" (C),        // 5
        "m" (incRowC),  // 6
        "m" (incColC)   // 7
    : // register clobber list
        "rax",  "rbx",  "rcx",  "rdx",    "rsi",   "rdi",
        "r8",   "r9",   "r10",  "r11",
        "xmm0", "xmm1", "xmm2",  "xmm3",  "xmm4",  "xmm5",  "xmm6",  "xmm7",
        "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15",
        "memory"
    );
}

} } // namespace avx, ulmBLAS


#endif // ULMBLAS_IMPL_LEVEL3_UKERNEL_AVX_UGEMM_TCC
