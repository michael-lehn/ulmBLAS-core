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

#ifndef ULMBLAS_IMPL_LEVEL3_UKERNEL_FMA_UGEMM_TCC
#define ULMBLAS_IMPL_LEVEL3_UKERNEL_FMA_UGEMM_TCC 1

#include <type_traits>
#include <ulmblas/impl/level3/ukernel/ref/ugemm.h>

namespace ulmBLAS { namespace fma {

template <typename Index>
typename std::enable_if<std::is_convertible<Index, std::int64_t>::value
                     && BlockSize<double>::MR==4
                     && BlockSize<double>::NR==12
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



    "vmovapd             0*32(%%rdx),        %%ymm1 \n\t"
    "vmovapd             1*32(%%rdx),        %%ymm2 \n\t"
    "vmovapd             2*32(%%rdx),        %%ymm3 \n\t"

    "vxorpd                  %%ymm4,         %%ymm4,          %%ymm4    \n\t"
    "vxorpd                  %%ymm5,         %%ymm5,          %%ymm5    \n\t"
    "vxorpd                  %%ymm6,         %%ymm6,          %%ymm6    \n\t"
    "vxorpd                  %%ymm7,         %%ymm7,          %%ymm7    \n\t"
    "vxorpd                  %%ymm8,         %%ymm8,          %%ymm8    \n\t"
    "vxorpd                  %%ymm9,         %%ymm9,          %%ymm9    \n\t"
    "vxorpd                  %%ymm10,        %%ymm10,         %%ymm10   \n\t"
    "vxorpd                  %%ymm11,        %%ymm11,         %%ymm11   \n\t"
    "vxorpd                  %%ymm12,        %%ymm12,         %%ymm12   \n\t"
    "vxorpd                  %%ymm13,        %%ymm13,         %%ymm13   \n\t"
    "vxorpd                  %%ymm14,        %%ymm14,         %%ymm14   \n\t"
    "vxorpd                  %%ymm15,        %%ymm15,         %%ymm15   \n\t"

    "movq                    $3*32,          %%r13                      \n\t"
    "movq                    $4* 8,          %%r12                      \n\t"

    "jmp                     check%=                                    \n\t"

    "loop%=:                                                            \n\t"

    "vbroadcastsd       0* 8(%%rsi),          %%ymm0                    \n\t"
    "addq                    %%r13,           %%rdx                     \n\t"
    "vfmadd231pd             %%ymm0,          %%ymm1,          %%ymm4   \n\t"
    "vfmadd231pd             %%ymm0,          %%ymm2,          %%ymm8   \n\t"
    "vfmadd231pd             %%ymm0,          %%ymm3,          %%ymm12  \n\t"

    "vbroadcastsd       1* 8(%%rsi),          %%ymm0                    \n\t"
    "decq                    %%rdi                                      \n\t"
    "vfmadd231pd             %%ymm0,          %%ymm1,          %%ymm5   \n\t"
    "vfmadd231pd             %%ymm0,          %%ymm2,          %%ymm9   \n\t"
    "vfmadd231pd             %%ymm0,          %%ymm3,          %%ymm13  \n\t"


    "vbroadcastsd       2* 8(%%rsi),          %%ymm0                    \n\t"
    "addq                    %%r12,           %%rsi                     \n\t"
    "vfmadd231pd             %%ymm0,          %%ymm1,          %%ymm6   \n\t"
    "vfmadd231pd             %%ymm0,          %%ymm2,          %%ymm10  \n\t"
    "vfmadd231pd             %%ymm0,          %%ymm3,          %%ymm14  \n\t"

    "vbroadcastsd      -1* 8(%%rsi),          %%ymm0                    \n\t"
    "vfmadd231pd             %%ymm0,          %%ymm1,          %%ymm7   \n\t"
    "vmovapd            0*32(%%rdx),          %%ymm1                    \n\t"
    "vfmadd231pd             %%ymm0,          %%ymm2,          %%ymm11  \n\t"
    "vmovapd            1*32(%%rdx),          %%ymm2                    \n\t"
    "vfmadd231pd            %%ymm0,           %%ymm3,          %%ymm15  \n\t"
    "vmovapd            2*32(%%rdx),          %%ymm3                    \n\t"

    "check%=:                                                           \n\t"
    "testq                   %%rdi,           %%rdi                     \n\t"
    "jg                      loop%=                                     \n\t"

    "movq      %3,           %%rdi                      \n\t"  // alpha
    "vbroadcastsd           (%%rdi),          %%ymm0    \n\t"

    "vmulpd                  %%ymm0,         %%ymm4,        %%ymm4      \n\t"
    "vmulpd                  %%ymm0,         %%ymm5,        %%ymm5      \n\t"
    "vmulpd                  %%ymm0,         %%ymm6,        %%ymm6      \n\t"
    "vmulpd                  %%ymm0,         %%ymm7,        %%ymm7      \n\t"
    "vmulpd                  %%ymm0,         %%ymm8,        %%ymm8      \n\t"
    "vmulpd                  %%ymm0,         %%ymm9,        %%ymm9      \n\t"
    "vmulpd                  %%ymm0,         %%ymm10,       %%ymm10     \n\t"
    "vmulpd                  %%ymm0,         %%ymm11,       %%ymm11     \n\t"
    "vmulpd                  %%ymm0,         %%ymm12,       %%ymm12     \n\t"
    "vmulpd                  %%ymm0,         %%ymm13,       %%ymm13     \n\t"
    "vmulpd                  %%ymm0,         %%ymm14,       %%ymm14     \n\t"
    "vmulpd                  %%ymm0,         %%ymm15,       %%ymm15     \n\t"


    "leaq                    (,%%r8,8),       %%r8                            \n\t"
    "leaq                    (,%%r9,8),       %%r9                            \n\t"

    "leaq                    (,%%r9,2),       %%r10        # 2*incColC        \n\t"
    "leaq                    (%%r10,%%r9),    %%r11        # 3*incColC       \n\t"
    "leaq                    (%%rcx,%%r10,2), %%rdx        # C + 4*incColC   \n\t"
    "leaq                    (%%rdx,%%r10,2), %%rax        # C + 8*incColC   \n\t"

    // check if beta == 0
    "movq      %4,           %%rdi                      \n\t"  // beta
    "vbroadcastsd           (%%rdi),          %%ymm0    \n\t"
    "vxorpd                  %%ymm1,          %%ymm1,          %%ymm1  \n\t"
    "vucomisd                %%xmm0,          %%xmm1                   \n\t"
    "je                      beta_zero%=                               \n\t"

    // case: beta != 0
    "#                                                                          \n\t"
    "#       Update C(0,0:3)                                                    \n\t"
    "#                                                                          \n\t"
    "vmovlpd                 (%%rcx),         %%xmm1,          %%xmm1           \n\t"
    "vmovhpd                 (%%rcx,%%r9),    %%xmm1,          %%xmm1          \n\t"
    "vmovlpd                 (%%rcx,%%r10),   %%xmm2,          %%xmm2          \n\t"
    "vmovhpd                 (%%rcx,%%r11),   %%xmm2,          %%xmm2          \n\t"

    "vextractf128            $1,              %%ymm4,          %%xmm3            \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm4,          %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rcx)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rcx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rcx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rcx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(0,4:7)                                                    \n\t"
    "#                                                                          \n\t"
    "vmovlpd                 (%%rdx),         %%xmm1,          %%xmm1           \n\t"
    "vmovhpd                 (%%rdx,%%r9),    %%xmm1,          %%xmm1           \n\t"
    "vmovlpd                 (%%rdx,%%r10),   %%xmm2,          %%xmm2           \n\t"
    "vmovhpd                 (%%rdx,%%r11),   %%xmm2,          %%xmm2           \n\t"

    "vextractf128            $1,              %%ymm8,          %%xmm3           \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm8,          %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rdx)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rdx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rdx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rdx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(0,8:11)                                                   \n\t"
    "#                                                                          \n\t"
    "vmovlpd                 (%%rax),         %%xmm1,          %%xmm1           \n\t"
    "vmovhpd                 (%%rax,%%r9),    %%xmm1,          %%xmm1           \n\t"
    "vmovlpd                 (%%rax,%%r10),   %%xmm2,          %%xmm2           \n\t"
    "vmovhpd                 (%%rax,%%r11),   %%xmm2,          %%xmm2           \n\t"

    "vextractf128            $1,              %%ymm12,         %%xmm3           \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm12,         %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rax)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rax,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rax,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rax,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(1,0:3)                                                    \n\t"
    "#                                                                          \n\t"
    "addq                    %%r8,            %%rcx                             \n\t"
    "addq                    %%r8,            %%rdx                             \n\t"
    "addq                    %%r8,            %%rax                             \n\t"

    "vmovlpd                 (%%rcx),         %%xmm1,          %%xmm1           \n\t"
    "vmovhpd                 (%%rcx,%%r9),    %%xmm1,          %%xmm1           \n\t"
    "vmovlpd                 (%%rcx,%%r10),   %%xmm2,          %%xmm2           \n\t"
    "vmovhpd                 (%%rcx,%%r11),   %%xmm2,          %%xmm2           \n\t"

    "vextractf128            $1,              %%ymm5,          %%xmm3           \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm5,          %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rcx)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rcx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rcx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rcx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(1,4:7)                                                    \n\t"
    "#                                                                          \n\t"
    "vmovlpd                 (%%rdx),          %%xmm1,          %%xmm1          \n\t"
    "vmovhpd                 (%%rdx,%%r9),     %%xmm1,          %%xmm1          \n\t"
    "vmovlpd                 (%%rdx,%%r10),    %%xmm2,          %%xmm2          \n\t"
    "vmovhpd                 (%%rdx,%%r11),    %%xmm2,          %%xmm2          \n\t"

    "vextractf128            $1,             %%ymm9,          %%xmm3            \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm9,          %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rdx)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rdx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rdx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rdx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(1,8:11)                                                   \n\t"
    "#                                                                          \n\t"
    "vmovlpd                 (%%rax),         %%xmm1,          %%xmm1           \n\t"
    "vmovhpd                 (%%rax,%%r9),    %%xmm1,          %%xmm1           \n\t"
    "vmovlpd                 (%%rax,%%r10),   %%xmm2,          %%xmm2           \n\t"
    "vmovhpd                 (%%rax,%%r11),   %%xmm2,          %%xmm2           \n\t"

    "vextractf128            $1,              %%ymm13,         %%xmm3           \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm13,         %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rax)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rax,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rax,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rax,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(2,0:3)                                                    \n\t"
    "#                                                                          \n\t"
    "addq                    %%r8,            %%rcx                             \n\t"
    "addq                    %%r8,            %%rdx                             \n\t"
    "addq                    %%r8,            %%rax                             \n\t"

    "vmovlpd                 (%%rcx),         %%xmm1,          %%xmm1           \n\t"
    "vmovhpd                 (%%rcx,%%r9),    %%xmm1,          %%xmm1           \n\t"
    "vmovlpd                 (%%rcx,%%r10),   %%xmm2,          %%xmm2           \n\t"
    "vmovhpd                 (%%rcx,%%r11),   %%xmm2,          %%xmm2           \n\t"

    "vextractf128            $1,              %%ymm6,          %%xmm3           \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm6,          %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rcx)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rcx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rcx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rcx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(2,4:7)                                                    \n\t"
    "#                                                                          \n\t"
    "vmovlpd                 (%%rdx),         %%xmm1,          %%xmm1           \n\t"
    "vmovhpd                 (%%rdx,%%r9),    %%xmm1,          %%xmm1           \n\t"
    "vmovlpd                 (%%rdx,%%r10),   %%xmm2,          %%xmm2           \n\t"
    "vmovhpd                 (%%rdx,%%r11),   %%xmm2,          %%xmm2           \n\t"

    "vextractf128            $1,              %%ymm10,         %%xmm3           \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm10,         %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rdx)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rdx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rdx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rdx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(2,8:11)                                                   \n\t"
    "#                                                                          \n\t"
    "vmovlpd                 (%%rax),         %%xmm1,          %%xmm1           \n\t"
    "vmovhpd                 (%%rax,%%r9),    %%xmm1,          %%xmm1           \n\t"
    "vmovlpd                 (%%rax,%%r10),   %%xmm2,          %%xmm2           \n\t"
    "vmovhpd                 (%%rax,%%r11),   %%xmm2,          %%xmm2           \n\t"

    "vextractf128            $1,              %%ymm14,         %%xmm3           \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm14,         %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rax)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rax,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rax,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rax,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(3,0:3)                                                    \n\t"
    "#                                                                          \n\t"
    "addq                    %%r8,            %%rcx                             \n\t"
    "addq                    %%r8,            %%rdx                             \n\t"
    "addq                    %%r8,            %%rax                             \n\t"

    "vmovlpd                 (%%rcx),         %%xmm1,          %%xmm1           \n\t"
    "vmovhpd                 (%%rcx,%%r9),    %%xmm1,          %%xmm1           \n\t"
    "vmovlpd                 (%%rcx,%%r10),   %%xmm2,          %%xmm2           \n\t"
    "vmovhpd                 (%%rcx,%%r11),   %%xmm2,          %%xmm2           \n\t"

    "vextractf128            $1,              %%ymm7,          %%xmm3           \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm7,          %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rcx)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rcx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rcx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rcx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(3,4:7)                                                    \n\t"
    "#                                                                          \n\t"
    "vmovlpd                 (%%rdx),         %%xmm1,          %%xmm1           \n\t"
    "vmovhpd                 (%%rdx,%%r9),    %%xmm1,          %%xmm1           \n\t"
    "vmovlpd                 (%%rdx,%%r10),   %%xmm2,          %%xmm2           \n\t"
    "vmovhpd                 (%%rdx,%%r11),   %%xmm2,          %%xmm2           \n\t"

    "vextractf128            $1,              %%ymm11,         %%xmm3           \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm11,         %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rdx)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rdx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rdx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rdx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(3,8:11)                                                   \n\t"
    "#                                                                          \n\t"
    "vmovlpd                 (%%rax),         %%xmm1,          %%xmm1           \n\t"
    "vmovhpd                 (%%rax,%%r9),    %%xmm1,          %%xmm1           \n\t"
    "vmovlpd                 (%%rax,%%r10),   %%xmm2,          %%xmm2           \n\t"
    "vmovhpd                 (%%rax,%%r11),   %%xmm2,          %%xmm2           \n\t"

    "vextractf128            $1,              %%ymm15,         %%xmm3           \n\t"

    "vmulpd                  %%xmm0,          %%xmm1,          %%xmm1           \n\t"
    "vaddpd                  %%xmm1,          %%xmm15,         %%xmm1           \n\t"
    "vmulpd                  %%xmm0,          %%xmm2,          %%xmm2           \n\t"
    "vaddpd                  %%xmm2,          %%xmm3,          %%xmm2           \n\t"

    "vmovlpd                 %%xmm1,          (%%rax)                           \n\t"
    "vmovhpd                 %%xmm1,          (%%rax,%%r9)                      \n\t"
    "vmovlpd                 %%xmm2,          (%%rax,%%r10)                     \n\t"
    "vmovhpd                 %%xmm2,          (%%rax,%%r11)                     \n\t"

    "jmp done%=                              \n\t"

    // case: beta == 0
    "beta_zero%=:                           \n\t"
    "#                                                                          \n\t"
    "#       Update C(0,0:3)                                                    \n\t"
    "#                                                                          \n\t"

    "vextractf128            $1,              %%ymm4,          %%xmm3            \n\t"

    "vmovlpd                 %%xmm4,          (%%rcx)                           \n\t"
    "vmovhpd                 %%xmm4,          (%%rcx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rcx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rcx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(0,4:7)                                                    \n\t"
    "#                                                                          \n\t"

    "vextractf128            $1,              %%ymm8,          %%xmm3           \n\t"

    "vmovlpd                 %%xmm8,          (%%rdx)                           \n\t"
    "vmovhpd                 %%xmm8,          (%%rdx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rdx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rdx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(0,8:11)                                                   \n\t"
    "#                                                                          \n\t"

    "vextractf128            $1,              %%ymm12,         %%xmm3           \n\t"

    "vmovlpd                 %%xmm12,         (%%rax)                           \n\t"
    "vmovhpd                 %%xmm12,         (%%rax,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rax,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rax,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(1,0:3)                                                    \n\t"
    "#                                                                          \n\t"
    "addq                    %%r8,            %%rcx                             \n\t"
    "addq                    %%r8,            %%rdx                             \n\t"
    "addq                    %%r8,            %%rax                             \n\t"

    "vextractf128            $1,              %%ymm5,          %%xmm3           \n\t"

    "vmovlpd                 %%xmm5,          (%%rcx)                           \n\t"
    "vmovhpd                 %%xmm5,          (%%rcx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rcx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rcx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(1,4:7)                                                    \n\t"
    "#                                                                          \n\t"

    "vextractf128            $1,             %%ymm9,          %%xmm3            \n\t"

    "vmovlpd                 %%xmm9,          (%%rdx)                           \n\t"
    "vmovhpd                 %%xmm9,          (%%rdx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rdx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rdx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(1,8:11)                                                   \n\t"
    "#                                                                          \n\t"

    "vextractf128            $1,              %%ymm13,         %%xmm3           \n\t"

    "vmovlpd                 %%xmm13,         (%%rax)                           \n\t"
    "vmovhpd                 %%xmm13,         (%%rax,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rax,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rax,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(2,0:3)                                                    \n\t"
    "#                                                                          \n\t"
    "addq                    %%r8,            %%rcx                             \n\t"
    "addq                    %%r8,            %%rdx                             \n\t"
    "addq                    %%r8,            %%rax                             \n\t"

    "vextractf128            $1,              %%ymm6,          %%xmm3           \n\t"

    "vmovlpd                 %%xmm6,          (%%rcx)                           \n\t"
    "vmovhpd                 %%xmm6,          (%%rcx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rcx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rcx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(2,4:7)                                                    \n\t"
    "#                                                                          \n\t"

    "vextractf128            $1,              %%ymm10,         %%xmm3           \n\t"

    "vmovlpd                 %%xmm10,         (%%rdx)                           \n\t"
    "vmovhpd                 %%xmm10,         (%%rdx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rdx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rdx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(2,8:11)                                                   \n\t"
    "#                                                                          \n\t"

    "vextractf128            $1,              %%ymm14,         %%xmm3           \n\t"

    "vmovlpd                 %%xmm14,         (%%rax)                           \n\t"
    "vmovhpd                 %%xmm14,         (%%rax,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rax,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rax,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(3,0:3)                                                    \n\t"
    "#                                                                          \n\t"
    "addq                    %%r8,            %%rcx                             \n\t"
    "addq                    %%r8,            %%rdx                             \n\t"
    "addq                    %%r8,            %%rax                             \n\t"

    "vextractf128            $1,              %%ymm7,          %%xmm3           \n\t"

    "vmovlpd                 %%xmm7,          (%%rcx)                           \n\t"
    "vmovhpd                 %%xmm7,          (%%rcx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rcx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rcx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(3,4:7)                                                    \n\t"
    "#                                                                          \n\t"

    "vextractf128            $1,              %%ymm11,         %%xmm3           \n\t"

    "vmovlpd                 %%xmm11,         (%%rdx)                           \n\t"
    "vmovhpd                 %%xmm11,         (%%rdx,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rdx,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rdx,%%r11)                     \n\t"

    "#                                                                          \n\t"
    "#       Update C(3,8:11)                                                   \n\t"
    "#                                                                          \n\t"

    "vextractf128            $1,              %%ymm15,         %%xmm3           \n\t"

    "vmovlpd                 %%xmm15,         (%%rax)                           \n\t"
    "vmovhpd                 %%xmm15,         (%%rax,%%r9)                      \n\t"
    "vmovlpd                 %%xmm3,          (%%rax,%%r10)                     \n\t"
    "vmovhpd                 %%xmm3,          (%%rax,%%r11)                     \n\t"

    "done%=:                           \n\t"
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
        "rax",  "rcx",  "rdx",    "rsi",   "rdi",
        "r8",   "r9",   "r10",  "r11", "r12", "r13",
        "xmm0", "xmm1", "xmm2",  "xmm3",  "xmm4",  "xmm5",  "xmm6",  "xmm7",
        "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15",
        "memory"
    );
}

#define DGEMM_INPUT_GS_BETA_NZ \
    "vmovlpd    (%%rcx        ),  %%xmm0,  %%xmm0  \n\t" \
    "vmovhpd    (%%rcx,%%rsi,1),  %%xmm0,  %%xmm0  \n\t" \
    "vmovlpd    (%%rcx,%%rsi,2),  %%xmm1,  %%xmm1  \n\t" \
    "vmovhpd    (%%rcx,%%r13  ),  %%xmm1,  %%xmm1  \n\t" \
    "vperm2f128 $0x20,   %%ymm1,  %%ymm0,  %%ymm0  \n\t" /*\
    "vmovlps    (%%rcx,%%rsi,4),  %%xmm2,  %%xmm2  \n\t" \
    "vmovhps    (%%rcx,%%r15  ),  %%xmm2,  %%xmm2  \n\t" \
    "vmovlps    (%%rcx,%%r13,2),  %%xmm1,  %%xmm1  \n\t" \
    "vmovhps    (%%rcx,%%r10  ),  %%xmm1,  %%xmm1  \n\t" \
    "vperm2f128 $0x20,   %%ymm1,  %%ymm2,  %%ymm2  \n\t"*/

#define DGEMM_OUTPUT_GS_BETA_NZ \
    "vextractf128  $1, %%ymm0,  %%xmm1           \n\t" \
    "vmovlpd           %%xmm0,  (%%rcx        )  \n\t" \
    "vmovhpd           %%xmm0,  (%%rcx,%%rsi  )  \n\t" \
    "vmovlpd           %%xmm1,  (%%rcx,%%rsi,2)  \n\t" \
    "vmovhpd           %%xmm1,  (%%rcx,%%r13  )  \n\t" /*\
    "vextractf128  $1, %%ymm2,  %%xmm1           \n\t" \
    "vmovlpd           %%xmm2,  (%%rcx,%%rsi,4)  \n\t" \
    "vmovhpd           %%xmm2,  (%%rcx,%%r15  )  \n\t" \
    "vmovlpd           %%xmm1,  (%%rcx,%%r13,2)  \n\t" \
    "vmovhpd           %%xmm1,  (%%rcx,%%r10  )  \n\t"*/

template <typename Index>
typename std::enable_if<std::is_convertible<Index, std::int64_t>::value
                     && BlockSize<double>::MR==8
                     && BlockSize<double>::NR==6
                     && BlockSize<double>::align==32,
         void>::type
ugemm(Index kc, double alpha,
      const double *A, const double *B,
      double beta,
      double *C, Index incRowC_, Index incColC_,
      const double *, const double *)
{
    int64_t k_iter  = kc / 4;
    int64_t k_left  = kc % 4;

    int64_t rs_c    = incRowC_;
    int64_t cs_c    = incColC_;

    double *pAlpha  = &alpha;
    double *pBeta   = &beta;

    __asm__ volatile
    (
    "                                            \n\t"
    "vzeroall                                    \n\t" // zero all xmm/ymm registers.
    "                                            \n\t"
    "                                            \n\t"
    "movq                %2, %%rax               \n\t" // load address of a.
    "movq                %3, %%rbx               \n\t" // load address of b.
    //"movq                %9, %%r15               \n\t" // load address of b_next.
    "                                            \n\t"
    "addq           $32 * 4, %%rax               \n\t"
    "                                            \n\t" // initialize loop by pre-loading
    "vmovaps           -4 * 32(%%rax), %%ymm0    \n\t"
    "vmovaps           -3 * 32(%%rax), %%ymm1    \n\t"
    "                                            \n\t"
    "movq                %6, %%rcx               \n\t" // load address of c
    "movq                %8, %%rdi               \n\t" // load cs_c
    "leaq        (,%%rdi,8), %%rdi               \n\t" // cs_c *= sizeof(double)
    "                                            \n\t"
    "leaq   (%%rdi,%%rdi,2), %%r13               \n\t" // r13 = 3*cs_c;
    "leaq   (%%rcx,%%r13,1), %%rdx               \n\t" // rdx = c + 3*cs_c;
    "prefetcht0   7 * 8(%%rcx)                   \n\t" // prefetch c + 0*cs_c
    "prefetcht0   7 * 8(%%rcx,%%rdi)             \n\t" // prefetch c + 1*cs_c
    "prefetcht0   7 * 8(%%rcx,%%rdi,2)           \n\t" // prefetch c + 2*cs_c
    "prefetcht0   7 * 8(%%rdx)                   \n\t" // prefetch c + 3*cs_c
    "prefetcht0   7 * 8(%%rdx,%%rdi)             \n\t" // prefetch c + 4*cs_c
    "prefetcht0   7 * 8(%%rdx,%%rdi,2)           \n\t" // prefetch c + 5*cs_c
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "movq      %0, %%rsi                         \n\t" // i = k_iter;
    "testq  %%rsi, %%rsi                         \n\t" // check i via logical AND.
    "je     .DCONSIDKLEFT%=                      \n\t" // if i == 0, jump to code that
    "                                            \n\t" // contains the k_left loop.
    "                                            \n\t"
    "                                            \n\t"
    ".DLOOPKITER%=:                              \n\t" // MAIN LOOP
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t" // iteration 0
    "prefetcht0  16 * 32(%%rax)                  \n\t"
    "                                            \n\t"
    "vbroadcastsd       0 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd       1 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm4    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm5    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm6    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm7    \n\t"
    "                                            \n\t"
    "vbroadcastsd       2 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd       3 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm8    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm9    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm10   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm11   \n\t"
    "                                            \n\t"
    "vbroadcastsd       4 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd       5 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm12   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm13   \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm14   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm15   \n\t"
    "                                            \n\t"
    "vmovaps           -2 * 32(%%rax), %%ymm0    \n\t"
    "vmovaps           -1 * 32(%%rax), %%ymm1    \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t" // iteration 1
    "vbroadcastsd       6 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd       7 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm4    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm5    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm6    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm7    \n\t"
    "                                            \n\t"
    "vbroadcastsd       8 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd       9 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm8    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm9    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm10   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm11   \n\t"
    "                                            \n\t"
    "vbroadcastsd      10 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd      11 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm12   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm13   \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm14   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm15   \n\t"
    "                                            \n\t"
    "vmovaps            0 * 32(%%rax), %%ymm0    \n\t"
    "vmovaps            1 * 32(%%rax), %%ymm1    \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t" // iteration 2
    "prefetcht0  20 * 32(%%rax)                  \n\t"
    "                                            \n\t"
    "vbroadcastsd      12 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd      13 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm4    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm5    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm6    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm7    \n\t"
    "                                            \n\t"
    "vbroadcastsd      14 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd      15 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm8    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm9    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm10   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm11   \n\t"
    "                                            \n\t"
    "vbroadcastsd      16 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd      17 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm12   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm13   \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm14   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm15   \n\t"
    "                                            \n\t"
    "vmovaps            2 * 32(%%rax), %%ymm0    \n\t"
    "vmovaps            3 * 32(%%rax), %%ymm1    \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t" // iteration 3
    "vbroadcastsd      18 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd      19 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm4    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm5    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm6    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm7    \n\t"
    "                                            \n\t"
    "vbroadcastsd      20 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd      21 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm8    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm9    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm10   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm11   \n\t"
    "                                            \n\t"
    "vbroadcastsd      22 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd      23 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm12   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm13   \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm14   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm15   \n\t"
    "                                            \n\t"
    "addq           $4 * 8 * 8, %%rax            \n\t" // a += 4*8 (unroll x mr)
    "addq           $4 * 6 * 8, %%rbx            \n\t" // b += 4*6 (unroll x nr)
    "                                            \n\t"
    "vmovaps           -4 * 32(%%rax), %%ymm0    \n\t"
    "vmovaps           -3 * 32(%%rax), %%ymm1    \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "decq   %%rsi                                \n\t" // i -= 1;
    "jne    .DLOOPKITER%=                        \n\t" // iterate again if i != 0.
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    ".DCONSIDKLEFT%=:                            \n\t"
    "                                            \n\t"
    "movq      %1, %%rsi                         \n\t" // i = k_left;
    "testq  %%rsi, %%rsi                         \n\t" // check i via logical AND.
    "je     .DPOSTACCUM%=                        \n\t" // if i == 0, we're done; jump to end.
    "                                            \n\t" // else, we prepare to enter k_left loop.
    "                                            \n\t"
    "                                            \n\t"
    ".DLOOPKLEFT%=:                                \n\t" // EDGE LOOP
    "                                            \n\t"
    "prefetcht0  16 * 32(%%rax)                  \n\t"
    "                                            \n\t"
    "vbroadcastsd       0 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd       1 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm4    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm5    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm6    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm7    \n\t"
    "                                            \n\t"
    "vbroadcastsd       2 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd       3 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm8    \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm9    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm10   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm11   \n\t"
    "                                            \n\t"
    "vbroadcastsd       4 *  8(%%rbx), %%ymm2    \n\t"
    "vbroadcastsd       5 *  8(%%rbx), %%ymm3    \n\t"
    "vfmadd231pd       %%ymm0, %%ymm2, %%ymm12   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm2, %%ymm13   \n\t"
    "vfmadd231pd       %%ymm0, %%ymm3, %%ymm14   \n\t"
    "vfmadd231pd       %%ymm1, %%ymm3, %%ymm15   \n\t"
    "                                            \n\t"
    "addq           $1 * 8 * 8, %%rax            \n\t" // a += 1*8 (unroll x mr)
    "addq           $1 * 6 * 8, %%rbx            \n\t" // b += 1*6 (unroll x nr)
    "                                            \n\t"
    "vmovaps           -4 * 32(%%rax), %%ymm0    \n\t"
    "vmovaps           -3 * 32(%%rax), %%ymm1    \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "decq   %%rsi                                \n\t" // i -= 1;
    "jne    .DLOOPKLEFT%=                        \n\t" // iterate again if i != 0.
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    ".DPOSTACCUM%=:                              \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "movq         %4, %%rax                      \n\t" // load address of alpha
    "movq         %5, %%rbx                      \n\t" // load address of beta 
    "vbroadcastsd    (%%rax), %%ymm0             \n\t" // load alpha and duplicate
    "vbroadcastsd    (%%rbx), %%ymm3             \n\t" // load beta and duplicate
    "                                            \n\t"
    "vmulpd           %%ymm0,  %%ymm4,  %%ymm4   \n\t" // scale by alpha
    "vmulpd           %%ymm0,  %%ymm5,  %%ymm5   \n\t"
    "vmulpd           %%ymm0,  %%ymm6,  %%ymm6   \n\t"
    "vmulpd           %%ymm0,  %%ymm7,  %%ymm7   \n\t"
    "vmulpd           %%ymm0,  %%ymm8,  %%ymm8   \n\t"
    "vmulpd           %%ymm0,  %%ymm9,  %%ymm9   \n\t"
    "vmulpd           %%ymm0,  %%ymm10, %%ymm10  \n\t"
    "vmulpd           %%ymm0,  %%ymm11, %%ymm11  \n\t"
    "vmulpd           %%ymm0,  %%ymm12, %%ymm12  \n\t"
    "vmulpd           %%ymm0,  %%ymm13, %%ymm13  \n\t"
    "vmulpd           %%ymm0,  %%ymm14, %%ymm14  \n\t"
    "vmulpd           %%ymm0,  %%ymm15, %%ymm15  \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "movq                %7, %%rsi               \n\t" // load rs_c
    "leaq        (,%%rsi,8), %%rsi               \n\t" // rsi = rs_c * sizeof(double)
    "                                            \n\t"
    "leaq   (%%rcx,%%rsi,4), %%rdx               \n\t" // load address of c +  4*rs_c;
    "                                            \n\t"
    "leaq   (%%rsi,%%rsi,2), %%r13               \n\t" // r13 = 3*rs_c;
    //"leaq   (%%rsi,%%rsi,4), %%r15               \n\t" // r15 = 5*rs_c;
    //"leaq   (%%r13,%%rsi,4), %%r10               \n\t" // r10 = 7*rs_c;
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t" // determine if
    "                                            \n\t" //    c    % 32 == 0, AND
    "                                            \n\t" //  8*cs_c % 32 == 0, AND
    "                                            \n\t" //    rs_c      == 1
    "                                            \n\t" // ie: aligned, ldim aligned, and
    "                                            \n\t" // column-stored
    "                                            \n\t"
    "cmpq       $8, %%rsi                        \n\t" // set ZF if (8*rs_c) == 8.
    "sete           %%bl                         \n\t" // bl = ( ZF == 1 ? 1 : 0 );
    "testq     $31, %%rcx                        \n\t" // set ZF if c & 32 is zero.
    "setz           %%bh                         \n\t" // bh = ( ZF == 0 ? 1 : 0 );
    "testq     $31, %%rdi                        \n\t" // set ZF if (8*cs_c) & 32 is zero.
    "setz           %%al                         \n\t" // al = ( ZF == 0 ? 1 : 0 );
    "                                            \n\t" // and(bl,bh) followed by
    "                                            \n\t" // and(bh,al) will reveal result
    "                                            \n\t"
    "                                            \n\t" // now avoid loading C if beta == 0
    "                                            \n\t"
    "vxorpd    %%ymm0,  %%ymm0,  %%ymm0          \n\t" // set ymm0 to zero.
    "vucomisd  %%xmm0,  %%xmm3                   \n\t" // set ZF if beta == 0.
    "je      .DBETAZERO%=                        \n\t" // if ZF = 1, jump to beta == 0 case
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t" // check if aligned/column-stored
    "andb     %%bl, %%bh                         \n\t" // set ZF if bl & bh == 1.
    "andb     %%bh, %%al                         \n\t" // set ZF if bh & al == 1.
    "jne     .DCOLSTORED%=                       \n\t" // jump to column storage case
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    ".DGENSTORED%=:                              \n\t"
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm4,  %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm6,  %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm8,  %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm10, %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm12, %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm14, %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "                                            \n\t"
    "                                            \n\t"
    "movq      %%rdx, %%rcx                      \n\t" // rcx = c + 4*rs_c
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm5,  %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm7,  %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm9,  %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm11, %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm13, %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    DGEMM_INPUT_GS_BETA_NZ
    "vfmadd213pd      %%ymm15, %%ymm3,  %%ymm0   \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "jmp    .DDONE%=                             \n\t" // jump to end.
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    ".DCOLSTORED%=:                              \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps    (%%rcx),       %%ymm0            \n\t"
    "vfmadd213pd      %%ymm4,  %%ymm3,  %%ymm0   \n\t"
    "vmovaps          %%ymm0,  (%%rcx)           \n\t"
    "addq      %%rdi, %%rcx                      \n\t"
    "vmovaps    (%%rdx),       %%ymm1            \n\t"
    "vfmadd213pd      %%ymm5,  %%ymm3,  %%ymm1   \n\t"
    "vmovaps          %%ymm1,  (%%rdx)           \n\t"
    "addq      %%rdi, %%rdx                      \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps    (%%rcx),       %%ymm0            \n\t"
    "vfmadd213pd      %%ymm6,  %%ymm3,  %%ymm0   \n\t"
    "vmovaps          %%ymm0,  (%%rcx)           \n\t"
    "addq      %%rdi, %%rcx                      \n\t"
    "vmovaps    (%%rdx),       %%ymm1            \n\t"
    "vfmadd213pd      %%ymm7,  %%ymm3,  %%ymm1   \n\t"
    "vmovaps          %%ymm1,  (%%rdx)           \n\t"
    "addq      %%rdi, %%rdx                      \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps    (%%rcx),       %%ymm0            \n\t"
    "vfmadd213pd      %%ymm8,  %%ymm3,  %%ymm0   \n\t"
    "vmovaps          %%ymm0,  (%%rcx)           \n\t"
    "addq      %%rdi, %%rcx                      \n\t"
    "vmovaps    (%%rdx),       %%ymm1            \n\t"
    "vfmadd213pd      %%ymm9,  %%ymm3,  %%ymm1   \n\t"
    "vmovaps          %%ymm1,  (%%rdx)           \n\t"
    "addq      %%rdi, %%rdx                      \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps    (%%rcx),       %%ymm0            \n\t"
    "vfmadd213pd      %%ymm10, %%ymm3,  %%ymm0   \n\t"
    "vmovaps          %%ymm0,  (%%rcx)           \n\t"
    "addq      %%rdi, %%rcx                      \n\t"
    "vmovaps    (%%rdx),       %%ymm1            \n\t"
    "vfmadd213pd      %%ymm11, %%ymm3,  %%ymm1   \n\t"
    "vmovaps          %%ymm1,  (%%rdx)           \n\t"
    "addq      %%rdi, %%rdx                      \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps    (%%rcx),       %%ymm0            \n\t"
    "vfmadd213pd      %%ymm12, %%ymm3,  %%ymm0   \n\t"
    "vmovaps          %%ymm0,  (%%rcx)           \n\t"
    "addq      %%rdi, %%rcx                      \n\t"
    "vmovaps    (%%rdx),       %%ymm1            \n\t"
    "vfmadd213pd      %%ymm13, %%ymm3,  %%ymm1   \n\t"
    "vmovaps          %%ymm1,  (%%rdx)           \n\t"
    "addq      %%rdi, %%rdx                      \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps    (%%rcx),       %%ymm0            \n\t"
    "vfmadd213pd      %%ymm14, %%ymm3,  %%ymm0   \n\t"
    "vmovaps          %%ymm0,  (%%rcx)           \n\t"
    //"addq      %%rdi, %%rcx                      \n\t"
    "vmovaps    (%%rdx),       %%ymm1            \n\t"
    "vfmadd213pd      %%ymm15, %%ymm3,  %%ymm1   \n\t"
    "vmovaps          %%ymm1,  (%%rdx)           \n\t"
    //"addq      %%rdi, %%rdx                      \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "jmp    .DDONE%=                             \n\t" // jump to end.
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    ".DBETAZERO%=:                               \n\t"
    "                                            \n\t" // check if aligned/column-stored
    "andb     %%bl, %%bh                         \n\t" // set ZF if bl & bh == 1.
    "andb     %%bh, %%al                         \n\t" // set ZF if bh & al == 1.
    "jne     .DCOLSTORBZ%=                       \n\t" // jump to column storage case
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    ".DGENSTORBZ%=:                              \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm4,  %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm6,  %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm8,  %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm10, %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm12, %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm14, %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "                                            \n\t"
    "                                            \n\t"
    "movq      %%rdx, %%rcx                      \n\t" // rcx = c + 4*rs_c
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm5,  %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm7,  %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm9,  %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm11, %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm13, %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "addq      %%rdi, %%rcx                      \n\t" // c += cs_c;
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps           %%ymm15, %%ymm0           \n\t"
    DGEMM_OUTPUT_GS_BETA_NZ
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "jmp    .DDONE%=                             \n\t" // jump to end.
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    ".DCOLSTORBZ%=:                              \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps          %%ymm4,  (%%rcx)           \n\t"
    "addq      %%rdi, %%rcx                      \n\t"
    "vmovaps          %%ymm5,  (%%rdx)           \n\t"
    "addq      %%rdi, %%rdx                      \n\t"
    "                                            \n\t"
    "vmovaps          %%ymm6,  (%%rcx)           \n\t"
    "addq      %%rdi, %%rcx                      \n\t"
    "vmovaps          %%ymm7,  (%%rdx)           \n\t"
    "addq      %%rdi, %%rdx                      \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps          %%ymm8,  (%%rcx)           \n\t"
    "addq      %%rdi, %%rcx                      \n\t"
    "vmovaps          %%ymm9,  (%%rdx)           \n\t"
    "addq      %%rdi, %%rdx                      \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps          %%ymm10, (%%rcx)           \n\t"
    "addq      %%rdi, %%rcx                      \n\t"
    "vmovaps          %%ymm11, (%%rdx)           \n\t"
    "addq      %%rdi, %%rdx                      \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps          %%ymm12, (%%rcx)           \n\t"
    "addq      %%rdi, %%rcx                      \n\t"
    "vmovaps          %%ymm13, (%%rdx)           \n\t"
    "addq      %%rdi, %%rdx                      \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "vmovaps          %%ymm14, (%%rcx)           \n\t"
    "                                            \n\t"
    "vmovaps          %%ymm15, (%%rdx)           \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    "                                            \n\t"
    ".DDONE%=:                                   \n\t"
    "                                            \n\t"

    : // output operands (none)
    : // input operands
      "m" (k_iter), // 0
      "m" (k_left), // 1
      "m" (A),      // 2
      "m" (A),      // 3
      "m" (pAlpha),  // 4
      "m" (pBeta),   // 5
      "m" (C),      // 6
      "m" (rs_c),   // 7
      "m" (cs_c)    // 8
    : // register clobber list
      "rax", "rbx", "rcx", "rdx", "rsi", "rdi", 
      "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15",
      "xmm0", "xmm1", "xmm2", "xmm3",
      "xmm4", "xmm5", "xmm6", "xmm7",
      "xmm8", "xmm9", "xmm10", "xmm11",
      "xmm12", "xmm13", "xmm14", "xmm15",
      "memory"
    );
}

} } // namespace fma, ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_UKERNEL_FMA_UGEMM_TCC
