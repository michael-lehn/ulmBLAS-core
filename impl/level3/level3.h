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

#ifndef ULMBLAS_IMPL_LEVEL3_LEVEL3_H
#define ULMBLAS_IMPL_LEVEL3_LEVEL3_H 1

#include <ulmblas/impl/level3/gemm.h>
#include <ulmblas/impl/level3/helmm.h>
#include <ulmblas/impl/level3/helr2k.h>
#include <ulmblas/impl/level3/helrk.h>
#include <ulmblas/impl/level3/heumm.h>
#include <ulmblas/impl/level3/heur2k.h>
#include <ulmblas/impl/level3/heurk.h>
#include <ulmblas/impl/level3/level3.h>
#include <ulmblas/impl/level3/mkernel/mgemm.h>
#include <ulmblas/impl/level3/mkernel/msylrk.h>
#include <ulmblas/impl/level3/mkernel/msyurk.h>
#include <ulmblas/impl/level3/mkernel/mtrlmm.h>
#include <ulmblas/impl/level3/mkernel/mtrlsm.h>
#include <ulmblas/impl/level3/mkernel/mtrumm.h>
#include <ulmblas/impl/level3/mkernel/mtrusm.h>
#include <ulmblas/impl/level3/pack/gepack.h>
#include <ulmblas/impl/level3/pack/helpack.h>
#include <ulmblas/impl/level3/pack/heupack.h>
#include <ulmblas/impl/level3/pack/sylpack.h>
#include <ulmblas/impl/level3/pack/syupack.h>
#include <ulmblas/impl/level3/pack/trlpack.h>
#include <ulmblas/impl/level3/pack/trlspack.h>
#include <ulmblas/impl/level3/pack/trupack.h>
#include <ulmblas/impl/level3/pack/truspack.h>
#include <ulmblas/impl/level3/sylmm.h>
#include <ulmblas/impl/level3/sylr2k.h>
#include <ulmblas/impl/level3/sylrk.h>
#include <ulmblas/impl/level3/syumm.h>
#include <ulmblas/impl/level3/syur2k.h>
#include <ulmblas/impl/level3/syurk.h>
#include <ulmblas/impl/level3/trlmm.h>
#include <ulmblas/impl/level3/trlsm.h>
#include <ulmblas/impl/level3/trumm.h>
#include <ulmblas/impl/level3/trusm.h>
#include <ulmblas/impl/level3/ukernel/ref/ugemm.h>
#include <ulmblas/impl/level3/ukernel/ref/utrlsm.h>
#include <ulmblas/impl/level3/ukernel/ref/utrusm.h>
#include <ulmblas/impl/level3/ukernel/sse/ugemm.h>
#include <ulmblas/impl/level3/ukernel/sse/utrlsm.h>
#include <ulmblas/impl/level3/ukernel/sse/utrusm.h>
#include <ulmblas/impl/level3/ukernel/ugemm.h>
#include <ulmblas/impl/level3/ukernel/usylrk.h>
#include <ulmblas/impl/level3/ukernel/usyurk.h>
#include <ulmblas/impl/level3/ukernel/utrlsm.h>
#include <ulmblas/impl/level3/ukernel/utrusm.h>

#endif // ULMBLAS_IMPL_LEVEL3_LEVEL3_H
