/*
 *   Copyright (c) 2007, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <cassert>
#include <extensions/flens/sparse_blas.h>
#include <iostream>

#ifdef MKL
#    ifdef MAC
#        include <Intel_MKL/mkl_spblas.h>
#    else
#        include <mkl_spblas.h>
#    endif
#endif

namespace flens {

//-- csr - compressed sparse row - (the Intel variant for crs) -----------------

#ifdef MKL_SCS
void
csrmv(Transpose Trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE k, float alpha, char *matdescra,
      float  *values, FLENS_DEFAULT_INDEXTYPE *columns,  FLENS_DEFAULT_INDEXTYPE *pointerB, FLENS_DEFAULT_INDEXTYPE *pointerE,
      float *x, float beta, float *y)
{
        char trans = (Trans==NoTrans) ? 'N' : 'T';

        mkl_scsrmv (&trans, &m, &k, &alpha, matdescra,
                    values, columns, pointerB, pointerE,
                    x, &beta, y);
}

void
csrmv(Transpose Trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE k, double alpha, char *matdescra,
      double  *values, FLENS_DEFAULT_INDEXTYPE *columns,  FLENS_DEFAULT_INDEXTYPE *pointerB, FLENS_DEFAULT_INDEXTYPE *pointerE,
      double *x, double beta, double *y)
{
        char trans = (Trans==NoTrans) ? 'N' : 'T';

        mkl_dcsrmv (&trans, &m, &k, &alpha, matdescra,
                    values, columns, pointerB, pointerE,
                    x, &beta, y);
}

void
csrmm(Transpose transA, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k, float alpha, char *matdescrA,
      const float *values, const FLENS_DEFAULT_INDEXTYPE *columns, const FLENS_DEFAULT_INDEXTYPE *pointerB,
      const FLENS_DEFAULT_INDEXTYPE *pointerE, const float *B, FLENS_DEFAULT_INDEXTYPE ldb,
      float beta, float *C, FLENS_DEFAULT_INDEXTYPE ldc)
{
        char trans = (transA==NoTrans) ? 'N' : 'T';

        // TODO: what about constness?
        mkl_scsrmm(&trans, &m, &n, &k, &alpha, matdescrA,
                   const_cast<float *>(values),
                   const_cast<FLENS_DEFAULT_INDEXTYPE *>(columns),
                   const_cast<FLENS_DEFAULT_INDEXTYPE *>(pointerB),
                   const_cast<FLENS_DEFAULT_INDEXTYPE *>(pointerE),
                   const_cast<float *>(B),
                   &ldb, &beta, C, &ldc);
}

void
csrmm(Transpose transA, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k, double alpha, char *matdescrA,
      const double *values, const FLENS_DEFAULT_INDEXTYPE *columns, const FLENS_DEFAULT_INDEXTYPE *pointerB,
      const FLENS_DEFAULT_INDEXTYPE *pointerE, const double *B, FLENS_DEFAULT_INDEXTYPE ldb,
      double beta, double *C, FLENS_DEFAULT_INDEXTYPE ldc)
{
    char trans = (transA==NoTrans) ? 'N' : 'T';

    // TODO: what about constness?
    mkl_dcsrmm(&trans, &m, &n, &k, &alpha, matdescrA,
               const_cast<double *>(values),
               const_cast<FLENS_DEFAULT_INDEXTYPE *>(columns),
               const_cast<FLENS_DEFAULT_INDEXTYPE *>(pointerB),
               const_cast<FLENS_DEFAULT_INDEXTYPE *>(pointerE),
               const_cast<double *>(B),
               &ldb, &beta, C, &ldc);
}
#endif // MKL

} // namespace flens

