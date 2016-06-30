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

#ifndef EXTENSIONS_FLENS_LAPACK_H
#define EXTENSIONS_FLENS_LAPACK_H 1

#include <complex>
#include <flens/vectortypes/impl/densevector.h>
#include <flens/matrixtypes/general/impl/gematrix.h>
#include <cxxblas/drivers/drivers.h>

namespace flens {

using std::complex;

FLENS_DEFAULT_INDEXTYPE
potrf(cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda);

FLENS_DEFAULT_INDEXTYPE
getrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv);

FLENS_DEFAULT_INDEXTYPE
getrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv);

FLENS_DEFAULT_INDEXTYPE
getrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv);

FLENS_DEFAULT_INDEXTYPE
getrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv);

FLENS_DEFAULT_INDEXTYPE
gbtrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, float *ab, FLENS_DEFAULT_INDEXTYPE ldab, FLENS_DEFAULT_INDEXTYPE *ipiv);

FLENS_DEFAULT_INDEXTYPE
gbtrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, double *ab, FLENS_DEFAULT_INDEXTYPE ldab, FLENS_DEFAULT_INDEXTYPE *ipiv);

FLENS_DEFAULT_INDEXTYPE
getri(FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda, const FLENS_DEFAULT_INDEXTYPE *ipiv,
      float *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
getri(FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, const FLENS_DEFAULT_INDEXTYPE *ipiv,
      double *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
getri(FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda, const FLENS_DEFAULT_INDEXTYPE *ipiv,
      complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
getri(FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda, const FLENS_DEFAULT_INDEXTYPE *ipiv,
      complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
getrs(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, const float *a, FLENS_DEFAULT_INDEXTYPE lda,
      const FLENS_DEFAULT_INDEXTYPE *ipiv, float *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
getrs(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, const double *a, FLENS_DEFAULT_INDEXTYPE lda,
      const FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
gbtrs(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs,
      const float *ab, FLENS_DEFAULT_INDEXTYPE ldab, const FLENS_DEFAULT_INDEXTYPE *ipiv, float *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
gbtrs(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs,
      const double *ab, FLENS_DEFAULT_INDEXTYPE ldab, const FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
gesv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, float *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv, float *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
gesv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, double *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
gesv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv,
     complex<float> *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
gesv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv,
     complex<double> *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
gbsv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs, float *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     FLENS_DEFAULT_INDEXTYPE *ipiv, float *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
gbsv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs, double *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
gbsv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs, complex<float> *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     FLENS_DEFAULT_INDEXTYPE *ipiv, complex<float> *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
gbsv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs, complex<double> *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     FLENS_DEFAULT_INDEXTYPE *ipiv, complex<double> *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
trtrs(cxxblas::StorageUpLo upLo, cxxblas::Transpose trans, cxxblas::Diag diag, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs,
      const float *a, FLENS_DEFAULT_INDEXTYPE lda, float *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
trtrs(cxxblas::StorageUpLo upLo, cxxblas::Transpose trans, cxxblas::Diag diag, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs,
      const double *a, FLENS_DEFAULT_INDEXTYPE lda, double *b, FLENS_DEFAULT_INDEXTYPE ldb);

FLENS_DEFAULT_INDEXTYPE
geqrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda, float *tau, float *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
geqrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, double *tau, double *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
geqp3(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *jpvt, float *tau, float *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
geqp3(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *jpvt, double *tau, double *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
orgqr(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k, float *a, FLENS_DEFAULT_INDEXTYPE lda, const float *tau,
      float *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
orgqr(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k, double *a, FLENS_DEFAULT_INDEXTYPE lda, const double *tau,
      double *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
ormqr(cxxblas::Side side, cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k,
      const float *a, FLENS_DEFAULT_INDEXTYPE lda, const float *tau, float *c, FLENS_DEFAULT_INDEXTYPE ldc,
      float *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
ormqr(cxxblas::Side side, cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k,
      const double *a, FLENS_DEFAULT_INDEXTYPE lda, const double *tau, double *c, FLENS_DEFAULT_INDEXTYPE ldc,
      double *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
gels(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, float *a, FLENS_DEFAULT_INDEXTYPE lda,
     float *b, FLENS_DEFAULT_INDEXTYPE ldb, float *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
gels(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, double *a, FLENS_DEFAULT_INDEXTYPE lda,
     double *b, FLENS_DEFAULT_INDEXTYPE ldb, double *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
gels(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda,
     complex<float> *b, FLENS_DEFAULT_INDEXTYPE ldb, complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
gels(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda,
     complex<double> *b, FLENS_DEFAULT_INDEXTYPE ldb, complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
gelss(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, float *a, FLENS_DEFAULT_INDEXTYPE lda, float *b, FLENS_DEFAULT_INDEXTYPE ldb,
     float *s, float rcond, FLENS_DEFAULT_INDEXTYPE rank, float *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
gelss(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, double *a, FLENS_DEFAULT_INDEXTYPE lda, double *b, FLENS_DEFAULT_INDEXTYPE ldb,
     double *s, double rcond, FLENS_DEFAULT_INDEXTYPE rank, double *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
gelss(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda,
      complex<float> *b, FLENS_DEFAULT_INDEXTYPE ldb, complex<float> *s,
      complex<float> rcond, FLENS_DEFAULT_INDEXTYPE rank,
      complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
gelss(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda,
      complex<double> *b, FLENS_DEFAULT_INDEXTYPE ldb, complex<double> *s,
      complex<double> rcond, FLENS_DEFAULT_INDEXTYPE rank,
      complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork);

typedef FLENS_DEFAULT_INDEXTYPE sgees_select(float *vr, float *vi);
typedef FLENS_DEFAULT_INDEXTYPE dgees_select(double *vr, double *vi);
typedef FLENS_DEFAULT_INDEXTYPE cgees_select(complex<float> *v);
typedef FLENS_DEFAULT_INDEXTYPE zgees_select(complex<double> *v);

FLENS_DEFAULT_INDEXTYPE
gees(bool jobvs, bool sort, sgees_select *select,
     FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE &sdim, float *wr, float *wi,
     float *vs, FLENS_DEFAULT_INDEXTYPE ldvs, float *work, FLENS_DEFAULT_INDEXTYPE lwork, FLENS_DEFAULT_INDEXTYPE *bwork);

FLENS_DEFAULT_INDEXTYPE
gees(bool jobvs, bool sort, dgees_select *select,
     FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE &sdim, double *wr, double *wi,
     double *vs, FLENS_DEFAULT_INDEXTYPE ldvs, double *work, FLENS_DEFAULT_INDEXTYPE lwork, FLENS_DEFAULT_INDEXTYPE *bwork);

FLENS_DEFAULT_INDEXTYPE
gees(bool jobvs, bool sort, cgees_select *select,
     FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE &sdim, complex<float> *w,
     complex<float> *vs, FLENS_DEFAULT_INDEXTYPE ldvs,
     complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork, float *rwork, FLENS_DEFAULT_INDEXTYPE *bwork);

FLENS_DEFAULT_INDEXTYPE
gees(bool jobvs, bool sort, zgees_select *select,
     FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE &sdim, complex<double> *w,
     complex<double> *vs, FLENS_DEFAULT_INDEXTYPE ldvs,
     complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork, double *rwork, FLENS_DEFAULT_INDEXTYPE *bwork);

FLENS_DEFAULT_INDEXTYPE
geev(bool jobvl, bool jobvr, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda,
     float *wr, float *wi, float *vl, FLENS_DEFAULT_INDEXTYPE ldvl, float *vr, FLENS_DEFAULT_INDEXTYPE ldvr,
     float *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
geev(bool jobvl, bool jobvr, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda,
     double *wr, double *wi, double *vl, FLENS_DEFAULT_INDEXTYPE ldvl, double *vr, FLENS_DEFAULT_INDEXTYPE ldvr,
     double *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
geev(bool jobvl, bool jobvr, FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda,
     complex<float> *w,
     complex<float> *vl, FLENS_DEFAULT_INDEXTYPE ldvl,
     complex<float> *vr, FLENS_DEFAULT_INDEXTYPE ldvr,
     complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork, float *rwork);

FLENS_DEFAULT_INDEXTYPE
geev(bool jobvl, bool jobvr, FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda,
     complex<double> *w,
     complex<double> *vl, FLENS_DEFAULT_INDEXTYPE ldvl,
     complex<double> *vr, FLENS_DEFAULT_INDEXTYPE ldvr,
     complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork, double *rwork);
     
//- ggev
FLENS_DEFAULT_INDEXTYPE
ggev(bool jobvl, bool jobvr, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, double *b, FLENS_DEFAULT_INDEXTYPE ldb,
     double *wr, double *wi, double *beta, double *vl, FLENS_DEFAULT_INDEXTYPE ldvl, double *vr, FLENS_DEFAULT_INDEXTYPE ldvr,
     double *work, FLENS_DEFAULT_INDEXTYPE lwork);

//- syev
FLENS_DEFAULT_INDEXTYPE
syev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda,
     float *w, float *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
syev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda,
     double *w, double *work, FLENS_DEFAULT_INDEXTYPE lwork);

//- sbev
FLENS_DEFAULT_INDEXTYPE
sbev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kd, float *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     float *w, float *z, FLENS_DEFAULT_INDEXTYPE ldz, float *work);

FLENS_DEFAULT_INDEXTYPE
sbev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kd, double *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     double *w, double *z, FLENS_DEFAULT_INDEXTYPE ldz, double *work);

//- spev
FLENS_DEFAULT_INDEXTYPE
spev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, float *ap, float *w,
     float *z, FLENS_DEFAULT_INDEXTYPE ldz, float *work);

FLENS_DEFAULT_INDEXTYPE
spev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, double *ap, double *w,
     double *z, FLENS_DEFAULT_INDEXTYPE ldz, double *work);

//- heev
FLENS_DEFAULT_INDEXTYPE
heev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda,
     float *w, complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork, float *rwork );

FLENS_DEFAULT_INDEXTYPE
heev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda,
     double *w, complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork, double *rwork );

//- hbev
FLENS_DEFAULT_INDEXTYPE
hbev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kd, complex<float> *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     float *w, complex<float> *Z, FLENS_DEFAULT_INDEXTYPE ldz,
     complex<float> *work, float *rwork);

FLENS_DEFAULT_INDEXTYPE
hbev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kd, complex<double> *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     double *w, complex<double> *Z, FLENS_DEFAULT_INDEXTYPE ldz,
     complex<double> *work, double *rwork);

//- hpev
FLENS_DEFAULT_INDEXTYPE
hpev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, complex<float> *ap, float *w,
     complex<float> *Z, FLENS_DEFAULT_INDEXTYPE ldz, complex<float> *work, float *rwork);

FLENS_DEFAULT_INDEXTYPE
hpev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, complex<double> *ap, double *w,
     complex<double> *Z, FLENS_DEFAULT_INDEXTYPE ldz, complex<double> *work, double *rwork);

enum SVectorsJob {
    All=0,       // A
    SmallDim=1,  // S
    Overwrite=2, // O
    None=3       // N
};

FLENS_DEFAULT_INDEXTYPE
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda,
      float *s,
      float *u, FLENS_DEFAULT_INDEXTYPE ldu,
      float *vt, FLENS_DEFAULT_INDEXTYPE ldvt,
      float *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda,
      double *s,
      double *u, FLENS_DEFAULT_INDEXTYPE ldu,
      double *vt, FLENS_DEFAULT_INDEXTYPE ldvt,
      double *work, FLENS_DEFAULT_INDEXTYPE lwork);

FLENS_DEFAULT_INDEXTYPE
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda,
      float *s,
      complex<float> *u, FLENS_DEFAULT_INDEXTYPE ldu,
      complex<float> *vt, FLENS_DEFAULT_INDEXTYPE ldvt,
      complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork, float *rwork);

FLENS_DEFAULT_INDEXTYPE
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda,
      double *s,
      complex<double> *u, FLENS_DEFAULT_INDEXTYPE ldu,
      complex<double> *vt, FLENS_DEFAULT_INDEXTYPE ldvt,
      complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork, double *rwork);

FLENS_DEFAULT_INDEXTYPE
gesdd(SVectorsJob jobz,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda,
      float *s,
      float *u, FLENS_DEFAULT_INDEXTYPE ldu,
      float *vt, FLENS_DEFAULT_INDEXTYPE ldvt,
      float *work, FLENS_DEFAULT_INDEXTYPE lwork, FLENS_DEFAULT_INDEXTYPE *iwork);

FLENS_DEFAULT_INDEXTYPE
gesdd(SVectorsJob jobz,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda,
      double *s,
      double *u, FLENS_DEFAULT_INDEXTYPE ldu,
      double *vt, FLENS_DEFAULT_INDEXTYPE ldvt,
      double *work, FLENS_DEFAULT_INDEXTYPE lwork, FLENS_DEFAULT_INDEXTYPE *iwork);

FLENS_DEFAULT_INDEXTYPE
gesdd(SVectorsJob jobz,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda,
      float *s,
      complex<float> *u, FLENS_DEFAULT_INDEXTYPE ldu,
      complex<float> *vt, FLENS_DEFAULT_INDEXTYPE ldvt,
      complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork, float *rwork, FLENS_DEFAULT_INDEXTYPE *iwork);

FLENS_DEFAULT_INDEXTYPE
gesdd(SVectorsJob jobz,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda,
      double *s,
      complex<double> *u, FLENS_DEFAULT_INDEXTYPE ldu,
      complex<double> *vt, FLENS_DEFAULT_INDEXTYPE ldvt,
      complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork, double *rwork, FLENS_DEFAULT_INDEXTYPE *iwork);

FLENS_DEFAULT_INDEXTYPE
gecon(char norm,
    FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda,
    double anorm,
    double *rcond,
    double *work, FLENS_DEFAULT_INDEXTYPE *iwork);


} // namespace flens

#endif // EXTENSIONS_FLENS_LAPACK_H

