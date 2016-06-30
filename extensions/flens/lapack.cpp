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
#include <complex>
#include <extensions/flens/lapack.h>

namespace flens {

extern "C" {

    void
    dpotrf_(char *uplo, FLENS_DEFAULT_INDEXTYPE *n, double *ab, FLENS_DEFAULT_INDEXTYPE *ldab, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgetrf_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, float *a, FLENS_DEFAULT_INDEXTYPE *lda, FLENS_DEFAULT_INDEXTYPE *ipiv, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgetrf_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, double *a, FLENS_DEFAULT_INDEXTYPE *lda, FLENS_DEFAULT_INDEXTYPE *ipiv, FLENS_DEFAULT_INDEXTYPE *info);

    void
    cgetrf_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE *lda, FLENS_DEFAULT_INDEXTYPE *ipiv, FLENS_DEFAULT_INDEXTYPE *info);

    void
    zgetrf_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE *lda, FLENS_DEFAULT_INDEXTYPE *ipiv, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgbtrf_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kl, FLENS_DEFAULT_INDEXTYPE *ku,
            float *ab, FLENS_DEFAULT_INDEXTYPE *ldab, FLENS_DEFAULT_INDEXTYPE *ipiv, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgbtrf_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kl, FLENS_DEFAULT_INDEXTYPE *ku,
            double *ab, FLENS_DEFAULT_INDEXTYPE *ldab, FLENS_DEFAULT_INDEXTYPE *ipiv, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgetri_(FLENS_DEFAULT_INDEXTYPE *n, float *a, FLENS_DEFAULT_INDEXTYPE *lda, const FLENS_DEFAULT_INDEXTYPE *ipiv, float *work,
            FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgetri_(FLENS_DEFAULT_INDEXTYPE *n, double *a, FLENS_DEFAULT_INDEXTYPE *lda, const FLENS_DEFAULT_INDEXTYPE *ipiv, double *work,
            FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    cgetri_(FLENS_DEFAULT_INDEXTYPE *n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE *lda, const FLENS_DEFAULT_INDEXTYPE *ipiv,
            complex<float> *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    zgetri_(FLENS_DEFAULT_INDEXTYPE *n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE *lda, const FLENS_DEFAULT_INDEXTYPE *ipiv,
            complex<double> *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgetrs_(char *trans, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs, const float *a, FLENS_DEFAULT_INDEXTYPE *lda,
            const FLENS_DEFAULT_INDEXTYPE *ipiv, float *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgetrs_(char *trans, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs, const double *a, FLENS_DEFAULT_INDEXTYPE *lda,
            const FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgbtrs_(char *trans, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kl, FLENS_DEFAULT_INDEXTYPE *ku, FLENS_DEFAULT_INDEXTYPE *nrhs, const float *ab,
            FLENS_DEFAULT_INDEXTYPE *ldab, const FLENS_DEFAULT_INDEXTYPE *ipiv, float *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgbtrs_(char *trans, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kl, FLENS_DEFAULT_INDEXTYPE *ku, FLENS_DEFAULT_INDEXTYPE *nrhs, const double *ab,
            FLENS_DEFAULT_INDEXTYPE *ldab, const FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgesv_(FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs, float *a, FLENS_DEFAULT_INDEXTYPE *lda,
           FLENS_DEFAULT_INDEXTYPE *ipiv, float *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgesv_(FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs, double *a, FLENS_DEFAULT_INDEXTYPE *lda,
           FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    cgesv_(FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs, complex<float> *a, FLENS_DEFAULT_INDEXTYPE *lda,
           FLENS_DEFAULT_INDEXTYPE *ipiv, complex<float> *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    zgesv_(FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs, complex<double> *a, FLENS_DEFAULT_INDEXTYPE *lda,
           FLENS_DEFAULT_INDEXTYPE *ipiv, complex<double> *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgbsv_(FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kl, FLENS_DEFAULT_INDEXTYPE *ku, FLENS_DEFAULT_INDEXTYPE *nrhs, float *ab, FLENS_DEFAULT_INDEXTYPE *ldab, FLENS_DEFAULT_INDEXTYPE *ipiv,
           float *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgbsv_(FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kl, FLENS_DEFAULT_INDEXTYPE *ku, FLENS_DEFAULT_INDEXTYPE *nrhs, double *ab, FLENS_DEFAULT_INDEXTYPE *ldab,
           FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    cgbsv_(FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kl, FLENS_DEFAULT_INDEXTYPE *ku, FLENS_DEFAULT_INDEXTYPE *nrhs, complex<float> *ab,
           FLENS_DEFAULT_INDEXTYPE *ldab, FLENS_DEFAULT_INDEXTYPE *ipiv, complex<float> *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    zgbsv_(FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kl, FLENS_DEFAULT_INDEXTYPE *ku, FLENS_DEFAULT_INDEXTYPE *nrhs, complex<double> *ab,
           FLENS_DEFAULT_INDEXTYPE *ldab, FLENS_DEFAULT_INDEXTYPE *ipiv, complex<double> *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    strtrs_(char *uplo, char *trans, char *diag, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs,
            const float *a, FLENS_DEFAULT_INDEXTYPE *lda, float *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dtrtrs_(char *uplo, char *trans, char *diag, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs,
            const double *a, FLENS_DEFAULT_INDEXTYPE *lda, double *b, FLENS_DEFAULT_INDEXTYPE *ldb, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgeqrf_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, float *a, FLENS_DEFAULT_INDEXTYPE *lda, float *tau,
            float *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgeqrf_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, double *a, FLENS_DEFAULT_INDEXTYPE *lda, double *tau,
            double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgeqp3_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, float *a, FLENS_DEFAULT_INDEXTYPE *lda, FLENS_DEFAULT_INDEXTYPE* jpvt, float *tau,
            float *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgeqp3_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, double *a, FLENS_DEFAULT_INDEXTYPE *lda, FLENS_DEFAULT_INDEXTYPE* jpvt, double *tau,
            double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sorgqr_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *k, float *a, FLENS_DEFAULT_INDEXTYPE *lda, const float *tau,
            float *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dorgqr_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *k, double *a, FLENS_DEFAULT_INDEXTYPE *lda, const double *tau,
            double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sormqr_(char *side, char *trans, FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *k,
            const float *a, FLENS_DEFAULT_INDEXTYPE *lda, const float *tau, float *c, FLENS_DEFAULT_INDEXTYPE *ldc,
            float *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dormqr_(char *side, char *trans, FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *k,
            const double *a, FLENS_DEFAULT_INDEXTYPE *lda, const double *tau, double *c, FLENS_DEFAULT_INDEXTYPE *ldc,
            double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgels_(char *trans, FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs, float *a, FLENS_DEFAULT_INDEXTYPE *lda,
           float *b, FLENS_DEFAULT_INDEXTYPE *ldb, float *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgels_(char *trans, FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs, double *a, FLENS_DEFAULT_INDEXTYPE *lda,
           double *b, FLENS_DEFAULT_INDEXTYPE *ldb, double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    cgels_(char *trans, FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs, complex<float> *a,
           FLENS_DEFAULT_INDEXTYPE *lda, complex<float> *b, FLENS_DEFAULT_INDEXTYPE *ldb,
           complex<float> *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    zgels_(char *trans, FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs, complex<double> *a,
           FLENS_DEFAULT_INDEXTYPE *lda, complex<double> *b, FLENS_DEFAULT_INDEXTYPE *ldb,
           complex<double> *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgelss_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs,
            float *a, FLENS_DEFAULT_INDEXTYPE *lda, float *b, FLENS_DEFAULT_INDEXTYPE *ldb,
            float *s, float *rcond, FLENS_DEFAULT_INDEXTYPE *rank,
            float *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgelss_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs,
            double *a, FLENS_DEFAULT_INDEXTYPE *lda, double *b, FLENS_DEFAULT_INDEXTYPE *ldb,
            double *s, double *rcond, FLENS_DEFAULT_INDEXTYPE *rank,
            double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    cgelss_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs,
            complex<float> *a, FLENS_DEFAULT_INDEXTYPE *lda,
            complex<float> *b, FLENS_DEFAULT_INDEXTYPE *ldb,
            complex<float> *s, complex<float> *rcond, FLENS_DEFAULT_INDEXTYPE *rank,
            complex<float> *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    zgelss_(FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *nrhs,
            complex<double> *a, FLENS_DEFAULT_INDEXTYPE *lda,
            complex<double> *b, FLENS_DEFAULT_INDEXTYPE *ldb,
            complex<double> *s, complex<double> *rcond, FLENS_DEFAULT_INDEXTYPE *rank,
            complex<double> *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgees_(char * jobvs, char * sort, sgees_select *select,
           FLENS_DEFAULT_INDEXTYPE *n, float *a, FLENS_DEFAULT_INDEXTYPE *lda, FLENS_DEFAULT_INDEXTYPE *sdim, float *wr, float *wi,
           float *vs, FLENS_DEFAULT_INDEXTYPE *ldvs, float *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *bwork,
           FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgees_(char * jobvs, char * sort, dgees_select *select,
           FLENS_DEFAULT_INDEXTYPE *n, double *a, FLENS_DEFAULT_INDEXTYPE *lda, FLENS_DEFAULT_INDEXTYPE *sdim, double *wr, double *wi,
           double *vs, FLENS_DEFAULT_INDEXTYPE *ldvs, double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *bwork,
           FLENS_DEFAULT_INDEXTYPE *info);

    void
    cgees_(char * jobvs, char * sort, cgees_select *select,
           FLENS_DEFAULT_INDEXTYPE *n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE *lda, FLENS_DEFAULT_INDEXTYPE *sdim, complex<float> *w,
           complex<float> *vs, FLENS_DEFAULT_INDEXTYPE *ldvs,
           complex<float> *work, FLENS_DEFAULT_INDEXTYPE *lwork, float *rwork, FLENS_DEFAULT_INDEXTYPE *bwork,
           FLENS_DEFAULT_INDEXTYPE *info);

    void
    zgees_(char * jobvs, char * sort, zgees_select *select,
           FLENS_DEFAULT_INDEXTYPE *n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE *lda, FLENS_DEFAULT_INDEXTYPE *sdim, complex<double> *w,
           complex<double> *vs, FLENS_DEFAULT_INDEXTYPE *ldvs,
           complex<double> *work, FLENS_DEFAULT_INDEXTYPE *lwork, double *rwork, FLENS_DEFAULT_INDEXTYPE *bwork,
           FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgeev_(char *jobvl, char *jobvr, FLENS_DEFAULT_INDEXTYPE *n, float *a, FLENS_DEFAULT_INDEXTYPE *lda,
           float *wr, float *wi,
           float *vl, FLENS_DEFAULT_INDEXTYPE *ldvl,
           float *vr, FLENS_DEFAULT_INDEXTYPE *ldvr,
           float *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgeev_(char *jobvl, char *jobvr, FLENS_DEFAULT_INDEXTYPE *n, double *a, FLENS_DEFAULT_INDEXTYPE *lda,
           double *wr, double *wi,
           double *vl, FLENS_DEFAULT_INDEXTYPE *ldvl,
           double *vr, FLENS_DEFAULT_INDEXTYPE *ldvr,
           double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    cgeev_(char *jobvl, char *jobvr, FLENS_DEFAULT_INDEXTYPE *n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE *lda,
           complex<float> *w,
           complex<float> *vl, FLENS_DEFAULT_INDEXTYPE *ldvl,
           complex<float> *vr, FLENS_DEFAULT_INDEXTYPE *ldvr,
           complex<float> *work,FLENS_DEFAULT_INDEXTYPE *lwork,float *rwork,FLENS_DEFAULT_INDEXTYPE *info);

    void
    zgeev_(char *jobvl, char *jobvr, FLENS_DEFAULT_INDEXTYPE *n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE *lda,
           complex<double> *w,
           complex<double> *vl, FLENS_DEFAULT_INDEXTYPE *ldvl,
           complex<double> *vr, FLENS_DEFAULT_INDEXTYPE *ldvr,
           complex<double> *work,FLENS_DEFAULT_INDEXTYPE *lwork,double *rwork,FLENS_DEFAULT_INDEXTYPE *info);
           
    void
    dggev_(char *jobvl, char *jobvr, FLENS_DEFAULT_INDEXTYPE *n, double *a, FLENS_DEFAULT_INDEXTYPE *lda, 
           double *b, FLENS_DEFAULT_INDEXTYPE *ldb,
           double *wr, double *wi, 
           double *beta, 
           double *vl, FLENS_DEFAULT_INDEXTYPE *ldvl, 
           double *vr, FLENS_DEFAULT_INDEXTYPE *ldvr,
           double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    ssyev_(char *jobz, char *uplo, FLENS_DEFAULT_INDEXTYPE *n, float *a, FLENS_DEFAULT_INDEXTYPE *lda,
           float *w, float *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dsyev_(char *jobz, char *uplo, FLENS_DEFAULT_INDEXTYPE *n, double *a, FLENS_DEFAULT_INDEXTYPE *lda,
           double *w, double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    ssbev_(char *jobz, char *upLo, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kd, float *ab, FLENS_DEFAULT_INDEXTYPE *ldab,
           float *w, float *z, FLENS_DEFAULT_INDEXTYPE *ldz, float *work, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dsbev_(char *jobz, char *upLo, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kd, double *ab, FLENS_DEFAULT_INDEXTYPE *ldab,
           double *w, double *z, FLENS_DEFAULT_INDEXTYPE *ldz, double *work, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sspev_(char *jobz, char *upLo, FLENS_DEFAULT_INDEXTYPE *n, float *ap, float *w,
           float *z, FLENS_DEFAULT_INDEXTYPE *ldz, float *work, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dspev_(char *jobz, char *upLo, FLENS_DEFAULT_INDEXTYPE *n, double *ap, double *w,
           double *z, FLENS_DEFAULT_INDEXTYPE *ldz, double *work, FLENS_DEFAULT_INDEXTYPE *info);

    void
    cheev_(char *jobz, char *uplo, FLENS_DEFAULT_INDEXTYPE *n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE *lda,
           float *w, complex<float> *work, FLENS_DEFAULT_INDEXTYPE *lwork, float *rwork,
           FLENS_DEFAULT_INDEXTYPE *info);

    void
    zheev_(char *jobz, char *uplo, FLENS_DEFAULT_INDEXTYPE *n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE *lda,
           double *w, complex<double> *work, FLENS_DEFAULT_INDEXTYPE *lwork, double *rwork,
           FLENS_DEFAULT_INDEXTYPE *info);

    FLENS_DEFAULT_INDEXTYPE
    chbev_(char *jobz, char *upLo, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kd, complex<float> *ab,
           FLENS_DEFAULT_INDEXTYPE *ldab, float *w, complex<float> *Z, FLENS_DEFAULT_INDEXTYPE *ldz,
           complex<float> *work, float *rwork, FLENS_DEFAULT_INDEXTYPE *info);

    FLENS_DEFAULT_INDEXTYPE
    zhbev_(char *jobz, char *upLo, FLENS_DEFAULT_INDEXTYPE *n, FLENS_DEFAULT_INDEXTYPE *kd, complex<double> *ab,
           FLENS_DEFAULT_INDEXTYPE *ldab, double *w, complex<double> *Z, FLENS_DEFAULT_INDEXTYPE *ldz,
           complex<double> *work, double *rwork, FLENS_DEFAULT_INDEXTYPE *info);

    FLENS_DEFAULT_INDEXTYPE
    chpev_(char *jobz, char *upLo, FLENS_DEFAULT_INDEXTYPE *n, complex<float> *ap,
           float *w, complex<float> *Z, FLENS_DEFAULT_INDEXTYPE *ldz,
           complex<float> *work, float *rwork, FLENS_DEFAULT_INDEXTYPE *info);

    FLENS_DEFAULT_INDEXTYPE
    zhpev_(char *jobz, char *upLo, FLENS_DEFAULT_INDEXTYPE *n, complex<double> *ap,
           double *w, complex<double> *Z, FLENS_DEFAULT_INDEXTYPE *ldz,
           complex<double> *work, double *rwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgesvd_(char *jobu, char *jobvt, FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, float *a, FLENS_DEFAULT_INDEXTYPE *lda,
            float *s, float *u, FLENS_DEFAULT_INDEXTYPE *ldu, float *vt, FLENS_DEFAULT_INDEXTYPE *ldvt,
            float *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgesvd_(char *jobu, char *jobvt, FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, double *a, FLENS_DEFAULT_INDEXTYPE *lda,
            double *s, double *u, FLENS_DEFAULT_INDEXTYPE *ldu, double *vt, FLENS_DEFAULT_INDEXTYPE *ldvt,
            double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    cgesvd_(char *jobu, char *jobvt,
            FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE *lda, float *s,
            complex<float> *u, FLENS_DEFAULT_INDEXTYPE *ldu,
            complex<float> *vt, FLENS_DEFAULT_INDEXTYPE *ldvt,
            complex<float> *work, FLENS_DEFAULT_INDEXTYPE *lwork, float *rwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    zgesvd_(char *jobu, char *jobvt,
            FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE *lda, double *s,
            complex<double> *u, FLENS_DEFAULT_INDEXTYPE *ldu,
            complex<double> *vt, FLENS_DEFAULT_INDEXTYPE *ldvt,
            complex<double> *work, FLENS_DEFAULT_INDEXTYPE *lwork, double *rwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    sgesdd_(char *jobz, FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, float *a, FLENS_DEFAULT_INDEXTYPE *lda,
            float *s, float *u, FLENS_DEFAULT_INDEXTYPE *ldu, float *vt, FLENS_DEFAULT_INDEXTYPE *ldvt,
            float *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *iwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgesdd_(char *jobz, FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, double *a, FLENS_DEFAULT_INDEXTYPE *lda,
            double *s, double *u, FLENS_DEFAULT_INDEXTYPE *ldu, double *vt, FLENS_DEFAULT_INDEXTYPE *ldvt,
            double *work, FLENS_DEFAULT_INDEXTYPE *lwork, FLENS_DEFAULT_INDEXTYPE *iwork, FLENS_DEFAULT_INDEXTYPE *info);

    void
    cgesdd_(char *jobz,
            FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE *lda, float *s,
            complex<float> *u, FLENS_DEFAULT_INDEXTYPE *ldu,
            complex<float> *vt, FLENS_DEFAULT_INDEXTYPE *ldvt,
            complex<float> *work, FLENS_DEFAULT_INDEXTYPE *lwork, float *rwork, FLENS_DEFAULT_INDEXTYPE *iwork,
            FLENS_DEFAULT_INDEXTYPE *info);

    void
    zgesdd_(char *jobzu,
            FLENS_DEFAULT_INDEXTYPE *m, FLENS_DEFAULT_INDEXTYPE *n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE *lda, double *s,
            complex<double> *u, FLENS_DEFAULT_INDEXTYPE *ldu,
            complex<double> *vt, FLENS_DEFAULT_INDEXTYPE *ldvt,
            complex<double> *work, FLENS_DEFAULT_INDEXTYPE *lwork, double *rwork, FLENS_DEFAULT_INDEXTYPE *iwork,
            FLENS_DEFAULT_INDEXTYPE *info);

    void
    dgecon_(char *norm, FLENS_DEFAULT_INDEXTYPE *n, double *a, FLENS_DEFAULT_INDEXTYPE *lda, double *anorm,
            double *rcond, double *work, FLENS_DEFAULT_INDEXTYPE *iwork, FLENS_DEFAULT_INDEXTYPE *info);
}

FLENS_DEFAULT_INDEXTYPE
potrf(cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';
    dpotrf_(&_upLo,&n,a,&lda,&info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
getrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv)
{
    FLENS_DEFAULT_INDEXTYPE info;
    sgetrf_(&m, &n, a, &lda, ipiv, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
getrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dgetrf_(&m, &n, a, &lda, ipiv, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
getrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv)
{
    FLENS_DEFAULT_INDEXTYPE info;
    cgetrf_(&m, &n, a, &lda, ipiv, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
getrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv)
{
    FLENS_DEFAULT_INDEXTYPE info;
    zgetrf_(&m, &n, a, &lda, ipiv, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gbtrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, float *ab, FLENS_DEFAULT_INDEXTYPE ldab, FLENS_DEFAULT_INDEXTYPE *ipiv)
{
    FLENS_DEFAULT_INDEXTYPE info;
    sgbtrf_(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gbtrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, double *ab, FLENS_DEFAULT_INDEXTYPE ldab, FLENS_DEFAULT_INDEXTYPE *ipiv)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dgbtrf_(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
getri(FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda, const FLENS_DEFAULT_INDEXTYPE *ipiv,
      float *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    sgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
getri(FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, const FLENS_DEFAULT_INDEXTYPE *ipiv,
      double *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
getri(FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda, const FLENS_DEFAULT_INDEXTYPE *ipiv,
      complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    cgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}


FLENS_DEFAULT_INDEXTYPE
getri(FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda, const FLENS_DEFAULT_INDEXTYPE *ipiv,
      complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    zgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
getrs(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, const float *a, FLENS_DEFAULT_INDEXTYPE lda,
      const FLENS_DEFAULT_INDEXTYPE *ipiv, float *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';

    sgetrs_(&_trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
getrs(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, const double *a, FLENS_DEFAULT_INDEXTYPE lda,
      const FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';

    dgetrs_(&_trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gbtrs(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs,
      const float *ab, FLENS_DEFAULT_INDEXTYPE ldab, const FLENS_DEFAULT_INDEXTYPE *ipiv, float *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';
    sgbtrs_(&_trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gbtrs(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs,
      const double *ab, FLENS_DEFAULT_INDEXTYPE ldab, const FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';
    dgbtrs_(&_trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gesv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, float *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv, float *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gesv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, double *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gesv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv,
     complex<float> *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    cgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gesv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *ipiv,
     complex<double> *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    zgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gbsv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs, float *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     FLENS_DEFAULT_INDEXTYPE *ipiv, float *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    sgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gbsv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs, double *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     FLENS_DEFAULT_INDEXTYPE *ipiv, double *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gbsv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs, complex<float> *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     FLENS_DEFAULT_INDEXTYPE *ipiv, complex<float> *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    cgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gbsv(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kl, FLENS_DEFAULT_INDEXTYPE ku, FLENS_DEFAULT_INDEXTYPE nrhs, complex<double> *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     FLENS_DEFAULT_INDEXTYPE *ipiv, complex<double> *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    zgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
trtrs(cxxblas::StorageUpLo upLo, cxxblas::Transpose trans, cxxblas::Diag diag, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs,
      const float *a, FLENS_DEFAULT_INDEXTYPE lda, float *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';
    char _diag = (diag==cxxblas::Unit) ? 'U' : 'N';

    strtrs_(&_upLo, &_trans, &_diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
trtrs(cxxblas::StorageUpLo upLo, cxxblas::Transpose trans, cxxblas::Diag diag, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs,
      const double *a, FLENS_DEFAULT_INDEXTYPE lda, double *b, FLENS_DEFAULT_INDEXTYPE ldb)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';
    char _diag = (diag==cxxblas::Unit) ? 'U' : 'N';

    dtrtrs_(&_upLo, &_trans, &_diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
geqrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda, float *tau, float *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    sgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
geqrf(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, double *tau, double *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
geqp3(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *jpvt, float *tau, float *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    sgeqp3_(&m, &n, a, &lda, jpvt, tau, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
geqp3(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE *jpvt, double *tau, double *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dgeqp3_(&m, &n, a, &lda, jpvt, tau, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
orgqr(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k, float *a, FLENS_DEFAULT_INDEXTYPE lda, const float *tau,
      float *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    sorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
orgqr(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k, double *a, FLENS_DEFAULT_INDEXTYPE lda, const double *tau,
      double *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
ormqr(cxxblas::Side side, cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k,
      const float *a, FLENS_DEFAULT_INDEXTYPE lda, const float *tau, float *c, FLENS_DEFAULT_INDEXTYPE ldc,
      float *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _side = (side==cxxblas::Left) ? 'L' : 'R';
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';

    sormqr_(&_side, &_trans, &m, &n, &k, a, &lda, tau,
            c, &ldc, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
ormqr(cxxblas::Side side, cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k,
      const double *a, FLENS_DEFAULT_INDEXTYPE lda, const double *tau, double *c, FLENS_DEFAULT_INDEXTYPE ldc,
      double *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _side = (side==cxxblas::Left) ? 'L' : 'R';
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';

    dormqr_(&_side, &_trans, &m, &n, &k, a, &lda, tau,
            c, &ldc, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gels(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, float *a, FLENS_DEFAULT_INDEXTYPE lda,
     float *b, FLENS_DEFAULT_INDEXTYPE ldb, float *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';

    sgels_(&_trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gels(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, double *a, FLENS_DEFAULT_INDEXTYPE lda,
     double *b, FLENS_DEFAULT_INDEXTYPE ldb, double *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';

    dgels_(&_trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gels(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda,
     complex<float> *b, FLENS_DEFAULT_INDEXTYPE ldb, complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';

    cgels_(&_trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gels(cxxblas::Transpose trans, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda,
     complex<double> *b, FLENS_DEFAULT_INDEXTYPE ldb, complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _trans = (trans==cxxblas::NoTrans) ? 'N' : 'T';

    zgels_(&_trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gelss(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, float *a, FLENS_DEFAULT_INDEXTYPE lda, float *b, FLENS_DEFAULT_INDEXTYPE ldb,
     float *s, float rcond, FLENS_DEFAULT_INDEXTYPE rank, float *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    sgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
            work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gelss(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, double *a, FLENS_DEFAULT_INDEXTYPE lda, double *b, FLENS_DEFAULT_INDEXTYPE ldb,
     double *s, double rcond, FLENS_DEFAULT_INDEXTYPE rank, double *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
            work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gelss(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda,
      complex<float> *b, FLENS_DEFAULT_INDEXTYPE ldb, complex<float> *s,
      complex<float> rcond, FLENS_DEFAULT_INDEXTYPE rank, complex<float> *work,
      FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    cgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
            work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gelss(FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE nrhs, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda,
      complex<double> *b, FLENS_DEFAULT_INDEXTYPE ldb, complex<double> *s,
      complex<double> rcond, FLENS_DEFAULT_INDEXTYPE rank, complex<double> *work,
      FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    zgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
            work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
geev(bool jobvl, bool jobvr, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda,
     float *wr, float *wi,
     float *vl, FLENS_DEFAULT_INDEXTYPE ldvl,
     float *vr, FLENS_DEFAULT_INDEXTYPE ldvr,
     float *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    sgeev_(&_jobvl, &_jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
           work, &lwork, &info);
    return info;
}

//- syev
FLENS_DEFAULT_INDEXTYPE
syev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda,
     float *w, float *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    ssyev_(&_jobz, &_upLo, &n, a, &lda, w, work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
syev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda,
     double *w, double *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    dsyev_(&_jobz, &_upLo, &n, a, &lda, w, work, &lwork, &info);
    return info;
}

//- sbev
FLENS_DEFAULT_INDEXTYPE
sbev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kd, float *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     float *w, float *z, FLENS_DEFAULT_INDEXTYPE ldz, float *work)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    ssbev_(&_jobz, &_upLo, &n, &kd, ab, &ldab, w, z, &ldz, work, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
sbev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kd, double *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     double *w, double *z, FLENS_DEFAULT_INDEXTYPE ldz, double *work)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    dsbev_(&_jobz, &_upLo, &n, &kd, ab, &ldab, w, z, &ldz, work, &info);
    return info;
}

//- spev
FLENS_DEFAULT_INDEXTYPE
spev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, float *ap, float *w,
     float *z, FLENS_DEFAULT_INDEXTYPE ldz, float *work)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    sspev_(&_jobz, &_upLo, &n, ap, w, z, &ldz, work, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
spev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, double *ap, double *w,
     double *z, FLENS_DEFAULT_INDEXTYPE ldz, double *work)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    dspev_(&_jobz, &_upLo, &n, ap, w, z, &ldz, work, &info);
    return info;
}

//- heev
FLENS_DEFAULT_INDEXTYPE
heev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda, float *w,
     complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork, float *rwork )
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    cheev_(&_jobz, &_upLo, &n, a, &lda, w, work, &lwork, rwork, &info);

    return info;
}

FLENS_DEFAULT_INDEXTYPE
heev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda, double *w,
     complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork, double *rwork )
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    zheev_(&_jobz, &_upLo, &n, a, &lda, w, work, &lwork, rwork, &info);

    return info;
}

//- hbev
FLENS_DEFAULT_INDEXTYPE
hbev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kd, complex<float> *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     float *w, complex<float> *Z, FLENS_DEFAULT_INDEXTYPE ldz,
     complex<float> *work, float *rwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    chbev_(&_jobz, &_upLo, &n, &kd, ab, &ldab, w, Z, &ldz, work, rwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
hbev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE kd, complex<double> *ab, FLENS_DEFAULT_INDEXTYPE ldab,
     double *w, complex<double> *Z, FLENS_DEFAULT_INDEXTYPE ldz,
     complex<double> *work, double *rwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    zhbev_(&_jobz, &_upLo, &n, &kd, ab, &ldab, w, Z, &ldz, work, rwork, &info);
    return info;
}

//- hpev
FLENS_DEFAULT_INDEXTYPE
hpev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, complex<float> *ap, float *w,
     complex<float> *Z, FLENS_DEFAULT_INDEXTYPE ldz, complex<float> *work, float *rwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    chpev_(&_jobz, &_upLo, &n, ap, w, Z, &ldz, work, rwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
hpev(bool jobz, cxxblas::StorageUpLo upLo, FLENS_DEFAULT_INDEXTYPE n, complex<double> *ap, double *w,
     complex<double> *Z, FLENS_DEFAULT_INDEXTYPE ldz, complex<double> *work, double *rwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==cxxblas::Upper) ? 'U' : 'L';

    zhpev_(&_jobz, &_upLo, &n, ap, w, Z, &ldz, work, rwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gees(bool jobvs, bool sort, sgees_select *select,
     FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE &sdim, float *wr, float *wi,
     float *vs, FLENS_DEFAULT_INDEXTYPE ldvs, float *work, FLENS_DEFAULT_INDEXTYPE lwork, FLENS_DEFAULT_INDEXTYPE *bwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobvs = (jobvs==true) ? 'V' : 'N';
    char _sort  = (sort==true) ? 'S' : 'N';

    sgees_(&_jobvs, &_sort, select, &n, a, &lda, &sdim,
           wr, wi, vs, &ldvs,
           work, &lwork, bwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gees(bool jobvs, bool sort, dgees_select *select,
     FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE &sdim, double *wr, double *wi,
     double *vs, FLENS_DEFAULT_INDEXTYPE ldvs, double *work, FLENS_DEFAULT_INDEXTYPE lwork, FLENS_DEFAULT_INDEXTYPE *bwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobvs = (jobvs==true) ? 'V' : 'N';
    char _sort  = (sort==true) ? 'S' : 'N';

    dgees_(&_jobvs, &_sort, select, &n, a, &lda, &sdim,
           wr, wi, vs, &ldvs,
           work, &lwork, bwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gees(bool jobvs, bool sort, cgees_select *select,
     FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE &sdim, complex<float> *w,
     complex<float> *vs, FLENS_DEFAULT_INDEXTYPE ldvs,
     complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork, float *rwork, FLENS_DEFAULT_INDEXTYPE *bwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobvs = (jobvs==true) ? 'V' : 'N';
    char _sort  = (sort==true) ? 'S' : 'N';

    cgees_(&_jobvs, &_sort, select, &n, a, &lda, &sdim,
           w, vs, &ldvs,
           work, &lwork, rwork, bwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gees(bool jobvs, bool sort, zgees_select *select,
     FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda, FLENS_DEFAULT_INDEXTYPE &sdim, complex<double> *w,
     complex<double> *vs, FLENS_DEFAULT_INDEXTYPE ldvs,
     complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork, double *rwork, FLENS_DEFAULT_INDEXTYPE *bwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobvs = (jobvs==true) ? 'V' : 'N';
    char _sort  = (sort==true) ? 'S' : 'N';

    zgees_(&_jobvs, &_sort, select, &n, a, &lda, &sdim,
           w, vs, &ldvs,
           work, &lwork, rwork, bwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
geev(bool jobvl, bool jobvr, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda,
     double *wr, double *wi,
     double *vl, FLENS_DEFAULT_INDEXTYPE ldvl,
     double *vr, FLENS_DEFAULT_INDEXTYPE ldvr,
     double *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    dgeev_(&_jobvl, &_jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
           work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
geev(bool jobvl, bool jobvr, FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda,
     complex<float> *w,
     complex<float> *vl, FLENS_DEFAULT_INDEXTYPE ldvl,
     complex<float> *vr, FLENS_DEFAULT_INDEXTYPE ldvr,
     complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork, float *rwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    cgeev_(&_jobvl, &_jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
           work, &lwork, rwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
geev(bool jobvl, bool jobvr, FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda,
     complex<double> *w,
     complex<double> *vl, FLENS_DEFAULT_INDEXTYPE ldvl,
     complex<double> *vr, FLENS_DEFAULT_INDEXTYPE ldvr,
     complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork, double *rwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    zgeev_(&_jobvl, &_jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
           work, &lwork, rwork, &info);
    return info;
}

//- ggev
FLENS_DEFAULT_INDEXTYPE
ggev(bool jobvl, bool jobvr, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda, 
     double *b, FLENS_DEFAULT_INDEXTYPE ldb,
     double *wr, double *wi, 
     double *beta, 
     double *vl, FLENS_DEFAULT_INDEXTYPE ldvl, 
     double *vr, FLENS_DEFAULT_INDEXTYPE ldvr,
     double *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';
    
    dggev_(&_jobvl, &_jobvr, &n, a, &lda, b, &ldb, wr, wi, beta, vl, &ldvl, vr, &ldvr,
           work, &lwork, &info);
    return info;
}

static char jobchar[4] = {'A','S','O','N'};

FLENS_DEFAULT_INDEXTYPE
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda,                // A
      float *s,                                      // singular values
      float *u, FLENS_DEFAULT_INDEXTYPE ldu,                         // left singular vectors
      float *vt, FLENS_DEFAULT_INDEXTYPE ldvt,                       // right singular vectors
      float *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    sgesvd_(&jobchar[jobu], &jobchar[jobvt],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda,                // A
      double *s,                                      // singular values
      double *u, FLENS_DEFAULT_INDEXTYPE ldu,                         // left singular vectors
      double *vt, FLENS_DEFAULT_INDEXTYPE ldvt,                       // right singular vectors
      double *work, FLENS_DEFAULT_INDEXTYPE lwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dgesvd_(&jobchar[jobu], &jobchar[jobvt],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda,       // A
      float *s,                                      // singular values
      complex<float> *u, FLENS_DEFAULT_INDEXTYPE ldu,                // left singular vectors
      complex<float> *vt, FLENS_DEFAULT_INDEXTYPE ldvt,              // right singular vectors
      complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork, float *rwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    cgesvd_(&jobchar[jobu], &jobchar[jobvt],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, rwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda,       // A
      double *s,                                      // singular values
      complex<double> *u, FLENS_DEFAULT_INDEXTYPE ldu,                // left singular vectors
      complex<double> *vt, FLENS_DEFAULT_INDEXTYPE ldvt,              // right singular vectors
      complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork, double *rwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    zgesvd_(&jobchar[jobu], &jobchar[jobvt],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, rwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gesdd(SVectorsJob jobz,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, float *a, FLENS_DEFAULT_INDEXTYPE lda,                // A
      float *s,                                      // singular values
      float *u, FLENS_DEFAULT_INDEXTYPE ldu,                         // left singular vectors
      float *vt, FLENS_DEFAULT_INDEXTYPE ldvt,                       // right singular vectors
      float *work, FLENS_DEFAULT_INDEXTYPE lwork, FLENS_DEFAULT_INDEXTYPE *iwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    sgesdd_(&jobchar[jobz],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, iwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gesdd(SVectorsJob jobz,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda,                // A
      double *s,                                      // singular values
      double *u, FLENS_DEFAULT_INDEXTYPE ldu,                         // left singular vectors
      double *vt, FLENS_DEFAULT_INDEXTYPE ldvt,                       // right singular vectors
      double *work, FLENS_DEFAULT_INDEXTYPE lwork, FLENS_DEFAULT_INDEXTYPE *iwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dgesdd_(&jobchar[jobz],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, iwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gesdd(SVectorsJob jobz,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<float> *a, FLENS_DEFAULT_INDEXTYPE lda,       // A
      float *s,                                      // singular values
      complex<float> *u, FLENS_DEFAULT_INDEXTYPE ldu,                // left singular vectors
      complex<float> *vt, FLENS_DEFAULT_INDEXTYPE ldvt,              // right singular vectors
      complex<float> *work, FLENS_DEFAULT_INDEXTYPE lwork, float *rwork, FLENS_DEFAULT_INDEXTYPE *iwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    cgesdd_(&jobchar[jobz],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, rwork, iwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gesdd(SVectorsJob jobz,
      FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, complex<double> *a, FLENS_DEFAULT_INDEXTYPE lda,       // A
      double *s,                                      // singular values
      complex<double> *u, FLENS_DEFAULT_INDEXTYPE ldu,                // left singular vectors
      complex<double> *vt, FLENS_DEFAULT_INDEXTYPE ldvt,              // right singular vectors
      complex<double> *work, FLENS_DEFAULT_INDEXTYPE lwork, double *rwork, FLENS_DEFAULT_INDEXTYPE *iwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    zgesdd_(&jobchar[jobz],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, rwork, iwork, &info);
    return info;
}

FLENS_DEFAULT_INDEXTYPE
gecon(char norm,
      FLENS_DEFAULT_INDEXTYPE n, double *a, FLENS_DEFAULT_INDEXTYPE lda,                // A
      double anorm,                             // the norm of A
      double *rcond,                            // reciprocal condition number
      double *work, FLENS_DEFAULT_INDEXTYPE *iwork)
{
    FLENS_DEFAULT_INDEXTYPE info;
    dgecon_(&norm, &n, a, &lda, &anorm, rcond, work, iwork, &info);
    return info;
}


} // namespace flens

