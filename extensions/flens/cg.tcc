/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <cassert>
#include <limits>
#include <iomanip>
#include <fstream>
#include <lawa/aux/densify.h>
#include <extensions/flens/sparse_blas_flens.h>

namespace lawa {

template <typename MA, typename VX, typename VB>
FLENS_DEFAULT_INDEXTYPE
cg(const MA &A, VX &x, const VB &b, typename _cg<VB>::T tol,
   FLENS_DEFAULT_INDEXTYPE maxIterations)
{
    typename _cg<VB>::T alpha, beta, rNormSquare, rNormSquarePrev;
    typename _cg<VB>::AuxVector Ap, r, p;

    if (x.length()!=A.numCols()) {
        x.engine().resize(A.numCols());
    }

    r = A*x - b;
    p = -1*r;
    rNormSquare = r*r;
    for (FLENS_DEFAULT_INDEXTYPE k=1; k<=maxIterations; k++) {
        #ifdef SOLVER_DEBUG
            std::cerr << "k = " << k << ", rho = " << sqrt(rNormSquare)
                << std::endl;
        #endif
        if (sqrt(rNormSquare)<=tol) {
            return k-1;
        }
        Ap = A*p;
        alpha = rNormSquare/(p * Ap);
        x += alpha*p;
        r += alpha*Ap;

        rNormSquarePrev = rNormSquare;
        rNormSquare = r*r;
        beta = rNormSquare/rNormSquarePrev;
        p = beta*p - r;
    }

    std::cerr << "Warning! Max iterations reached! rho = " << sqrt(rNormSquare)
              << std::endl;
    #ifdef SOLVER_DEBUG
        std::cerr << "A =\n" << A << std::endl;
    #endif
    return maxIterations;
}

template <typename MA, typename VX, typename VB>
FLENS_DEFAULT_INDEXTYPE
cgls(const MA &A, VX &x, const VB &b, typename _cg<VB>::T tol,
     FLENS_DEFAULT_INDEXTYPE maxIterations)
{
    using namespace flens;
    using namespace flens::blas;
    typename _cg<VB>::T alpha, beta, gammaPrev, gamma, b_norm;
    typename _cg<VB>::AuxVector r, q, s, p;

    assert(b.length()==A.numRows());

    std::cerr << "   cgls: tol = " << tol << std::endl;
    if (x.length()!=A.numCols()) {
        x.engine().resize(A.numCols());
    }
    /*
    for (FLENS_DEFAULT_INDEXTYPE i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        x(i) = 0;
    }
    */

    b_norm = b*b;
    if (std::sqrt(b_norm) < 1e-15) {
        for (FLENS_DEFAULT_INDEXTYPE i=x.firstIndex(); i<=x.lastIndex(); ++i) {
            x(i) = 0;
        }
        return 0;
    }

    r = A*x;
    r += (-1.)*b;
    r *= -1.;
    mv(cxxblas::Trans, typename _cg<VB>::T(1), A, r, typename _cg<VB>::T(0), s);
    p = s;
    gammaPrev = s*s;
#ifdef SOLVER_DEBUG
    std::ofstream gammafile("CGLS_Convergence.txt");
    gammafile << "# Norm(A'*r)^2  Norm(b-Ax)^2" << std::endl; 
#endif
    for (FLENS_DEFAULT_INDEXTYPE k=1; k<=maxIterations; k++) {
        #ifdef SOLVER_DEBUG
            typename _cg<VB>::T res = r*r;
            gammafile << sqrt(gammaPrev) << " " << sqrt(res) << std::endl;
            std::cerr << "k = " << k << ", gamma = " << sqrt(gammaPrev)
            << std::endl;
        #endif
        q = A*p;
        alpha = gammaPrev/(q*q);
        x +=   alpha *p;
        r += (-alpha)*q;
        mv(cxxblas::Trans, typename _cg<VB>::T(1), A, r, typename _cg<VB>::T(0), s);
        gamma = s*s;
        if (sqrt(gamma)<=tol) {
            return k-1;
        }
        beta  = gamma/gammaPrev;
        p *= beta;
        p += s;
        gammaPrev = gamma;
    }
#ifdef SOLVER_DEBUG
    gammafile.close();
#endif
    return maxIterations;
}


// Algorithm 9.2, Y. Saad: Iterative Methods for Sparse Linear Systems
// for solving Ax=b with P^T A P u = P^T b, u=P^{-1} x
// Note the role of P and P^T is switched.
// This algorithm uses split preconditioning defined as in
// (6.2), K. Urban: Wavelet Methods for Elliptic PDEs
// Also non-symmetric preconditioning is possible as mentioned in remark 6.2
template <typename Prec, typename MA, typename VX, typename VB>
FLENS_DEFAULT_INDEXTYPE
pcg(const Prec &P, const MA &A, VX &x, const VB &b,
    typename _cg<VB>::T tol, FLENS_DEFAULT_INDEXTYPE maxIterations)
{
    typename _cg<VB>::T pNormSquare, alpha, beta, rHatq, rHatqPrev;
    typename _cg<VB>::AuxVector r, rHat, p, Ap;

    if (x.length()!=A.numCols()) {
        x.engine().resize(A.numCols());
    }
    r = b - A*x;
    if(r*r == 0){
        return 0;
    }    
    rHat = transpose(P)*r;
    p = P*rHat;
    rHatq = rHat*rHat;
    for (FLENS_DEFAULT_INDEXTYPE k=1; k<=maxIterations; k++) {
        Ap = A*p;
        alpha = rHatq/(Ap*p);
        x += alpha*p;
        rHat = rHat - alpha * transpose(P)*Ap;
        rHatqPrev = rHatq;
        rHatq =rHat*rHat;
        beta = rHatq/rHatqPrev;
        p *= beta;
        p += P*rHat;
        pNormSquare = p*p;
        #ifdef SOLVER_DEBUG
            std::cerr << "k = " << k << ", rho = " <<  sqrt(pNormSquare) 
                << std::endl;
        #endif
        if (sqrt(pNormSquare)<=tol) {
            return k-1;
        }
    }
    return maxIterations;
}

/*
template <typename Prec, typename MA, typename VX, typename VB>
FLENS_DEFAULT_INDEXTYPE
pcg(const Prec &B, const MA &A, VX &x, const VB &b,
    typename _cg<VB>::T tol, FLENS_DEFAULT_INDEXTYPE maxIterations)
{
    typename _cg<VB>::T pNormSquare, alpha, beta, rq, rqPrev;
    typename _cg<VB>::AuxVector r, q, p, Ap;

    r = A*x - b;
    q = B*r;

    p = q;
    // TODO: next line results in an error with T = FLENS_DEFAULT_INDEXTYPE double. WHY???
    // p = -q;
    p *= -1;
    rq = r*q;

    for (FLENS_DEFAULT_INDEXTYPE k=1; k<=maxIterations; k++) {
        pNormSquare = p*p;
        if (sqrt(pNormSquare)<tol) {
            return k-1;
        }
        Ap = A*p;
        alpha = rq/(p*Ap);
        x += alpha*p;

        r += alpha*Ap;
        q = B*r;

        rqPrev = rq;
        rq = r*q;
        beta = rq/rqPrev;
        p = beta*p - q;
    }
    return maxIterations;
}
*/
} // namespace lawa

