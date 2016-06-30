/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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
 #include <iomanip>
 #include <lawa/righthandsides/righthandsides.h>
 #include <lawa/operators/operators.h>
 #include <lawa/functiontypes/separablefunction2d.h>
 #include <lawa/methods/adaptive/datastructures/coefficients.h>
 #include <lawa/methods/adaptive/datastructures/index.h>
 #include <lawa/methods/adaptive/righthandsides/rhs.h>
 #include <lawa/settings/enum.h>
 #include <cmath>
 #include <flens/flens.cxx>

 
namespace lawa {

template <typename T, typename Index, typename MA, typename RHS>
T
estimateError_Au_M_f(MA &A, RHS &F, const Coefficients<Lexicographical,T,Index> & u,
                     const IndexSet<Index> &LambdaRow)
{
    Coefficients<Lexicographical,T,Index> Au(u.d,u.d_), f(u.d,u.d_), res(u.d,u.d_);
    Au = mv_sparse(LambdaRow,A,u);
    f  = F(LambdaRow);
    res = Au-f;
    return res.norm(2.)/f.norm(2.);
}

template <typename T, typename Index, typename MA, typename RHS>
T
computeErrorInH1Norm(MA &A_H, RHS &F_H, const Coefficients<Lexicographical,T,Index> & u,
                     T HNormOfExactSolution, bool optimized)
{
    Coefficients<Lexicographical,T,Index> Au, f_H;
    IndexSet<Index> LambdaRow;
    LambdaRow = supp(u);
    if (optimized) {
        Au = A_H.mv(LambdaRow,u);
    }
    else {
        Au = mv_sparse(LambdaRow,A_H,u);
    }
    //Au = mv(supp(u),A_H,u);
    //std::cerr << "Attention: Slow matrix vector product..." << std::endl;
    T uAu = u*Au;
    f_H   = F_H(supp(u));
    T fu  = f_H*u;
    std::cerr << "   Estimated Energy norm squared: " << uAu << std::endl;

    return std::sqrt(fabs(std::pow(HNormOfExactSolution,(T)2.)- 2*fu + uAu));

}

template <typename T, typename Index, typename SOLVER, typename MA_H, typename RHS_H>
void
postprocessing_H1(SOLVER& Solver, MA_H &A_H1, RHS_H &F_H1, T H1norm, const char* filename)
{

    Coefficients<Lexicographical,T,Index> u;
    std::ofstream file(filename);

    for (FLENS_DEFAULT_INDEXTYPE i=0; i<(FLENS_DEFAULT_INDEXTYPE)(Solver.solutions.size()); ++i) {
        u = Solver.solutions[i];
        T ErrorH1Norm = computeErrorInH1Norm(A_H1, F_H1, u, H1norm);
        file      << supp(u).size() << " " << Solver.linsolve_iterations[i] << " "
                  << Solver.times[i] << " " << Solver.residuals[i] << " "
                  << ErrorH1Norm << std::endl;
        std::cerr << supp(u).size() << " " << Solver.linsolve_iterations[i] << " "
                  << Solver.times[i] << " " << Solver.residuals[i] << " "
                  << ErrorH1Norm << std::endl;
    }
    file.close();
}

template<typename T, typename Preconditioner>
T
estimate_SpaceTimeError_L0T_L2(Coefficients<Lexicographical,T,Index2D> & u, 
                               Coefficients<Lexicographical,T,Index2D> & u_exact,
                               const Preconditioner &P)
{
    T error_est = 0.0;
    IndexSet<Index2D> Lambda = supp(u);
    IndexSet<Index2D> ExpandedLambda = supp(u_exact);
    
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
    
    for (const_set_it it=ExpandedLambda.begin(); it!=ExpandedLambda.end(); ++it) {
        T prec = P((*it));
        if( Lambda.count(*it) > 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += prec * (u[*it] - u_exact[*it]) * prec * (u[*it] - u_exact[*it]);
        }
        else{
            //std::cout << (*it) << " u_exact = " << u_exact[*it] << std::endl;
            error_est += prec * u_exact[*it] * prec * u_exact[*it];
        }
    }
    for (const_set_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        T prec = P((*it));
        if( ExpandedLambda.count(*it) == 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += prec * u[*it] * prec * u[*it];
        }
    }
    return std::sqrt(error_est);
}

template<typename T, typename Preconditioner>
T
estimate_SpaceTimeError_L0T_H1(Coefficients<Lexicographical,T,Index2D> & u, 
                               Coefficients<Lexicographical,T,Index2D> & u_exact,
                               const Preconditioner &P)
{
    T error_est = 0.0;
    IndexSet<Index2D> Lambda = supp(u);
    IndexSet<Index2D> ExpandedLambda = supp(u_exact);
    
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
    
    for (const_set_it it=ExpandedLambda.begin(); it!=ExpandedLambda.end(); ++it) {
        T prec = P((*it));
        if( Lambda.count(*it) > 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += pow2i<T>(2*(*it).index2.j) * prec * (u[*it] - u_exact[*it]) * prec * (u[*it] - u_exact[*it]);
        }
        else{
            //std::cout << (*it) << " u_exact = " << u_exact[*it] << std::endl;
            error_est += pow2i<T>(2*(*it).index2.j) * prec * u_exact[*it] * prec * u_exact[*it];
        }
    }
    for (const_set_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        T prec = P((*it));
        if( ExpandedLambda.count(*it) == 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += pow2i<T>(2*(*it).index2.j) * prec * u[*it] * prec * u[*it];
        }
    }
    return std::sqrt(error_est);
}

template<typename T, typename Preconditioner>
T
estimate_SpaceTimeError_W0T(Coefficients<Lexicographical,T,Index2D> & u, 
                            Coefficients<Lexicographical,T,Index2D> & u_exact,
                            const Preconditioner &P)
{
    T error_est = 0.0;
    IndexSet<Index2D> Lambda = supp(u);
    IndexSet<Index2D> ExpandedLambda = supp(u_exact);
    
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
    
    for (const_set_it it=ExpandedLambda.begin(); it!=ExpandedLambda.end(); ++it) {
        T prec = P((*it));
        if( Lambda.count(*it) > 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += (pow2i<T>(2*(*it).index1.j - 2*(*it).index2.j) + pow2i<T>(2*(*it).index2.j)) 
                            * prec * (u[*it] - u_exact[*it]) * prec * (u[*it] - u_exact[*it]);
        }
        else{
            //std::cout << (*it) << " u_exact = " << u_exact[*it] << std::endl;
            error_est += (pow2i<T>(2*(*it).index1.j - 2*(*it).index2.j) + pow2i<T>(2*(*it).index2.j)) 
                            * prec * u_exact[*it] * prec * u_exact[*it];
        }
    }
    for (const_set_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        T prec = P((*it));
        if( ExpandedLambda.count(*it) == 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += pow2i<T>(2*(*it).index2.j) * prec * u[*it] * prec * u[*it];
        }
    }
    
    return std::sqrt(error_est);
}


template <typename T, typename Basis2D>
T
h1Error2D(Basis2D& basis2d,  SeparableFunction2D<T>&  u,
                             SeparableFunction2D<T>&  ux,
                             SeparableFunction2D<T>&  uy,
          Coefficients<Lexicographical, T, Index2D>&  uh,
          T tolu,
          T tolux,
          T toluy)
{
    typedef     SeparableRHS2D<T, Basis2D>                   INTEGRAL;
    typedef     NoPreconditioner<T, Index2D>                 NoPrec;
    typedef     RHS<T, Index2D, INTEGRAL, NoPrec>            _RHS;
    typedef     Integral<Gauss,
                typename Basis2D::FirstBasisType,
                typename Basis2D::FirstBasisType>            WAVX;
    typedef     Integral<Gauss,
                typename Basis2D::SecondBasisType,
                typename Basis2D::SecondBasisType>           WAVY;
    typedef     Coefficients<Lexicographical, T, Index2D>    Coefficients;

    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor>>   nodeltas;

    INTEGRAL      uint   (basis2d, u,  nodeltas, nodeltas, 100);
    INTEGRAL      uxint  (basis2d, ux, nodeltas, nodeltas, 100);
    INTEGRAL      uyint  (basis2d, uy, nodeltas, nodeltas, 100);
    INTEGRAL      uxintdp(basis2d, ux, nodeltas, nodeltas, 100, 1, 0);
    INTEGRAL      uyintdp(basis2d, uy, nodeltas, nodeltas, 100, 0, 1);

    NoPrec        noprec;

    _RHS          uref   (uint,    noprec);
    _RHS          uxref  (uxint,   noprec);
    _RHS          uyref  (uyint,   noprec);
    _RHS          uxdpref(uxintdp, noprec);
    _RHS          uydpref(uyintdp, noprec);

    WAVX          psix(basis2d.first,  basis2d.first);
    WAVY          psiy(basis2d.second, basis2d.second);

    IndexSet<Index2D>       Lambda;
    T gamma = 0.2;
    getSparseGridIndexSet(basis2d, Lambda, 1, 0, gamma);

    Coefficients    errorvec;
    T               error = 0.;

    // L2 error
    Coefficients    Uref;
    sample_f(basis2d, Lambda, uref, Uref, tolu, true);
    errorvec = Uref - uh;
    error   += std::pow(errorvec.norm(2.), 2.);
    Uref.clear();
    errorvec.clear();

    // H1 seminorm error
    // x part
    sample_f(basis2d, Lambda, uxref, Uref, tolux, true);
    error += std::pow(Uref.norm(2.), 2.);
    Uref.clear();

    T ux_uhx = 0.;
    for (auto& lambda : uh) {
        ux_uhx += lambda.second*uxdpref(lambda.first);
    }
    error -= 2.*ux_uhx;

    T uhx_uhx = 0.;
    for (auto& lambda1 : uh) {
        for (auto& lambda2 : uh) {
            if (lambda2.first.index2.xtype == lambda1.first.index2.xtype &&
                lambda2.first.index2.j     == lambda1.first.index2.j     &&
                lambda2.first.index2.k     == lambda1.first.index2.k) {

                uhx_uhx += lambda1.second*lambda2.second*
                           psix(lambda1.first.index1.j, lambda1.first.index1.k,
                                                 lambda1.first.index1.xtype, 1,
                               lambda2.first.index1.j, lambda2.first.index1.k,
                                                 lambda2.first.index1.xtype, 1);
            }
        }
    }
    error += uhx_uhx;

    // y part
    sample_f(basis2d, Lambda, uyref, Uref, toluy, true);
    error += std::pow(Uref.norm(2.), 2.);
    Uref.clear();

    T uy_uhy = 0.;
    for (auto& lambda : uh) {
        uy_uhy += lambda.second*uydpref(lambda.first);
    }
    error -= 2.*uy_uhy;

    T uhy_uhy = 0.;
    for (auto& lambda1 : uh) {
        for (auto& lambda2 : uh) {
            if (lambda2.first.index1.xtype == lambda1.first.index1.xtype &&
                lambda2.first.index1.j     == lambda1.first.index1.j     &&
                lambda2.first.index1.k     == lambda1.first.index1.k) {

                uhy_uhy += lambda1.second*lambda2.second*
                           psiy(lambda1.first.index2.j, lambda1.first.index2.k,
                                                 lambda1.first.index2.xtype, 1,
                               lambda2.first.index2.j, lambda2.first.index2.k,
                                                 lambda2.first.index2.xtype, 1);
            }
        }
    }
    error += uhy_uhy;

    return std::sqrt(error);
}

} // namespace lawa

