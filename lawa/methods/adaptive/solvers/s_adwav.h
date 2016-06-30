/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2014  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

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
#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_H 1

#include <iostream>
#include <vector>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/matrixoperations.h>
#include <lawa/methods/adaptive/algorithms/linearsystemsolvers.h>
#include <lawa/methods/adaptive/postprocessing/postprocessing.h>

namespace lawa {

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
class S_ADWAV {
    public:
        S_ADWAV(const Basis &basis, MA &A, RHS &F, T contraction, T start_threshTol,
                T _linTol=1e-6, T _resTol=1e-4, FLENS_DEFAULT_INDEXTYPE _NumOfIterations=10, FLENS_DEFAULT_INDEXTYPE _MaxItsPerThreshTol=5,
                T _eps=1e-2, FLENS_DEFAULT_INDEXTYPE _MaxSizeLambda = 400, T _resStopTol=0.1,
                std::vector<FLENS_DEFAULT_INDEXTYPE> _Jmaxvec = std::vector<FLENS_DEFAULT_INDEXTYPE>(0));

        //solver for symmetric elliptic problems
        void solve(const IndexSet<Index> &Initial_Lambda, const char *linsolvertype,
                   const char *filename, FLENS_DEFAULT_INDEXTYPE assemble_matrix=2, T H1norm=0.);
        //solver for symmetric elliptic problems
        void solve_cg(const IndexSet<Index> &Initial_Lambda, FLENS_DEFAULT_INDEXTYPE assemble_matrix=1, T H1norm=0.);
        //solver for symmetric elliptic problems without B-Splines
        void solve_cg_WO_XBSpline(const IndexSet<Index> &Initial_Lambda, FLENS_DEFAULT_INDEXTYPE assemble_matrix=1, T H1norm=0.);
        //solver for elliptic problems
        void solve_gmres(const IndexSet<Index> &Initial_Lambda, FLENS_DEFAULT_INDEXTYPE assemble_matrix=1);
        void solve_gmresm(const IndexSet<Index> &Initial_Lambda, FLENS_DEFAULT_INDEXTYPE assemble_matrix=1);
        //solver for indefinite problems
        void solve_cgls(const IndexSet<Index> &Initial_Lambda, FLENS_DEFAULT_INDEXTYPE assemble_matrix=1);
        
        void
        set_parameters(T _contraction, T _threshTol, T _linTol=1e-6, T _resTol=1e-4, 
                       FLENS_DEFAULT_INDEXTYPE _NumOfIterations=10, FLENS_DEFAULT_INDEXTYPE _MaxItsPerThreshTol=5, T _eps=1e-2, 
                       FLENS_DEFAULT_INDEXTYPE _MaxSizeLambda = 400, T _resStopTol=0.1, 
                       std::vector<FLENS_DEFAULT_INDEXTYPE> _Jmaxvec = std::vector<FLENS_DEFAULT_INDEXTYPE>(0));
        void
        get_parameters(T& _contraction, T& _threshTol, T& _linTol, T& _resTol, 
                       FLENS_DEFAULT_INDEXTYPE& _NumOfIterations, FLENS_DEFAULT_INDEXTYPE& _MaxItsPerThreshTol, T& _eps, 
                       FLENS_DEFAULT_INDEXTYPE& _MaxSizeLambda, T& _resStopTol,
                       std::vector<FLENS_DEFAULT_INDEXTYPE>& _Jmaxvec);
    
        std::vector<Coefficients<Lexicographical,T,Index> > solutions;
        std::vector<T>               residuals;
        std::vector<T>               times;
        std::vector<T>               linsolve_iterations;
        std::vector<T>               toliters;

        const Basis &basis;
    
        bool relative_thresh;

    private:
        MA &A;
        RHS &F;
        T contraction, threshTol, linTol, resTol;
        FLENS_DEFAULT_INDEXTYPE NumOfIterations; 
        FLENS_DEFAULT_INDEXTYPE MaxItsPerThreshTol;
        T eps;
        FLENS_DEFAULT_INDEXTYPE MaxSizeLambda;
        T resStopTol;
        std::vector<FLENS_DEFAULT_INDEXTYPE> Jmaxvec;

};

} //namespace lawa

#include <lawa/methods/adaptive/solvers/s_adwav.tcc>

#endif //LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_H

