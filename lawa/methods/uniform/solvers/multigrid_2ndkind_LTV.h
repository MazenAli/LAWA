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
 
#ifndef LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRID_2NDKIND_LTV_H
#define LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRID_2NDKIND_LTV_H 1

#include <lawa/methods/uniform/datastructures/datastructures.h>
#include <lawa/methods/uniform/solvers/thetascheme1d_LTV.h>
#include <lawa/righthandsides/homogeneousrhs.h>

namespace lawa{
    
/* MultiGrid of the 2nd kind
 *      This class solves periodic problems using the MultiGrid Algorithm of the second kind
 *      (see W. Hackbusch "Fast Numerical Solution Of Time-Periodic Parabolic Problems by a 
 *       Multigrid Method", SIAM J. Sci. Stat. Comp. 2 (2), pp. 198-206 (1981))
 */    
template<typename T, typename PrimalBasis, typename DualBasis, typename BilinearForm, typename RHSIntegral>
class MultiGrid_2ndKind_LTV{
    
    public: 
        
        MultiGrid_2ndKind_LTV(PrimalBasis& _b, DualBasis& _b_, BilinearForm& a, RHSIntegral& rhs, 
                          T theta, T deltaT, int timesteps, int minLevel);
        
        flens::DenseVector<flens::Array<T> > 
        solve(flens::DenseVector<flens::Array<T> >& u0, int maxLevel);
    
    private:
        
        PrimalBasis& b;
        DualBasis& b_;
        typedef ThetaScheme1D_LTV<T, PrimalBasis, BilinearForm, RHSIntegral>        FullThetaScheme;
        typedef TimeStepping<T, FullThetaScheme>                                    FullTimeStepMethod;    
        typedef ThetaScheme1D_LTV<T, PrimalBasis, BilinearForm, HomogeneousRHS<T> > HomThetaScheme;
        typedef TimeStepping<T, HomThetaScheme>                                     HomTimeStepMethod;
        typedef FixedPointSolver<T, HomTimeStepMethod>                              ThetaFPSolver;
        
        FullThetaScheme         full_theta;
        FullTimeStepMethod      full_ts;
        HomogeneousRHS<T>       hom_rhs;
        HomThetaScheme          hom_theta;
        HomTimeStepMethod       hom_ts;
        ThetaFPSolver           fp;


        class MG_2ndKind_Smoother{
            private:
                MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>* mg_ptr;
        
            public: 
                MG_2ndKind_Smoother(MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>* ref);
                
                flens::DenseVector<flens::Array<T> > 
                solve(flens::DenseVector<flens::Array<T> > u, flens::DenseVector<flens::Array<T> > f);  
                
                void setLevel(int level){ mg_ptr->hom_ts.setLevel(level);}
                
                flens::DenseVector<flens::Array<T> >
                getResiduum(flens::DenseVector<flens::Array<T> >& u,flens::DenseVector<flens::Array<T> >& f)
                { return mg_ptr->hom_ts.getResiduum(u) - f;}
        };
    
        class MG_2ndKind_Solver{
            private:
                MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>* mg_ptr;
        
            public: 
                MG_2ndKind_Solver(MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>* ref);
                
                flens::DenseVector<flens::Array<T> > 
                solve(flens::DenseVector<flens::Array<T> > u, flens::DenseVector<flens::Array<T> > f);  
                
                void setLevel(int level){ mg_ptr->fp.setLevel(level);}
        };
        
        friend class MG_2ndKind_Smoother;
        friend class MG_2ndKind_Solver;
        
        typedef MultigridSolver<T, PrimalBasis, DualBasis, MG_2ndKind_Smoother, MG_2ndKind_Solver> Multigrid;
        MG_2ndKind_Smoother mg_smoother;
        MG_2ndKind_Solver mg_solver;
        Multigrid mg;

};



} // namespace lawa

#include <lawa/methods/uniform/solvers/multigrid_2ndkind_LTV.tcc>

#endif // LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRID_2NDKIND_LTV_H

