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


#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_WEIGHTEDADAPTIVEIDENTITYOPERATOR1D_PG_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_WEIGHTEDADAPTIVEIDENTITYOPERATOR1D_PG_H 1

#include <lawa/settings/enum.h>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/operators/pdeoperators1d/weightedpdeoperator1d_pg.h>
#include <lawa/preconditioners/preconditioners.h>
#include <lawa/methods/adaptive/datastructures/mapmatrix.h>

#ifdef TRONE
    #include <tr1/unordered_map>
#elif BOOST
    #include <boost/unordered_map.hpp>
#elif CONEONE
	#include <unordered_map>
#else
    #include <ext/hash_map>
#endif

namespace lawa {

/* Adaptive PDE VaryingCoefficients OPERATOR: Petrov Galerkin version
 *
 *    a(v,u) =  diffusion(x) * Integral(v_x * u_x) +  convection_f(x) * Integral(v * u_x)
 *              + reaction_f(x) * Integral(v * u)
 *
 */
template <typename T, typename TrialBasis, typename TestBasis>
struct AdaptiveWeightedPDEOperator1D_PG
{
    typedef flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> >        SparseMatrixT;
    typedef NoCompression<T,Index1D,TrialBasis>                             NoCompression1D;
    typedef NoPreconditioner<T,Index1D>                                     NoPreconditioner1D;
    typedef WeightedPDEOperator1D_PG<T,TrialBasis,TestBasis>                WeightedPDEOp1D;
//    typedef MapMatrix<T,Index1D,WeightedPDEOp1D,
//                      NoCompression1D,NoPreconditioner1D>                   DataWeightedPDEOp1D;
#ifdef TRONE
    typedef typename std::tr1::unordered_map<Entry<Index1D>, T, entry_hashfunction<Index1D>,
                                                             entry_eqfunction<Index1D> > DataWeightedPDEOp1D;
#elif BOOST
    typedef typename boost::unordered_map<Entry<Index1D>, T, entry_hashfunction<Index1D>,
                                                                 entry_eqfunction<Index1D> > DataWeightedPDEOp1D;
#elif CONEONE
    typedef typename std::unordered_map<Entry<Index1D>, T, entry_hashfunction<Index1D>,
                                                                 entry_eqfunction<Index1D> > DataWeightedPDEOp1D;
#else
    typedef typename __gnu_cxx::hash_map<Entry<Index1D>, T, entry_hashfunction<Index1D>,
                                                             entry_eqfunction<Index1D> > DataWeightedPDEOp1D;
#endif


    AdaptiveWeightedPDEOperator1D_PG(const TrialBasis& _trialbasis1d, const TestBasis& _testbasis1d, 
								  Function<T> &_reaction_f, Function<T> &_convection_f, Function<T>& _diffusion_f,
                                  int order=10,
                                  bool reactionIsZero=false, bool convectionIsZero=false,
                                  bool diffusionIsZero=false);

    T
    operator()(const Index1D &row_index, const Index1D &col_index);

    T
    operator()(XType xtype_row, int j_row, long k_row, XType xtype_col, int j_col, long k_col);

    void
    clear();

    const TrialBasis&           trialbasis1d;
    const TestBasis&            testbasis1d;
    NoCompression1D             compression1d;
    const WeightedPDEOp1D       weightedpdeop1d;
    NoPreconditioner1D          prec1d;
    DataWeightedPDEOp1D         data;

};

}   //namespace lawa

#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptiveweightedpdeoperator1d_pg.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_WEIGHTEDADAPTIVEIDENTITYOPERATOR1D_PG_H

