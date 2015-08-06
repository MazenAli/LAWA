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

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_RB_SYSTEM_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_RB_SYSTEM_H_

#include <flens/flens.cxx>
#include <array>
#include <vector>
#include <lawa/methods/rb/datastructures/thetastructure.h>
#include <lawa/methods/rb/datastructures/rb_parameters.h>


namespace lawa {

/* RB_System:
 * 		This class contains all N-dependent data and functions needed
 * 		for the computation of a RB approximation.
 */
template <typename T, typename ParamType>
class RB_System {

public:
    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
    typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

	RB_System(ThetaStructure<ParamType>& _thetas_a,
			  ThetaStructure<ParamType>& _thetas_f);

	DenseVectorT
	get_rb_solution(std::size_t N, ParamType& mu);

	DenseVectorT
	get_rb_solution(std::vector<std::size_t> indices, ParamType& mu);
	
	void
    get_rb_LGS(std::size_t N, ParamType& mu, FullColMatrixT& A, DenseVectorT& F);

	virtual T
	get_errorbound(const DenseVectorT& u_N, ParamType& mu);
	
	virtual T
	get_errorbound(std::vector<std::size_t> indices, const DenseVectorT& u_N, ParamType& mu);

	T
	get_errorbound_accuracy(const DenseVectorT& u_N, ParamType& mu, std::vector<T> eps_f, std::vector<std::vector<T> > eps_a);

	T
	residual_dual_norm(const DenseVectorT& u_N, ParamType& mu);

	T
	residual_dual_norm(std::vector<std::size_t> indices, const DenseVectorT& u_N, ParamType& mu);

	virtual T
	alpha_LB(ParamType& mu);

	std::size_t
	Q_f();

	std::size_t
	Q_a();

	void
	write_rb_data(const std::string& directory_name = "offline_data");

	void
	read_rb_data(const std::string& directory_name = "offline_data", int nb = -1);

	void
	remove_basisfunction(std::size_t nb);

	RB_Parameters<ParamType> 					rb_params;

    ThetaStructure<ParamType>& 					thetas_a;
    ThetaStructure<ParamType>& 					thetas_f;

    std::vector<FullColMatrixT>     			RB_A_matrices;
    std::vector<DenseVectorT>       			RB_F_vectors;
    FullColMatrixT  							RB_inner_product;

    FullColMatrixT                              F_F_representor_norms; // Speicherbedarf kann verringert werden..
    std::vector<FullColMatrixT>                 A_F_representor_norms;
    std::vector<std::vector<FullColMatrixT> >   A_A_representor_norms; //.. Ausnutzen der Symmetrie (Matrix als Vektor)

};


} // namespace lawa

#include <lawa/methods/rb/datastructures/rb_system.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_RB_SYSTEM_H_ */
