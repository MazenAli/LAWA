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

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_RB_BASE_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_RB_BASE_H_

#include <vector>
#include <lawa/methods/rb/datastructures/rb_system.h>
#include <lawa/methods/rb/datastructures/thetastructure.h>
#include <lawa/methods/rb/datastructures/rb_parameters.h>

namespace lawa {

/* RB_Base:
 * 		General offline computations for the construction of
 * 		a reduced basis model.
 */
template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
class RB_Base{

	typedef typename DataType::ValueType 						T;
	typedef typename RB_Greedy_Parameters<ParamType>::intArray 	intArrayType;

public:
	RB_Base(RB_Model& _rb_system, TruthModel& _rb_truth);

	std::size_t
	n_bf();

	void
	train_Greedy(std::size_t N = 0);
	
	void
	train_Greedy(std::vector<ParamType>& Xi_train, std::size_t N = 0);

	/*
	void
	train_strong_Greedy();
	*/

    DataType
    reconstruct_u_N(typename RB_Model::DenseVectorT u, std::size_t N);
    
    DataType
    reconstruct_u_N(typename RB_Model::DenseVectorT u, std::vector<std::size_t> indices);
    
    DataType
    reconstruct_res_repr_N(typename RB_Model::DenseVectorT u, std::size_t N, ParamType& mu);

	void
	write_basisfunctions(const std::string& directory_name = "offline_data/bf", int nb = -1);

    void
	write_snapshot(DataType& u, const std::string& directory_name = "offline_data/snap",  int nb = -1);

	void
	read_basisfunctions(const std::string& directory_name = "offline_data/bf", int nb = -1);

	void
	write_rieszrepresentors(const std::string& directory_name = "offline_data/representors", int nb = -1);

	void
	read_rieszrepresentors(const std::string& directory_name = "offline_data/representors", int nb = -1);

	DataType&
	get_bf(std::size_t i);

	void
	read_greedy_info(const char* filename, int nb = -1);
	
	void
    read_repr_accuracies(const char* filename, int Nmax);

	virtual T
	get_direct_errorbound(const typename RB_Model::DenseVectorT& u_N, ParamType& mu, DataType& res_repr);

	virtual T
	update_direct_errorbound(const typename RB_Model::DenseVectorT& u_N, ParamType& mu, DataType& res_repr);

	std::vector<ParamType>
	generate_uniform_paramset(ParamType min_param, ParamType max_param, intArrayType param_nb, intArrayType log_scaling);

	void
	print_paramset(std::vector<ParamType> PSet);

	void
	remove_basisfunction(std::size_t nb, bool from_greedy_info);

	RB_Greedy_Parameters<ParamType> 	greedy_params;

	std::vector<DataType> 				rb_basisfunctions;

public:

	void
	add_to_basis(const DataType& u);

	void
	add_to_RB_structures(const DataType& bf);

    virtual
	void
	calculate_Riesz_RHS_information(bool update = false);

    virtual
    void
	calculate_Riesz_LHS_information(const DataType& bf, bool update = false	);

	void
	recalculate_A_F_norms();

	T
	find_max_errest(std::size_t N, std::vector<ParamType>& Xi_train, ParamType& current_param, std::map<ParamType, DataType>& truth_sols);

	RB_Model& 		rb_system;
	TruthModel& 	rb_truth;

	std::vector<DataType>				F_representors;  // Dim: 1 x Q_f
	std::vector<std::vector<DataType> > A_representors;  // Dim: n x Q_a

	std::vector<T>						eps_F_representors;  // Dim: 1 x Q_f
	std::vector<std::vector<T> > 		eps_A_representors;  // Dim: n x Q_a

	RB_Greedy_Information<ParamType>	greedy_info;

};

} // namespace lawa

#include <lawa/methods/rb/datastructures/rb_base.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_RB_BASE_H_ */
