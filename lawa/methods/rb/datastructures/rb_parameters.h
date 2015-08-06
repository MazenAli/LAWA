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

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_RB_PARAMETERS_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_RB_PARAMETERS_H_

#include <array>
#include <lawa/settings/enum.h>

namespace lawa {

enum TrainingType {weak, strong, strong_adaptive, weak_direct};
enum SnapshotTolReductionCrit {repeated_param, conv_rate_degradation};

/* Helper struct to get size of an array
 */
template<typename ParamType>
struct ParamInfo {};

template<typename T,std::size_t N>
struct ParamInfo<std::array<T,N> >
{
	typedef T ValueType;
	static const std::size_t dim = N;

	static void
	print(std::array<T,N> param);
};

/* Parameters for an RB Training
 *
 */
template<typename ParamType>
struct RB_Greedy_Parameters{

    VariationalFormulationType  problem_type;
	TrainingType 				training_type;
	SnapshotTolReductionCrit 	snapshot_tol_red_crit;

	typedef typename std::array<std::size_t,ParamInfo<ParamType>::dim> intArray;
	double 		tol;
    double      deter_rate;
	std::size_t 		Nmax;
	ParamType 	min_param;
	ParamType 	max_param;
	intArray	nb_training_params;
	intArray	log_scaling;

	bool 		print_info;
	std::string print_file;
	bool 		verbose;
	bool 		write_during_training;
	std::string trainingdata_folder;
	bool 		print_paramset;
	bool 		erase_snapshot_params;
	bool		orthonormalize_bfs;
	bool		tighten_tol;
	bool		tighten_tol_rieszA;
	bool		tighten_tol_rieszF;
	double		tighten_tol_reduction;
	bool		update_snapshot;			// Do not recompute snapshots at same param from scratch
	bool		update_rieszF;
	bool		update_rieszA;
	bool 		coarsen_rieszA_for_update;
	bool		test_estimator_equivalence;
	bool		tighten_estimator_accuracy;
	double 		riesz_constant_X;
	double 		riesz_constant_Y;
	bool 		write_direct_representors;
	double 		min_error_reduction;		// if snapshot_tol_red_crit = conv_rate_degradation
	double		refSolution_tol_factor; 	// factor between accuracy of "exact" sol and that of "truth"
    bool        read_truth_sols;            // for strong Greedy
    std::size_t nb_existing_truth_sols;     // for strong Greedy
	
	
	RB_Greedy_Parameters(
                        VariationalFormulationType  _problem_type = Galerkin,
						TrainingType _training_type = weak,
						double _tol = 1e-2,
                        double _deter_rate = 0.5,
						std::size_t _Nmax = 20,
						ParamType _min_param = ParamType(),
						ParamType _max_param = ParamType(),
						intArray  _training_params_per_dim = intArray(),
						intArray _log_scaling = intArray(),
						bool _print_info = true,
						std::string _print_file = "greedy_info.txt",
						bool _verbose = true,
						bool _write_during_training = true,
						std::string _trainingdata_folder = "training_data",
						bool _print_paramset = false,
						bool _erase_snapshot_params = false,
						bool _orthonormalize_bfs = true,
						bool _tighten_tol	= false,
						SnapshotTolReductionCrit _snapshot_tol_red_crit = repeated_param,
						bool _tighten_tol_rieszA = false,
						bool _tighten_tol_rieszF = false,
						double _tighten_tol_reduction = 0.1,
						bool _update_snapshot = false,
						bool _update_rieszF = false,
						bool _update_rieszA = false,
						bool _coarsen_rieszA_for_update = false,
						bool _test_estimator_equivalence = false,
						bool _tighten_estimator_accuracy = false,
						double _riesz_constant_X = 1.,
						double _riesz_constant_Y = 1.,
						bool _write_direct_representors = false,
						double _min_error_reduction = 0.5,
						double	_refSolution_tol_factor = 0.1,
						bool _read_truth_sols = false,
						std::size_t _nb_existing_truth_sols = 0);

	void print();
};

template<typename ParamType>
struct RB_Greedy_Information{
	std::vector<double> 	 			   		greedy_errors;
	std::vector<ParamType>	 			   		snapshot_params;
	std::vector<std::size_t>		 	   		u_size;					// Dim 1 x n_bf
	std::vector<std::vector<std::size_t> >  	repr_f_size;		    // Dim 1 x Q_f in each iteration (if refined)
	std::vector<std::vector<std::size_t> > 		repr_a_size;			// Dim n_bf x Q_a
	std::vector<std::vector<std::size_t> >		repr_r_size;			// Dim n_bf x n_train (only for weak direct Greedy!)

	std::vector<std::vector<double> >			eps_res_bound;
	std::vector<std::vector<double> >			eps_aff;

	std::vector<std::vector<double> >					accuracies_f_reprs;			// Dim N x Q_f
	std::vector<std::vector<std::vector<double> > >		accuracies_a_reprs;			// Dim N x n_bf x Q_a

	void print(const char* filename = "greedy_info.txt");

	void read(const char* filename, std::size_t Qf, std::size_t Qa, int nb = -1);
	
    void read_repr_accuracies(const char* filename, std::size_t Qf, std::size_t Qa, int Nmax);

	void print_bound_info(const char* filename = "greedy_info_eps_errest");
};

/* Parameters for a RB solution
 */
template<typename ParamType>
struct RB_Parameters {

	SolverCall  call;
	ParamType 	ref_param;

	bool verbose;

	RB_Parameters(SolverCall _call = call_cg,
				  ParamType _ref_param = ParamType(),
				  bool _verbose = true);

	void print();
};


} // namespace lawa

#include <lawa/methods/rb/datastructures/rb_parameters.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_RB_PARAMETERS_H_ */
