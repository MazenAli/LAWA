#include <iostream>
#include <fstream>
#include <array>

namespace lawa {

template<typename T,std::size_t N>
void
ParamInfo<std::array<T,N> >::print(std::array<T,N> mu)
{
	std::cout << "[";
	for(auto& el : mu){
		std::cout << " " << std::fixed << std::right << std::setw(18) << el;
	}
	std::cout << "] ";
}

template<typename ParamType>
RB_Greedy_Parameters<ParamType>::RB_Greedy_Parameters(
        VariationalFormulationType _problem_type,
		TrainingType _training_type, double _tol, std::size_t _Nmax,
		ParamType _min_param, ParamType _max_param,
		intArray _training_params_per_dim,
		intArray _log_scaling, bool _print_info,
		std::string _print_file, bool _verbose,
		bool _write_during_training, std::string _trainingdata_folder,
		bool _print_paramset, bool _erase_snapshot_params,
		bool _orthonormalize_bfs, bool _tighten_tol,
		SnapshotTolReductionCrit _snapshot_tol_red_crit,
		bool _tighten_tol_rieszA, bool _tighten_tol_rieszF, double _tighten_tol_reduction,
		bool _update_snapshot, bool _update_rieszF, bool _update_rieszA, bool _coarsen_rieszA_for_update,
		bool _test_estimator_equivalence, bool _tighten_estimator_accuracy,
		double _riesz_constant_X, double _riesz_constant_Y,
		bool _write_direct_representors, double _min_error_reduction,
		double	_refSolution_tol_factor,
		bool _read_truth_sols, std::size_t _nb_existing_truth_sols)
 : problem_type(_problem_type),
   training_type(_training_type), snapshot_tol_red_crit(_snapshot_tol_red_crit),
   tol(_tol), Nmax(_Nmax), min_param(_min_param), max_param(_max_param),
   nb_training_params(_training_params_per_dim), log_scaling(_log_scaling),
   print_info(_print_info), print_file(_print_file), verbose(_verbose),
   write_during_training(_write_during_training), trainingdata_folder(_trainingdata_folder),
   print_paramset(_print_paramset),
   erase_snapshot_params(_erase_snapshot_params),
   orthonormalize_bfs(_orthonormalize_bfs),
   tighten_tol(_tighten_tol),
   tighten_tol_rieszA(_tighten_tol_rieszA),
   tighten_tol_rieszF(_tighten_tol_rieszF),
   tighten_tol_reduction(_tighten_tol_reduction),
   update_snapshot(_update_snapshot),
   update_rieszF(_update_rieszF),
   update_rieszA(_update_rieszA),
   coarsen_rieszA_for_update(_coarsen_rieszA_for_update),
   test_estimator_equivalence(_test_estimator_equivalence),
   tighten_estimator_accuracy(_tighten_estimator_accuracy),
   riesz_constant_X(_riesz_constant_X),
   riesz_constant_Y(_riesz_constant_Y),
   write_direct_representors(_write_direct_representors),
   min_error_reduction(_min_error_reduction),
   refSolution_tol_factor(_refSolution_tol_factor),
   read_truth_sols(_read_truth_sols),
   nb_existing_truth_sols(_nb_existing_truth_sols)
{}

template<typename ParamType>
void
RB_Greedy_Parameters<ParamType>::print()
{
	std::cout << "###### RB Training Parameters ######################" << std::endl;
	std::cout << std::left << std::setw(24) << "# Training Type:" << std::setw(20) <<  ((training_type==strong || training_type==strong_adaptive)?"strong Greedy":"weak Greedy")
			  << (training_type==weak_direct?" (direct)":"") << (training_type==strong_adaptive?" (w.r.t. reference solution)":"") << std::endl;
	std::cout << std::left << std::setw(24) << "# tol:" << std::setw(20) <<  tol << std::endl;
	std::cout << std::left << std::setw(24) << "# Nmax:" << std::setw(20) <<  Nmax << std::endl;
	std::cout << std::left << std::setw(24) << "# Min_Param:" << std::setw(2) << "[ ";
		for(auto& el : min_param){
			std::cout << std::right << std::setw(4) << el << " ";
		}
		std::cout << std::left << " ]" << std::endl;
	std::cout << std::left << std::setw(24) << "# Max_Param:" << std::setw(2) << "[ ";
		for(auto& el : max_param){
			std::cout << std::right << std::setw(4) << el << " ";
		}
		std::cout << std::left << " ]" << std::endl;
	std::cout << std::left << std::setw(24) << "# Nb of train. params:" << std::setw(2) << "[ ";
		for(auto& el : nb_training_params){
			std::cout << std::right << std::setw(4) << el << " ";
		}
		std::cout << std::left << " ]" << std::endl;
	std::cout << std::left << std::setw(24) << "# Log Scaling:" << std::setw(2) << "[ ";
		for(auto& el : log_scaling){
			std::cout << std::right << std::setw(4) << el << " ";
		}
		std::cout << std::left << " ]" << std::endl;
	std::cout << std::left << std::setw(32) << "# print_info:" << std::setw(20) <<  (print_info?"true":"false") << std::endl;
	std::cout << std::left << std::setw(32) << "# print in folder:" << std::setw(20) <<  trainingdata_folder << std::endl;
	std::cout << std::left << std::setw(32) << "# verbose:" << std::setw(20) <<  (verbose?"true":"false") << std::endl;
	std::cout << std::left << std::setw(32) << "# write during training:" << std::setw(20) <<  (write_during_training?"true":"false") << std::endl;
	std::cout << std::left << std::setw(32) << "# print_paramset:" << std::setw(20) <<  (print_paramset?"true":"false") << std::endl;
	std::cout << std::left << std::setw(32) << "# erase snapshot params:" << std::setw(20) <<  (erase_snapshot_params?"true":"false") << std::endl;
	std::cout << std::left << std::setw(32) << "# orthonormalize bfs:" << std::setw(20) <<  (orthonormalize_bfs?"true":"false") << std::endl;
	std::cout << std::left << std::setw(32) << "# tighten tolerance:" << std::setw(20) <<  (tighten_tol?"true":"false")
			 << ((tighten_tol && snapshot_tol_red_crit == repeated_param)? " (at repeated param choice)": " (at Greedy conv. degradation)") << std::endl;
	std::cout << std::left << std::setw(32) << "# tighten tolerance Riesz A:" << std::setw(20) <<  (tighten_tol_rieszA?"true":"false") << std::endl;
	std::cout << std::left << std::setw(32) << "# tighten tolerance Riesz F:" << std::setw(20) <<  (tighten_tol_rieszF ?"true":"false") << std::endl;
	std::cout << std::left << std::setw(32) << "# tighten tolerance red.:" << std::setw(20) <<  tighten_tol_reduction << std::endl;
	std::cout << std::left << std::setw(32) << "# update snapshots:" << std::setw(20) << (update_snapshot ?"true":"false") << std::endl;
	std::cout << std::left << std::setw(32) << "# update Riesz Repr F:" << std::setw(20) << (update_rieszF ?"true":"false") << std::endl;
	std::cout << std::left << std::setw(32) << "# update Riesz Repr A:" << std::setw(10) << (update_rieszA ?"true":"false") << "  " << (coarsen_rieszA_for_update?"(coarsened)":"") << std::endl;
	std::cout << std::left << std::setw(32) << "# test estimator equivalence:" << std::setw(20) << (test_estimator_equivalence?"true":"false")
																			   << " " << (tighten_estimator_accuracy?"(with adjustment)":"(no adjustment)") << std::endl;
	std::cout << std::left << std::setw(32) << "#  riesz_constant_X:" << std::setw(20) << riesz_constant_X << std::endl;
	std::cout << std::left << std::setw(32) << "#  riesz_constant_Y:" << std::setw(20) << riesz_constant_Y << std::endl;
	std::cout << std::left << std::setw(32) << "# write direct representors:" << std::setw(20) <<  (write_direct_representors?"true":"false") << std::endl;
	std::cout << std::left << std::setw(32) << "# min_error_reduction (if conv_rate_degrad.):" << std::setw(20) <<  min_error_reduction << std::endl;
	std::cout << std::left << std::setw(32) << "# reference sol accuracy (w.r.t snapshot acc):" << std::setw(20) <<  refSolution_tol_factor << std::endl;
	std::cout << std::left << std::setw(32) << "# read in truth sols? :" << std::setw(20) <<  read_truth_sols << std::endl;
	std::cout << std::left << std::setw(32) << "#       yes: how many?:" << std::setw(20) <<  nb_existing_truth_sols << std::endl;
	std::cout << "####################################################" << std::endl << std::endl;

}

template<typename ParamType>
void
RB_Greedy_Information<ParamType>::RB_Greedy_Information::
print(const char* filename)
{
	std::size_t pdim = ParamInfo<ParamType>::dim;
    std::ofstream infofile(filename);
    if(infofile.is_open()){
    	infofile << "# N Error Mu Size(u) Size(ReprF_qf) ... Size(Repr_A_bf_qa)" << std::endl;
    	for(std::size_t i=0; i < greedy_errors.size(); ++i){
    		// Write Error
    		infofile << i << " " << greedy_errors[i] << " ";
    		// Write Parameter
    		for(std::size_t d = 0; d < pdim; ++d){
    			infofile << snapshot_params[i][d] << " ";
    		}
    		// Write Size of Solution
    		infofile << u_size[i] << " ";
    		// Write Size of Riesz Reprs of F
    		for(auto& el : repr_f_size[i]){
    			infofile << el << " ";
    		}
    		// Write Size of Riesz Reprs of A
			for(auto& el : repr_a_size[i]){
				infofile << el << " ";
			}
    		infofile << std::endl;
    	}
        infofile.close();
    }
    else{
    	std::cerr << "Error opening file " << filename << " for writing! " << std::endl;
    }

    if(repr_r_size.size() > 0){
    	std::string file(filename);
    	file += "_repr_res.txt";
    	std::ofstream infofile_r(file.c_str());
    	if(infofile_r.is_open()){
    		for(std::size_t n = 0; n < repr_r_size[0].size(); ++n){
        		for(std::size_t i = 0; i < repr_r_size.size(); ++i){
        			infofile_r << repr_r_size[i][n] << " ";
        		}
        		infofile_r << std::endl;
    		}

    		infofile_r.close();
    	}
        else{
        	std::cerr << "Error opening file " << file << " for writing! " << std::endl;
        }
    }

    std::string accfilename(filename);
    accfilename += "_repr_accuracies ";
    std::ofstream accfile(accfilename);
    if(accfile.is_open()){
    	accfile << "# N Eps_F_1 ... Eps_F_Qf Eps_A_1_1 .. Eps_A_Qa_1 Eps_A_1_2 .. Eps_A_Qa_2 .. " << std::endl;
    	for(std::size_t N = 0; N < accuracies_f_reprs.size(); ++N){
    		accfile << N;
    		for(std::size_t qf = 0; qf < accuracies_f_reprs[N].size(); ++qf){
    			accfile << " " << accuracies_f_reprs[N][qf];
    		}

    		for(std::size_t n = 0; n < accuracies_a_reprs[N].size(); ++n){
    			for(std::size_t qa = 0; qa < accuracies_a_reprs[N][n].size(); ++qa){
    				accfile << " " << accuracies_a_reprs[N][n][qa] << " ";
    			}
    		}
    		accfile << std::endl;
    	}

    	accfile.close();
    }
    else{
    	std::cerr << "Error opening file " << accfilename << " for writing! " << std::endl;
    }
}

template<typename ParamType>
void
RB_Greedy_Information<ParamType>::RB_Greedy_Information::
read(const char* filename, std::size_t Qf, std::size_t Qa, int nb)
{
	std::size_t pdim = ParamInfo<ParamType>::dim;
    std::ifstream infofile(filename);
    if(infofile.is_open()){
    	std::string header;
    	std::getline(infofile, header);

    	int i;
    	typename ParamInfo<ParamType>::ValueType err, size;

    	int countmax = (nb < 0) ? 10000 : nb;
    	int count = 0;
    	while(!infofile.eof() && count < countmax){
    		// Read Error
    		infofile >> i;
        	std::cout << "i = " << i << std::endl;
        	infofile >> err;
        	std::cout << "err = " << err << std::endl;
        	greedy_errors.push_back(err);

    		// Read Parameter
        	ParamType mu;
    		for(std::size_t d = 0; d < pdim; ++d){
    			infofile >> mu[d];
            	std::cout << "mu_" << d <<" = " << mu[d] << std::endl;
    		}
    		snapshot_params.push_back(mu);

    		// Read Size of Solution
    		infofile >> size;
    		std::cout << "size u = " << size << std::endl;
    		u_size.push_back(size);

			// Read Size of Riesz Reprs of F
    		std::vector<std::size_t> sizes_f;
    		for(std::size_t d = 0; d < Qf; ++d){
				infofile >> size;
				std::cout << "size f_" << d << " = " << size << std::endl;
				sizes_f.push_back(size);
			}
    		repr_f_size.push_back(sizes_f);

    		// Read Size of Riesz Reprs of A
    		std::vector<std::size_t> sizes;
    		for(std::size_t d = 0; d < Qa; ++d){
    			infofile >> size;
        		std::cout << "size a_" << d << " = " << size << std::endl;
        		sizes.push_back(size);
    		}
    		repr_a_size.push_back(sizes);

			count++;
    	}

        infofile.close();

        // Test
        print("test.txt");

    }
    else{
    	std::cerr << "Error opening file " << filename << " for reading! " << std::endl;
    }
}

template<typename ParamType>
void
RB_Greedy_Information<ParamType>::RB_Greedy_Information::
read_repr_accuracies(const char* filename, std::size_t Qf, std::size_t Qa, int Nmax)
{
    std::ifstream accfile(filename);
    if(accfile.is_open()){
    	std::string header;
    	std::getline(accfile, header);
    	
        int N;
    	for(std::size_t n = 0; n < Nmax; ++n){
    		accfile >> N;
    		
            std::vector<double> eps_f_n(Qf);
    		for(std::size_t qf = 0; qf < Qf; ++qf){
    			accfile >> eps_f_n[qf];
    		}
            accuracies_f_reprs.push_back(eps_f_n);

            std::vector<double>                 eps_a_n_q(Qa);
            std::vector<std::vector<double> >   eps_a_n;
    		for(std::size_t m = 0; m <= n; ++m){
    			for(std::size_t qa = 0; qa < Qa; ++qa){
    				accfile >> eps_a_n_q[qa];
    			}
                eps_a_n.push_back(eps_a_n_q);
    		}
    		accuracies_a_reprs.push_back(eps_a_n);
    	}
    	accfile.close();
    }
    else{
    	std::cerr << "Error opening file " << filename << " for writing! " << std::endl;
    }
}

template<typename ParamType>
void
RB_Greedy_Information<ParamType>::RB_Greedy_Information::
print_bound_info(const char* filename)
{
	for(std::size_t i = 0; i < eps_res_bound.size(); ++i){
		std::stringstream filename_boundN;
		filename_boundN << filename << "_bound_N_" << i << ".txt";
		std::ofstream file(filename_boundN.str().c_str());
		if(file.is_open()){
			for(auto& entry : eps_res_bound[i]){
				file << entry << std::endl;
			}
			file.close();
		}
	    else{
	    	std::cerr << "Error opening file " << filename_boundN.str() << " for writing! " << std::endl;
	    }
	}

	for(std::size_t i = 0; i < eps_aff.size(); ++i){
		std::stringstream filename_epsN;
		filename_epsN << filename << "_aff_N_" << i << ".txt";
		std::ofstream file(filename_epsN.str().c_str());
		if(file.is_open()){
			for(auto& entry : eps_aff[i]){
				file << entry << std::endl;
			}
			file.close();
		}
	    else{
	    	std::cerr << "Error opening file " << filename_epsN.str() << " for writing! " << std::endl;
	    }
	}
}

template<typename ParamType>
RB_Parameters<ParamType>::RB_Parameters(SolverCall _call, ParamType _ref_param, bool _verbose)
 : call(_call), ref_param(_ref_param), verbose(_verbose)
{}

template<typename ParamType>
void
RB_Parameters<ParamType>::print()
{
	std::cout << "###### RB Parameters #################" << std::endl;
	std::cout << std::left << std::setw(24) << "# SolverType:" << std::setw(20);
		switch(call){
			case call_cg: 	std::cout << "cg" << std::endl;
							break;
			case call_gmres:std::cout << "gmres" << std::endl;
							break;
			default:		std::cout << "<unrecognized>" << std::endl;
							break;
		}
	std::cout << std::left << std::setw(24) << "# Ref_Param:" << std::setw(2) << "[ ";
		for(auto& el : ref_param){
			std::cout << std::right << std::setw(4) << el << " ";
		}
		std::cout << std::left << " ]" << std::endl;
	std::cout << std::left << std::setw(24) << "# verbose:" << std::setw(20) <<  (verbose?"true":"false") << std::endl;
	std::cout << "#########################################" << std::endl << std::endl;
}

} // namespace lawa
