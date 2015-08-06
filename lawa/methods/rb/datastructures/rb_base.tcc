#include <algorithm>
#include <iomanip>
#include <sys/stat.h>
#include <tuple>

namespace lawa {

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
RB_Base<RB_Model,TruthModel, DataType, ParamType>::RB_Base(RB_Model& _rb_system, TruthModel& _rb_truth)
 : rb_system(_rb_system), rb_truth(_rb_truth)
{}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
std::size_t
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
n_bf()
{
	return rb_basisfunctions.size();
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
DataType&
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
get_bf(std::size_t i)
{
	assert(i < rb_basisfunctions.size());
	return rb_basisfunctions[i];
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
train_Greedy(std::size_t N){
    //------------------------------------------------//
	//      Parameter training set generation
	//------------------------------------------------//
	std::vector<ParamType> Xi_train = generate_uniform_paramset(greedy_params.min_param,
																greedy_params.max_param,
																greedy_params.nb_training_params,
																greedy_params.log_scaling);
																
    train_Greedy(Xi_train, N);
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
train_Greedy(std::vector<ParamType>& Xi_train, std::size_t N)
{

	//================================================//
	//      Initial Computations
	//================================================//

	if(greedy_params.verbose){
		std::cout.precision(12);
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=========      OFFLINE TRAINING  " << (greedy_params.training_type==strong?"(STRONG)":"(WEAK)")
	    												   << (greedy_params.training_type==weak_direct?" (direct)":"") << "            ================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl << std::endl;

	    std::cout << "Starting with N = " << N << std::endl << std::endl;

	}

	if(greedy_params.print_paramset){
		std::cout << "Training Parameters: " << std::endl;
		print_paramset(Xi_train);
	}


	//------------------------------------------------//
	//      Strong Greedy snapshot calculations
	//------------------------------------------------//
	// Preparations for strong Greedy (Output Folder / Calculation of Truth Sols)
	std::map<ParamType, DataType> truth_sols;
	if(greedy_params.training_type == strong || greedy_params.training_type == strong_adaptive){
		std::string truthdatafolder;
		std::string paraminfo_filename;
		std::ofstream paraminfo_file;
		if(greedy_params.write_during_training || greedy_params.read_truth_sols){
			truthdatafolder = greedy_params.trainingdata_folder + "/truthsolutions";

            if(greedy_params.write_during_training){
    			// Make a directory to store all the basisfunction files
    			if(mkdir(greedy_params.trainingdata_folder.c_str(), 0777) == -1)
    			{
    				if(rb_system.rb_params.verbose){
    					  std::cerr << "         [In RB_Base::train_strong_Greedy: Directory "
    							    << greedy_params.trainingdata_folder << " already exists, overwriting contents.]" << std::endl;
    				}
    			}
    			if(mkdir(truthdatafolder.c_str(), 0777) == -1)
    			{
    				if(rb_system.rb_params.verbose){
    					  std::cerr << "         [In RB_Base::train_strong_Greedy: Directory "
    							    << truthdatafolder << " already exists, overwriting contents.]" << std::endl;
    				}
    			}

    			// Open file for parameter information
    			paraminfo_filename = truthdatafolder + "/paraminfo.txt";
    		}
		}

		// Calculate all truth solutions
		std::size_t count = 0;

		// If we want to compute the error w.r.t. a reference solution (!= truth),
		// lower the AWGM tolerance (needed in adaptive RBM)
		T snapshot_acc;
		if(greedy_params.training_type == strong_adaptive){
			snapshot_acc = rb_truth.access_solver().access_params().tol;
			rb_truth.access_solver().access_params().tol *= greedy_params.refSolution_tol_factor;
		}
	    for (auto& mu : Xi_train) {
            DataType u;
	        if(greedy_params.read_truth_sols && greedy_params.nb_existing_truth_sols >= count){
	            std::stringstream filename;
	    		filename << truthdatafolder << "/truthsol_" << count << ".txt";
                readCoeffVector2D(u, filename.str().c_str());
	        }
	        else{
    			u = rb_truth.get_truth_solution(mu);	            
	        }
	    	truth_sols.insert(std::make_pair(mu, u));
	    	count++;

	    	if(greedy_params.write_during_training){
	    		if(!(greedy_params.read_truth_sols && greedy_params.nb_existing_truth_sols > count)){
    	    		std::stringstream filename;
    	    		filename << truthdatafolder << "/truthsol_" << count << ".txt";
    			    saveCoeffVector2D(u, rb_truth.get_trialbasis(), filename.str().c_str());

    			    paraminfo_file.open(paraminfo_filename.c_str(), std::fstream::app);
    			    paraminfo_file << count;
    			    for(auto& d : mu){
    			    	paraminfo_file << " " << d;
    			    }
    			    paraminfo_file << std::endl;
    			    paraminfo_file.close();
    			}
	    	}
	    }
	    // Reset snapshot tolerance
		if(greedy_params.training_type == strong_adaptive){
			rb_truth.access_solver().access_params().tol = snapshot_acc;
		}
	}

	//------------------------------------------------//
	//      F Riesz representors
	//------------------------------------------------//
	if(N == 0){
		// In order to be able to calculate empty error bounds,
		// we have to calculate the Riesz Representors for F
		calculate_Riesz_RHS_information();
	}

	//================================================//
	//      Greedy Iterations
	//================================================//

	ParamType current_param = Xi_train[0];
	T max_error = 0;
	bool repeated_iteration = false;
	long double init_tol = rb_truth.access_solver().access_params().tol;
    std::vector<std::pair<ParamType, unsigned int>> repeated_snaps;
	do {

		if(greedy_params.verbose){
		    std::cout << "||---------------------------------------------------------------------||" << std::endl;
		    std::cout << "||-------------- Greedy Search for new parameter  ---------------------||" << std::endl;
		    std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}

		if(greedy_params.test_estimator_equivalence){
			std::vector<double> bounds, eps_aff;
			greedy_info.eps_res_bound.push_back(bounds);
			greedy_info.eps_aff.push_back(eps_aff);
		}

		if(greedy_params.training_type == weak_direct){
			std::vector<std::size_t> res_repr_sizes;
			greedy_info.repr_r_size.push_back(res_repr_sizes);
		}

		//------------------------------------------------//
		//  Find Parameter with maximum error (estimator)
		//------------------------------------------------//
		T error_est;
		max_error = find_max_errest(N, Xi_train, current_param, truth_sols);
        if (max_error <= greedy_params.tol) {
            std::cout << "Greedy tolerance reached\n";
            std::cout << "Greedy error     : " << max_error << std::endl;
            std::cout << "Greedy tolerance : " << greedy_params.tol << std::endl;
            break;
        }

		//------------------------------------------------//
		//      Check conditions for tolerance tightening
		//------------------------------------------------//
		bool update_snapshot = false;
		DataType u;
		if(greedy_params.tighten_tol){

			bool tol_red_crit = false;

			switch(greedy_params.snapshot_tol_red_crit){
			case repeated_param:
			{
				// find (first occurence of) current parameter
				for(auto& mu : greedy_info.snapshot_params){
					bool is_current_param = true;
					for(std::size_t i = 0; i < mu.size(); ++i){
						//if(mu[i]!=current_param[i]){
						if(std::abs(mu[i] - current_param[i]) > 1e-05){
							is_current_param = false;
							break;
						}
                        // has parameter repeated before?
                        if (i == mu.size()-1) {
                            std::cout << "Repeating parameter is ["
                                      << current_param[0]
                                      << ", "
                                      << current_param[1]
                                      << "]\n";
                            bool repeated_before = false;
                            for (auto& pair : repeated_snaps) {
                                std::cout << "Checking repeated_snaps at ["
                                          << pair.first[0]
                                          << ", "
                                          << pair.first[1]
                                          << "]\n";
                                repeated_before = true;
                                ParamType&  mu_rep = pair.first;
                                for (std::size_t j = 0; j < mu_rep.size();
                                                        ++j) {
                                    if(std::abs(mu_rep[j] - current_param[j]) > 1e-05) {
                                        repeated_before = false;
                                        break;
                                    }
                                }
                                // yes
                                if (repeated_before) {
                                    std::cout << "Yes, parameter has repeated before\n";
                                    std::cout << "Increasing repetitions number to "
                                              << ++pair.second
                                              << std::endl;
                                    break;
                                }
                            }
                            // no
                            if (!repeated_before) {
                                std::cout << "No, first time\n";
                                std::pair<ParamType, unsigned int>
                                             rep_snap(current_param, 1);
                                repeated_snaps.push_back(rep_snap);
                            }
                        }
					}
					if(is_current_param){
						tol_red_crit = true;
						update_snapshot = true;
						break;
					}
				}
				break;
			}
			case conv_rate_degradation:
			{
				// First check criterion
				if(repeated_iteration == false && N >= 2 && log10(greedy_info.greedy_errors[N-1]/max_error)/log10(greedy_info.greedy_errors[N-2]/greedy_info.greedy_errors[N-1]) < greedy_params.min_error_reduction){
					if(greedy_params.verbose){
						std::cout << std::endl<< "----> Error reduction: " << std::endl;
						std::cout << "      Greedy Errors: N-1 = " <<greedy_info.greedy_errors[N-2] << ", N = " << greedy_info.greedy_errors[N-1] << ", error est = " << max_error << std::endl;
						std::cout << "          with log10: N-1 = " << log10(greedy_info.greedy_errors[N-2]) << ", N = " << log10(greedy_info.greedy_errors[N-1]) << ", error est = " << log10(max_error) << std::endl;
						std::cout << "          -> reduction factor " << log10(greedy_info.greedy_errors[N-1]/max_error)/log10(greedy_info.greedy_errors[N-2]/greedy_info.greedy_errors[N-1])
								  << " < " << greedy_params.min_error_reduction << std::endl;
						std::cout << "=====> Recalculate Snapshot Nb. " << N << std::endl << std::endl;

					}
					tol_red_crit = true;
					repeated_iteration = true;

					// Reset training to last iteration
					u = rb_basisfunctions[rb_basisfunctions.size()-1];

					max_error = greedy_info.greedy_errors[greedy_info.greedy_errors.size()-1];
					current_param = greedy_info.snapshot_params[greedy_info.snapshot_params.size()-1];

					remove_basisfunction(N, true);
					if(greedy_params.verbose){
						std::cout << "       Removed old snapshot and corresponding information" << std::endl;
					}

					N = N-1;

					update_snapshot=true;
				}
				else{
					repeated_iteration = false;
				}

				break;
			}
			default:
                std::cerr << "Snapshot Tolerance Adaptation Method not implemented yet!" << std::endl;
                exit(1);
			}


			if(tol_red_crit){
                unsigned int max_reps = 0;
                // current reduction level
                for (auto& pair : repeated_snaps) {
                    if (pair.second > max_reps) {
                        max_reps = pair.second;
                    }
                }
                std::cout << "====> Current parameter repetition level is "
                          << max_reps << std::endl;
				rb_truth.access_solver().access_params().tol = init_tol*
                    std::pow(greedy_params.tighten_tol_reduction,
                            (double) max_reps);
				if(greedy_params.verbose){
					std::cout << std::endl<< "====> -------------------------------------------------------------- <=====" << std::endl;
					std::cout << "====>     Reduced snapshot tolerance to " << std::scientific << rb_truth.access_solver().access_params().tol << std::endl;
					std::cout << "====> -------------------------------------------------------------- <=====" << std::endl;

				}

				if(greedy_params.tighten_tol_rieszA){
					rb_truth.access_RieszSolver_A().access_params().tol *= greedy_params.tighten_tol_reduction;
					if(greedy_params.verbose){
						std::cout << std::endl<< "====> -------------------------------------------------------------- <=====" << std::endl;
						std::cout << "====>     Reduced Riesz Repr A tolerance to " << std::scientific << rb_truth.access_RieszSolver_A().access_params().tol << std::endl;
						std::cout << "====> -------------------------------------------------------------- <=====" << std::endl;
					}
				}
				if(greedy_params.tighten_tol_rieszF){
					rb_truth.access_RieszSolver_F().access_params().tol *= greedy_params.tighten_tol_reduction;
					if(greedy_params.verbose){
						std::cout << std::endl<< "====> -------------------------------------------------------------- <=====" << std::endl;
						std::cout << "====>     Reduced Riesz Repr F tolerance to " << std::scientific << rb_truth.access_RieszSolver_F().access_params().tol << std::endl;
						std::cout << "====> -------------------------------------------------------------- <=====" << std::endl;
					}
					// Save old Riesz Representor
					if(greedy_params.write_during_training){
						if(mkdir(greedy_params.trainingdata_folder.c_str(), 0777) == -1)
						{
							if(greedy_params.verbose){
								  std::cerr << "         [In RB_Base::train_Greedy: Directory "
											<< greedy_params.trainingdata_folder << " already exists, overwriting contents.]" << std::endl;
							}
						}
						std::string repr_folder = greedy_params.trainingdata_folder + "/representors";
						if(mkdir(repr_folder.c_str(), 0777) == -1)
						{
							if(greedy_params.verbose){
								  std::cerr << "         [In RB_Base::train_Greedy: Directory "
											<< repr_folder << " already exists, overwriting contents.]" << std::endl;
							}
						}
						std::stringstream old_repr_folder;
						old_repr_folder << greedy_params.trainingdata_folder << "/representors/F_repr_before_it_" << N+1;
						write_rieszrepresentors(old_repr_folder.str(), 0);
					}
					// Recalculate F Representors
					bool update = false;
					if(greedy_params.update_rieszF){
						update = true;
					}
					calculate_Riesz_RHS_information(update);
					recalculate_A_F_norms();
				}
			}
		}

		greedy_info.greedy_errors.push_back(max_error);
		greedy_info.snapshot_params.push_back(current_param);


		//------------------------------------------------//
		//      Update parameter training set
		//------------------------------------------------//
		// Remove parameter from indexset
		if(greedy_params.erase_snapshot_params){
			// find current parameter
			for(auto& mu : Xi_train){
				bool is_current_param = true;
				for(std::size_t i = 0; i < mu.size(); ++i){
					//if(mu[i]!=current_param[i]){
					if(std::abs(mu[i] - current_param[i]) > 1e-05){
						is_current_param = false;
						break;
					}
				}
				if(is_current_param){
					auto mu_it = std::find(Xi_train.begin(), Xi_train.end(), mu);
					Xi_train.erase(mu_it);
					break;
				}
			}
		}


		//------------------------------------------------//
		//      Compute next snapshot
		//------------------------------------------------//
		if(greedy_params.verbose){
			std::cout << std::endl << "Greedy Error = " << std::scientific << max_error << std::endl << std::endl;
		}

		if(greedy_params.verbose){
		    std::cout << "||=====================================================================||" << std::endl;
		    std::cout << "||       SNAPSHOT  at ";  ParamInfo<ParamType>::print(current_param); std::cout << std::endl;
		    std::cout << "||=====================================================================||" << std::endl << std::endl;
		}

        exit(0);

		if(rb_truth.access_solver().access_params().print_info){
			std::stringstream filename;
			filename << greedy_params.print_file << "_bf" << N+1 << ".txt";
			rb_truth.access_solver().access_params().info_filename = filename.str();
		}

		if(false/*update_snapshot && greedy_params.update_snapshot*/){

			if(greedy_params.snapshot_tol_red_crit == conv_rate_degradation){
				if(greedy_params.verbose){
					std::cout << " Repeating truth computation for snapshot nb " << N+1 << ", starting with old solution." << std::endl << std::endl;
				}
			}
			else{
				// Find last corresponding snapshot:
				// Find (last occurence of) current parameter
				int last_snapshot_index;
				for(auto& mu : greedy_info.snapshot_params){
					bool is_current_param = true;
					for(std::size_t i = 0; i < mu.size(); ++i){
						//if(mu[i]!=current_param[i]){
						if(std::abs(mu[i] - current_param[i]) > 1e-05){
							is_current_param = false;
							break;
						}
					}
					if(is_current_param){
						if(&mu - &greedy_info.snapshot_params[0] < (int)greedy_info.snapshot_params.size()-1){
							last_snapshot_index = &mu - &greedy_info.snapshot_params[0];
						}
					}
				}

				if(greedy_params.verbose){
					std::cout << " Starting truth computation with snapshot nb " << last_snapshot_index+1 << std::endl << std::endl;
				}

				u = rb_basisfunctions[last_snapshot_index];
			}


			// Start solver with copy of old snapshot
			rb_truth.get_truth_solution(current_param, u);

			if(greedy_params.training_type == strong){
				truth_sols.erase(current_param);
		    	truth_sols.insert(std::make_pair(current_param, u));
			}
		}
		else{
			if(greedy_params.training_type == strong){
				for(auto& entry : truth_sols){
					bool is_current_mu = true;
					for(std::size_t i = 0; i < current_param.size(); ++i){
						//if(entry.first[i]!=current_param[i]){
						if(std::abs(entry.first[i] - current_param[i]) > 1e-05){
							is_current_mu = false;
							break;
						}
					}
					if(is_current_mu){
						u = entry.second;
					}
				}
			}
			else{
                double delta    = 1./2.;
                double beta_mu  = rb_system.alpha_LB(current_param);
				//rb_truth.access_solver().access_params().tol_primal =
                //    greedy_params.tol*delta*beta_mu;
				rb_truth.access_solver().access_params().tol =
                    greedy_params.tol*delta*beta_mu*beta_mu;
                std::cout << "beta(mu) = " << rb_system.alpha_LB(current_param)
                          << std::endl;
                std::cout << "delta    = " << delta << std::endl;
                std::cout << "tol      = " << greedy_params.tol << std::endl;
                std::cout << "and thus...\n";
                std::cout << "Tolerance for snapshot set to\n"
                          << "Dual   : "
                          << std::scientific
                          << rb_truth.access_solver().access_params().tol
                //          << std::endl
                //          << "Primal : "
                //          << std::scientific
                //          << rb_truth.access_solver().access_params().tol_primal
                          << std::endl;
				u = rb_truth.get_truth_solution(current_param);
			}
		}

		//------------------------------------------------//
		//      Update RB basis and infos
		//------------------------------------------------//

        std::string snap_folder = greedy_params.trainingdata_folder + "/snap";
        mkdir(greedy_params.trainingdata_folder.c_str(), 0777);
        write_snapshot(u, snap_folder, (int)N+1);
		add_to_basis(u);
		N++;

		greedy_info.u_size.push_back(u.size());
		std::vector<std::size_t> f_sizes;
		for(auto& el : F_representors){
			f_sizes.push_back(el.size());
		}
		greedy_info.repr_f_size.push_back(f_sizes);
		std::vector<std::size_t> a_sizes;
		for(auto& el : A_representors[n_bf()-1]){
			a_sizes.push_back(el.size());
		}
		greedy_info.repr_a_size.push_back(a_sizes);
		greedy_info.accuracies_f_reprs.push_back(eps_F_representors);
		greedy_info.accuracies_a_reprs.push_back(eps_A_representors);

		if(greedy_params.write_during_training){

			if(greedy_params.verbose){
				std::cout << "=====>  Writing RB Greedy Training information to file " << std::endl << std::endl;
			}
			if(mkdir(greedy_params.trainingdata_folder.c_str(), 0777) == -1)
			{
				if(greedy_params.verbose){
					  std::cerr << "         [In RB_Base::train_Greedy: Directory "
							    << greedy_params.trainingdata_folder << " already exists, overwriting contents.]" << std::endl;
				}
			}

			// Write Basis Functions
			std::string bf_folder = greedy_params.trainingdata_folder + "/bf";
			write_basisfunctions(bf_folder, (int)N);

			// Write Riesz Representors
			std::string repr_folder = greedy_params.trainingdata_folder + "/representors";
			write_rieszrepresentors(repr_folder);

			// Write RB Data
			rb_system.write_rb_data(greedy_params.trainingdata_folder);

			// Write Training Information
			greedy_info.print(greedy_params.print_file.c_str());

			// If gathered, write Error Estimator Bound Information
			if(greedy_params.test_estimator_equivalence){
				std::string eps_file = greedy_params.print_file + "_eps";
				greedy_info.print_bound_info(eps_file.c_str());
			}
		}

	} while( (N < greedy_params.Nmax) && (max_error > greedy_params.tol) );

	if(greedy_params.print_info){
		greedy_info.print();
	}
}

/*template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
train_strong_Greedy()
{
	if(greedy_params.verbose){
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=========      OFFLINE TRAINING                      ================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl << std::endl;

	}

	std::vector<ParamType> Xi_train = generate_uniform_paramset(greedy_params.min_param,
																greedy_params.max_param,
																greedy_params.nb_training_params);
	if(greedy_params.print_paramset){
		std::cout << "Training Parameters: " << std::endl;
		print_paramset(Xi_train);
	}

	std::string truthdatafolder;
	std::string paraminfo_filename;
	std::ofstream paraminfo_file;
	if(greedy_params.write_during_training){
		truthdatafolder = greedy_params.trainingdata_folder + "/truthsolutions";

		// Make a directory to store all the basisfunction files
		if(mkdir(greedy_params.trainingdata_folder.c_str(), 0777) == -1)
		{
			if(rb_system.rb_params.verbose){
				  std::cerr << "         [In RB_Base::train_strong_Greedy: Directory "
						    << greedy_params.trainingdata_folder << " already exists, overwriting contents.]" << std::endl;
			}
		}
		if(mkdir(truthdatafolder.c_str(), 0777) == -1)
		{
			if(rb_system.rb_params.verbose){
				  std::cerr << "         [In RB_Base::train_strong_Greedy: Directory "
						    << truthdatafolder << " already exists, overwriting contents.]" << std::endl;
			}
		}

		// Open file for parameter information
		paraminfo_filename = truthdatafolder + "/paraminfo.txt";
	}

	// Calculate all truth solutions
	std::map<ParamType, DataType> truth_sols;
	std::size_t count = 0;
    for (auto& mu : Xi_train) {
		DataType u = rb_truth.get_truth_solution(mu);
    	truth_sols.insert(std::make_pair(mu, u));
    	count++;

    	if(greedy_params.write_during_training){
    		std::stringstream filename;
    		filename << truthdatafolder << "/truthsol_" << count << ".txt";
		    saveCoeffVector2D(u, rb_truth.get_trialbasis(), filename.str().c_str());

		    paraminfo_file.open(paraminfo_filename.c_str(), std::fstream::app);
		    paraminfo_file << count;
		    for(auto& d : mu){
		    	paraminfo_file << " " << d;
		    }
		    paraminfo_file << std::endl;
		    paraminfo_file.close();
    	}
    }

	// In order to be able to calculate empty error bounds,
	// we have to calculate the Riesz Representors for F
	calculate_Riesz_RHS_information();

	for(auto& el : F_representors){
		greedy_info.repr_f_size.push_back(el.size());
	}

	ParamType current_param;
	std::size_t N = 0;
	T max_error = 0;
	do {

		if(greedy_params.verbose){
		    std::cout << "||---------------------------------------------------------------------||" << std::endl;
		    std::cout << "||-------------- Greedy Search for new parameter  ---------------------||" << std::endl;
		    std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}

		T error_est;
		max_error = 0;
		for(auto& mu : Xi_train){
            typename RB_Model::DenseVectorT u_N = rb_system.get_rb_solution(N, mu);

            //error_est = rb_system.get_errorbound(u_N,mu);
            DataType diff;
            if(N > 0){
            	diff = reconstruct_u_N(u_N,N);
            }
            diff -= (*truth_sols.find(mu)).second;
            error_est = std::sqrt(rb_truth.innprod_Y_u_u(diff,diff));

            if(greedy_params.verbose){
            	std::cout << "    u_N = " << u_N;
            	std::cout << "    Mu = ";
            	ParamInfo<ParamType>::print(mu);
            	std::cout << " : Error = " << std::scientific << std::setw(15) << error_est << std::endl;
            }

            if(error_est > max_error){
            	max_error = error_est;
            	current_param = mu;
            }
		}

		greedy_info.greedy_errors.push_back(max_error);
		greedy_info.snapshot_params.push_back(current_param);

		// Remove parameter from indexset
		if(greedy_params.erase_snapshot_params){
			for(auto& mu : Xi_train){
				bool is_current_param = true;
				for(std::size_t i = 0; i < mu.size(); ++i){
					if(mu[i]!=current_param[i]){
						is_current_param = false;
						break;
					}
				}
				if(is_current_param){
					auto mu_it = std::find(Xi_train.begin(), Xi_train.end(), mu);
					Xi_train.erase(mu_it);
					break;
				}
			}
		}


		if(greedy_params.verbose){
			std::cout << std::endl << "Greedy Error = " << std::scientific << max_error << std::endl << std::endl;
		}

		if(greedy_params.verbose){
		    std::cout << "||=====================================================================||" << std::endl;
		    std::cout << "||       SNAPSHOT  at ";  ParamInfo<ParamType>::print(current_param); std::cout << std::endl;
		    std::cout << "||=====================================================================||" << std::endl << std::endl;
		}

		if(rb_truth.access_solver().access_params().print_info){
			std::stringstream filename;
			filename << greedy_params.print_file << "_bf" << N+1 << ".txt";
			rb_truth.access_solver().access_params().info_filename = filename.str();
		}

		//DataType u = rb_truth.get_truth_solution(current_param);
		DataType u;
		for(auto& entry : truth_sols){
			bool is_current_mu = true;
			for(std::size_t i = 0; i < current_param.size(); ++i){
				if(entry.first[i]!=current_param[i]){
					is_current_mu = false;
					break;
				}
			}
			if(is_current_mu){
				u = entry.second;
			}
		}
		add_to_basis(u);
		N++;

		greedy_info.u_size.push_back(u.size());
		std::vector<std::size_t> a_sizes;
		for(auto& el : A_representors[n_bf()-1]){
			a_sizes.push_back(el.size());
		}
		greedy_info.repr_a_size.push_back(a_sizes);

		if(greedy_params.write_during_training){

			if(greedy_params.verbose){
				std::cout << "=====>  Writing RB Greedy Training information to file " << std::endl << std::endl;
			}
			if(mkdir(greedy_params.trainingdata_folder.c_str(), 0777) == -1)
			{
				if(greedy_params.verbose){
					  std::cerr << "         [In RB_Base::train_Greedy: Directory "
							    << greedy_params.trainingdata_folder << " already exists, overwriting contents.]" << std::endl;
				}
			}

			// Write Basis Functions
			std::string bf_folder = greedy_params.trainingdata_folder + "/bf";
			write_basisfunctions(bf_folder, (int)N);

			// Write Riesz Representors
			std::string repr_folder = greedy_params.trainingdata_folder + "/representors";
			write_rieszrepresentors(repr_folder);
			if(N == 1){
				write_rieszrepresentors(repr_folder);
			}

			// Write RB Data
			rb_system.write_rb_data(greedy_params.trainingdata_folder);

			// Write Training Information
			greedy_info.print(greedy_params.print_file.c_str());
		}

	} while( (N < greedy_params.Nmax) && (max_error > greedy_params.tol) );

	if(greedy_params.print_info){
		greedy_info.print();
	}
}*/

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
DataType
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
reconstruct_u_N(typename RB_Model::DenseVectorT u, std::size_t N)
{
	assert(N <= n_bf());
	assert(u.length() > 0);
	assert(u.length() >= (int)N);

	DataType u_full;
	for (unsigned int i = 1; i <= N; ++i) {
		u_full +=  u(i) * rb_basisfunctions[i-1];
	}

	return u_full;
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
DataType
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
reconstruct_u_N(typename RB_Model::DenseVectorT u, std::vector<std::size_t> indices)
{
    size_t N = u.length();
	assert(N <= n_bf());
	assert(u.length() > 0);
	assert(u.length() >= (int)N);

	DataType u_full;
	for (unsigned int i = 1; i <= N; ++i) {
		//std::cout << "i = " << u(i) << ", bf index = " << indices[i-1]  << " with size " << rb_basisfunctions[indices[i-1]].size() << std::endl;
		u_full +=  u(i) * rb_basisfunctions[indices[i-1]];
	}

	return u_full;
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
DataType
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
reconstruct_res_repr_N(typename RB_Model::DenseVectorT u, std::size_t N, ParamType& mu)
{
	assert(N <= n_bf());
	assert(u.length() >= (int)N);


	DataType res_full;
	for(std::size_t i = 0; i < rb_system.Q_f(); ++i){
		res_full += rb_system.thetas_f.eval(i,mu) * F_representors[i];
	}
	for(std::size_t n = 1; n <= N; ++n){
		for(std::size_t i = 0; i < rb_system.Q_a(); ++i){
			res_full += u(n) * rb_system.thetas_a.eval(i,mu) * A_representors[n-1][i];
		}
	}

	return res_full;
}




template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
typename DataType::ValueType
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
get_direct_errorbound(const typename RB_Model::DenseVectorT& u_N, ParamType& mu, DataType& res_repr)
{

	DataType u_approx = reconstruct_u_N(u_N, u_N.length());
	if(greedy_params.verbose){
		std::cout << "||---------------------------------------------------------------------||" << std::endl;
		std::cout << "||-------------- Calculate Direct Riesz Representors for Residual -----||" << std::endl;
		std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
	}
	res_repr = rb_truth.get_riesz_representor_res(u_approx, mu);
    return  res_repr.norm(2.)*std::sqrt(greedy_params.riesz_constant_Y) / rb_system.alpha_LB(mu);
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
typename DataType::ValueType
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
update_direct_errorbound(const typename RB_Model::DenseVectorT& u_N, ParamType& mu, DataType& res_repr)
{

	DataType u_approx = reconstruct_u_N(u_N, u_N.length());
	if(greedy_params.verbose){
		std::cout << "||-------------- Update Direct Riesz Representors for Residual -----||" << std::endl << std::endl;
	}
	rb_truth.get_riesz_representor_res(u_approx, mu, res_repr);
    return  res_repr.norm(2.)*std::sqrt(greedy_params.riesz_constant_Y) / rb_system.alpha_LB(mu);
}


template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
std::vector<ParamType>
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
generate_uniform_paramset(ParamType min_param, ParamType max_param, intArrayType param_nb, intArrayType log_scaling)
{
	std::size_t pdim = ParamInfo<ParamType>::dim;

	// Calculate step lengths in each parameter dimension
    std::vector<T> h(pdim);
    for(std::size_t d = 0; d < pdim; ++d) {
    	if(param_nb[d] > 1){
    		if(log_scaling[d] == 0){
                h[d] = (max_param[d] - min_param[d]) / ((T)param_nb[d]-1.);
    		}
    		else{
    			h[d] = std::log10(max_param[d] / min_param[d]) / ((T)param_nb[d]-1.);
    		}
    	}
    	else{
    		h[d] = 0;
    	}
    }

    std::vector<ParamType> Xi_train;

    std::vector<std::size_t>	   index(pdim,0);

    while(true){
        ParamType new_mu;
        for(std::size_t i = 0; i < pdim; ++i){
        	if(log_scaling[i]==0){
            	new_mu[i] = std::min(min_param[i] + index[i]*h[i], max_param[i]);
        	}
        	else{
            	new_mu[i] = std::min( std::pow(10., log10(min_param[i]) + index[i]*h[i]), max_param[i]);
        	}
        }
        Xi_train.push_back(new_mu);

        for(int i = pdim-1; ; --i){
        	if(i<0){
        		return Xi_train;
        	}
        	index[i]++;
        	if(index[i]==param_nb[i]){
        		index[i] = 0;
        	}
        	else{
        		break;
        	}
        }
    }


    /*
    // Generate Parameter Grid
    if(pdim == 1) {
        for (std::size_t i = 0; i < greedy_params.nb_training_params[0]; ++i){
            ParamType new_mu;
            new_mu[0]  = std::min(greedy_params.min_param[0] + i*h[0], greedy_params.max_param[0]);
            Xi_train.push_back(new_mu);
        }
    }
    else{
        if (pdim == 2) {
            for (std::size_t i = 0; i < greedy_params.nb_training_params[0]; ++i){
                for (std::size_t j = 0; j < greedy_params.nb_training_params[1]; ++j){
                	ParamType new_mu;
                    new_mu[0]  = std::min(greedy_params.min_param[0] + i*h[0], greedy_params.max_param[0]);
                    new_mu[1]  = std::min(greedy_params.min_param[1] + j*h[1], greedy_params.max_param[1]);
                    Xi_train.push_back(new_mu);
                }
            }
        }
        else {
            std::cerr << "Generate Trainingsset for dim = " << pdim << " : Not implemented yet " << std::endl;
        }
    }
    */

    return Xi_train;
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
print_paramset(std::vector<ParamType> PSet)
{
	for(auto& el : PSet){
		ParamInfo<ParamType>::print(el);
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
add_to_basis(const DataType& u)
{
	if(greedy_params.verbose){
	    std::cout << "||---------------------------------------------------------------------||" << std::endl;
	    std::cout << "||-------------- Adding Snapshot to Basis   ---------------------------||" << std::endl;
	    std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
	}

	// =========== Orthogonalization ======================= //

	DataType new_bf = u;

	if(greedy_params.orthonormalize_bfs){
		for(auto& bf : rb_basisfunctions){
			new_bf = new_bf - bf * rb_truth.innprod_Y_u_u(bf, u);
		}
	}

	T new_bf_norm_sq = rb_truth.innprod_Y_u_u(new_bf, new_bf);
	new_bf.scale(1./std::sqrt(new_bf_norm_sq));


	rb_basisfunctions.push_back(new_bf);

	if(greedy_params.verbose){
	    std::cout << "||------- GRAM-SCHMIDT Norm: " << std::setw(10) << std::sqrt(new_bf_norm_sq) << " --------------------------||" << std::endl;

	    std::cout << std::endl << "||------- UPDATE RB Structures ----------------------------------------||" << std::endl;
	}


	add_to_RB_structures(new_bf);

	calculate_Riesz_LHS_information(new_bf);

}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
add_to_RB_structures(const DataType& bf)
{
	std::size_t Qa = rb_system.Q_a();
	std::size_t Qf = rb_system.Q_f();

	// ===== Update RB_A ====== //

	// We need some matrix to start with...
	bool first_bf = false;
    if (rb_system.RB_A_matrices.size() < Qa) {
        for (std::size_t q_a = 0; q_a < Qa; ++q_a) {
        	typename RB_Model::FullColMatrixT A(1,1);
            rb_system.RB_A_matrices.push_back(A);
        }
        first_bf = true;
    }

    for (std::size_t q_a = 0; q_a < Qa; ++q_a) {
    	// "Pad" RB_A_matrices with zeros to new size
    	std::size_t n;
		if (first_bf == true) {
			rb_system.RB_A_matrices[q_a].engine().resize(1, 1);
			n = 1;
		}
		else {
			typename RB_Model::FullColMatrixT tmp(rb_system.RB_A_matrices[q_a]);
			rb_system.RB_A_matrices[q_a].engine().resize(tmp.numRows()+1, tmp.numCols()+1);
			rb_system.RB_A_matrices[q_a](tmp.rows(), tmp.cols()) = tmp;

			n = rb_system.RB_A_matrices[q_a].numRows();
			for(unsigned int i = 1; i <= n; ++i) {
				rb_system.RB_A_matrices[q_a](i, n) = 0.;
				rb_system.RB_A_matrices[q_a](n, i) = 0.;
			}

		}

		// Compute new entries (one column, one row)
        for (unsigned int i = 1; i <= n; ++i) {
        	rb_system.RB_A_matrices[q_a](n,i) = rb_truth.lhs_u_u(q_a, bf, rb_basisfunctions[i-1]);
        }
        for (unsigned int i = 1; i < n; ++i) {
        	rb_system.RB_A_matrices[q_a](i,n) = rb_truth.lhs_u_u(q_a, rb_basisfunctions[i-1], bf);
        }

        if(greedy_params.verbose){
    		std::cout << std::endl << "||------- RB_A (" << q_a << ")  -----------------------------------------||" << std::endl;
    		std::cout << rb_system.RB_A_matrices[q_a] << std::endl;
        }
    }

	// ===== Update RB_F ====== //

	// We need some matrix to start with...
    if (rb_system.RB_F_vectors.size() < Qf) {
    	rb_system.RB_F_vectors.resize(Qf);
    	first_bf = true;
    }
    for (std::size_t q_f = 0; q_f < Qf; ++q_f) {
    	// Add one zero entry to RB_F_vectors
    	std::size_t n;
    	if (first_bf == true) {
    		rb_system.RB_F_vectors[q_f].engine().resize(1);
    		n = 1;
    	}
    	else {
    		typename RB_Model::DenseVectorT tmp(rb_system.RB_F_vectors[q_f]);
    		rb_system.RB_F_vectors[q_f].engine().resize(tmp.length()+1);
    		n = rb_system.RB_F_vectors[q_f].length();
    		rb_system.RB_F_vectors[q_f](tmp.range()) = tmp;
    		rb_system.RB_F_vectors[q_f](n) = 0.;

    	}

    	rb_system.RB_F_vectors[q_f](n) = rb_truth.rhs_u(q_f, bf);

    	if(greedy_params.verbose){
    		std::cout << std::endl << "||------- RB_F (" << q_f << ")  -----------------------------------------||" << std::endl;
    		std::cout << rb_system.RB_F_vectors[q_f] << std::endl;
    	}
    }

	// ===== Update RB_InnerProduct ====== //

    std::size_t n;
    if (rb_system.RB_inner_product.numRows()==0) {
        rb_system.RB_inner_product.engine().resize(1, 1);
        n=1;
    }
    else {
    	typename RB_Model::FullColMatrixT tmp(rb_system.RB_inner_product);
    	rb_system.RB_inner_product.engine().resize(tmp.numRows()+1, tmp.numCols()+1);
    	rb_system. RB_inner_product(tmp.rows(), tmp.cols()) = tmp;
    	n = rb_system. RB_inner_product.numRows();
        for(unsigned int i = 1; i <= n; ++i) {
        	rb_system.RB_inner_product(i, n) = 0.;
        	rb_system.RB_inner_product(n, i) = 0.;
        }
    }

    for (unsigned int i = 1; i <= n; ++i) {
    	rb_system.RB_inner_product(n, i) = rb_truth.innprod_Y_u_u(bf, rb_basisfunctions[i-1]);
		rb_system.RB_inner_product(i, n) = rb_system.RB_inner_product(n,i);
    }

	if(greedy_params.verbose){
		std::cout << std::endl << "||------- RB_InnerProduct  ----------------------------------||" << std::endl;
		std::cout << rb_system.RB_inner_product << std::endl;
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
calculate_Riesz_RHS_information(bool update)
{

	if(!update){
		F_representors.clear();
		eps_F_representors.resize(rb_system.Q_f());

		if(greedy_params.verbose){
			std::cout << "||---------------------------------------------------------------------||" << std::endl;
			std::cout << "||-------------- Calculate Riesz Representors for F -------------------||" << std::endl;
			std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}
	}
	else{

		if(greedy_params.verbose){
			std::cout << "||---------------------------------------------------------------------||" << std::endl;
			std::cout << "||-------------- Update Riesz Representors for F --------   -----------||" << std::endl;
			std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}
	}
	// Calculate the Riesz Representors for F
	std::size_t Qf = rb_system.Q_f();
	for (unsigned int i = 0; i < Qf; ++i) {

		if(greedy_params.verbose){
			std::cout << "||------- Component Nb " << std::setw(3) << i << "  -------------------------------------------||" << std::endl << std::endl;
		}

		T eps;
		if(update){
			eps = rb_truth.get_riesz_representor_f(i, F_representors[i]);
		}
		else{
			DataType c;
			eps = rb_truth.get_riesz_representor_f(i, c);
			if(std::fabs(eps-1.) < 1e-8){
				eps = eps_F_representors[i] / std::sqrt(greedy_params.riesz_constant_Y);
			}
			F_representors.push_back(c);
		}
		eps_F_representors[i] = eps * std::sqrt(greedy_params.riesz_constant_Y);
	}

	// Update the Riesz Representor Norms
	rb_system.F_F_representor_norms.engine().resize((int)rb_system.Q_f(), (int)rb_system.Q_f());
	DataType vec1, vec2;
	for(std::size_t qf1 = 1; qf1 <= Qf; ++qf1) {
		for (std::size_t qf2 = qf1; qf2 <= Qf; ++qf2) {
			vec1 = F_representors[qf1-1];
			vec2 = F_representors[qf2-1];

			rb_system.F_F_representor_norms(qf1,qf2) = rb_truth.innprod_Y_v_v(vec1, vec2);

			if(qf1 != qf2) {
				rb_system.F_F_representor_norms(qf2,qf1) = rb_system.F_F_representor_norms(qf1,qf2);
			}
		}
	}

	if(greedy_params.verbose){
		std::cout << std::endl << "||------- Representor Norms F x F  ------------------------------------||" << std::endl;
		std::cout << rb_system.F_F_representor_norms << std::endl;
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
calculate_Riesz_LHS_information(const DataType& bf, bool update)
{

	// Check if we only do an update and search for "old" bf index
	int last_snapshot_index = -1;
	if(update || greedy_params.update_rieszA){
	    if(update){
            last_snapshot_index = greedy_info.snapshot_params.size()-1;
	    }
	    else{
    		for(auto& mu : greedy_info.snapshot_params){
    			bool is_last_param = true;
    			for(std::size_t i = 0; i < mu.size(); ++i){
    				if(mu[i]!=greedy_info.snapshot_params[greedy_info.snapshot_params.size()-1][i]){
    					is_last_param = false;
    					break;
    				}
    			}
    			if(is_last_param){
    				if(&mu - &greedy_info.snapshot_params[0] < (int)greedy_info.snapshot_params.size()-1){
    					last_snapshot_index = &mu - &greedy_info.snapshot_params[0];
    				}
    			}
    		}	        
	    }

		
		if(last_snapshot_index >= 0){

			if(greedy_params.verbose){
				std::cout << "||---------------------------------------------------------------------||" << std::endl;
				std::cout << "||-------------- Update Riesz Representors for A -------------------||" << std::endl;
				std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;

				std::cout << " Starting representor computation with snapshot nb " << last_snapshot_index+1 << std::endl << std::endl;
			}
		}
		else{
			if(greedy_params.verbose){
				std::cout << "||---------------------------------------------------------------------||" << std::endl;
				std::cout << "||-------------- Calculate Riesz Representors for A -------------------||" << std::endl;
				std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
			}
		}
	}
	else{
		if(greedy_params.verbose){
			std::cout << "||---------------------------------------------------------------------||" << std::endl;
			std::cout << "||-------------- Calculate Riesz Representors for A -------------------||" << std::endl;
			std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}
	}


	// Calculate the Riesz Representors for A
	std::size_t Qa = rb_system.Q_a();
	std::vector<DataType> 	new_A_reprs(Qa);
	std::vector<T>			new_A_eps(Qa);
	for (unsigned int i = 0; i < Qa; ++i) {

		if(greedy_params.verbose){
			std::cout << "||------- Component Nb " << std::setw(3) << i << "  -------------------------------------------||" << std::endl << std::endl;
		}

		DataType c;
		T eps;
		if(update || greedy_params.update_rieszA){
			// Start solver with copy of old Riesz representor
			c = A_representors[last_snapshot_index][i];
			std::cout << "Old Representor: " << c.size() << " indizes" << std::endl;
			eps = rb_truth.get_riesz_representor_a(i, bf, c, greedy_params.coarsen_rieszA_for_update);
			if(std::fabs(eps-1.) < 1e-8){
				eps = eps_A_representors[last_snapshot_index][i] / std::sqrt(greedy_params.riesz_constant_Y);
			}
		}
		else{
			eps = rb_truth.get_riesz_representor_a(i, bf, c);
		}
		new_A_reprs[i] = c;
		new_A_eps[i] = eps * std::sqrt(greedy_params.riesz_constant_Y);
	}

	if(update){
		A_representors[A_representors.size()-1] = new_A_reprs;
		eps_A_representors[eps_A_representors.size()-1] = new_A_eps;
	}
	else{
		A_representors.push_back(new_A_reprs);
		eps_A_representors.push_back(new_A_eps);
	}


	// Update the Riesz Representor Norms A x A
	int N = rb_basisfunctions.size();
	DataType vec1, vec2;
    for(int n1 = 0; n1 < N; ++n1) {
    	typename RB_Model::FullColMatrixT A_n1_N(Qa, Qa);
        for(std::size_t qa1 = 1; qa1 <= Qa; ++qa1) {
        	for(std::size_t qa2 = qa1; qa2 <= Qa; ++qa2) {
        		vec1 = A_representors[n1][qa1-1];
        		vec2 = A_representors[N-1][qa2-1];
        		A_n1_N(qa1, qa2) = rb_truth.innprod_Y_v_v(vec1, vec2);

        		if(qa1 != qa2){
        			if(n1 == N-1){
        				A_n1_N(qa2, qa1) = A_n1_N(qa1, qa2);
        			}
        			else{
        				vec1 = A_representors[n1][qa2-1];
        				vec2 = A_representors[N-1][qa1-1];
        				A_n1_N(qa2, qa1) = rb_truth.innprod_Y_v_v(vec1, vec2);
        			}
        		}
        	}
        }
        if(update){
            rb_system.A_A_representor_norms[n1][rb_system.A_A_representor_norms[n1].size()-1] = A_n1_N;
        }
        else{
            if(n1 == N-1){
            	std::vector<typename RB_Model::FullColMatrixT> newvec;
            	newvec.push_back(A_n1_N);
            	rb_system.A_A_representor_norms.push_back(newvec);
            }
            else{
            	rb_system.A_A_representor_norms[n1].push_back(A_n1_N);
            }
        }


    	if(greedy_params.verbose){
    		std::cout << std::endl << "||------- Representor Norms A["<< n1 <<"] x A["<< N-1 << "]  -----------------------------||" << std::endl;
    		std::cout << rb_system.A_A_representor_norms[n1][N-1-n1] << std::endl;
    	}
    }

	// Update the Riesz Representor Norms A x F
    std::size_t Qf = rb_system.Q_f();
    typename RB_Model::FullColMatrixT A_F(Qa, Qf);
    for(std::size_t qa = 1; qa <= Qa; ++qa) {
        for(std::size_t qf = 1; qf <= Qf; ++qf) {
            vec1 = A_representors[N-1][qa-1];
            vec2 = F_representors[qf-1];
            A_F(qa, qf) = rb_truth.innprod_Y_v_v(vec1, vec2);
        }
    }
    if(update){
        rb_system.A_F_representor_norms[rb_system.A_F_representor_norms.size()-1] = A_F;
    }
    else{
        rb_system.A_F_representor_norms.push_back(A_F);
    }
	if(greedy_params.verbose){
		std::cout << std::endl << "||------- Representor Norms A["<< N-1 <<"] x F  --------------------------------||" << std::endl;
		std::cout << rb_system.A_F_representor_norms[N-1] << std::endl;
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
recalculate_A_F_norms()
{

	rb_system.A_F_representor_norms.clear();

	for(std::size_t N = 1; N <= n_bf(); ++N){
		// Update the Riesz Representor Norms A x F
	    std::size_t Qa = rb_system.Q_a();
	    std::size_t Qf = rb_system.Q_f();
	    typename RB_Model::FullColMatrixT A_F(Qa, Qf);
	    for(std::size_t qa = 1; qa <= Qa; ++qa) {
	        for(std::size_t qf = 1; qf <= Qf; ++qf) {
	            DataType vec1 = A_representors[N-1][qa-1];
	            DataType vec2 = F_representors[qf-1];
	            A_F(qa, qf) = rb_truth.innprod_Y_v_v(vec1, vec2);
	        }
	    }
	    rb_system.A_F_representor_norms.push_back(A_F);
		if(greedy_params.verbose){
			std::cout << std::endl << "||------- Representor Norms A["<< N-1 <<"] x F  --------------------------------||" << std::endl;
			std::cout << rb_system.A_F_representor_norms[N-1] << std::endl;
		}
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
typename DataType::ValueType
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
find_max_errest(std::size_t N, std::vector<ParamType>& Xi_train, ParamType& current_param, std::map<ParamType, DataType>& truth_sols)
{

	T max_error;

	bool do_restart = true;
	while(do_restart){
		max_error = 0;
		T error_est = 0;

		do_restart = false;
		for(auto& mu : Xi_train){

			typename RB_Model::DenseVectorT u_N = rb_system.get_rb_solution(N, mu);

			switch(greedy_params.training_type){
			case weak:
			{
				error_est = rb_system.get_errorbound(u_N,mu);

				if(greedy_params.test_estimator_equivalence){
					//T eps_F = rb_truth.access_RieszSolver_F().access_params().tol * greedy_params.equivalence_tol_factor;
					//T eps_A = rb_truth.access_RieszSolver_A().access_params().tol * greedy_params.equivalence_tol_factor;

					T eps_aff = rb_system.get_errorbound_accuracy(u_N, mu, eps_F_representors, eps_A_representors);
					T bound = rb_system.alpha_LB(mu)*error_est*std::sqrt(greedy_params.riesz_constant_Y);

					if(eps_aff > bound && greedy_params.tighten_estimator_accuracy){
						std::cout << " !$!$! Equivalence Criterion violated: Eps_Aff = " << eps_aff << " > " << bound << "  !$!$!" << std::endl;
						std::cout << " !$!$! " << std::endl;
						// Can we do sth about it in this iteration?
						auto old_eps_a = eps_A_representors;
						std::vector<T> empty_eps(eps_F_representors.size());
						for(auto& el : old_eps_a[old_eps_a.size()-1]){
							el = 0;
						}

						T eps_zero = rb_system.get_errorbound_accuracy(u_N, mu, empty_eps, old_eps_a);
						if(eps_zero > bound){
							std::cerr << "!$!$!    No Chance: We have to recalculate old representors!! " << std::endl;
						}
						else{
							// Find accuracies that would be ok (hopefully)
							auto test_eps_f = eps_F_representors;
							auto test_eps_a = eps_A_representors;

							T old_eps_aff = eps_aff;
							// Reduce first the accuracies for a(N)
							int count = 0;
							int count_A = 0;
							int count_F = 0;

							// Reduce accuracies of A and F equally
							// (but only if it decreases eps_aff)
							// (Alternative: test reducing eps_a as far as possible, then eps_f)
							do{
								if(count_A == count){
									for(auto& el : test_eps_a[test_eps_a.size()-1]){
										el *= 0.5;
									}
									eps_aff = rb_system.get_errorbound_accuracy(u_N, mu, test_eps_f, test_eps_a);
									if(eps_aff < old_eps_aff) count_A++;
									std::cout << "!$!$!    Reducing eps_A by 0.5: Eps_Aff = " << eps_aff << std::endl;
								}
								if(count_A > count) old_eps_aff = eps_aff;

								if(eps_aff > bound && count_F == count){
									for(auto& el : test_eps_f){
										el *= 0.5;
									}
									eps_aff = rb_system.get_errorbound_accuracy(u_N, mu, test_eps_f, test_eps_a);
									if(eps_aff < old_eps_aff) count_F++;
									std::cout << "!$!$!    Reducing eps_F by 0.5: Eps_Aff = " << eps_aff << std::endl;
								}
								if(count_F > count) old_eps_aff = eps_aff;

								count++;

							}while(eps_aff > bound || (count > count_A && count > count_F));

							// Set new solver tolerances
							for(int i = 0; i < count_A; ++i){
								rb_truth.access_RieszSolver_A().access_params().tol *= 0.5;
							}
							for(int i = 0; i < count_F; ++i){
								rb_truth.access_RieszSolver_F().access_params().tol *= 0.5;
							}

							// Compute updates
							if(count_F > 0){
								calculate_Riesz_RHS_information(true);
								for(int i = 0; i < rb_system.Q_f(); ++i){
									greedy_info.repr_f_size[greedy_info.repr_f_size.size()-1][i] = F_representors[i].size();
								}
							}
							if(count_A > 0){
								calculate_Riesz_LHS_information(rb_basisfunctions[rb_basisfunctions.size()-1], true);
								for(int i = 0; i < rb_system.Q_a(); ++i){
									greedy_info.repr_a_size[greedy_info.repr_a_size.size()-1][i] = A_representors[A_representors.size()-1][i].size();
								}
							}

							greedy_info.eps_aff[greedy_info.eps_aff.size()-1].clear();
							greedy_info.eps_res_bound[greedy_info.eps_res_bound.size()-1].clear();

							do_restart = true;

							break;
						}
					}
					greedy_info.eps_aff[greedy_info.eps_aff.size()-1].push_back(eps_aff);
					greedy_info.eps_res_bound[greedy_info.eps_res_bound.size()-1].push_back(bound);
				}
				break;
			}
			case strong_adaptive:
			case strong:
			{
				DataType diff;
				if(N > 0){
					diff = reconstruct_u_N(u_N,N);
				}
				diff -= (*truth_sols.find(mu)).second;
				error_est = std::sqrt(rb_truth.innprod_Y_u_u(diff,diff));
				break;
			}
			case weak_direct:
				DataType res_repr;
				error_est = get_direct_errorbound(u_N,mu, res_repr);

				if(greedy_params.test_estimator_equivalence){

					// Initial Check if we have to update the representor
					T bound = rb_system.alpha_LB(mu)*error_est;
					T effective_eps = rb_truth.access_RieszSolver_Res().access_params().tol*std::sqrt(greedy_params.riesz_constant_Y);
					std::cout << "Eps = " << effective_eps << ", Bound: " << bound << std::endl;

					// Save old values, so they can be restored
					T initial_tol = rb_truth.access_RieszSolver_Res().access_params().tol;
					std::size_t max_its = rb_truth.access_RieszSolver_Res().access_params().max_its;
					// Set maximum iteration number to 1, so that we can check after each iteration
					rb_truth.access_RieszSolver_Res().access_params().max_its = 0;
					// Set tolerance small, so that we alway do 1 iteration  (including extension of index set)
					rb_truth.access_RieszSolver_Res().access_params().tol = 1e-14;

					DataType res_repr_old = res_repr;
					T res_repr_resnorm = 1.;
					T old_res_repr_resnorm;
					while(effective_eps > bound){

						// Save old representor
						res_repr_old = res_repr;
						old_res_repr_resnorm = res_repr_resnorm;

						// Get new riesz representor
						DataType u_approx = reconstruct_u_N(u_N, u_N.length());
						if(greedy_params.verbose){
							std::cout << "||-------------- Update Direct Riesz Representors for Residual -----||" << std::endl << std::endl;
						}

						// Attention: We return estimated residual of solution,
						// but in res_repr is the already extended index set
						// (which doesn't make any difference for the norm, as the extension coefficients are zero)!
						res_repr_resnorm = rb_truth.get_riesz_representor_res(u_approx, mu, res_repr, old_res_repr_resnorm);
						effective_eps = res_repr_resnorm*std::sqrt(greedy_params.riesz_constant_Y);
						error_est = res_repr.norm(2.)*std::sqrt(greedy_params.riesz_constant_Y) / rb_system.alpha_LB(mu);
						bound = res_repr.norm(2.)*std::sqrt(greedy_params.riesz_constant_Y);
						std::cout << "Eps = " << effective_eps << ", Bound: " << bound << std::endl;

					}

					greedy_info.eps_res_bound[greedy_info.eps_res_bound.size()-1].push_back(bound);
					greedy_info.eps_aff[greedy_info.eps_aff.size()-1].push_back(effective_eps);
					// Reset tolerance to initial value, so that it is not reduced for the next parameters
					rb_truth.access_RieszSolver_Res().access_params().tol = initial_tol;
					rb_truth.access_RieszSolver_Res().access_params().max_its = max_its;

					// Get res_repr without extension
					for(auto it = res_repr.begin(); it != res_repr.end();){
						if(res_repr_old.find((*it).first) == res_repr_old.end()){
							it = res_repr.erase(it);
						}
						else{
							++it;
						}
					}

				}

				greedy_info.repr_r_size[greedy_info.repr_r_size.size()-1].push_back(res_repr.size());

				if(greedy_params.write_direct_representors){

					if(greedy_params.verbose){
						std::cout << "=====>  Writing Direct Riesz Representor to file " << std::endl << std::endl;
					}
					mkdir(greedy_params.trainingdata_folder.c_str(), 0777);

					// Write Riesz Representors
					std::string repr_folder = greedy_params.trainingdata_folder + "/direct_representors";
					mkdir(repr_folder.c_str(), 0777);

					std::stringstream filename;
					filename << repr_folder << "/Res_representor_N_" <<  N << "_Mu";
					for(auto& el : mu){
						filename << "_"<< el;
					}
					filename  << ".txt";
					saveCoeffVector2D(res_repr, rb_truth.get_testbasis(), filename.str().c_str());
				}

			}

			if(do_restart) break;

			if(greedy_params.verbose){
				std::cout << "    u_N = " << u_N;
				std::cout << "    Mu = ";
				ParamInfo<ParamType>::print(mu);
				std::cout << " : Error = " << std::scientific << std::setw(15) << error_est << std::endl;
			}

			if(error_est > max_error){
				max_error = error_est;
				current_param = mu;
			}
		}
	}


	return max_error;
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
remove_basisfunction(std::size_t nb, bool from_greedy_info){

	assert(nb <= n_bf());


	if(greedy_params.orthonormalize_bfs == false){

		// Remove calN-data
		rb_basisfunctions.erase(rb_basisfunctions.begin()+(nb-1));
		A_representors.erase(A_representors.begin() + (nb-1));

		// Remove N-data
		rb_system.remove_basisfunction(nb);

		// Remove meta-data
		if(from_greedy_info){
			greedy_info.greedy_errors.erase(greedy_info.greedy_errors.begin()+nb-1);
			greedy_info.snapshot_params.erase(greedy_info.snapshot_params.begin()+nb-1);
			greedy_info.u_size.erase(greedy_info.u_size.begin()+nb-1);
			greedy_info.repr_a_size.erase(greedy_info.repr_a_size.begin()+nb-1);

			if(greedy_info.repr_r_size.size() >= nb){
				greedy_info.repr_r_size.erase(greedy_info.repr_r_size.begin()+nb-1);
			}
		}


	}
	else{
		// Re-Orthogonalization
		std::cerr << "$$$ !!!!!!!!!!!!!!!!!!!! $$$" << std::endl;
		std::cerr << "     Orthogonalized Snapshots! Didn't know what to do! Didn't erase anything" << std::endl;
		std::cerr << "$$$ !!!!!!!!!!!!!!!!!!!! $$$" << std::endl;
	}


}


template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
write_basisfunctions(const std::string& directory_name, int nb){

	if(rb_system.rb_params.verbose){
		std::cout << "=====>  Writing RB BasisFunctions to file " << std::endl << std::endl;
	}

	// Make a directory to store all the data files
	if(mkdir(directory_name.c_str(), 0777) == -1)
	{
		if(rb_system.rb_params.verbose){
			  std::cerr << "         [In RB_Base::write_basisfunctions: Directory "
					    << directory_name << " already exists, overwriting contents.]" << std::endl;
		}
	}

	std::string n_bf_filename = directory_name + "/n_bf.txt";
	std::ofstream n_bf_file(n_bf_filename.c_str());
	n_bf_file <<  rb_basisfunctions.size() << std::endl;
	n_bf_file.close();

	if(nb < 0){
		for(std::size_t i = 0; i < rb_basisfunctions.size(); ++i){
		    std::stringstream filename;
		    filename << directory_name << "/bf_" << i+1 << ".txt";
		    saveCoeffVector2D(rb_basisfunctions[i], rb_truth.get_trialbasis(), filename.str().c_str());
		}
	}
	else{
		assert((std::size_t)nb-1 < rb_basisfunctions.size());
	    std::stringstream filename;
	    filename << directory_name << "/bf_" << nb << ".txt";
	    saveCoeffVector2D(rb_basisfunctions[nb-1], rb_truth.get_trialbasis(), filename.str().c_str());
	}


}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
write_snapshot(DataType& u, const std::string& directory_name, int nb){

	if(rb_system.rb_params.verbose){
		std::cout << "=====>  Writing RB snapshots to file " << std::endl << std::endl;
	}

	// Make a directory to store all the data files
	if(mkdir(directory_name.c_str(), 0777) == -1)
	{
		if(rb_system.rb_params.verbose){
			  std::cerr << "         [In RB_Base::write_snapshot: Directory "
					    << directory_name << " already exists, overwriting contents.]" << std::endl;
		}
	}

	std::string n_bf_filename = directory_name + "/n_snap.txt";
	std::ofstream n_bf_file(n_bf_filename.c_str());
	n_bf_file <<  rb_basisfunctions.size() + 1 << std::endl;
	n_bf_file.close();

	if(nb < 0){
		for(std::size_t i = 0; i < rb_basisfunctions.size() + 1; ++i){
		    std::stringstream filename;
		    filename << directory_name << "/snap_" << i+1 << ".txt";
		    saveCoeffVector2D(rb_basisfunctions[i], rb_truth.get_trialbasis(), filename.str().c_str());
		}
	}
	else{
		assert((std::size_t)nb-1 < rb_basisfunctions.size()+1);
	    std::stringstream filename;
	    filename << directory_name << "/snap_" << nb << ".txt";
	    saveCoeffVector2D(u, rb_truth.get_trialbasis(), filename.str().c_str());
	}
}


template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
read_basisfunctions(const std::string& directory_name, int nb){

	unsigned int n_bf;
	if(nb < 0){
		std::string n_bf_filename = directory_name + "/n_bf.txt";

		std::ifstream n_bf_file(n_bf_filename.c_str());
		if(n_bf_file.is_open()){
			n_bf_file >> n_bf;
			n_bf_file.close();
		}
		else{
			std::cerr << "Unable to read number of basis functions: " << strerror(errno) << std::endl;
			exit(1);
		}
	}
	else{
		n_bf = nb;
	}


	rb_basisfunctions.clear();
	for(unsigned int i = 1; i <= n_bf; ++i){
		std::stringstream filename;
		filename << directory_name << "/bf_" << i << ".txt";

		DataType bf_coeffs;
		readCoeffVector2D(bf_coeffs, filename.str().c_str(),false);
		rb_basisfunctions.push_back(bf_coeffs);

		if(rb_system.rb_params.verbose){
			std::cout << " Read " << filename.str() << std::endl;
		}
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
write_rieszrepresentors(const std::string& directory_name, int nb)
{
	// Make a directory to store all the data files
	if(mkdir(directory_name.c_str(), 0777) == -1)
	{
		if(rb_system.rb_params.verbose){
			  std::cerr << "         [In RB_Base::write_rieszrepresentors: Directory "
					    << directory_name << " already exists, overwriting contents.]" << std::endl;
		}
	}

	if(nb <= 0){
		for(std::size_t i = 0; i < rb_system.Q_f(); ++i){
			std::stringstream filename;
			filename << directory_name << "/F_representor_" << i+1 << ".txt";
			saveCoeffVector2D(F_representors[i], rb_truth.get_testbasis(), filename.str().c_str());
		}
		if(nb < 0){
			for(std::size_t n = 0; n < rb_basisfunctions.size(); ++n){
				for(std::size_t i = 0; i < rb_system.Q_a(); ++i){
					std::stringstream filename;
					filename << directory_name << "/A_representor_" << i+1 << "_" << n+1 << ".txt";
					saveCoeffVector2D(A_representors[n][i], rb_truth.get_testbasis(), filename.str().c_str());
				}
			}
		}
	}
	else{
		assert((std::size_t)nb < rb_basisfunctions.size());

		for(std::size_t i = 0; i < rb_system.Q_a(); ++i){
			std::stringstream filename;
			filename << directory_name << "/A_representor_" << i+1 << "_" << nb+1 << ".txt";
			saveCoeffVector2D(A_representors[nb][i], rb_truth.get_testbasis(), filename.str().c_str());
		}
	}

}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
read_rieszrepresentors(const std::string& directory_name, int nb)
{

	F_representors.clear();
	for(std::size_t i = 0; i < rb_system.Q_f(); ++i){
		std::stringstream filename;
		filename << directory_name << "/F_representor_" << i+1 << ".txt";

		DataType rieszf_coeffs;
		readCoeffVector2D(rieszf_coeffs, filename.str().c_str(),false);
		F_representors.push_back(rieszf_coeffs);
        eps_F_representors.push_back(rb_truth.access_RieszSolver_F().access_params().tol);

		if(rb_system.rb_params.verbose){
			std::cout << " Read " << filename.str() << std::endl;
            std::cout << "      Set accuracy to " << *(eps_F_representors.end()-1) << std::endl;
		}
	}

	A_representors.clear();
	std::size_t n_bf = (nb < 0) ? rb_basisfunctions.size() : nb;

	for(std::size_t n = 0; n < n_bf; ++n){
		std::vector<DataType> A_reprs_n;
		std::vector<T> eps_A_n;
        
		for(std::size_t i = 0; i < rb_system.Q_a(); ++i){
			std::stringstream filename;
			filename << directory_name << "/A_representor_" << i+1 << "_" << n+1 << ".txt";

			DataType riesza_coeffs;
			readCoeffVector2D(riesza_coeffs, filename.str().c_str(),false);
			A_reprs_n.push_back(riesza_coeffs);
			
			eps_A_n.push_back(rb_truth.access_RieszSolver_A().access_params().tol);
            
			if(rb_system.rb_params.verbose){
				std::cout << " Read " << filename.str() << std::endl;
				std::cout << "      Set accuracy to " << *(eps_A_n.end()-1) << std::endl;
			}
		}
		A_representors.push_back(A_reprs_n);
		eps_A_representors.push_back(eps_A_n);
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
read_greedy_info(const char* filename, int nb){
	greedy_info.read(filename, rb_system.Q_f(), rb_system.Q_a(), nb);
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
read_repr_accuracies(const char* filename, int Nmax){
	greedy_info.read_repr_accuracies(filename, rb_system.Q_f(), rb_system.Q_a(), Nmax);
    eps_F_representors = *(greedy_info.accuracies_f_reprs.end()-1);
    eps_A_representors = *(greedy_info.accuracies_a_reprs.end()-1);
    if(rb_system.rb_params.verbose){
        std::cout << "     Set repr_F accuracies to: [";
        for(auto& a : eps_F_representors){
            std::cout << " " << a;
        }
        std::cout << "]" << std::endl;
        
        std::cout << "     Set repr_A accuracies to: ";
        for(auto& a : eps_A_representors){
            std::cout << " [";
            for(auto& el : a){
                std::cout << " " << el;
            }
            std::cout << "]" << std::endl;
        }
    }
}

} // namespace lawa
