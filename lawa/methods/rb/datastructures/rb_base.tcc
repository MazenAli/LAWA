#include <algorithm>
#include <iomanip>
#include <sys/stat.h>
#include <tuple>
#include <cassert>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/algorithms/sample.h>
#include <chrono>

namespace lawa {

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType,
          typename _IndexSet>
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::RB_Base(RB_Model& _rb_system,
        TruthModel& _rb_truth, const _IndexSet& _Lambda)
 : rb_system(_rb_system), rb_truth(_rb_truth), Lambda(_Lambda)
{}

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
std::size_t
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
n_bf()
{
	return rb_basisfunctions.size();
}


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
          void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
set_tol(double tol_)
{
    tol = tol_;
}


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
          void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
set_mw(bool _IsMW)
{
    IsMW = _IsMW;
}


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
DataType&
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
get_bf(std::size_t i)
{
	assert(i < rb_basisfunctions.size());
	return rb_basisfunctions[i];
}

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
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

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
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
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double>                      elapsed;
	if(N == 0){
		// In order to be able to calculate empty error bounds,
		// we have to calculate the Riesz Representors for F
        start   = std::chrono::system_clock::now();
		calculate_Riesz_RHS_information();
        end     = std::chrono::system_clock::now();
        elapsed = end-start;
        std::cout << "Time for computing RHS information " << elapsed.count()
                  << "s\n";
	}

	//================================================//
	//      Greedy Iterations
	//================================================//

	ParamType current_param = Xi_train[0];
	T max_error = 0;
    std::vector<std::pair<ParamType, unsigned FLENS_DEFAULT_INDEXTYPE>> repeated_snaps;
	do {

		if(greedy_params.verbose){
		    std::cout << "||---------------------------------------------------------------------||" << std::endl;
		    std::cout << "||-------------- Greedy Search for new parameter  ---------------------||" << std::endl;
		    std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}

		if(greedy_params.training_type == weak_direct){
			std::vector<std::size_t> res_repr_sizes;
			greedy_info.repr_r_size.push_back(res_repr_sizes);
		}

		//------------------------------------------------//
		//  Find Parameter with maximum error (estimator)
		//------------------------------------------------//
		T error_est;

        start     = std::chrono::system_clock::now();
        max_error = find_max_errest(N, Xi_train, current_param, truth_sols);
        end       = std::chrono::system_clock::now();
        elapsed   = end-start;
        std::cout << "Time for computing next mu " << elapsed.count()
                  << "s\n";
        if (max_error <= greedy_params.tol) {
            std::cout << "Greedy tolerance reached\n";
            std::cout << "Greedy error     : " << max_error << std::endl;
            std::cout << "Greedy tolerance : " << greedy_params.tol << std::endl;
            break;
        }

		greedy_info.greedy_errors.push_back(max_error);
		greedy_info.snapshot_params.push_back(current_param);


		//------------------------------------------------//
		//      Compute next snapshot
		//------------------------------------------------//
        DataType u;
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

        if(greedy_params.training_type == strong){
            for(auto& entry : truth_sols){
                bool is_current_mu = true;
                for(std::size_t i = 0; i < current_param.size(); ++i){
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
            double delta    = greedy_params.deter_rate;
            double beta_mu  = rb_system.alpha_LB(current_param);
            double factor   = 0.;

            if (greedy_params.problem_type == PetrovGalerkin) {
                rb_truth.access_solver().access_params().tol_primal =
                greedy_params.tol*delta*beta_mu;
                factor = beta_mu*beta_mu*delta;
            } else {
                factor = beta_mu*delta;
            }

            rb_truth.access_solver().access_params().tol =
                greedy_params.tol*factor;
            if (greedy_params.verbose) {
                std::cout << "beta(mu) = " << rb_system.alpha_LB(current_param)
                          << std::endl;
                std::cout << "delta    = " << delta << std::endl;
                std::cout << "tol      = " << greedy_params.tol << std::endl;
                std::cout << "and thus...\n";
                std::cout << "Tolerance for snapshot set to\n"
                          << std::scientific
                          << rb_truth.access_solver().access_params().tol
                          << std::endl;

                if (greedy_params.problem_type == PetrovGalerkin) {
                    std::cout << "Primal : "
                              << std::scientific
                              << rb_truth.access_solver().access_params().tol_primal
                              << std::endl;
                }
            }

            start   = std::chrono::system_clock::now();
            u = rb_truth.get_truth_solution(current_param);
            end     = std::chrono::system_clock::now();
            elapsed = end-start;
            std::cout << "Time for computing snapshot " << elapsed.count()
                      << "s\n";
        }

		//------------------------------------------------//
		//      Update RB basis and infos
		//------------------------------------------------//

        std::string snap_folder = greedy_params.trainingdata_folder + "/snap";
        mkdir(greedy_params.trainingdata_folder.c_str(), 0777);
        write_snapshot(u, snap_folder, (FLENS_DEFAULT_INDEXTYPE)N+1);
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
			write_basisfunctions(bf_folder, (FLENS_DEFAULT_INDEXTYPE)N);

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


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
DataType
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
reconstruct_u_N(typename RB_Model::DenseVectorT u, std::size_t N)
{
	assert(N <= n_bf());
	assert(u.length() > 0);
	assert(u.length() >= (FLENS_DEFAULT_INDEXTYPE)N);

	DataType u_full;
	for (unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= N; ++i) {
		u_full +=  u(i) * rb_basisfunctions[i-1];
	}

	return u_full;
}

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
DataType
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
reconstruct_u_N(typename RB_Model::DenseVectorT u, std::vector<std::size_t> indices)
{
    size_t N = u.length();
	assert(N <= n_bf());
	assert(u.length() > 0);
	assert(u.length() >= (FLENS_DEFAULT_INDEXTYPE)N);

	DataType u_full;
	for (unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= N; ++i) {
		u_full +=  u(i) * rb_basisfunctions[indices[i-1]];
	}

	return u_full;
}

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
DataType
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
reconstruct_res_repr_N(typename RB_Model::DenseVectorT u, std::size_t N, ParamType& mu)
{
	assert(N <= n_bf());
	assert(u.length() >= (FLENS_DEFAULT_INDEXTYPE)N);


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


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
typename DataType::ValueType
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
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

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
typename DataType::ValueType
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
update_direct_errorbound(const typename RB_Model::DenseVectorT& u_N, ParamType& mu, DataType& res_repr)
{

	DataType u_approx = reconstruct_u_N(u_N, u_N.length());
	if(greedy_params.verbose){
		std::cout << "||-------------- Update Direct Riesz Representors for Residual -----||" << std::endl << std::endl;
	}
	rb_truth.get_riesz_representor_res(u_approx, mu, res_repr);
    return  res_repr.norm(2.)*std::sqrt(greedy_params.riesz_constant_Y) / rb_system.alpha_LB(mu);
}


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
std::vector<ParamType>
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
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

        for(FLENS_DEFAULT_INDEXTYPE i = pdim-1; ; --i){
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

    return Xi_train;
}

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
print_paramset(std::vector<ParamType> PSet)
{
	for(auto& el : PSet){
		ParamInfo<ParamType>::print(el);
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
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
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed;
    start   = std::chrono::system_clock::now();
    calculate_Riesz_LHS_information(new_bf);
    end     = std::chrono::system_clock::now();
    elapsed = end-start;
    std::cout << "Time for computing LHS information " << elapsed.count()
                  << "s\n";
}

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
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
			rb_system.RB_A_matrices[q_a].engine().resize((FLENS_DEFAULT_INDEXTYPE)1, (FLENS_DEFAULT_INDEXTYPE)1);
			n = 1;
		}
		else {
			typename RB_Model::FullColMatrixT tmp(rb_system.RB_A_matrices[q_a]);
			rb_system.RB_A_matrices[q_a].engine().resize(tmp.numRows()+1, tmp.numCols()+1);
			rb_system.RB_A_matrices[q_a](tmp.rows(), tmp.cols()) = tmp;

			n = rb_system.RB_A_matrices[q_a].numRows();
			for(unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= n; ++i) {
				rb_system.RB_A_matrices[q_a](i, n) = 0.;
				rb_system.RB_A_matrices[q_a](n, i) = 0.;
			}

		}

		// Compute new entries (one column, one row)
        for (unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= n; ++i) {
        	rb_system.RB_A_matrices[q_a](n,i) = rb_truth.lhs_u_u(q_a, bf, rb_basisfunctions[i-1]);
        }
        for (unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i < n; ++i) {
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
    		rb_system.RB_F_vectors[q_f].engine().resize((FLENS_DEFAULT_INDEXTYPE)1);
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
        rb_system.RB_inner_product.engine().resize((FLENS_DEFAULT_INDEXTYPE)1, (FLENS_DEFAULT_INDEXTYPE)1);
        n=1;
    }
    else {
    	typename RB_Model::FullColMatrixT tmp(rb_system.RB_inner_product);
    	rb_system.RB_inner_product.engine().resize(tmp.numRows()+1, tmp.numCols()+1);
    	rb_system. RB_inner_product(tmp.rows(), tmp.cols()) = tmp;
    	n = rb_system. RB_inner_product.numRows();
        for(unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= n; ++i) {
        	rb_system.RB_inner_product(i, n) = 0.;
        	rb_system.RB_inner_product(n, i) = 0.;
        }
    }

    for (unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= n; ++i) {
    	rb_system.RB_inner_product(n, i) = rb_truth.innprod_Y_u_u(bf, rb_basisfunctions[i-1]);
		rb_system.RB_inner_product(i, n) = rb_system.RB_inner_product(n,i);
    }

	if(greedy_params.verbose){
		std::cout << std::endl << "||------- RB_InnerProduct  ----------------------------------||" << std::endl;
		std::cout << rb_system.RB_inner_product << std::endl;
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType,
          typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
calculate_Riesz_RHS_information()
{
    calc_wav_RHS();

    std::size_t Qf = this->rb_system.Q_f();
    this->rb_system.F_F_representor_norms.engine().resize((FLENS_DEFAULT_INDEXTYPE) Qf,
                                                           (FLENS_DEFAULT_INDEXTYPE) Qf);
    for(std::size_t qf1 = 1; qf1 <= Qf; ++qf1) {
        for (std::size_t qf2 = qf1; qf2 <= Qf; ++qf2) {

            this->rb_system.F_F_representor_norms(qf1,qf2) =
                  dotProduct(this->F_representors[qf1-1],
                             this->F_representors[qf2-1]);

            if(qf1 != qf2) {
                this->rb_system.F_F_representor_norms(qf2,qf1) = this->rb_system.F_F_representor_norms(qf1,qf2);
            }
        }
    }

}


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model, TruthModel, DataType, ParamType, _IndexSet>::
calc_wav_RHS()
{
    this->F_representors.clear();

    auto& F      = this->rb_truth.access_RieszSolver_F().get_rhs();
    auto& P      = this->rb_truth.access_solver().get_testprec();
    auto& basis  = this->rb_truth.access_solver().get_testbasis();

    std::size_t Qf = this->rb_system.Q_f();
    assert(Qf>0);
    this->F_representors.resize(Qf);

    for (std::size_t i=0; i<Qf; ++i) {
        F.set_active_comp(i);
        sample_f(basis, Lambda, F, P,
                 this->F_representors[i],
                 tol, IsMW);
        this->total += this->F_representors[i].size();
    }

    if (greedy_params.verbose) {
        std::cout << "Current total size of F is " << total << std::endl;
    }
}


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
calculate_Riesz_LHS_information(const DataType& bf)
{
    FLENS_DEFAULT_INDEXTYPE N = this->rb_basisfunctions.size();
    if (greedy_params.verbose) {
        std::cout << "Current N=" << N << std::endl;
    }
    calc_wav_LHS(bf);

    // Update the Riesz Representor Norms A x A
    std::size_t Qa = this->rb_system.Q_a();
    for(FLENS_DEFAULT_INDEXTYPE n1 = 0; n1 < N; ++n1) {
        typename RB_Model::FullColMatrixT A_n1_N(Qa, Qa);
        for(std::size_t qa1 = 1; qa1 <= Qa; ++qa1) {
            for(std::size_t qa2 = qa1; qa2 <= Qa; ++qa2) {
                A_n1_N(qa1, qa2) =
                dotProduct(this->A_representors[n1][qa1-1],
                           this->A_representors[N-1][qa2-1]);

                if(qa1 != qa2){
                    if(n1 == N-1){
                        A_n1_N(qa2, qa1) = A_n1_N(qa1, qa2);
                    }
                    else{
                        A_n1_N(qa2, qa1) =
                        dotProduct(this->A_representors[n1][qa2-1],
                                   this->A_representors[N-1][qa1-1]);
                    }
                }
            }
        }
        if(n1 == N-1){
            std::vector<typename RB_Model::FullColMatrixT> newvec;
            newvec.push_back(A_n1_N);
            this->rb_system.A_A_representor_norms.push_back(newvec);
        }
        else{
            this->rb_system.A_A_representor_norms[n1].push_back(A_n1_N);
        }
    }

    // Update the Riesz Representor Norms A x F
    std::size_t Qf = this->rb_system.Q_f();
    typename RB_Model::FullColMatrixT A_F(Qa, Qf);
    for(std::size_t qa = 1; qa <= Qa; ++qa) {
        for(std::size_t qf = 1; qf <= Qf; ++qf) {
            A_F(qa, qf) =
            dotProduct(this->A_representors[N-1][qa-1],
                       this->F_representors[qf-1]);
        }
    }
    this->rb_system.A_F_representor_norms.push_back(A_F);
}


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
calc_wav_LHS(const DataType& bf)
{
    std::size_t Qa = this->rb_system.Q_a();
    assert(Qa>0);
    std::vector<DataType>   temp(Qa);
    auto& A     = this->rb_truth.access_solver().get_lhs().get_localops();
    auto& P     = this->rb_truth.access_solver().get_testprec();
    auto& basis = this->rb_truth.access_solver().get_testbasis();

    auto _Lambda = supp(bf);

    for (std::size_t i=0; i<Qa; ++i) {
        sample_Au(basis, _Lambda, (*A[i]), P, bf, temp[i],
                  tol, IsMW);
        // Because Kristina
        for (auto& lambda : temp[i]) {
            lambda.second *= -1.;
        }
        total += temp[i].size();
    }

    this->A_representors.push_back(temp);
    if (greedy_params.verbose) {
        std::cout << "Current total size of A and F is " << total << std::endl;
    }
}


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
typename DataType::ValueType
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
find_max_errest(std::size_t N, std::vector<ParamType>& Xi_train, ParamType& current_param, std::map<ParamType, DataType>& truth_sols)
{

    T max_error = 0;
    T error_est = 0;

    for(auto& mu : Xi_train){

        typename RB_Model::DenseVectorT u_N = rb_system.get_rb_solution(N, mu);

        switch(greedy_params.training_type){
        case weak:
        {
            error_est = rb_system.get_errorbound(u_N,mu);
            break;
        }
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
        {
            DataType res_repr;
            error_est = get_direct_errorbound(u_N,mu, res_repr);

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
            break;
        }
        case strong_adaptive:
        {
            /* Nothing, implementation artifact */
            break;
        }
        }

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

    return max_error;
}

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
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


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
write_basisfunctions(const std::string& directory_name, FLENS_DEFAULT_INDEXTYPE nb){

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

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
write_snapshot(DataType& u, const std::string& directory_name, FLENS_DEFAULT_INDEXTYPE nb){

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


template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
read_basisfunctions(const std::string& directory_name, FLENS_DEFAULT_INDEXTYPE nb){

	unsigned FLENS_DEFAULT_INDEXTYPE n_bf;
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
	for(unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= n_bf; ++i){
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

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
write_rieszrepresentors(const std::string& directory_name, FLENS_DEFAULT_INDEXTYPE nb)
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

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
read_rieszrepresentors(const std::string& directory_name, FLENS_DEFAULT_INDEXTYPE nb)
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

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
read_greedy_info(const char* filename, FLENS_DEFAULT_INDEXTYPE nb){
	greedy_info.read(filename, rb_system.Q_f(), rb_system.Q_a(), nb);
}

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
read_repr_accuracies(const char* filename, FLENS_DEFAULT_INDEXTYPE Nmax){
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


// Something still wrong with this procedure...
template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType, _IndexSet>::
calc_rb_data()
{
    // F_F_norms
    std::size_t Qf = this->rb_system.Q_f();
    this->rb_system.F_F_representor_norms.resize(Qf, Qf);
    for(std::size_t qf1 = 1; qf1 <= Qf; ++qf1) {
        for (std::size_t qf2 = qf1; qf2 <= Qf; ++qf2) {

            this->rb_system.F_F_representor_norms(qf1,qf2) =
                  dotProduct(this->F_representors[qf1-1],
                             this->F_representors[qf2-1]);

            if(qf1 != qf2) {
                this->rb_system.F_F_representor_norms(qf2,qf1) = this->rb_system.F_F_representor_norms(qf1,qf2);
            }
        }
    }

    // A_A_norms
    std::size_t Qa      = this->rb_system.Q_a();
    std::size_t n_bf    = this->n_bf();
    this->rb_system.A_A_representor_norms.resize(n_bf);
    for(std::size_t n1 = 0; n1 < n_bf; ++n1) {
        this->rb_system.A_A_representor_norms[n1].clear();

        for (std::size_t n2 = n1; n2 < n_bf; ++n2) {
            typename RB_Model::FullColMatrixT A_n1_n2(Qa, Qa);
            for(std::size_t qa1 = 1; qa1 <= Qa; ++qa1) {
                for(std::size_t qa2 = qa1; qa2 <= Qa; ++qa2) {
                    A_n1_n2(qa1, qa2) =
                    dotProduct(this->A_representors[n1][qa1-1],
                               this->A_representors[n2][qa2-1]);

                    if(qa1 != qa2){
                        if(n1 == n2){
                            A_n1_n2(qa2, qa1) = A_n1_n2(qa1, qa2);
                        }
                        else{
                            A_n1_n2(qa2, qa1) =
                            dotProduct(this->A_representors[n1][qa2-1],
                                       this->A_representors[n2][qa1-1]);
                        }
                    }
                }
            }

            this->rb_system.A_A_representor_norms[n1].push_back(A_n1_n2);
        }
    }


    // A_F_norms
    this->rb_system.A_F_representor_norms.resize(n_bf);
    for(std::size_t n = 0; n < n_bf; ++n) {
        typename RB_Model::FullColMatrixT A_F(Qa, Qf);
        for(std::size_t qa = 1; qa <= Qa; ++qa) {
            for(std::size_t qf = 1; qf <= Qf; ++qf) {
                A_F(qa, qf) =
                dotProduct(this->A_representors[n][qa-1],
                           this->F_representors[qf-1]);
            }
        }
        this->rb_system.A_F_representor_norms[n].resize(Qa, Qf);
        this->rb_system.A_F_representor_norms[n] = A_F;
    }

    // RB_A
    this->rb_system.RB_A_matrices.resize(Qa);
    for (std::size_t q_a = 0; q_a < Qa; ++q_a) {
        this->rb_system.RB_A_matrices[q_a].resize(n_bf, n_bf);
        for (std::size_t i = 1; i <= n_bf; ++i) {
            for (std::size_t j = 1; j <= n_bf; ++j) {
                this->rb_system.RB_A_matrices[q_a](j,i) = this->rb_truth.lhs_u_u(q_a, this->rb_basisfunctions[j-1],
                                                                                this->rb_basisfunctions[i-1]);
            }
        }
        std::cout << std::endl << "||------- RB_A (" << q_a << ")  -----------------------------------------||" << std::endl;
        std::cout << this->rb_system.RB_A_matrices[q_a] << std::endl;
    }

    // RB_F
    this->rb_system.RB_F_vectors.resize(Qf);
    for (std::size_t q_f = 0; q_f < Qf; ++q_f) {
        this->rb_system.RB_F_vectors[q_f].resize(n_bf);
        for (std::size_t n = 1; n <= n_bf; ++n) {
            this->rb_system.RB_F_vectors[q_f](n) = this->rb_truth.rhs_u(q_f, this->rb_basisfunctions[n-1]);
        }
        std::cout << std::endl << "||------- RB_F (" << q_f+1 << ")  -----------------------------------------||" << std::endl;
        std::cout << this->rb_system.RB_F_vectors[q_f] << std::endl;
    }

    // RB_InnerProduct
    this->rb_system.RB_inner_product.resize(n_bf, n_bf);
    for (unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= n_bf; ++i) {
        for (unsigned FLENS_DEFAULT_INDEXTYPE j = i; j <= n_bf; ++j) {
            this->rb_system.RB_inner_product(i, j) = this->rb_truth.innprod_Y_u_u(this->rb_basisfunctions[i-1],
                                                                            this->rb_basisfunctions[j-1]);
            this->rb_system.RB_inner_product(j, i) = this->rb_system.RB_inner_product(i,j);
        }
    }

    std::cout << std::endl << "||------- RB_InnerProduct  ----------------------------------||" << std::endl;
    std::cout << this->rb_system.RB_inner_product << std::endl;
}


} // namespace lawa
