#include <iostream>
#include <fstream>
#include <iomanip>

namespace lawa {

AWGM_PG_Parameters::AWGM_PG_Parameters(double _tol, double _tol_primal, double _alpha, std::size_t _max_its,
		std::size_t _max_basissize, bool _reset_res,
        StableExpansionVersion _stable_exp_u, StableExpansionVersion _stable_exp_res,
        ResidualConstruction _res_construction,
        bool _print_info, bool _verbose, bool _plot_solution, bool _verbose_extra,
		std::size_t _hashmapsize_trial, std::size_t _hashmapsize_test,
		std::string _info_filename,	std::string _plot_filename,
		bool _write_intermediary_solutions, std::string _intermediary_solutions_filename,
		bool _clear_solver)
: tol(_tol), tol_primal(_tol_primal),
  alpha(_alpha), max_its(_max_its), max_basissize(_max_basissize),
  reset_res(_reset_res), stable_exp_u(_stable_exp_u), stable_exp_res(_stable_exp_res),
  res_construction(_res_construction), print_info(_print_info),
  verbose(_verbose), plot_solution(_plot_solution), verbose_extra(_verbose_extra),
  hashmapsize_trial(_hashmapsize_trial), hashmapsize_test(_hashmapsize_test),
  info_filename(_info_filename), plot_filename(_plot_filename),
  write_intermediary_solutions(_write_intermediary_solutions),
  intermediary_solutions_filename(_intermediary_solutions_filename),
  clear_solver(_clear_solver)
  {}

void
AWGM_PG_Parameters::print()
{
	std::cout << "###### AWGM Parameters #################" << std::endl;
	std::cout << std::left << std::setw(24) << "# tol:" << std::setw(20) <<  tol << std::endl;
	std::cout << std::left << std::setw(24) << "# alpha:" << std::setw(20) <<  alpha << std::endl;
	std::cout << std::left << std::setw(24) << "# max_its:" << std::setw(20) <<  max_its << std::endl;
	std::cout << std::left << std::setw(24) << "# max_basissize:" << std::setw(20) <<  max_basissize << std::endl;
	std::cout << std::left << std::setw(24) << "# reset_res:" << std::setw(20) <<  (reset_res?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# stable_exp_u:" << std::setw(20) <<  (stable_exp_u == FullExpansion?"Full": (stable_exp_u == WoMixedHWExpansion? "WoMixedHW" : "OnlyTemporalHW")) << std::endl;
	std::cout << std::left << std::setw(24) << "# stable_exp_res:" << std::setw(20) <<  (stable_exp_res == FullExpansion?"Full": (stable_exp_res == WoMixedHWExpansion? "WoMixedHW" : "OnlyTemporalHW")) << std::endl;
	std::cout << std::left << std::setw(24) << "# res_construction:" << std::setw(20) <<  (res_construction == SimpleStableExpansion?"Simple":
																						  (res_construction == DoubleStableExpansion? "DoubleStable" :
																						  (res_construction == DoubleMTExtension ? "DoubleMTCone" : "Parallel"))) << std::endl;
	std::cout << std::left << std::setw(24) << "# print_info:" << std::setw(20) <<  (print_info?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# verbose:" << std::setw(20) <<  (verbose?"true":"false") << (verbose_extra?" (extra)":"") << std::endl;
	std::cout << std::left << std::setw(24) << "# plot_solution:" << std::setw(20) <<  (plot_solution?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# hashmapsize trial:" << std::setw(20) <<  hashmapsize_trial << std::endl;
	std::cout << std::left << std::setw(24) << "# hashmapsize test:" << std::setw(20) <<  hashmapsize_test << std::endl;
	std::cout << std::left << std::setw(24) << "# info_filename:" << std::setw(20) <<  info_filename << std::endl;
	std::cout << std::left << std::setw(24) << "# plot_filename:" << std::setw(20) <<  plot_filename << std::endl;
	std::cout << std::left << std::setw(24) << "# write intermed sols:" << std::setw(20) <<  (write_intermediary_solutions?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# intermed sols file:" << std::setw(20) <<  intermediary_solutions_filename << std::endl;
	std::cout << std::left << std::setw(24) << "# clear solver:" << std::setw(20) <<  clear_solver << std::endl;
	std::cout << "#########################################" << std::endl << std::endl;
}


IS_Parameters::IS_Parameters(bool _adaptive_tol, std::size_t _max_its, double _init_tol,
			  double _res_reduction, double _absolute_tol, bool _verbose)
: adaptive_tol(_adaptive_tol), max_its(_max_its), init_tol(_init_tol),
  res_reduction(_res_reduction), absolute_tol(_absolute_tol), tol(adaptive_tol?_init_tol:_absolute_tol),verbose(_verbose)
{}

void
IS_Parameters::print()
{
	std::cout << "###### Inner Solver Parameters ##########" << std::endl;
	std::cout << std::left << std::setw(18) << "# adaptive_tol:" << std::setw(16) <<  (adaptive_tol?"true":"false") << std::endl;
	std::cout << std::left << std::setw(18) << "# max_its:" << std::setw(16) <<  max_its << std::endl;
	std::cout << std::left << std::setw(18) << "# init_tol:" << std::setw(16) <<  init_tol << std::endl;
	std::cout << std::left << std::setw(18) << "# res_reduction:" << std::setw(16) <<  res_reduction << std::endl;
	std::cout << std::left << std::setw(18) << "# absolute_tol:" << std::setw(16) <<  absolute_tol << std::endl;
	std::cout << std::left << std::setw(18) << "# verbose:" << std::setw(16) <<  (verbose?"true":"false") << std::endl;
	std::cout << "#########################################" << std::endl << std::endl;
}

AWGM_Parameters::AWGM_Parameters(double _tol, double _tol_primal, double _alpha, std::size_t _max_its,
		std::size_t _max_basissize, bool _print_info,
		bool _verbose, bool _plot_solution, bool _verbose_extra,
		std::size_t _hashmapsize, std::string _info_filename, std::string _plot_filename,
		bool _clear_solver)
: tol(_tol), tol_primal(_tol_primal), alpha(_alpha), max_its(_max_its), max_basissize(_max_basissize),
  print_info(_print_info), verbose(_verbose), plot_solution(_plot_solution),
  verbose_extra(_verbose_extra), hashmapsize(_hashmapsize),
  info_filename(_info_filename), plot_filename(_plot_filename),
  clear_solver(_clear_solver)
{}

void
AWGM_Parameters::print()
{
	std::cout << "###### AWGM Parameters #################" << std::endl;
	std::cout << std::left << std::setw(24) << "# tol:" << std::setw(20) <<  tol << std::endl;
	std::cout << std::left << std::setw(24) << "# alpha:" << std::setw(20) <<  alpha << std::endl;
	std::cout << std::left << std::setw(24) << "# max_its:" << std::setw(20) <<  max_its << std::endl;
	std::cout << std::left << std::setw(24) << "# max_basissize:" << std::setw(20) <<  max_basissize << std::endl;
	std::cout << std::left << std::setw(24) << "# print_info:" << std::setw(20) <<  (print_info?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# verbose:" << std::setw(20) <<  (verbose?"true":"false") << (verbose_extra?" (extra)":"") << std::endl;
	std::cout << std::left << std::setw(24) << "# plot_solution:" << std::setw(20) <<  (plot_solution?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# hashmapsize:" << std::setw(20) <<  hashmapsize << std::endl;
	std::cout << std::left << std::setw(24) << "# info_filename:" << std::setw(20) <<  info_filename << std::endl;
	std::cout << std::left << std::setw(24) << "# plot_filename:" << std::setw(20) <<  plot_filename << std::endl;
	std::cout << std::left << std::setw(24) << "# clear solver:" << std::setw(20) <<  clear_solver << std::endl;
	std::cout << "#########################################" << std::endl << std::endl;
}

ISWGM_Parameters::ISWGM_Parameters(bool _print_info,
									bool _verbose,
									bool _plot_solution,
									std::string _info_filename,
									std::string _plot_filename)
 : print_info(_print_info), verbose(_verbose), plot_solution(_plot_solution),
   info_filename(_info_filename), plot_filename(_plot_filename)
{}

void
ISWGM_Parameters::print()
{
	std::cout << "###### ISWGM Parameters #################" << std::endl;
	std::cout << std::left << std::setw(24) << "# print_info:" << std::setw(20) <<  (print_info?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# verbose:" << std::setw(20) <<  (verbose?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# plot_solution:" << std::setw(20) <<  (plot_solution?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# info_filename:" << std::setw(20) <<  info_filename << std::endl;
	std::cout << std::left << std::setw(24) << "# plot_filename:" << std::setw(20) <<  plot_filename << std::endl;
	std::cout << "#########################################" << std::endl << std::endl;
}

void
AWGM_PG_Information::print(const char* filename)
{
    std::ofstream infofile(filename);
    if(infofile.is_open()){
    	infofile << "# It Res Res_NE SizeTrial SizeTest SizeTestResNE SizeTestRes CGLS_Its" << std::endl;
    	for(std::size_t i=0; i < awgm_res.size(); ++i){
    		infofile << i << " " << awgm_res[i] << " " << awgm_resNE[i] << " "
    				<< sizeLambdaTrial[i] << " " << sizeLambdaTest[i] << " "
    				<< sizeLambdaResNE[i] << " " << sizeLambdaRes[i] << " "
    				<< cgls_its[i] << std::endl;
    	}
        infofile.close();
    }
    else{
    	std::cerr << "Error opening file " << filename << " for writing! " << std::endl;
    }
}

void
AWGM_PG_Information::reset()
{
	awgm_res.clear();
	awgm_resNE.clear();
	sizeLambdaTrial.clear();
	sizeLambdaTest.clear();
	sizeLambdaResNE.clear();
	sizeLambdaRes.clear();
	cgls_its.clear();
}

void
AWGM_Information::print(const char* filename)
{
    std::ofstream infofile(filename);
    if(infofile.is_open()){
    	infofile << "# It Res SizeLambda SizeTestRes CG_Its" << std::endl;
    	for(std::size_t i=0; i < awgm_res.size(); ++i){
    		infofile << i << " " << awgm_res[i] << " " << sizeLambda[i] << " "
    		        << sizeLambdaRes[i] << " " << cg_its[i] << std::endl;
    	}
        infofile.close();
    }
    else{
    	std::cerr << "Error opening file " << filename << " for writing! " << std::endl;
    }
}

void
AWGM_Information::reset()
{
	awgm_res.clear();
	sizeLambda.clear();
	sizeLambdaRes.clear();
	cg_its.clear();
}

void
ISWGM_Information::print(const char* filename)
{
    std::ofstream infofile(filename);
    if(infofile.is_open()){
    	infofile << "# It Res" << std::endl;
    	for(std::size_t i=0; i < is_res.size(); ++i){
    		infofile << i << " " << is_res[i] << std::endl;
    	}
        infofile.close();
    }
    else{
    	std::cerr << "Error opening file " << filename << " for writing! " << std::endl;
    }
}

void
ISWGM_Information::reset()
{
	is_res.clear();
}

void
ISWGM_PG_Information::print(const char* filename)
{
    std::ofstream infofile(filename);
    if(infofile.is_open()){
    	infofile << "# It Res ResNE" << std::endl;
    	for(std::size_t i=0; i < is_res.size(); ++i){
    		infofile << i << " " << is_res[i] << " " << is_resNE[i] << std::endl;
    	}
        infofile.close();
    }
    else{
    	std::cerr << "Error opening file " << filename << " for writing! " << std::endl;
    }
}

void
ISWGM_PG_Information::reset()
{
	is_res.clear();
	is_resNE.clear();
}


std::ostream& operator<<(std::ostream& s,
                         const HTAWGM_Params& params)
{
    s << "uzero      = " << params.uzero      << std::endl;
    s << "tol_awgm   = " << params.tol_awgm   << std::endl;
    s << "gamma0     = " << params.gamma0     << std::endl;
    s << "gamma1     = " << params.gamma1     << std::endl;
    s << "gammait    = " << params.gammait    << std::endl;
    s << "maxit_awgm = " << params.maxit_awgm << std::endl;
    s << "maxit_pcg  = " << params.maxit_pcg  << std::endl;
    s << "delta1_pcg = " << params.delta1_pcg << std::endl;
    s << "delta2_pcg = " << params.delta2_pcg << std::endl;
    s << "delta3_pcg = " << params.delta3_pcg << std::endl;
    s << "alpha      = " << params.alpha      << std::endl;
    s << "nrmA       = " << params.nrmA      << std::endl;
    s << "omega      = " << params.omega      << std::endl;

    return s;
}


std::ostream& operator<<(std::ostream& s,
                         const HTRICH_Params& params)
{
    s << "uzero       = " << params.uzero       << std::endl;
    s << "nrmA        = " << params.nrmA        << std::endl;
    s << "cA          = " << params.cA          << std::endl;
    s << "omega       = " << params.omega       << std::endl;
    s << "beta1       = " << params.beta1       << std::endl;
    s << "beta2       = " << params.beta2       << std::endl;
    s << "kappa1      = " << params.kappa1      << std::endl;
    s << "kappa2      = " << params.kappa2      << std::endl;
    s << "kappa3      = " << params.kappa3      << std::endl;
    s << "eps0        = " << params.eps0        << std::endl;
    s << "rho         = " << params.rho         << std::endl;
    s << "tol_rich    = " << params.tol_rich    << std::endl;
    s << "trunc       = " << params.trunc       << std::endl;
    s << "maxit_rich  = " << params.maxit_rich  << std::endl;
    s << "maxit_inner = " << params.maxit_inner << std::endl;

    return s;
}


std::ostream& operator<<(std::ostream& s,
                         const Rank1UP_Params& params)
{
    s << "update     = " << params.update    << std::endl;
    s << "check_res  = " << params.check_res << std::endl;
    s << "orthog     = " << params.orthog    << std::endl;
    s << "sw         = " << params.sw        << std::endl;
    s << "balance    = " << params.balance   << std::endl;
    s << "tol_als    = " << params.tol_als   << std::endl;
    s << "tol_cg     = " << params.tol_cg    << std::endl;
    s << "max_sweep  = " << params.max_sweep << std::endl;
    s << "maxit_cg   = " << params.maxit_cg  << std::endl;
    s << "verbose    = " << params.verbose   << std::endl;

    return s;
}

std::ostream& operator<<(std::ostream& s,
                         const OptTTCoreParams& params)
{
    s << "maxit     = " << params.maxit << std::endl;
    s << "tol       = " << params.tol   << std::endl;
    s << "stag      = " << params.stag  << std::endl;

    return s;
}

std::ostream& operator<<(std::ostream& s,
                         const GreedyALSParams& params)
{
    s << "maxit     = " << params.maxit << std::endl;
    s << "tol       = " << params.tol   << std::endl;
    s << "stag      = " << params.stag  << std::endl;

    return s;
}

std::ostream& operator<<(std::ostream& s,
                         const AgALSParams& params)
{
    s << "tol       = " << params.tol       << std::endl;
    s << "maxit     = " << params.maxit     << std::endl;
    s << "gamma     = " << params.gamma     << std::endl;
    s << "bulk      = " << params.bulk      << std::endl;
    s << "rndinit   = " << params.rndinit   << std::endl;
    s << "r1update  = " << params.r1update  << std::endl;
    s << "coreopt   = " << params.coreopt   << std::endl;
    s << "greedyals = " << params.greedyals << std::endl;

    return s;
}

} // namespace lawa
