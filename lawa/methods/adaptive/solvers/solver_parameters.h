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

#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_SOLVER_PARAMETERS_H_
#define LAWA_METHODS_ADAPTIVE_SOLVERS_SOLVER_PARAMETERS_H_

#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>

namespace lawa {

enum StableExpansionVersion {
	FullExpansion,				// include all combinations of higher wavelet neighbours (i.e. (lambda_t + 1, lambda_x + 1)
	WoMixedHWExpansion,			// include only higher wavelets in one coordinate direction (i.e. (lambda_t + 1, lambda_x), (lambda_t, lambda_x + 1)
	OnlyTemporalHWExpansion		// include only higher wavelets in time
};

enum ResidualConstruction {
	SimpleStableExpansion,		// LambdaHat --(MT-Cone)--> XiHat --(StableExp)--> XiCheck
	DoubleStableExpansion,		// LambdaHat --(MT-Cone)--(StableExp)--> XiCheck --(StableExp)--> XiHat
	ParallelMTCones,			// LambdaHat --(MT-Cone)--> XiHat, LambdaCheck --(MT-Cone)--> XiCheck
	DoubleMTExtension,			// LambdaHat --(MT-Cone)--> XiHat_tmp --(StableExp)--> XiCheck, XiHat_tmp --(MT-Cone)--> XiHat
};

/**
 * Parameters for the Adaptive Wavelet-Galerkin Method
 *  (= the outer solver), Petrov-Galerkin version
 */
struct AWGM_PG_Parameters{

    double        			tol;
    double                  tol_primal;
    double        			alpha;
    std::size_t    			max_its;
    std::size_t    			max_basissize;
    bool           			reset_res;
    StableExpansionVersion 	stable_exp_u;
    StableExpansionVersion 	stable_exp_res;
    ResidualConstruction	res_construction;


    bool        	print_info;
    bool         	verbose;
    bool        	plot_solution;
    bool        	verbose_extra;
    std::size_t     hashmapsize_trial;
    std::size_t     hashmapsize_test;
    std::string 	info_filename;
    std::string 	plot_filename;
    bool			write_intermediary_solutions;
    std::string		intermediary_solutions_filename;
    bool            clear_solver;

    AWGM_PG_Parameters(double _tol = 5e-03,
                    double _tol_primal = 0.,
                    double _alpha = 0.95,
                    std::size_t _max_its = 100,
                    std::size_t _max_basissize = 400000,
                    bool _reset_res = false,
                    StableExpansionVersion _stable_exp_u = OnlyTemporalHWExpansion,
                    StableExpansionVersion _stable_exp_res = OnlyTemporalHWExpansion,
                    ResidualConstruction _res_construction = SimpleStableExpansion,
                    bool _print_info = true,
                    bool _verbose = true,
                    bool _plot_solution = false,
                    bool _verbose_extra = false,
                    std::size_t _hashmapsize_trial = 10,
                    std::size_t _hashmapsize_test = 10,
                    std::string _info_filename = "awgm_cgls_conv_info.txt",
                    std::string _plot_filename = "awgm_cgls_u_plot",
                    bool _write_intermediary_solutions = false,
                    std::string _intermediary_solutions_filename = "awgm_cgls_u",
                    bool _clear_solver = false);

    void print();
};

/**
 * Parameters for the Adaptive Wavelet-Galerkin Method
 *  (= the outer solver), Galerkin version
 */
struct AWGM_Parameters{

    double     tol;
    double     tol_primal; // dummy variable
    double     alpha;
    std::size_t     max_its;
    std::size_t     max_basissize;

    bool    print_info;
    bool     verbose;
    bool    plot_solution;
    bool    verbose_extra;
    std::size_t     hashmapsize;
    std::string info_filename;
    std::string plot_filename;
    bool            clear_solver;


    AWGM_Parameters(double _tol = 5e-03,
                    double _tol_primal = 0.,
                    double _alpha = 0.95,
                    std::size_t _max_its = 100,
                    std::size_t _max_basissize = 400000,
                    bool _print_info = true,
                    bool _verbose = true,
                    bool _plot_solution = false,
                    bool _verbose_extra = false,
                    std::size_t _hashmapsize = 10,
                    std::string _info_filename = "awgm_cg_conv_info.txt",
                    std::string _plot_filename = "awgm_cg_u_plot",
                    bool _clear_solver = false);

    void print();
};

/**
 * Parameters for the inner solver (cg/cgls)
 */
struct IS_Parameters{
    bool            adaptive_tol;
    std::size_t     max_its;
    double          init_tol;
    double          res_reduction;
    double          absolute_tol;
    
    double          tol;

    bool            verbose;

    IS_Parameters(bool _adaptive_tol = true,
                  std::size_t _max_its = 100,
                  double _init_tol = 0.001,
                  double _res_reduction = 0.01,
                  double _absolute_tol = 1e-8,
                  bool _verbose = true);

    void print();
};

/**
 * Parameters for a indexset-based Wavelet-Galerkin Method
 *  (= the outer solver)
 */
struct ISWGM_Parameters{

    bool    		print_info;
    bool    		verbose;
    bool    		plot_solution;
    std::string 	info_filename;
    std::string 	plot_filename;
    
    double          tol;
    std::size_t     max_its;

    ISWGM_Parameters(bool _print_info = true,
                    bool _verbose = true,
                    bool _plot_solution = false,
                    std::string _info_filename = "iswgm_conv_info.txt",
                    std::string _plot_filename = "iswgm_u_plot");

    void print();
};

/**
 * Gathers information that is interesting during a awgm-solve
 * and which can later be printed out (Petrov-Galerkin version)
 */
struct AWGM_PG_Information{
    std::vector<double>         awgm_res, awgm_resNE;
    std::vector<std::size_t>    sizeLambdaTrial, sizeLambdaTest,
                                sizeLambdaResNE, sizeLambdaRes,
                                cgls_its;

    void print(const char* filename = "awgm_cgls_conv_info.txt");

    void reset();
};

/**
 * Gathers information that is interesting during a awgm-solve
 * and which can later be printed out (alerkin version)
 */
struct AWGM_Information{
    std::vector<double>         awgm_res;
    std::vector<std::size_t>    sizeLambda, sizeLambdaRes,
                                cg_its;

    void print(const char* filename = "awgm_cg_conv_info.txt");

    void reset();
};

/**
 * Gathers information that is interesting during a iswgm-solve
 * and which can later be printed out (Galerkin version)
 */
struct ISWGM_Information{
    std::vector<double>         is_res;

    void print(const char* filename = "iswgm_is_conv_info.txt");

    void reset();
};

/**
 * Gathers information that is interesting during a iswgm-solve
 * and which can later be printed out (Petrov-Galerkin version)
 */
struct ISWGM_PG_Information{
    std::vector<double>         is_res, is_resNE;

    void print(const char* filename = "iswgm_pg_is_conv_info.txt");

    void reset();
};


/* HTAWGM solver parameters */
struct HTAWGM_Params
{
    bool     uzero      = true;  /* u_0=0? */
    double   tol_awgm   = 1e-08; /* AWGM solve accuracy */
    double   gamma0     = 1e-02; /* galerkin_pcg adaptive accuracy */
    double   gamma1     = 1e-01; /* galerkin_pcg adaptive accuracy */
    unsigned gammait    = 10;    /* galerkin_pcg adaptive accuracy */
    unsigned maxit_awgm = 5e+01; /* Maxit for AWGM */
    unsigned maxit_pcg  = 5e+01; /* Maxit for galerkin_pcg */
    double   delta1_pcg = 1e-01; /* Adaptive trunc tolerance search dir */
    double   delta2_pcg = 1e-01; /* Adaptive trunc tolerance search dir*/
    double   delta3_pcg = 1e-01; /* Adaptive trunc tolerance solution */
    double   dres_pcg   = 2.;    /* Scaling adjustment for res trunc */
    double   alpha      = 0.5;  /* Bulk chasing parameter */
    double   nrmA       = 50;    /* Estimate for operator norm */
    double   omega      = 0.1;   /* Estimate for relative res eval precision */
};

std::ostream& operator<<(std::ostream& s,
                         const HTAWGM_Params& params);

/* HTRICH solver parameters */
struct HTRICH_Params
{
    bool     uzero       = true;
    double   nrmA        = 50;
    double   cA          = 10;
    double   omega       = 0.25;
    double   beta1       = 0.5;
    double   beta2       = 0.5;
    double   kappa1      = 0.25;
    double   kappa2      = 0.25;
    double   kappa3      = 0.25;
    double   eps0        = 100;
    double   rho         = 0.5;
    double   tol_rich    = 1e-06;
    double   trunc       = 1e-01;
    unsigned maxit_rich  = 1e+02;
    unsigned maxit_inner = 1e+02;
};

std::ostream& operator<<(std::ostream& s,
                         const HTRICH_Params& params);

/* Rank 1 symmetric update parameters */
struct Rank1UP_Params
{
    bool     update     = false;
    bool     check_res  = false;
    bool     orthog     = false;
    bool     sw         = true;
    double   balance    = 50.;
    double   tol_als    = 1e-02;
    double   tol_cg     = 1e-08;
    unsigned max_sweep  = 1e+02;
    unsigned maxit_cg   = 3e+02;
    bool     verbose    = false;

};

std::ostream& operator<<(std::ostream& s,
                         const Rank1UP_Params& params);

struct OptTTCoreParams
{
    unsigned maxit   = 10;
    double   tol     = 1e-08;
    double   stag    = 1e-06;
};

std::ostream& operator<<(std::ostream& s,
                         const OptTTCoreParams& params);

struct GreedyALSParams
{
    unsigned maxit   = 10;
    double   tol     = 1e-04;
    double   stag    = 1e-04;
};

std::ostream& operator<<(std::ostream& s,
                         const GreedyALSParams& params);

struct AgALSParams
{
    double          tol     = 1e-08;
    unsigned        maxit   = 3e+01;
    double          gamma   = 1e-01;
    double          bulk    = .9;
    double          rndinit = 1.;
    Rank1UP_Params  r1update;
    OptTTCoreParams coreopt;
    GreedyALSParams greedyals;
};

std::ostream& operator<<(std::ostream& s,
                         const AgALSParams& params);

struct AdaptiveLeafParams
{
    double   tol          = 1e-06;
    unsigned maxit        = 100;
    double   gamma        = 0.1;
    unsigned cg_maxit     = 100;
    bool     cg_verbose   = false;
    bool     verbose      = false;
    double   alpha        = 0.5;
    bool     bulk_verbose = false;
};

struct Rank1AdaptiveAlsParams
{
    unsigned           max_sweep    = 10;
    double             stag         = 1e-02;
    bool               verbose      = false;
    AdaptiveLeafParams adaptiveLeaf;
};

struct AdaptiveGreedyParams
{
    bool                   verbose = false;
    double                 tol     = 1e-06;
    unsigned               maxit   = 30;
    Rank1AdaptiveAlsParams r1Als;
    OptTTCoreParams        dmrgTTcore;
};

} // namespace lawa

#include <lawa/methods/adaptive/solvers/solver_parameters.tcc>

#endif /* LAWA_METHODS_ADAPTIVE_SOLVERS_SOLVER_PARAMETERS_H_ */
