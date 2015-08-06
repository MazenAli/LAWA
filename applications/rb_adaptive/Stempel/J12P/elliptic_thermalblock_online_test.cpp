#include "elliptic_thermalblock_problem_mw.h"
#include <fstream>
#include <sstream>

#define J       12
#define GAMMA   0.

void
read_paramfile(string paramfilename, vector<ParamType>& v){
    T mu1, mu2;
    ifstream paramfile(paramfilename);
    if(paramfile.is_open()){
    	while(!paramfile.eof()){
        	paramfile >> mu1 >> mu2;
        	ParamType mu = {{mu1, mu2}};
        	v.push_back(mu);
    	}
    	paramfile.close();
    }
    else{
    	cerr << "Couldn't open " << paramfilename << " for reading!" << endl;
    }
}

int main (int argc, char* argv[]) {

	//===============================================================//
	//========= PROBLEM SETUP  =======================//
	//===============================================================//
	
    int d   = 2;
    int d_  = 2;
    int j0  = 2;

    if(argc < 5 || argc > 7){
    	cerr << "Usage: " << argv[0] << " offline_data_folder output_name readSols(0/1) updateSols(0/1) (InputSolutionFolder) (OutputSolutionFolder)" << endl;
    	exit(1);
     }

    string offline_folder = argv[1];
    string output = argv[2];
    int readSols = atoi(argv[3]);
    int updateSols = atoi(argv[4]);
    string solfolder_in, solfolder_out;
    if(argc >= 6){
        solfolder_in = argv[5];
    }
    if(argc == 7){
        solfolder_out = argv[6];
    }

    /// Basis initialization
    Basis_X     basis_x(d,d_,j0);
    Basis_Y     basis_y(d,0);
    //Basis_Y     basis_y(d,d_,j0);
    //basis_x.enforceBoundaryCondition<DirichletBC>();
    basis_y.enforceBoundaryCondition<DirichletBC>();

    Basis2D basis2d(basis_x,basis_y);
    
    /// Initialization of operators
    DenseVectorT no_singPts;
    Function<T> zero_Fct(zero_fct,no_singPts);
    Function<T> chi_omega_0_Fct(chi_omega_0,no_singPts);
    Function<T> chi_omega_1_Fct(chi_omega_1,no_singPts);

    // Bilinear Forms
    Identity1D_X 		    IdentityBil_x(basis_x);
    Identity1D_Y 	        IdentityBil_y(basis_y);
    Laplace1D_X 	        LaplaceBil_x(basis_x);
    Laplace1D_Y 	        LaplaceBil_y(basis_y);
    WeightedLaplace1D_X     WLaplaceBil_x_0(basis_x, zero_Fct, zero_Fct, chi_omega_0_Fct, 10,true, true, false);
    WeightedLaplace1D_X     WLaplaceBil_x_1(basis_x, zero_Fct, zero_Fct, chi_omega_1_Fct, 10,true, true, false);
    WeightedIdentity1D_X    WIdentityBil_x_0(basis_x, chi_omega_0_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    WeightedIdentity1D_X    WIdentityBil_x_1(basis_x, chi_omega_1_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    
    RefIdentity1D_X 		RefIdentityBil_x(basis_x.refinementbasis);
    RefIdentity1D_Y 	    RefIdentityBil_y(basis_y.refinementbasis);
    RefLaplace1D_X 	        RefLaplaceBil_x(basis_x.refinementbasis);
    RefLaplace1D_Y 	        RefLaplaceBil_y(basis_y.refinementbasis);
    RefWeightedLaplace1D_X  RefWLaplaceBil_x_0(basis_x.refinementbasis, zero_Fct, zero_Fct, chi_omega_0_Fct, 10,true, true, false);
    RefWeightedLaplace1D_X  RefWLaplaceBil_x_1(basis_x.refinementbasis, zero_Fct, zero_Fct, chi_omega_1_Fct, 10,true, true, false);
    RefWeightedIdentity1D_X RefWIdentityBil_x_0(basis_x.refinementbasis, chi_omega_0_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    RefWeightedIdentity1D_X RefWIdentityBil_x_1(basis_x.refinementbasis, chi_omega_1_Fct, zero_Fct, zero_Fct, 10, false, true, true);

    /// Initialization of local operator
    LOp_Id1D_X              lOp_Id1D_x      (basis_x, basis_x, RefIdentityBil_x,   IdentityBil_x);
    LOp_Id1D_Y              lOp_Id1D_y      (basis_y, basis_y, RefIdentityBil_y,   IdentityBil_y);
    LOp_Lapl1D_X            lOp_Lapl1D_x    (basis_x, basis_x, RefLaplaceBil_x,    LaplaceBil_x);
    LOp_Lapl1D_Y            lOp_Lapl1D_y    (basis_y, basis_y, RefLaplaceBil_y,    LaplaceBil_y);
    LOp_WLapl1D_X           lOp_WLapl1D_x_0 (basis_x, basis_x, RefWLaplaceBil_x_0, WLaplaceBil_x_0);
    LOp_WLapl1D_X           lOp_WLapl1D_x_1 (basis_x, basis_x, RefWLaplaceBil_x_1, WLaplaceBil_x_1);
    LOp_WId1D_X             lOp_WId1D_x_0   (basis_x, basis_x, RefWIdentityBil_x_0, WIdentityBil_x_0);
    LOp_WId1D_X             lOp_WId1D_x_1   (basis_x, basis_x, RefWIdentityBil_x_1, WIdentityBil_x_1);


    IndexSet<Index2D> LambdaExact;
    getSparseGridIndexSet(basis2d, LambdaExact, J, 0, GAMMA);

    LOp_Id_Id_2D			localIdentityIdentityOp2D		(lOp_Id1D_x, 		lOp_Id1D_y);
    LOp_Id_Lapl_2D			localIdentityLaplaceOp2D		(lOp_Id1D_x, 		lOp_Lapl1D_y);
    LOp_Lapl_Id_2D			localLaplaceIdentityOp2D		(lOp_Lapl1D_x, 		lOp_Id1D_y);
    LOp_WLapl_Id_2D			localWeightLaplaceIdentityOp2D_0(lOp_WLapl1D_x_0, 	lOp_Id1D_y, LambdaExact);
    LOp_WLapl_Id_2D			localWeightLaplaceIdentityOp2D_1(lOp_WLapl1D_x_1, 	lOp_Id1D_y, LambdaExact);
    LOp_WId_Lapl_2D			localWeightIdentityLaplaceOp2D_0(lOp_WId1D_x_0, 	lOp_Lapl1D_y, LambdaExact);
    LOp_WId_Lapl_2D			localWeightIdentityLaplaceOp2D_1(lOp_WId1D_x_1, 	lOp_Lapl1D_y, LambdaExact);

    localIdentityIdentityOp2D.setJ(9);
    localIdentityLaplaceOp2D.setJ(9);
    localLaplaceIdentityOp2D.setJ(9);
    localWeightLaplaceIdentityOp2D_0.setJ(9);
    localWeightLaplaceIdentityOp2D_1.setJ(9);
    localWeightIdentityLaplaceOp2D_0.setJ(9);
    localWeightIdentityLaplaceOp2D_1.setJ(9);

    /// Initialization of preconditioner
    Prec2D prec(basis2d);
    NoPrec2D noPrec;

    /// Initialization of rhs

    /// Right Hand Side:
    ///     No Singular Supports in both dimensions
    DenseVectorT sing_support;
    ///      Forcing Functions
    SeparableFunction2D<T> F1_Fct(f_1_x, sing_support, f_1_y, sing_support);
    SeparableFunction2D<T> F2_Fct(f_1_x, sing_support, f_2_y, sing_support);
    SeparableFunction2D<T> F3_Fct(f_1_x, sing_support, f_3_y, sing_support);
    SeparableFunction2D<T> F4_Fct(f_2_x, sing_support, f_1_y, sing_support);
    SeparableFunction2D<T> F5_Fct(f_2_x, sing_support, f_2_y, sing_support);
    SeparableFunction2D<T> F6_Fct(f_2_x, sing_support, f_3_y, sing_support);
    SeparableFunction2D<T> F7_Fct(f_3_x, sing_support, f_1_y, sing_support);
    SeparableFunction2D<T> F8_Fct(f_3_x, sing_support, f_2_y, sing_support);
    SeparableFunction2D<T> F9_Fct(f_3_x, sing_support, f_3_y, sing_support);
    ///     Peaks: points and corresponding coefficients
    ///             (heights of jumps in derivatives)
    FullColMatrixT nodeltas;
    SeparableRhsIntegral2D			rhs_1(basis2d, F1_Fct, nodeltas, nodeltas, 20);
    SeparableRhsIntegral2D			rhs_2(basis2d, F2_Fct, nodeltas, nodeltas, 20);
    SeparableRhsIntegral2D			rhs_3(basis2d, F3_Fct, nodeltas, nodeltas, 20);
    SeparableRhsIntegral2D			rhs_4(basis2d, F4_Fct, nodeltas, nodeltas, 20);
    SeparableRhsIntegral2D			rhs_5(basis2d, F5_Fct, nodeltas, nodeltas, 20);
    SeparableRhsIntegral2D			rhs_6(basis2d, F6_Fct, nodeltas, nodeltas, 20);
    SeparableRhsIntegral2D			rhs_7(basis2d, F7_Fct, nodeltas, nodeltas, 20);
    SeparableRhsIntegral2D			rhs_8(basis2d, F8_Fct, nodeltas, nodeltas, 20);
    SeparableRhsIntegral2D			rhs_9(basis2d, F9_Fct, nodeltas, nodeltas, 20);
    SeparableRhs           			F1(rhs_1,noPrec, LambdaExact);
    SeparableRhs           			F2(rhs_2,noPrec, LambdaExact);
    SeparableRhs           			F3(rhs_3,noPrec, LambdaExact);
    SeparableRhs           			F4(rhs_4,noPrec, LambdaExact);
    SeparableRhs           			F5(rhs_5,noPrec, LambdaExact);
    SeparableRhs           			F6(rhs_6,noPrec, LambdaExact);
    SeparableRhs           			F7(rhs_7,noPrec, LambdaExact);
    SeparableRhs           			F8(rhs_8,noPrec, LambdaExact);
    SeparableRhs           			F9(rhs_9,noPrec, LambdaExact);


	//===============================================================//
	//===============  RB SETUP =====================================//
	//===============================================================//

    // Affine Decompositions:
    // 	Left Hand Side
    vector<ThetaStructure<ParamType>::ThetaFct>	lhs_theta_fcts;
    lhs_theta_fcts.push_back(no_theta);
    lhs_theta_fcts.push_back(no_theta);
    lhs_theta_fcts.push_back(theta_1);
    lhs_theta_fcts.push_back(theta_1);
    ThetaStructure<ParamType> lhs_theta(lhs_theta_fcts);

    vector<AbstractLocalOperator2D<T>* > lhs_ops;
    lhs_ops.push_back(&localWeightLaplaceIdentityOp2D_0);
    lhs_ops.push_back(&localWeightIdentityLaplaceOp2D_0);
    lhs_ops.push_back(&localWeightLaplaceIdentityOp2D_1);
    lhs_ops.push_back(&localWeightIdentityLaplaceOp2D_1);

    Affine_Op_2D affine_lhs(lhs_theta, lhs_ops);

    H1_InnProd_2D innprod(localIdentityIdentityOp2D, localLaplaceIdentityOp2D, localIdentityLaplaceOp2D);

    // Right Hand Side
    vector<ThetaStructure<ParamType>::ThetaFct> rhs_theta_fcts;
    rhs_theta_fcts.push_back(theta_chi_1);
    rhs_theta_fcts.push_back(theta_chi_2);
    rhs_theta_fcts.push_back(theta_chi_3);
    rhs_theta_fcts.push_back(theta_chi_4);
    rhs_theta_fcts.push_back(theta_chi_5);
    rhs_theta_fcts.push_back(theta_chi_6);
    rhs_theta_fcts.push_back(theta_chi_7);
    rhs_theta_fcts.push_back(theta_chi_8);
    rhs_theta_fcts.push_back(theta_chi_9);
    ThetaStructure<ParamType> rhs_theta(rhs_theta_fcts);
    vector<SeparableRhs*> rhs_fcts;
    rhs_fcts.push_back(&F1);
    rhs_fcts.push_back(&F2);
    rhs_fcts.push_back(&F3);
    rhs_fcts.push_back(&F4);
    rhs_fcts.push_back(&F5);
    rhs_fcts.push_back(&F6);
    rhs_fcts.push_back(&F7);
    rhs_fcts.push_back(&F8);
    rhs_fcts.push_back(&F9);

    Affine_Rhs_2D affine_rhs(rhs_theta, rhs_fcts);
    RieszF_Rhs_2D rieszF_rhs(rhs_fcts);
    RieszA_Rhs_2D rieszA_rhs(lhs_ops);


	//===============================================================//
	//===============  AWGM =========================================//
	//===============================================================//

    /* IS Parameters Default Values
	bool adaptive_tol = true;
	size_t max_its = 100;
	double init_tol = 0.001;
	double res_reduction = 0.01;
	double absolute_tol = 1e-8;
	bool verbose = true;
	*/

    /* AWGM Parameters Default Values
    double tol = 5e-03;
	double alpha = 0.7;
	size_t max_its = 100;
	size_t max_basissize = 400000;
	bool print_info = true;
	bool verbose = true;
	bool plot_solution = false;
	bool verbose_extra = false; //(print added wavelet indizes)
	size_t hashmapsize = 10;
	std::string info_filename = "awgm_cg_conv_info.txt";
	std::string plot_filename = "awgm_cg_u_plot";
	*/

    IS_Parameters is_parameters;
    AWGM_Parameters awgm_truth_parameters, awgm_riesz_f_parameters, awgm_riesz_a_parameters;
    awgm_truth_parameters.max_its = 1000;
    
    is_parameters.adaptive_tol = true;
    is_parameters.init_tol = 0.0001;
    is_parameters.res_reduction = 0.01;

    //----------- Solver ---------------- //

    T gamma = 0.2;
    IndexSet<Index2D> Lambda;
    getSparseGridIndexSet(basis2d,Lambda,10,0,gamma);

    MT_AWGM_Truth awgm_u(basis2d, affine_lhs, affine_rhs, prec, awgm_truth_parameters, is_parameters);
    awgm_u.set_sol(dummy);
    awgm_u.awgm_params.tol = 1e-08;
    awgm_u.awgm_params.max_basissize = 3e+06;
    awgm_u.set_initial_indexset(Lambda);

    MT_AWGM_Riesz_F awgm_rieszF(basis2d, innprod, rieszF_rhs, prec, awgm_riesz_f_parameters, is_parameters);
    awgm_rieszF.set_sol(dummy);
    awgm_rieszF.set_initial_indexset(Lambda);
    awgm_rieszF.awgm_params.tol = 5e-03;
    awgm_rieszF.awgm_params.info_filename = "stempel_rieszF_conv_info.txt";

    MT_AWGM_Riesz_A awgm_rieszA(basis2d, innprod, rieszA_rhs, prec, awgm_riesz_a_parameters, is_parameters);
    awgm_rieszA.set_sol(dummy);
    awgm_rieszA.set_initial_indexset(Lambda);
    awgm_rieszA.awgm_params.tol = 5e-03;
    awgm_rieszA.awgm_params.info_filename = "stempel_rieszA_conv_info.txt";

    MTTruthSolver rb_truth(awgm_u, awgm_rieszF, awgm_rieszA, innprod, affine_lhs, rieszF_rhs);

    //----------- RB System ---------------- //

    LB_Base<ParamType, MTTruthSolver> lb_base(rb_truth, lhs_theta, false);
    IndexSet<Index2D> Lambda_Alpha_sparse, Lambda_Alpha_full ;

/*
    //getSparseGridIndexSet(basis2d_trial,LambdaTrial_Alpha_sparse,3,0,gamma);
    getFullIndexSet(basis2d, Lambda_Alpha_full, 4,4,0);

    cout << "++ Assembling Matrices for Alpha Computation ... " << endl << endl;
    lb_base.assemble_matrices_for_alpha_computation(Lambda_Alpha_full);
    cout << "++ ... done " << endl;
*/
    RB_Model rb_system(lhs_theta, rhs_theta, lb_base);
    rb_system.read_alpha("alphas.txt");

    rb_system.rb_params.ref_param = {{1., 1.}};
    rb_system.rb_params.call = call_cg;

    //----------- RB Base ---------------- //

    RB_BaseModel rb_base(rb_system, rb_truth, LambdaExact);

    /* RB Greedy Parameters Default Values
		TrainingType _training_type = weak,
		double _tol = 1e-2,
		std::size_t _Nmax = 20,
		ParamType _min_param = ParamType(),
		ParamType _max_param = ParamType(),
		intArray  _training_params_per_dim = intArray(),
		intArray  _log_scaling = intArray(),
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
		bool _equivalence_tol_factor = 1.,
		bool _write_direct_representors = false,
		T _min_error_reduction = 0.5);
     */


    /* RB Parameters Default Values
      		SolverCall call = call_cg,
			ParamType ref_param = ParamType(),
			bool verbose = true
     */

    ParamType mu_min = {{0.01, 0}};
    ParamType mu_max = {{10, 1}};

    rb_base.greedy_params.tol = 1e-04;
    rb_base.greedy_params.min_param = mu_min;
    rb_base.greedy_params.max_param = mu_max;
    rb_base.greedy_params.Nmax = 	20;
    rb_base.greedy_params.nb_training_params = {{20, 9}};
    rb_base.greedy_params.log_scaling = {{1, 0}};
    rb_base.greedy_params.print_file = "stempel_greedy_info.txt";
    rb_base.greedy_params.trainingdata_folder = "training_data_stempel";
    rb_base.greedy_params.print_paramset = true;
    rb_base.greedy_params.erase_snapshot_params = false;
    rb_base.greedy_params.orthonormalize_bfs = false;
    rb_base.greedy_params.tighten_tol = false;
    rb_base.greedy_params.tighten_tol_rieszA = true;
    rb_base.greedy_params.tighten_tol_rieszF = true;
    rb_base.greedy_params.update_snapshot = true;
    rb_base.greedy_params.update_rieszF = true;
    rb_base.greedy_params.update_rieszA = false;
    rb_base.greedy_params.test_estimator_equivalence = true;
	
    cout << "Parameters Training: " << std::endl << std::endl;
    rb_base.greedy_params.print();
    rb_system.rb_params.print();

    cout << "Parameters Truth Solver: " << std::endl << std::endl;
    awgm_u.awgm_params.print();

    awgm_u.is_params.print();

    cout << "Parameters Riesz Solver F : " << std::endl << std::endl;
    awgm_rieszF.awgm_params.print();

    cout << "Parameters Riesz Solver A : " << std::endl << std::endl;
    awgm_rieszA.awgm_params.print();

    // Read RB Data
    rb_system.read_rb_data("training_data_stempel");
    //rb_base.read_basisfunctions("training_data_stempel/bf");
    //rb_base.read_rieszrepresentors("training_data_stempel/representors");

    // Read Test Parameters
    std::vector<ParamType> Xi_test;
    read_paramfile("Xitest.txt", Xi_test);


    cout << "===============================" << endl;
    cout << "Test Parameters: " << endl;
    for(auto& mu : Xi_test){
    	cout << mu[0] << " " << mu[1] << endl;
    }
    cout << "===============================" << endl;


    // Test Reduced Solutions
    std::map<ParamType, std::vector<T> > errs, errbounds;
    int Nmax = rb_system.RB_inner_product.numRows();

    DenseVectorT av_err(Nmax);
    DenseVectorT av_errest(Nmax);

	// For all parameters:
    for(auto& mu : Xi_test){
        
        cout << "Mu = [" << mu[0] << " " << mu[1] << "]" << endl;
        
    	stringstream ufilename_in, ufilename_out;
    	if(argc >= 6){
    		ufilename_in << solfolder_in << "/";
    	}
    	if(argc == 7){
            ufilename_out << solfolder_out << "/";
    	}
    	ufilename_in << "u_" << mu[0] << "_" << mu[1] << ".txt";
    	ufilename_out << "u_" << mu[0] << "_" << mu[1] << ".txt";

    	/*DataType u;
    	if(!readSols){
        	// Compute truth solution
        	u = rb_truth.get_truth_solution(mu);
        	saveCoeffVector2D(u, basis2d, ufilename_out.str().c_str());
    	}
    	else{
    		readCoeffVector2D(u, ufilename_in.str().c_str());
    		
    		if(updateSols){
    		    rb_truth.get_truth_solution(mu, u);
    		    saveCoeffVector2D(u, basis2d, ufilename_out.str().c_str());
    		}
    	}*/

    	std::vector<T> errs_mu, errbounds_mu;

    	T err, errbound;
		for(size_t n = 1; n <= Nmax; ++n){

    		// Compute RB solution
    		DenseVectorT u_N = rb_system.get_rb_solution(n, mu);
            cout << "u_" << n << " = " << u_N << endl;
    		// Compute real error
    		//DataType u_approx = rb_base.reconstruct_u_N(u_N, n);
    		//u_approx -= u;
    		err = 0.;//u_approx.norm(2.);
            cout << "Err = " << err << endl;
    		errs_mu.push_back(err);
    		
    		av_err(n) += err;
    		
    		// Compute error estimator
    		errbound = rb_system.get_errorbound(u_N, mu);
    		errbounds_mu.push_back(errbound);
    		av_errest(n) += errbound;
    		
		}

    	errs.insert(make_pair(mu, errs_mu));
    	errbounds.insert(make_pair(mu, errbounds_mu));
    }

    av_err *= 1./Xi_test.size();
	av_errest *= 1./Xi_test.size();

    ofstream errfile(output + "_OnlineTest_Errors.txt");
    if(errfile.is_open()){
    	errfile << "Mu1 Mu2 Err(N=1) Err(N=2) ...." << endl;
    	for(auto& entry : errs){
    		errfile << entry.first[0] << " " << entry.first[1];
    		for(auto& e : entry.second){
    			errfile << " " << e;
    		}
    		errfile << endl;
    	}
    	errfile.close();
    }
    else{
    	cerr << "Couldn't open Test_Errors.txt for writing!" << endl;
    }

    ofstream errestfile(output + "_OnlineTest_ErrorEstimators.txt");
    if(errestfile.is_open()){
    	errestfile << "Mu1 Mu2 ErrEst(N=1) ErrEst(N=2) ...." << endl;
    	for(auto& entry : errbounds){
    		errestfile << entry.first[0] << " " << entry.first[1];
    		for(auto& e : entry.second){
    			errestfile << " " << e;
    		}
    		errestfile << endl;
    	}
    	errestfile.close();
    }
    else{
    	cerr << "Couldn't open Test_Errors.txt for writing!" << endl;
    }

    
    ofstream av_file(output + "_OnlineTest_Averages.txt");
    if(av_file.is_open()){
        av_file << "# N avErr avErrEst" << endl;
        for(int i = 1; i <= Nmax; ++i){
            av_file << i << " " << i << " " << av_err(i) << " " << av_errest(i) << endl;
        }
        
        av_file.close();        
    }
    else{
        cerr << "Couldn't open file " << output << "_OnlineTest_Averages.txt" << " for writing! " << endl;
    }

    return 0;
}

