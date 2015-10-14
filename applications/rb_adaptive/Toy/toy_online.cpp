#include "toy.h"
#include <fstream>
#include <sstream>


struct precon
{
    T
    operator() (const Index2D&)
    {
        return 1.;
    }
};

T
zero(T, T)
{
    return 0.;
}


void
read_paramfile(string paramfilename, vector<ParamType>& v){
    T mu1;
    ifstream paramfile(paramfilename);
    if(paramfile.is_open()){
    	while(!paramfile.eof()){
        	paramfile >> mu1;
        	ParamType mu = {{mu1}};
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
	
    int d   = 4;
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
    _Basis     basis(d,j0);
    Basis2D basis2d(basis,basis);
    
    /// Initialization of operators
    DenseVectorT no_singPts;
    Function<T> zero_Fct(zero_fct,no_singPts);

    // Bilinear Forms
    Identity1D  		    IdentityBil(basis);
    Laplace1D   	        LaplaceBil(basis);
    RefIdentity1D    		RefIdentityBil(basis.refinementbasis);
    RefLaplace1D 	        RefLaplaceBil(basis.refinementbasis);

    /// Initialization of local operator
    LOp_Id1D              lOp_Id1D      (basis, basis, RefIdentityBil,   IdentityBil);
    LOp_Lapl1D            lOp_Lapl1D    (basis, basis, RefLaplaceBil,    LaplaceBil);

    LOp_Id_Id_2D			localIdentityIdentityOp2D		(lOp_Id1D, 		lOp_Id1D);
    LOp_Id_Lapl_2D			localIdentityLaplaceOp2D		(lOp_Id1D, 		lOp_Lapl1D);
    LOp_Lapl_Id_2D			localLaplaceIdentityOp2D		(lOp_Lapl1D, 		lOp_Id1D);

    localIdentityIdentityOp2D.setJ(9);
    localIdentityLaplaceOp2D.setJ(9);
    localLaplaceIdentityOp2D.setJ(9);

    /// Initialization of preconditioner
    Prec2D prec(basis2d);
    NoPrec2D noPrec;

    /// Initialization of rhs

    /// Right Hand Side:
    ///     No Singular Supports in both dimensions
    DenseVectorT sing_support;
    ///      Forcing Functions
    Function2D<T> F1_Fct(f1, sing_support, sing_support);
    Function2D<T> F2_Fct(f2, sing_support, sing_support);
    Function2D<T> F3_Fct(f3, sing_support, sing_support);
    Function2D<T> F4_Fct(f4, sing_support, sing_support);
    Function2D<T> F5_Fct(f5, sing_support, sing_support);
    Function2D<T> F6_Fct(f6, sing_support, sing_support);
    Function2D<T> F7_Fct(f7, sing_support, sing_support);
    Function2D<T> F8_Fct(f8, sing_support, sing_support);
    Function2D<T> F9_Fct(f9, sing_support, sing_support);
    ///     Peaks: points and corresponding coefficients
    ///             (heights of jumps in derivatives)
    FullColMatrixT nodeltas;
    RhsIntegral2D			rhs_1(basis2d, F1_Fct, 100);
    RhsIntegral2D			rhs_2(basis2d, F2_Fct, 100);
    RhsIntegral2D			rhs_3(basis2d, F3_Fct, 100);
    RhsIntegral2D			rhs_4(basis2d, F4_Fct, 100);
    RhsIntegral2D			rhs_5(basis2d, F5_Fct, 100);
    RhsIntegral2D			rhs_6(basis2d, F6_Fct, 100);
    RhsIntegral2D			rhs_7(basis2d, F7_Fct, 100);
    RhsIntegral2D			rhs_8(basis2d, F8_Fct, 100);
    RhsIntegral2D			rhs_9(basis2d, F9_Fct, 100);
    Rhs           			F1(rhs_1,noPrec);
    Rhs           			F2(rhs_2,noPrec);
    Rhs           			F3(rhs_3,noPrec);
    Rhs           			F4(rhs_4,noPrec);
    Rhs           			F5(rhs_5,noPrec);
    Rhs           			F6(rhs_6,noPrec);
    Rhs           			F7(rhs_7,noPrec);
    Rhs           			F8(rhs_8,noPrec);
    Rhs           			F9(rhs_9,noPrec);


	//===============================================================//
	//===============  RB SETUP =====================================//
	//===============================================================//

    // Affine Decompositions:
    // 	Left Hand Side
    vector<ThetaStructure<ParamType>::ThetaFct>	lhs_theta_fcts;
    lhs_theta_fcts.push_back(no_theta);
    ThetaStructure<ParamType> lhs_theta(lhs_theta_fcts);

    vector<AbstractLocalOperator2D<T>* > lhs_ops;
    lhs_ops.push_back(&localLaplaceIdentityOp2D);
    lhs_ops.push_back(&localIdentityLaplaceOp2D);
    Flex_COp_2D                       A(lhs_ops);
    vector<Flex_COp_2D*>                 ops_vec;
    ops_vec.push_back(&A);

    Affine_Op_2D affine_lhs(lhs_theta, ops_vec);

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
    vector<Rhs*> rhs_fcts;
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


    IS_Parameters is_parameters;
    AWGM_Parameters awgm_truth_parameters, awgm_riesz_f_parameters, awgm_riesz_a_parameters;
    awgm_truth_parameters.max_its = 1000;
    
    is_parameters.adaptive_tol = true;
    is_parameters.init_tol = 0.0001;
    is_parameters.res_reduction = 0.01;

    //----------- Solver ---------------- //

    T gamma = 0.2;
    IndexSet<Index2D> Lambda;
    getSparseGridIndexSet(basis2d, Lambda, 1,0,gamma);

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

    RB_Model rb_system(lhs_theta, rhs_theta);

    rb_system.rb_params.ref_param = {{1.}};
    rb_system.rb_params.call = call_cg;

    //----------- RB Base ---------------- //

    RB_BaseModel rb_base(rb_system, rb_truth, Lambda);


    ParamType mu_min = {{0}};
    ParamType mu_max = {{1}};

    rb_base.greedy_params.tol = 1e-05;
    rb_base.greedy_params.min_param = mu_min;
    rb_base.greedy_params.max_param = mu_max;
    rb_base.greedy_params.Nmax = 	6;
    rb_base.greedy_params.nb_training_params = {{nb_stempel}};
    rb_base.greedy_params.log_scaling = {{0}};
    rb_base.greedy_params.print_file = "toy_greedy_info.txt";
    rb_base.greedy_params.trainingdata_folder = "training_data_toy";
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
    rb_system.read_rb_data("training_data_toy");
    rb_base.read_basisfunctions("training_data_toy/bf");

    // Read Test Parameters
    std::vector<ParamType> Xi_test;
    read_paramfile("Xitest.txt", Xi_test);


    cout << "===============================" << endl;
    cout << "Test Parameters: " << endl;
    for(auto& mu : Xi_test){
    	cout << mu[0] << endl;
    }
    cout << "===============================" << endl;


    // Test Reduced Solutions
    std::map<ParamType, std::vector<T> > errs, errbounds;
    int Nmax = rb_system.RB_inner_product.numRows();

    DenseVectorT av_err(Nmax);
    DenseVectorT av_errest(Nmax);

	// For all parameters:
    for(auto& mu : Xi_test){
        
        cout << "Mu = [" << mu[0]  << "]" << endl;
        
    	stringstream ufilename_in, ufilename_out;
    	if(argc >= 6){
    		ufilename_in << solfolder_in << "/";
    	}
    	if(argc == 7){
            ufilename_out << solfolder_out << "/";
    	}
    	ufilename_in << "u_" << mu[0]  << ".txt";
    	ufilename_out << "u_" << mu[0] << ".txt";

    	std::vector<T> errs_mu, errbounds_mu;

    	T err, errbound;
		for(size_t n = 1; n <= 1/*Nmax*/; ++n){

    		// Compute RB solution
    		DenseVectorT u_N = rb_system.get_rb_solution(n, mu);
            cout << "u_" << n << " = " << u_N << endl;
    		// Compute real error
    		DataType u_approx = rb_base.reconstruct_u_N(u_N, n);
    		//u_approx -= u;
    		err = 0.;//u_approx.norm(2.);
            cout << "Err = " << err << endl;
    		errs_mu.push_back(err);
    		
    		av_err(n) += err;

            // Plot
    		for (int i = 1; i<= 1; ++i) {
                std::string     plot;
                plot     = "plot_reduced_";
                plot     += std::to_string(i);

                std::cout << "Printing plot to file " << plot << std::endl;
                precon _p;
                plot2D(basis2d, u_approx, _p,
                                         zero,
                                         0., 1.,
                                         0., 1.,
                                         1e-03,
                                         plot.c_str());
            }
            exit(0);

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

