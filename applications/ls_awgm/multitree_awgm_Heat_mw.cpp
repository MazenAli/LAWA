#include "problem_Heat_mw.h"


int main () {

	//===============================================================//
	//========= PROBLEM SETUP  =======================//
	//===============================================================//
	
    int d   = 2;
    int d_  = 2;
    int j0  = 2;

    //getchar();

    /// Basis initialization
    TrialBasis_Time      basis_per(d,d_,j0);
    TestBasis_Time       basis_int(d,d_,j0);
    Basis_Space 		 basis_intbc(d,0);
    basis_intbc.enforceBoundaryCondition<DirichletBC>();

    Basis2D_Trial basis2d_trial(basis_per,basis_intbc);
    Basis2D_Test  basis2d_test(basis_int,basis_intbc);

    /// Initialization of operator

    // Bilinear Forms
    Convection1D_Time			ConvectionBil_t(basis_per, basis_int);
    Identity1D_Time 		    IdentityBil_t(basis_per, basis_int);
    Identity1D_Space 	        IdentityBil_x(basis_intbc, basis_intbc);
    Laplace1D_Space 	        LaplaceBil_x(basis_intbc, basis_intbc);
    
    RefConvection1D_Time 		RefConvectionBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Time 		    RefIdentityBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Space 	    RefIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefLaplace1D_Space 	        RefLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);

    // Transposed Bilinear Forms
    TranspConvection1D_Time 	TranspConvectionBil_t(basis_per, basis_int);
    TranspIdentity1D_Time 		TranspIdentityBil_t(basis_per, basis_int);
    TranspIdentity1D_Space 	    TranspIdentityBil_x(basis_intbc, basis_intbc);
    TranspLaplace1D_Space 	    TranspLaplaceBil_x(basis_intbc, basis_intbc);
    
    RefTranspConvection1D_Time 	RefTranspConvectionBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefTranspIdentity1D_Time 	RefTranspIdentityBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefTranspIdentity1D_Space 	RefTranspIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefTranspLaplace1D_Space 	RefTranspLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);

    /// Initialization of local operator
    LOp_Conv1D_Time		lOp_Conv1D_t(basis_int, basis_per, RefConvectionBil_t, ConvectionBil_t);
    LOp_Id1D_Time		lOp_Id1D_t  (basis_int, basis_per, RefIdentityBil_t, IdentityBil_t);
    LOp_Id1D_Space		lOp_Id1D_x  (basis_intbc, basis_intbc, RefIdentityBil_x, IdentityBil_x);
    LOp_Lapl1D_Space	lOp_Lapl1D_x(basis_intbc, basis_intbc, RefLaplaceBil_x, LaplaceBil_x);
    
    LOpT_Conv1D_Time	lOpT_Conv1D_t(basis_per, basis_int, RefTranspConvectionBil_t, TranspConvectionBil_t);
    LOpT_Id1D_Time		lOpT_Id1D_t  (basis_per, basis_int, RefTranspIdentityBil_t, TranspIdentityBil_t);
    LOpT_Id1D_Space		lOpT_Id1D_x  (basis_intbc, basis_intbc, RefTranspIdentityBil_x, TranspIdentityBil_x);
    LOpT_Lapl1D_Space	lOpT_Lapl1D_x(basis_intbc, basis_intbc, RefTranspLaplaceBil_x, TranspLaplaceBil_x);

    LOp_Conv_Id_2D		localConvectionIdentityOp2D(lOp_Conv1D_t, lOp_Id1D_x);
    LOp_Id_Lapl_2D		localIdentityLaplaceOp2D(lOp_Id1D_t, lOp_Lapl1D_x);
    
    LOpT_Conv_Id_2D		transpLocalConvectionIdentityOp2D(lOpT_Conv1D_t, lOpT_Id1D_x);
    LOpT_Id_Lapl_2D		transpLocalIdentityLaplaceOp2D(lOpT_Id1D_t, lOpT_Lapl1D_x);

    localConvectionIdentityOp2D.setJ(9);
    localIdentityLaplaceOp2D.setJ(9);
    transpLocalConvectionIdentityOp2D.setJ(9);
    transpLocalIdentityLaplaceOp2D.setJ(9);

    // Use CompoundLocalOperator2D
    COp_Heat            localOperator2D(localConvectionIdentityOp2D,localIdentityLaplaceOp2D);
    COpT_Heat           transpLocalOperator2D(transpLocalConvectionIdentityOp2D,transpLocalIdentityLaplaceOp2D);

    // Use FlexibleCompoundLocalOperator2D
//    vector<AbstractLocalOperator2D<T>* > localOperatorVec, transpLocalOperatorVec;
//    localOperatorVec.push_back(&localConvectionIdentityOp2D);
//    localOperatorVec.push_back(&localIdentityLaplaceOp2D);
//    transpLocalOperatorVec.push_back(&transpLocalConvectionIdentityOp2D);
//    transpLocalOperatorVec.push_back(&transpLocalIdentityLaplaceOp2D);
//    FlexibleCompoundLocalOperator2D       localOperator2D(localOperatorVec);
//    FlexibleCompoundLocalOperator2D  	    transpLocalOperator2D(transpLocalOperatorVec);

    /// Initialization of preconditioner
    LeftPrec2D leftPrec(basis2d_test);
    RightPrec2D rightPrec(basis2d_trial);

    NoPrec2D noPrec;

    /// Initialization of rhs

    /// Right Hand Side:
    ///     No Singular Supports in both dimensions
    DenseVectorT sing_support_x;
    DenseVectorT sing_support_t(n+1);
    for(size_t i = 0; i <= n; ++i){
    	sing_support_t(i+1) = i*l;
    }
    ///      Forcing Functions
    SeparableFunction2D<T> F_fct(f_t, sing_support_t, f_x, sing_support_x);
    ///     Peaks: points and corresponding coefficients
    ///             (heights of jumps in derivatives)
    FullColMatrixT nodeltas;
    SeparableRhsIntegral2D			rhs(basis2d_test, F_fct, nodeltas, nodeltas, 20);
    SeparableRhs           			F(rhs,noPrec);

	//===============================================================//
	//===============  AWGM =========================================//
	//===============================================================//


    /* AWGM PG Parameters Default Values
    double tol = 5e-03;
	double alpha = 0.7;
	size_t max_its = 100;
	size_t max_basissize = 400000;
	bool reset_res = false;
	StableExpansionVersion stable_exp_u = FullExpansion;
	StableExpansionVersion stable_exp_res = FullExpansion;
	bool print_info = true;
	bool verbose = true;
	bool plot_solution = false;
	bool verbose_extra = false; //(print added wavelet indizes)
	size_t hashmapsize_trial = 10;
	size_t hashmapsize_test = 10;
	std::string info_filename = "awgm_cgls_conv_info.txt";
	std::string plot_filename = "awgm_cgls_u_plot";
	bool write_intermediary_solutions = false;
    std::string intermediary_solutions_filename = "awgm_cgls_u";
	*/

    /* IS Parameters Default Values
	bool adaptive_tol = true;
	size_t max_its = 100;
	double init_tol = 0.001;
	double res_reduction = 0.01;
	double absolute_tol = 1e-8;
	bool verbose = true;
	*/

    // MultitreeAWGM with default values
    //MT_AWGM multitree_awgm(basis2d_trial, basis2d_test, localOperator2D, transLocalOperator2D,
    //    						F, rightPrec, leftPrec);


    // If you want other parameters
    AWGM_PG_Parameters awgm_parameters;
    IS_Parameters cgls_parameters;
    // .... set them here:
    awgm_parameters.tol = 0.0001;
    //awgm_parameters.stable_exp_u = OnlyTemporalHWExpansion;
    //awgm_parameters.stable_exp_res = OnlyTemporalHWExpansion;
    awgm_parameters.res_construction = DoubleStableExpansion; // ParallelMTCones
    awgm_parameters.plot_solution = false;
    awgm_parameters.verbose_extra = true;
    awgm_parameters.info_filename = "awgm_ExSaw_mw_DoubleExp_conv_info.txt";
    awgm_parameters.plot_filename = "awgm_ExSaw_mw_DoubleExp_u_plot";
    awgm_parameters.write_intermediary_solutions = true;
    awgm_parameters.intermediary_solutions_filename = "awgm_ExSaw_mw_DoubleExp_u_Iteration";
    
    cgls_parameters.adaptive_tol = true;
    cgls_parameters.init_tol = 1e-4;
    cgls_parameters.res_reduction = 0.01;
    cgls_parameters.max_its = 700;

    MT_AWGM multitree_awgm(basis2d_trial, basis2d_test, localOperator2D, transpLocalOperator2D,
    						F, rightPrec, leftPrec, awgm_parameters, cgls_parameters);


    multitree_awgm.awgm_params.print();
    multitree_awgm.is_params.print();

    multitree_awgm.set_sol(dummy);

    /// Initialization of solution vector and initial index sets
    Coefficients<Lexicographical,T,Index2D> u;

    T gamma = 0.2;
    IndexSet<Index2D> LambdaTrial, LambdaTest;
    getSparseGridIndexSet(basis2d_trial,LambdaTrial,2,0,gamma);
    getSparseGridIndexSet(basis2d_test ,LambdaTest ,2,1,gamma);

    Timer time;
    time.start();
    multitree_awgm.solve(u, LambdaTrial, LambdaTest);
    time.stop();
    cout << "Solution took " << time.elapsed() << " seconds" << endl;

    saveCoeffVector2D(u, basis2d_trial, "awgm_ExSaw_mw_DoubleExp_u.txt");
    return 0;
}

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, string filename)
{
    std::ifstream infile (filename.c_str());
    if (infile.is_open()) {
        cerr << "   Indexset file is open." << endl;
    }
    else {
        cerr << "   Indexset file " << filename.c_str()  << " is not open." << endl;
    }

	int t1,t2;
    int j1,j2;
    long k1,k2;
    T coeff;

    while(!infile.eof()) {

    	infile >> t1 >> j1 >> k1 >> t2 >> j2 >> k2 >> coeff;

        if (t1 == 1 && t2 == 1) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (t1 == 1 && t2 == 0) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (t1 == 0 && t2 == 1) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (t1 == 0 && t2 == 0) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else {
            std::cerr << "Could not read file." << std::endl;
            exit(1); return;
        }
    }
}

