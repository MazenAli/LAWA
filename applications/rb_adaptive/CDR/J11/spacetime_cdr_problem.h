#include <iostream>
#include <iomanip>
#include <utility>
#include <array>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

#include <lawa/lawa.h>

/*
#include <lawa/methods/rb/datastructures/thetastructure.h>
#include <lawa/methods/rb/operators/affinelocaloperator.h>
#include <lawa/methods/rb/righthandsides/affinerhs.h>
#include <lawa/methods/adaptive/solvers/multitreeawgm2.h>
#include <lawa/methods/rb/righthandsides/flexiblebilformrhs.h>
#include <lawa/methods/rb/datastructures/mt_truth.h>
#include <lawa/methods/rb/datastructures/rb_system.h>
#include <lawa/methods/rb/datastructures/rb_base.h>
#include <lawa/methods/rb/datastructures/rb_parameters.h>
#include <lawa/methods/rb/datastructures/lb_base.h>
*/


using namespace std;
using namespace lawa;

// quick and dirty finite wavelet expansion implementation
template <typename T, typename Index, typename RHSINTEGRAL,
          typename Preconditioner>
class SimpleRHS : public RHS<T, Index, RHSINTEGRAL, Preconditioner>
{
    private:
        const IndexSet<Index>&      LambdaTest;
    public:
        SimpleRHS(const RHSINTEGRAL& rhsintegral, Preconditioner& P,
                  IndexSet<Index>& _LambdaTest):
                  RHS<T, Index, RHSINTEGRAL, Preconditioner>(rhsintegral, P),
                  LambdaTest(_LambdaTest){}

        T
        operator()(const Index& lambda) override
        {
            typename IndexSet<Index>::const_iterator
                                        found = LambdaTest.find(lambda);

            if (found == LambdaTest.end()) {
                return 0.;
            } else {
                return RHS<T, Index, RHSINTEGRAL, Preconditioner>
                           ::operator()(lambda);
            }
        }

        Coefficients<Lexicographical,T,Index>
        operator()(const IndexSet<Index> &Lambda) override
        {
            typedef typename IndexSet<Index>::const_iterator const_set_it;
            Coefficients<Lexicographical,T,Index> ret;
            for (const_set_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
                T tmp = SimpleRHS<T,Index,RHSINTEGRAL,Preconditioner>::operator()(*lambda);
                ret[*lambda] = tmp;
            }
            return ret;
        }
};

template <typename LocalOperator1, typename LocalOperator2>
struct SimpleOp2D : LocalOperator2D<LocalOperator1, LocalOperator2>
{
    typedef typename LocalOperator1::T      T;
    const IndexSet<Index2D>&                LambdaTrial; // must be a multi-tree!
    const IndexSet<Index2D>&                LambdaTest; // must be a multi-tree!

    SimpleOp2D(LocalOperator1& _localop1, LocalOperator2& _localop2,
               IndexSet<Index2D>& _LambdaTrial,
               IndexSet<Index2D>& _LambdaTest):
               LocalOperator2D<LocalOperator1, LocalOperator2>
               (_localop1, _localop2), LambdaTrial(_LambdaTrial),
               LambdaTest(_LambdaTest){}

    void
    eval(const Coefficients<Lexicographical, T, Index2D>& v,
               Coefficients<Lexicographical, T, Index2D>& AAv) override
    {
        Coefficients<Lexicographical, T, Index2D>   temp = v;

        P(LambdaTrial, temp);
        P(LambdaTest, AAv);
        LocalOperator2D<LocalOperator1, LocalOperator2>::eval(temp, AAv);
    }
};

// quick and dirty implementation of new error estimator

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
class RB_Base_wavest : public RB_Base<RB_Model, TruthModel, DataType, ParamType>
{
    private:
        // sparse grid and thresh parameter
        const _IndexSet&    Lambda;
        double              ETA     = 1e-16;
        bool                recalcF = false;
        long                total   = 0;

    public:
        RB_Base_wavest(RB_Model& _rb_system, TruthModel& _rb_truth,
                       const _IndexSet& _Lambda)
                      :RB_Base<RB_Model, TruthModel, DataType, ParamType>
                       (_rb_system, _rb_truth),
                       Lambda(_Lambda){}


        void
        init()
        {
            this->greedy_params.tighten_tol_rieszA              = false;
            this->greedy_params.tighten_tol_rieszF              = false;
            this->greedy_params.update_rieszF                   = false;
            this->greedy_params.update_rieszA                   = false;
            this->greedy_params.test_estimator_equivalence      = false;
            this->greedy_params.tighten_estimator_accuracy      = false;
            std::cout << "Exact LambdaTest size is " << Lambda.size() << std::endl;
        }

    private:
        void
        calculate_Riesz_RHS_information(bool) override
        {
            std::cout << "Correct RHS being called\n";
            if (recalcF) {
                std::cerr << "Warning: recalculating RHS! This should not happen!\n";
            }
            calc_wav_RHS();
            recalcF = true;

            std::size_t Qf = this->rb_system.Q_f();
            this->rb_system.F_F_representor_norms.engine().resize((int) Qf, (int) Qf);
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


        void
        calculate_Riesz_LHS_information(const DataType& bf, bool) override
        {
            std::cout << "Correct LHS called\n";
            int N = this->rb_basisfunctions.size();
            std::cout << "Current N=" << N << std::endl;
            calc_wav_LHS(bf);

            // Update the Riesz Representor Norms A x A
            std::size_t Qa = this->rb_system.Q_a();
            for(int n1 = 0; n1 < N; ++n1) {
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


        void
        calc_wav_RHS()
        {
            assert(this->F_representors.size()==0);

            auto& F = this->rb_truth.access_RieszSolver_F().get_rhs();
            auto& P = this->rb_truth.access_RieszSolver_F().get_testprec();

            std::size_t Qf = this->rb_system.Q_f();
            assert(Qf>0);
            this->F_representors.resize(Qf);

            for (std::size_t i=0; i<Qf; ++i) {
                F.set_active_comp(i);
                this->F_representors[i] = F(Lambda, P); // precondtioned
                THRESH_NoCopy(this->F_representors[i], ETA);
                this->total += this->F_representors[i].size();
            }
            std::cout << "Current total size of F is " << total << std::endl;
        }


        void
        calc_wav_LHS(const DataType& bf)
        {
            std::size_t             Qa = this->rb_system.Q_a();
            assert(Qa>0);
            std::vector<DataType>   temp(Qa);
            auto& A = this->rb_truth.access_RieszSolver_A().get_rhs().get_bilformvec();
            auto& P = this->rb_truth.access_RieszSolver_A().get_testprec();


            for (std::size_t i=0; i<Qa; ++i) {
                FillWithZeros(Lambda, temp[i]);
                (*A[i]).eval(bf, temp[i]);

                // left preconditioning
                for (auto& el : temp[i]) {
                    el.second *= P(el.first)*(-1.);
                }

                THRESH_NoCopy(temp[i], ETA);
                total += temp[i].size();
            }

            this->A_representors.push_back(temp);
            std::cout << "Current total size of A and F is " << total << std::endl;
        }
};


//===============================================================//
//========= Simple RB System	  ===============================//
//===============================================================//

template <typename _T, typename _ParamType, typename LB_System>
class Simple_RB_System: public RB_System<_T,_ParamType> {

    typedef std::map<_ParamType, _T>   StabConstVec;

public:
	Simple_RB_System(ThetaStructure<_ParamType>& _thetas_a,
			  	  ThetaStructure<_ParamType>& _thetas_f,
			  	  LB_System& _lb_system)
	 : RB_System<_T,_ParamType>(_thetas_a, _thetas_f),
	   lb_system(_lb_system){}

    _T
	alpha_LB(_ParamType& mu)
    {
    	// If it is precalculated, return value
    	for(auto& el: alphas){
    	    bool is_mu = true;
			for(std::size_t i = 0; i < mu.size(); ++i){
				if( (mu[i]-(el.first)[i]) > 1e-4 ){
					is_mu = false;
					break;
				}
			}
			if(is_mu){
				return el.second;
			}
    	}
    	// Else we have to calculate it
    	_T alpha = lb_system.calculate_alpha(mu);
    	alphas.insert(std::make_pair(mu, alpha));
    	return alpha;
    }

    void
    read_alpha(const char* filename){
        ifstream alphaFile(filename);
        while(alphaFile.good()){
             _ParamType mu;
             _T alpha;
             size_t pdim = ParamInfo<_ParamType>::dim;
             for(size_t p = 0; p < pdim; ++p){
            	 alphaFile >> mu[p];
             }
             alphaFile >> alpha;
             alphas.insert(std::make_pair(mu, alpha));
        }
    }

    void
    write_alpha(const char* filename){
        size_t pdim = ParamInfo<_ParamType>::dim;
        ofstream alphaFile(filename);
        if(alphaFile.is_open()){
        	for(auto& a : alphas){
                for(size_t p = 0; p < pdim; ++p){
               	 alphaFile << a.first[p] << " ";
                }
                alphaFile << a.second << std::endl;
        	}
        }
        else{
        	std::cerr << "Couldn't open file " << filename << " for writing." << std::endl;
        	return;
        }
    }

private:
    StabConstVec	alphas;
    LB_System&		lb_system;
};


//===============================================================//
//========= TYPEDEFS  =======================//
//===============================================================//

//==== General ====//
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef Coefficients<Lexicographical,T,Index2D>						CoeffVector;

//==== Basis 1D & 2D ====//
typedef Basis<T, Primal, Periodic, CDF>	                            TrialBasis_Time;
//typedef Basis<T, Primal, Periodic, CDF>								TestBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>						    TestBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>							Basis_Space;
typedef TrialBasis_Time::RefinementBasis                            TrialRefBasis_Time;
typedef TestBasis_Time::RefinementBasis                           	TestRefBasis_Time;
typedef Basis_Space::RefinementBasis                                RefBasis_Space;

// !!! ATTENTION: Use typedefs in space for definition of space-time inner product !!!
//    => Assumption that test basis in time = space basis (+/- boundary conditions)

typedef TensorBasis2D<Adaptive,TrialBasis_Time,Basis_Space>         Basis2D_Trial;
typedef TensorBasis2D<Adaptive,TestBasis_Time,Basis_Space>          Basis2D_Test;

//==== Adaptive Operators ====//
typedef AdaptiveConvectionOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>                   Convection1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>                     Identity1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,Basis_Space,Basis_Space>                            Identity1D_Space;
typedef AdaptiveLaplaceOperator1D_PG<T,Basis_Space,Basis_Space>                             Laplace1D_Space;
typedef AdaptiveWeightedPDEOperator1D_PG<T,Basis_Space,Basis_Space>                         Convection1D_Space;
typedef AdaptiveIdentityOperator1D_PG<T,TrialBasis_Time,TrialBasis_Time>                    Identity1D_Time_Trial;
typedef AdaptiveConvectionOperator1D_PG<T,TrialBasis_Time,TrialBasis_Time>                  Convection1D_Time_Trial;

typedef TransposedAdaptiveConvectionOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>         TranspConvection1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>           TranspIdentity1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,Basis_Space,Basis_Space>                  TranspIdentity1D_Space;
typedef TransposedAdaptiveLaplaceOperator1D_PG<T,Basis_Space,Basis_Space>                   TranspLaplace1D_Space;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,Basis_Space,Basis_Space>               TranspConvection1D_Space;

typedef AdaptiveConvectionOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>             RefConvection1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>               RefIdentity1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,RefBasis_Space,RefBasis_Space>                      RefIdentity1D_Space;
typedef AdaptiveLaplaceOperator1D_PG<T,RefBasis_Space,RefBasis_Space>                       RefLaplace1D_Space;
typedef AdaptiveWeightedPDEOperator1D_PG<T,RefBasis_Space,RefBasis_Space>                   RefConvection1D_Space;
typedef AdaptiveIdentityOperator1D_PG<T,TrialRefBasis_Time,TrialRefBasis_Time>              RefIdentity1D_Time_Trial;
typedef AdaptiveConvectionOperator1D_PG<T,TrialRefBasis_Time,TrialRefBasis_Time>            RefConvection1D_Time_Trial;

typedef TransposedAdaptiveConvectionOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>   RefTranspConvection1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>     RefTranspIdentity1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,RefBasis_Space,RefBasis_Space>            RefTranspIdentity1D_Space;
typedef TransposedAdaptiveLaplaceOperator1D_PG<T,RefBasis_Space,RefBasis_Space>             RefTranspLaplace1D_Space;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,RefBasis_Space,RefBasis_Space>         RefTranspConvection1D_Space;

//==== LocalOperators ====//
typedef LocalOperator1D<TestBasis_Time,TrialBasis_Time,
                        RefConvection1D_Time,Convection1D_Time>                 LOp_Conv1D_Time;
typedef LocalOperator1D<TestBasis_Time,TrialBasis_Time,
                        RefIdentity1D_Time,Identity1D_Time>                     LOp_Id1D_Time;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefIdentity1D_Space,Identity1D_Space>	                LOp_Id1D_Space;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefLaplace1D_Space,Laplace1D_Space>	                    LOp_Lapl1D_Space;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefConvection1D_Space,Convection1D_Space>	            LOp_Conv1D_Space;
typedef LocalOperator1D<TrialBasis_Time,TrialBasis_Time,
                        RefIdentity1D_Time_Trial,Identity1D_Time_Trial>         LOp_Id1D_Time_Trial;
typedef LocalOperator1D<TrialBasis_Time,TrialBasis_Time,
                        RefConvection1D_Time_Trial,Convection1D_Time_Trial>     LOp_Conv1D_Time_Trial;

typedef LocalOperator1D<TrialBasis_Time,TestBasis_Time,
                        RefTranspConvection1D_Time,TranspConvection1D_Time>     LOpT_Conv1D_Time;
typedef LocalOperator1D<TrialBasis_Time,TestBasis_Time,
                        RefTranspIdentity1D_Time,TranspIdentity1D_Time>         LOpT_Id1D_Time;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefTranspIdentity1D_Space,TranspIdentity1D_Space>	    LOpT_Id1D_Space;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefTranspLaplace1D_Space,TranspLaplace1D_Space>	        LOpT_Lapl1D_Space;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefTranspConvection1D_Space,TranspConvection1D_Space>	LOpT_Conv1D_Space;

typedef SimpleOp2D<LOp_Conv1D_Time, LOp_Id1D_Space>                        LOp_Conv_Id_2D;
typedef SimpleOp2D<LOp_Id1D_Time, LOp_Id1D_Space>                          LOp_Id_Id_2D;
typedef SimpleOp2D<LOp_Id1D_Time, LOp_Lapl1D_Space>                        LOp_Id_Lapl_2D;
typedef SimpleOp2D<LOp_Id1D_Time, LOp_Conv1D_Space>                        LOp_Id_Conv_2D;
typedef SimpleOp2D<LOpT_Conv1D_Time, LOpT_Id1D_Space>                      LOpT_Conv_Id_2D;
typedef SimpleOp2D<LOpT_Id1D_Time, LOpT_Id1D_Space>                        LOpT_Id_Id_2D;
typedef SimpleOp2D<LOpT_Id1D_Time, LOpT_Lapl1D_Space>                      LOpT_Id_Lapl_2D;
typedef SimpleOp2D<LOpT_Id1D_Time, LOpT_Conv1D_Space>                      LOpT_Id_Conv_2D;


typedef LocalOperator2D<LOp_Id1D_Space,LOp_Id1D_Space>							LOp_Test_Id_Id_2D;
typedef LocalOperator2D<LOp_Id1D_Space,LOp_Lapl1D_Space>						LOp_Test_Id_Lapl_2D;
typedef LocalOperator2D<LOp_Id1D_Time_Trial,LOp_Id1D_Space>						LOp_Trial_Id_Id_2D;
typedef LocalOperator2D<LOp_Id1D_Time_Trial,LOp_Conv1D_Space>					LOp_Trial_Id_Conv_2D;
typedef LocalOperator2D<LOp_Id1D_Time_Trial,LOp_Lapl1D_Space>					LOp_Trial_Id_Lapl_2D;
typedef LocalOperator2D<LOp_Conv1D_Time_Trial,LOp_Id1D_Space>					LOp_Trial_Conv_Id_2D;

//==== CompoundOperators ====//
typedef FlexibleCompoundLocalOperator<Index2D,AbstractLocalOperator2D<T> > 		Flex_COp_2D;
typedef CompoundLocalOperator<Index2D,LOp_Conv_Id_2D,LOp_Id_Lapl_2D>    		COp_Heat;
typedef CompoundLocalOperator<Index2D,LOpT_Conv_Id_2D,LOpT_Id_Lapl_2D>    	    COpT_Heat;
typedef CompoundLocalOperator<Index2D,LOp_Test_Id_Id_2D, LOp_Test_Id_Lapl_2D>	COp_TestInnProd;
typedef CompoundLocalOperator<Index2D,LOp_Trial_Id_Id_2D, LOp_Trial_Id_Lapl_2D>	TrialInnProdY;

//==== Preconditioners ====//
typedef AdaptiveLeftNormPreconditioner2D<T,Basis2D_Test>            LeftPrec2D;
typedef AdaptiveRightNormPreconditioner2D_c<T,Basis2D_Trial>        RightPrec2D;
typedef NoPreconditioner<T, Index2D>								NoPrec2D;

//==== RightHandSides ====//
typedef SeparableRHS2D<T,Basis2D_Test>                              SeparableRhsIntegral2D;
typedef SimpleRHS<T,Index2D,SeparableRhsIntegral2D,
            NoPrec2D>                                         		SeparableRhs;
typedef SeparableRHS2D<T,Basis2D_Trial>                             SeparableRhsIntegral2D_Trial;
typedef RHS<T,Index2D,SeparableRhsIntegral2D_Trial,
            NoPrec2D>                                         		SeparableRhs_Trial;


//==== RB Stuff ====//
const size_t PDim = 2;

typedef CoeffVector 								 						DataType;
typedef array<T,PDim>	 													ParamType;

typedef AffineLocalOperator<Index2D,AbstractLocalOperator2D<T>,ParamType>	Affine_Op_2D;
typedef AffineRhs<T,Index2D,SeparableRhs,ParamType>							Affine_Rhs_2D;
typedef FlexibleCompoundRhs<T,Index2D,SeparableRhs>							RieszF_Rhs_2D;
typedef FlexibleBilformRhs<Index2D,AbstractLocalOperator2D<T> >				RieszA_Rhs_2D;
typedef FlexibleCompoundRhs<T,Index2D,SeparableRhs_Trial>					Flex_Rhs_2D;

typedef AffineBilformRhs<Index2D, AbstractLocalOperator2D<T>, ParamType>			AffineA_Rhs_2D;
typedef ResidualRhs<Index2D, AffineA_Rhs_2D, Affine_Rhs_2D, ParamType, DataType>	RieszRes_Rhs_2D;

typedef MultiTreeAWGM_PG<Index2D,Basis2D_Trial, Basis2D_Test,Affine_Op_2D,
		Affine_Op_2D,Affine_Rhs_2D,RightPrec2D,LeftPrec2D>					MT_AWGM_Truth;
typedef MultiTreeAWGM2<Index2D,Basis2D_Test, COp_TestInnProd,
		RieszF_Rhs_2D,LeftPrec2D>											MT_AWGM_Riesz_F;
typedef MultiTreeAWGM2<Index2D,Basis2D_Test, COp_TestInnProd,
		RieszA_Rhs_2D,LeftPrec2D>											MT_AWGM_Riesz_A;
typedef MultiTreeAWGM2<Index2D,Basis2D_Test, COp_TestInnProd,
		RieszRes_Rhs_2D,LeftPrec2D>											MT_AWGM_Riesz_Res;

typedef MT_Truth<DataType,ParamType,MT_AWGM_Truth,
				 MT_AWGM_Riesz_F,MT_AWGM_Riesz_A,TrialInnProdY,
				 Flex_COp_2D, Flex_Rhs_2D, MT_AWGM_Riesz_Res>				MTTruthSolver;

typedef LB_Base<ParamType, MTTruthSolver> 									LB_Type;
typedef Simple_RB_System<T,ParamType, LB_Type>								RB_Model;
typedef RB_Base_wavest<RB_Model,MTTruthSolver,DataType,ParamType,IndexSet<Index2D>>			
                                                                            RB_BaseModel;

T f_t(T t)
{
    return cos(2.L*M_PI*t);
}

/*T f_x(T x)
{
    return -4*(x-0.5)*(x-0.5) + 1;
}
*/
T f_x(T)
{
	return 1.;
}

T dummy(T, T)
{
    return 0;
}

T no_theta(const std::array<T,PDim>& /*mu*/)
{
	return 1.;
}

T theta_conv(const std::array<T,PDim>& mu)
{
	return mu[0];
}

T theta_reac(const std::array<T,PDim>& mu)
{
	return mu[1];
}

T
weight_convection(T x)
{
  return 0.5 - x;
}

T zero_fct(T /*x*/){
	return 0;
}
