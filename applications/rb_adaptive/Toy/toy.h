#include <iostream>
#include <iomanip>
#include <utility>
#include <array> 
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

// quick and dirty implementation of new error estimator

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
class RB_Base_wavest : public RB_Base<RB_Model, TruthModel, DataType, ParamType>
{
    private:
        // Sample parameters
        const _IndexSet&    Lambda;
        double              tol     = 1e-8;
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
            std::cout << "Init Lambda size is " << Lambda.size() << std::endl;
        }

    public:
        void
        calculate_Riesz_RHS_information(bool) override
        {
            std::cout << "Correct RHS called\n";
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
            this->F_representors.clear();
            assert(this->F_representors.size()==0);

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
                         tol, true);
                this->total += this->F_representors[i].size();
            }
            std::cout << "Current total size of F is " << total << std::endl;
        }


        void
        calc_wav_LHS(const DataType& bf)
        {
            std::size_t Qa = this->rb_system.Q_a();
            assert(Qa>0);
            std::vector<DataType>   temp(Qa);
            auto& A     = this->rb_truth.access_solver().get_lhs().get_localops();
            auto& P     = this->rb_truth.access_solver().get_testprec();
            auto& basis = this->rb_truth.access_solver().get_testbasis();


            for (std::size_t i=0; i<Qa; ++i) {
                sample_Au(basis, Lambda, (*A[i]), P, bf, temp[i],
                          tol, true);
                // Because Kristina
                for (auto& lambda : temp[i]) {
                    lambda.second *= -1.;
                }
                total += temp[i].size();
            }

            this->A_representors.push_back(temp);
            std::cout << "Current total size of A and F is " << total << std::endl;
        }

    public:
        // Something still wrong with this procedure...
        void
        calc_rb_data()
        {
            std::cout << "Calculating RB data\n";
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
            /*this->rb_system.A_A_representor_norms.resize(n_bf);
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
            */

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

            /*/ RB_A
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
            */

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

            /*/ RB_InnerProduct
            this->rb_system.RB_inner_product.resize(n_bf, n_bf);
            for (unsigned int i = 1; i <= n_bf; ++i) {
                for (unsigned int j = i; j <= n_bf; ++j) {
                    this->rb_system.RB_inner_product(i, j) = this->rb_truth.innprod_Y_u_u(this->rb_basisfunctions[i-1],
                                                                                    this->rb_basisfunctions[j-1]);
                    this->rb_system.RB_inner_product(j, i) = this->rb_system.RB_inner_product(i,j);
                }
            }

            std::cout << std::endl << "||------- RB_InnerProduct  ----------------------------------||" << std::endl;
            std::cout << this->rb_system.RB_inner_product << std::endl;
            */
        }
};

        //===============================================================//
        //========= Simple RB System	  ===============================//
        //===============================================================//

template <typename _T, typename _ParamType, typename LB_System>
class Simple_RB_System: public RB_System<_T,_ParamType> {

    public:
        Simple_RB_System(ThetaStructure<_ParamType>& _thetas_a,
                      ThetaStructure<_ParamType>& _thetas_f)
         : RB_System<_T,_ParamType>(_thetas_a, _thetas_f)
           {}

        _T
        alpha_LB(_ParamType& /*mu*/)
        {
            return 0.5;
        }

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
const DomainType      _DomainType    = Interval;
const FunctionSide    FctSide       = Orthogonal;
const Construction    Constr        = Multi;
const Construction    RefConstr     = MultiRefinement;
typedef Basis<T,FctSide, _DomainType, Constr>					_Basis;
typedef _Basis::RefinementBasis                                  RefBasis;

typedef TensorBasis2D<Adaptive,_Basis,_Basis>                     Basis2D;

//==== Adaptive Operators ====//
typedef AdaptiveIdentityOperator1D<T,FctSide, _DomainType, Constr>                Identity1D;
typedef AdaptiveLaplaceOperator1D<T,FctSide, _DomainType, Constr>                 Laplace1D;

typedef AdaptiveIdentityOperator1D<T,FctSide, _DomainType, RefConstr>             RefIdentity1D;
typedef AdaptiveLaplaceOperator1D<T,FctSide, _DomainType, RefConstr>              RefLaplace1D;

//==== LocalOperators ====//
typedef LocalOperator1D<_Basis,_Basis,RefIdentity1D,Identity1D>                   LOp_Id1D;
typedef LocalOperator1D<_Basis,_Basis,RefLaplace1D,Laplace1D>                     LOp_Lapl1D;

typedef LocalOperator2D<LOp_Id1D, LOp_Id1D>                        	                LOp_Id_Id_2D;
typedef LocalOperator2D<LOp_Lapl1D, LOp_Id1D>                                       LOp_Lapl_Id_2D;
typedef LocalOperator2D<LOp_Id1D, LOp_Lapl1D>                                       LOp_Id_Lapl_2D;

//==== CompoundOperators ====//
typedef FlexibleCompoundLocalOperator<Index2D,AbstractLocalOperator2D<T> > 		        Flex_COp_2D;
typedef CompoundLocalOperator<Index2D,LOp_Id_Id_2D,LOp_Lapl_Id_2D,LOp_Id_Lapl_2D>       H1_InnProd_2D;

//==== Preconditioners ====//
typedef H1NormPreconditioner2D<T, Basis2D>                          Prec2D;
typedef NoPreconditioner<T, Index2D>								NoPrec2D;

//==== RightHandSides ====//
typedef SmoothRHSWithAlignedSing2D<T, Basis2D, SparseGridGP>          RhsIntegral2D;
typedef RHS<T, Index2D, RhsIntegral2D, NoPrec2D>                    Rhs;


//==== RB Stuff ====//
const size_t PDim = 1;

typedef CoeffVector 								 						DataType;
typedef array<T,PDim>	 													ParamType;


typedef AffineLocalOperator<Index2D,Flex_COp_2D,ParamType>	        Affine_Op_2D;
typedef AffineRhs<T,Index2D,Rhs,ParamType>							Affine_Rhs_2D;
typedef FlexibleCompoundRhs<T,Index2D,Rhs>							RieszF_Rhs_2D;
typedef FlexibleBilformRhs<Index2D, Flex_COp_2D>				    RieszA_Rhs_2D;

typedef FlexibleCompoundRhs<T,Index2D,Rhs>     					Flex_Rhs_2D;

typedef AffineBilformRhs<Index2D, AbstractLocalOperator2D<T>, ParamType>			AffineA_Rhs_2D;
typedef ResidualRhs<Index2D, AffineA_Rhs_2D, Affine_Rhs_2D, ParamType, DataType>	RieszRes_Rhs_2D;


typedef MultiTreeAWGM2<Index2D,Basis2D, Affine_Op_2D,
		Affine_Rhs_2D,Prec2D>                           					MT_AWGM_Truth;
typedef MultiTreeAWGM2<Index2D,Basis2D, H1_InnProd_2D,
		RieszF_Rhs_2D,Prec2D>											    MT_AWGM_Riesz_F;
typedef MultiTreeAWGM2<Index2D,Basis2D, H1_InnProd_2D,
		RieszA_Rhs_2D,Prec2D>											    MT_AWGM_Riesz_A;
typedef MultiTreeAWGM2<Index2D,Basis2D, H1_InnProd_2D,
		RieszRes_Rhs_2D,Prec2D>												MT_AWGM_Riesz_Res;

typedef MT_Truth<DataType,ParamType,MT_AWGM_Truth,
				 MT_AWGM_Riesz_F,MT_AWGM_Riesz_A,H1_InnProd_2D,
				 Affine_Op_2D, Flex_Rhs_2D, MT_AWGM_Riesz_Res>				MTTruthSolver;

typedef LB_Base<ParamType, MTTruthSolver> 									        LB_Type;
typedef Simple_RB_System<T,ParamType, LB_Type>								        RB_Model;
typedef RB_Base_wavest<RB_Model,MTTruthSolver,DataType,ParamType,IndexSet<Index2D>> RB_BaseModel;

#define mx1 0.2
#define mx2 0.5
#define mx3 0.8
#define my1 0.2
#define my2 0.5
#define my3 0.8
#define a   1./35.
#define nb_stempel 9

#define mx4 0.3
#define mx5 0.4
#define mx6 0.6
#define mx7 0.7
#define my4 0.3
#define my5 0.4
#define my6 0.6
#define my7 0.7

#define addon 40
#define offset 1



double
ux(double x, double mx)
{
    return std::exp(-1*std::pow(x-mx, 2.)/(a*a));
}


double
xux(double x, double mx)
{
    return ux(x, mx)*(-1.*(x-mx)*2./(a*a));
}


double
u1(double x)
{
    return ux(x, mx1);
}


double
ux1(double x)
{
    return xux(x, mx1);
}


double
u2(double x)
{
    return ux(x, mx2);
}


double
ux2(double x)
{
    return xux(x, mx2);
}


double
u3(double x)
{
    return ux(x, mx3);
}


double
ux3(double x)
{
    return xux(x, mx3);
}


double
u4(double x)
{
    return ux(x, mx4);
}


double
ux4(double x)
{
    return xux(x, mx4);
}


double
u5(double x)
{
    return ux(x, mx5);
}


double
ux5(double x)
{
    return xux(x, mx5);
}


double
u6(double x)
{
    return ux(x, mx6);
}


double
ux6(double x)
{
    return xux(x, mx6);
}


double
u7(double x)
{
    return ux(x, mx7);
}


double
ux7(double x)
{
    return xux(x, mx7);
}


double
uxx(double x, double mx)
{
    return ux(x, mx)*(std::pow(x-mx,2.)*4./(a*a*a*a)
           -2./(a*a));
}


double
f(double x, double y, double m1, double m2)
{
    return -1.*(uxx(x, m1)*ux(y, m2)+uxx(y, m2)*ux(x, m1));
}


double
f1(double x, double y)
{
    return -1.*(uxx(x, mx1)*ux(y, my1)+uxx(y, my1)*ux(x, mx1));
}


double
f2(double x, double y)
{
    return -1.*(uxx(x, mx2)*ux(y, my1)+uxx(y, my1)*ux(x, mx2));
}


double
f3(double x, double y)
{
    return -1.*(uxx(x, mx3)*ux(y, my1)+uxx(y, my1)*ux(x, mx3));
}


double
f4(double x, double y)
{
    return -1.*(uxx(x, mx1)*ux(y, my2)+uxx(y, my2)*ux(x, mx1));
}


double
f5(double x, double y)
{
    return -1.*(uxx(x, mx2)*ux(y, my2)+uxx(y, my2)*ux(x, mx2));
}


double
f6(double x, double y)
{
    return -1.*(uxx(x, mx3)*ux(y, my2)+uxx(y, my2)*ux(x, mx3));
}


double
f7(double x, double y)
{
    return -1.*(uxx(x, mx1)*ux(y, my3)+uxx(y, my3)*ux(x, mx1));
}


double
f8(double x, double y)
{
    return -1.*(uxx(x, mx2)*ux(y, my3)+uxx(y, my3)*ux(x, mx2));
}


double
f9(double x, double y)
{
    return -1.*(uxx(x, mx3)*ux(y, my3)+uxx(y, my3)*ux(x, mx3));
}


double
f10(double x, double y)
{
    return f(x, y, mx1, my4);
}


double
f11(double x, double y)
{
    return f(x, y, mx1, my5);
}


double
f12(double x, double y)
{
    return f(x, y, mx1, my6);
}


double
f13(double x, double y)
{
    return f(x, y, mx1, my7);
}


double
f14(double x, double y)
{
    return f(x, y, mx2, my4);
}


double
f15(double x, double y)
{
    return f(x, y, mx2, my5);
}


double
f16(double x, double y)
{
    return f(x, y, mx2, my6);
}


double
f17(double x, double y)
{
    return f(x, y, mx2, my7);
}


double
f18(double x, double y)
{
    return f(x, y, mx3, my4);
}


double
f19(double x, double y)
{
    return f(x, y, mx3, my5);
}


double
f20(double x, double y)
{
    return f(x, y, mx3, my6);
}


double
f21(double x, double y)
{
    return f(x, y, mx3, my7);
}


double
f22(double x, double y)
{
    return f(x, y, mx4, my1);
}


double
f23(double x, double y)
{
    return f(x, y, mx4, my2);
}


double
f24(double x, double y)
{
    return f(x, y, mx4, my3);
}


double
f25(double x, double y)
{
    return f(x, y, mx4, my4);
}


double
f26(double x, double y)
{
    return f(x, y, mx4, my5);
}


double
f27(double x, double y)
{
    return f(x, y, mx4, my6);
}


double
f28(double x, double y)
{
    return f(x, y, mx4, my7);
}


double
f29(double x, double y)
{
    return f(x, y, mx5, my1);
}


double
f30(double x, double y)
{
    return f(x, y, mx5, my2);
}


double
f31(double x, double y)
{
    return f(x, y, mx5, my3);
}


double
f32(double x, double y)
{
    return f(x, y, mx5, my4);
}


double
f33(double x, double y)
{
    return f(x, y, mx5, my5);
}


double
f34(double x, double y)
{
    return f(x, y, mx5, my6);
}


double
f35(double x, double y)
{
    return f(x, y, mx5, my7);
}


double
f36(double x, double y)
{
    return f(x, y, mx6, my1);
}


double
f37(double x, double y)
{
    return f(x, y, mx6, my2);
}


double
f38(double x, double y)
{
    return f(x, y, mx6, my3);
}


double
f39(double x, double y)
{
    return f(x, y, mx6, my4);
}


double
f40(double x, double y)
{
    return f(x, y, mx6, my5);
}


double
f41(double x, double y)
{
    return f(x, y, mx6, my6);
}


double
f42(double x, double y)
{
    return f(x, y, mx6, my7);
}


double
f43(double x, double y)
{
    return f(x, y, mx7, my1);
}


double
f44(double x, double y)
{
    return f(x, y, mx7, my2);
}


double
f45(double x, double y)
{
    return f(x, y, mx7, my3);
}


double
f46(double x, double y)
{
    return f(x, y, mx7, my4);
}


double
f47(double x, double y)
{
    return f(x, y, mx7, my5);
}


double
f48(double x, double y)
{
    return f(x, y, mx7, my6);
}


double
f49(double x, double y)
{
    return f(x, y, mx7, my7);
}


T dummy(T, T)
{
    return 0;
}


T no_theta(const std::array<T,PDim>& /*mu*/)
{
	return 1.;
}


T theta_chi_1(const std::array<T,PDim>& mu)
{
    if(mu[0] <= 1./nb_stempel){
       return 1;
    }
	return 0;
}


T theta_chi_2(const std::array<T,PDim>& mu)
{
    if(mu[0] > 1./nb_stempel && mu[0] <= 2./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_3(const std::array<T,PDim>& mu)
{
    if(mu[0] > 2./nb_stempel && mu[0] <= 3./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_4(const std::array<T,PDim>& mu)
{
    if(mu[0] > 3./nb_stempel && mu[0] <= 4./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_5(const std::array<T,PDim>& mu)
{
    if(mu[0] > 4./nb_stempel && mu[0] <= 5./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_6(const std::array<T,PDim>& mu)
{
    if(mu[0] > 5./nb_stempel && mu[0] <= 6./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_7(const std::array<T,PDim>& mu)
{
    if(mu[0] > 6./nb_stempel && mu[0] <= 7./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_8(const std::array<T,PDim>& mu)
{
    if(mu[0] > 7./nb_stempel && mu[0] <= 8./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_9(const std::array<T,PDim>& mu)
{
    if(mu[0] > 8./nb_stempel && mu[0] <= 9./nb_stempel){
       return 1;
    }
	return 0;
}


T theta_chi_10(const std::array<T,PDim>& mu)
{
    if(mu[0] > offset && mu[0] <= (offset+1./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_11(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+1./addon) && mu[0] <= (offset+2./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_12(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+2./addon) && mu[0] <= (offset+3./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_13(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+3./addon) && mu[0] <= (offset+4./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_14(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+4./addon) && mu[0] <= (offset+5./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_15(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+5./addon) && mu[0] <= (offset+6./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_16(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+6./addon) && mu[0] <= (offset+7./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_17(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+7./addon) && mu[0] <= (offset+8./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_18(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+8./addon) && mu[0] <= (offset+9./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_19(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+9./addon) && mu[0] <= (offset+10./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_20(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+10./addon) && mu[0] <= (offset+11./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_21(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+11./addon) && mu[0] <= (offset+12./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_22(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+12./addon) && mu[0] <= (offset+13./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_23(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+13./addon) && mu[0] <= (offset+14./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_24(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+14./addon) && mu[0] <= (offset+15./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_25(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+15./addon) && mu[0] <= (offset+16./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_26(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+16./addon) && mu[0] <= (offset+17./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_27(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+17./addon) && mu[0] <= (offset+18./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_28(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+18./addon) && mu[0] <= (offset+19./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_29(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+19./addon) && mu[0] <= (offset+20./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_30(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+20./addon) && mu[0] <= (offset+21./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_31(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+21./addon) && mu[0] <= (offset+22./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_32(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+22./addon) && mu[0] <= (offset+23./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_33(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+23./addon) && mu[0] <= (offset+24./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_34(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+24./addon) && mu[0] <= (offset+25./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_35(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+25./addon) && mu[0] <= (offset+26./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_36(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+26./addon) && mu[0] <= (offset+27./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_37(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+27./addon) && mu[0] <= (offset+28./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_38(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+28./addon) && mu[0] <= (offset+29./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_39(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+29./addon) && mu[0] <= (offset+30./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_40(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+30./addon) && mu[0] <= (offset+31./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_41(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+31./addon) && mu[0] <= (offset+32./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_42(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+32./addon) && mu[0] <= (offset+33./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_43(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+33./addon) && mu[0] <= (offset+34./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_44(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+34./addon) && mu[0] <= (offset+35./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_45(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+35./addon) && mu[0] <= (offset+36./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_46(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+36./addon) && mu[0] <= (offset+37./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_47(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+37./addon) && mu[0] <= (offset+38./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_48(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+38./addon) && mu[0] <= (offset+39./addon)){
       return 1;
    }
	return 0;
}


T theta_chi_49(const std::array<T,PDim>& mu)
{
    if(mu[0] > (offset+39./addon) && mu[0] <= (offset+40./addon)){
       return 1;
    }
	return 0;
}


T zero_fct(T /*x*/){
	return 0;
}
