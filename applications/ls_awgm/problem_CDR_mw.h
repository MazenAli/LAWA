#include <iostream>
#include <iomanip>
#include <utility>
#include <cmath>
#include <array>

#include <lawa/lawa.h>


using namespace std;
using namespace lawa;

//===============================================================//
//========= TYPEDEFS  =======================//
//===============================================================//

//==== General ====//
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

//==== Basis 1D & 2D ====//
typedef Basis<T, Primal, Periodic, CDF>	                            TrialBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>						    TestBasis_Time;
typedef Basis<T,Orthogonal,Interval,Multi>							Basis_Space;
typedef TrialBasis_Time::RefinementBasis                            TrialRefBasis_Time;
typedef TestBasis_Time::RefinementBasis                           	TestRefBasis_Time;
typedef Basis_Space::RefinementBasis                                RefBasis_Space;

typedef TensorBasis2D<Adaptive,TrialBasis_Time,Basis_Space>         Basis2D_Trial;
typedef TensorBasis2D<Adaptive,TestBasis_Time,Basis_Space>          Basis2D_Test;

//==== Adaptive Operators ====//
typedef AdaptiveConvectionOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>                   Convection1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>                     Identity1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,Basis_Space,Basis_Space>                            Identity1D_Space;
typedef AdaptiveLaplaceOperator1D_PG<T,Basis_Space,Basis_Space>                             Laplace1D_Space;
typedef AdaptiveConvectionOperator1D_PG<T,Basis_Space,Basis_Space>                         	Convection1D_Space;

typedef TransposedAdaptiveConvectionOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>         TranspConvection1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>           TranspIdentity1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,Basis_Space,Basis_Space>                  TranspIdentity1D_Space;
typedef TransposedAdaptiveLaplaceOperator1D_PG<T,Basis_Space,Basis_Space>                   TranspLaplace1D_Space;
typedef TransposedAdaptiveConvectionOperator1D_PG<T,Basis_Space,Basis_Space>                TranspConvection1D_Space;

typedef AdaptiveConvectionOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>             RefConvection1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>               RefIdentity1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,RefBasis_Space,RefBasis_Space>                      RefIdentity1D_Space;
typedef AdaptiveLaplaceOperator1D_PG<T,RefBasis_Space,RefBasis_Space>                       RefLaplace1D_Space;
typedef AdaptiveConvectionOperator1D_PG<T,RefBasis_Space,RefBasis_Space>                   	RefConvection1D_Space;

typedef TransposedAdaptiveConvectionOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>   RefTranspConvection1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>     RefTranspIdentity1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,RefBasis_Space,RefBasis_Space>            RefTranspIdentity1D_Space;
typedef TransposedAdaptiveLaplaceOperator1D_PG<T,RefBasis_Space,RefBasis_Space>             RefTranspLaplace1D_Space;
typedef TransposedAdaptiveConvectionOperator1D_PG<T,RefBasis_Space,RefBasis_Space>          RefTranspConvection1D_Space;

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

typedef LocalOperator2D<LOp_Conv1D_Time, LOp_Id1D_Space>                        LOp_Conv_Id_2D;
typedef LocalOperator2D<LOp_Id1D_Time, LOp_Id1D_Space>                        	LOp_Id_Id_2D;
typedef LocalOperator2D<LOp_Id1D_Time, LOp_Lapl1D_Space>                        LOp_Id_Lapl_2D;
typedef LocalOperator2D<LOp_Id1D_Time, LOp_Conv1D_Space>                        LOp_Id_Conv_2D;
typedef LocalOperator2D<LOpT_Conv1D_Time, LOpT_Id1D_Space>                      LOpT_Conv_Id_2D;
typedef LocalOperator2D<LOpT_Id1D_Time, LOpT_Id1D_Space>                      	LOpT_Id_Id_2D;
typedef LocalOperator2D<LOpT_Id1D_Time, LOpT_Lapl1D_Space>                      LOpT_Id_Lapl_2D;
typedef LocalOperator2D<LOpT_Id1D_Time, LOpT_Conv1D_Space>                      LOpT_Id_Conv_2D;


//==== CompoundOperators ====//
typedef FlexibleCompoundLocalOperator<Index2D,AbstractLocalOperator2D<T> > 		Flex_COp_2D;
typedef CompoundLocalOperator<Index2D,LOp_Conv_Id_2D,LOp_Id_Lapl_2D,
                                        LOp_Id_Conv_2D, LOp_Id_Id_2D>    		Comp_COp_2D;

const size_t PDim = 2;
typedef array<T,PDim>	 													ParamType;
typedef AffineLocalOperator<Index2D,AbstractLocalOperator2D<T>,ParamType>	Affine_Op_2D;

//==== Preconditioners ====//
typedef AdaptiveLeftNormPreconditioner2D<T,Basis2D_Test>            LeftPrec2D;
typedef AdaptiveRightNormPreconditioner2D_c<T,Basis2D_Trial>        RightPrec2D;
typedef NoPreconditioner<T, Index2D>								NoPrec2D;

//==== RightHandSides ====//
typedef SmoothRHSWithAlignedSing2D<T, Basis2D_Test, FullGridGL>		RhsIntegral2D;
typedef RHS<T,Index2D,RhsIntegral2D, NoPrec2D>						AdaptiveRhs2D;


//==== Solver ====//

typedef MultiTreeAWGM_PG<Index2D,Basis2D_Trial, Basis2D_Test,Affine_Op_2D,
		Affine_Op_2D,AdaptiveRhs2D,RightPrec2D,LeftPrec2D>					MT_AWGM;
//typedef MultiTreeAWGM_PG<Index2D,Basis2D_Trial, Basis2D_Test,Flex_COp_2D,
//		Flex_COp_2D,SumOfSeparableRhs,RightPrec2D,LeftPrec2D>				MT_AWGM;

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, string filename);

void
readSolutionFromFile(Coefficients<Lexicographical,T,Index2D>& u, string filename);


// Curve parameters
const T 	displ = 1.;
const T 	center = 0.5;
const T 	ampl = 0.25;
const T 	damping = 1000.;

const T 	mu_lapl = 1.;
const T 	mu_conv = 1.;
const T 	mu_reac = 1.;

T c(T t)
{
    return center + ampl*std::sin(2.*M_PI*t);
}

T curve(T t, T x)
{
	T tmp1 = (x - displ*c(t));
	T tmp2 = tmp1 * tmp1;
    return std::exp(- damping * tmp2) * (tmp2 * (- 4. * mu_lapl*damping*damping)
    									+ tmp1 * (4.*M_PI*displ*ampl*damping*std::cos(2.*M_PI*t) - 2.*damping*mu_conv)
    									+ 2.*mu_lapl*damping + mu_reac);
}

T dummy(T, T)
{
    return 0;
}

T sol(T t, T x){
	return std::exp(- damping * (x - displ*c(t)) * (x - displ*c(t)));
}

T no_theta(const std::array<T,PDim>& /*mu*/)
{
	return 1.;
}

T theta_lapl(const std::array<T,PDim>& /*mu*/)
{
	return mu_lapl;
}

T theta_conv(const std::array<T,PDim>& /*mu*/)
{
	return mu_conv;
}

T theta_reac(const std::array<T,PDim>& /*mu*/)
{
	return mu_reac;
}
