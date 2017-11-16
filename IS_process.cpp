//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	IS_process.cpp
//  By : 1222781
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "IS_process.h"

#include "OptionBase.h"
#include "Input.h"
#include "rv_library.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

IS_process::IS_process(Input & inp, OptionBase & opt)
:   S_0_(inp.GetS0())
,   sig_(inp.Getsig())
,   r_(inp.Getr())
,	T_(inp.GetT())
,	X_(inp.GetX())
{
	mu_T_ = std::log(S_0_) + (r_ - 0.5*sig_*sig_)*T_;
    sig_T_ = sig_*std::sqrt(T_);
    opt_ = &opt;
    g_alpha_ = 2;
    
    double In_mode = Find_mode();
	
	double Mode_difference = (In_mode - X_) * opt.IdentityN();
    g_beta_ = (g_alpha_ - 1.0)/Mode_difference;           
    //g_beta_ = 0.02;    
    //g_beta_ = 0.2;   
}


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	Next_P() 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double IS_process::Next_P(double random) 
{
	double S = X_ + random * opt_->IdentityN();  // generates random gamma variate
    return hf(S)/g_S(S);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	mode_target()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double IS_process::mode_target(double S)
{
    return (S - X_)*(std::log(S) - mu_T_) - X_*sig_T_*sig_T_; 
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	Find_mode()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double IS_process::Find_mode()
{
	double O_IdentityN = opt_->IdentityN();

	static const long MAX_ITERATIONS = 100;
	static const double ROOT_TOLERANCE = 0.000000001;

    double first_S = X_;

    double target = 0.0;
    
    if(mode_target(first_S) >= 0) throw std::runtime_error("Bad root finding");
    
    double second_S = X_;        
    do
    {
        second_S += X_ * O_IdentityN;
    }
    while(mode_target(second_S) <= 0); 
        
    long i = 1;
          
    while (( (second_S - first_S) * O_IdentityN > ROOT_TOLERANCE) && (i < MAX_ITERATIONS))
    {    
        double root = 0.5*(second_S + first_S);
        
        double this_value = mode_target(root) ;
        
        if (this_value == target) return root;
        
        if (this_value < target) first_S = root;
        if (this_value > target) second_S = root;
        
        ++i;
	}
    
    return 0.5*(second_S + first_S);    
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	f_S() - Probability density function
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double IS_process::f_S(double S) const
{
    double nS = (std::log(S) - mu_T_)/sig_T_; 
    double w = S*sig_T_*rv::ROOT_TWO_PI;
    
    return std::exp(-0.5*nS*nS)/w;
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	hf() - 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double IS_process::hf(double S) const
{
	double p = f_S(S);						//probability of S
	double payoff = opt_ -> ComputePO(S);
    return p*payoff;
}
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	g_s() - option based
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double IS_process::g_S(double S) const
{
	double s = (S - X_) * opt_ -> IdentityN();
    return g_beta_*std::pow(g_beta_*s, g_alpha_ - 1)*std::exp(-g_beta_*s)/std::exp(rv::gammln(g_alpha_));
}
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
