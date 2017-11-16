//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	IS_method.cpp
//  By : 1222781
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "IS_method.h"

#include "Input.h"
#include "Output.h"

#include "Accumulator.h"

#include "IS_process.h"
#include "StopWatch.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

IS_method::IS_method(Input & inp, Output & out)
:   M_(inp.GetM())
{
    acc_ = new Accumulator(inp);
    
    out_ = &out;
    
    inp_ = &inp;
    
    run_status_ = false;
    
    t_ = 0;
}


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	getters
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double IS_method::GetOptionValue() const 
{
	if (run_status_ == false) throw std::runtime_error("IS_method:  Option Value requested before running");
    return acc_->GetOptionValue();
}

double IS_method::GetOptionSE() const 
{
	if (run_status_ == false) throw std::runtime_error("IS_method:  SE requested before running");
    return acc_->GetSE();
}

double IS_method::GetTime() const 
{
	if (run_status_ == false) throw std::runtime_error("IS_method:  Time requested before running");
    return t_;
}
                           
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	run()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void IS_method::run(OptionBase & opt) const
{
    if (run_status_ == true) throw std::runtime_error("IS_method:  invalid run request, method has been ran and is matured");
    
	StopWatch stw;
    stw.StartStopWatch();
    
	IS_process prc(* inp_, opt);
    
    std::default_random_engine uni_gen;
	
	double alpha = prc.Get_g_alpha();
	double beta = prc.Get_g_beta();
	
    std::gamma_distribution<double> gamma(alpha, 1.0/beta);
    
    
    
    for(long j = 1; j <= M_; ++j)
    {
		
		double random_gamma = gamma(uni_gen);
               
        double payoff = prc.Next_P(random_gamma);
        
        acc_->AddValue(payoff);       
    }
    
    t_ = stw.GetTime();
    run_status_ = true;
}
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
