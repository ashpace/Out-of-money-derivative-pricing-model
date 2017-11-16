//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	MC_method.cpp
//  By : 1222781
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "MC_method.h"

#include "Input.h"
#include "Output.h"
#include "Accumulator.h"
#include "OptionBase.h"
#include "GBM_process.h"
#include "StopWatch.h"

#include <vector>
#include <stdexcept>

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

MC_method::MC_method(Input & inp, Output & out)
:   M_(inp.GetM())
{
    acc_ = new Accumulator(inp);
    
    out_ = &out;
    
    inp_ = &inp;
    
    run_status_ = false;
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	getters
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double MC_method::GetOptionValue() const 
{
	if (run_status_ == false) throw std::runtime_error("MC_method:  value requested before running");
    return acc_->GetOptionValue();
}

double MC_method::GetOptionSE() const 
{
	if (run_status_ == false) throw std::runtime_error("MC_method:  SE requested before running");
    return acc_->GetSE();
}

double MC_method::GetTime() const 
{
	if (run_status_ == false) throw std::runtime_error("MC_method:  Time requested before running");
    return t_;
}
                           
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	run()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void MC_method::run(OptionBase & opt) const
{
    if (run_status_ == true) throw std::runtime_error("MC_method:  invalid run request, method has been ran and is matured");
    
    StopWatch stw;
    stw.StartStopWatch();
    
    GBM_process prc(* inp_);
    
    for(long j = 1; j <= M_; ++j)
    {
        
        double Next_S = prc.Next_S();
                
        double payoff = opt.ComputePO(Next_S);
        
        acc_->AddValue(payoff);       
    }
    
    t_ = stw.GetTime();
    
    run_status_ = true;
}
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
