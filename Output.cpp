//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	Output.cpp
//  By : 1222781
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "Output.h"

#include "Valuation.h"
#include "StopWatch.h"

#include "utility.h"

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	interface
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void Output::SetOutput(Valuation & val)
{
	expl_ = val.Get_expl();       //Explicit Solution
	
    is_se_ = val.Get_IS_SE();
    is_val_ = val.Get_IS_Value();
    is_t_ = val.Get_IS_Time();
    
    pl_se_ = val.Get_PL_SE();
    pl_val_ = val.Get_PL_Value();
    pl_t_ = val.Get_PL_Time();
    
    speed_up_ = val.Get_Speed_up();
}
 
void Output::DoOutput()
{
	ut::OutputLine("Explicit Solution = ", expl_);
	
	ut::OutputLine("");
	ut::OutputLine("IS Method Output");
    ut::OutputLine("Option value  ", is_val_);
    ut::OutputLine("se            ", is_se_);
    ut::OutputLine("Time taken    ", is_t_);
	ut::OutputLine("");
	
	ut::OutputLine("Plain Method Output");
    ut::OutputLine("Option value  ", pl_val_);
    ut::OutputLine("se            ", pl_se_);
    ut::OutputLine("Time taken    ", pl_t_);
    ut::OutputLine("");
    
    ut::OutputLine("speed up is   ", speed_up_);
}

void Output::OutputBanner(std::string strg)
{
    ut::OutputLine(strg);
}

void Output::OutputCounter(long j, long M, long interval)
{
    ut::OutputCounter(j, M, interval);
}
  
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
