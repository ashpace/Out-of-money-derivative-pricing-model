//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	Valuation.cpp
//  By : 1222781
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "Valuation.h"

#include "Input.h"
#include "Output.h"

#include "OptionBase.h"
#include "EuroCall.h"
#include "EuroPut.h"

#include "IS_process.h"

#include "Method_Base.h"
#include "IS_method.h"
#include "MC_method.h"

#include "AnalyticalFormulae.h"
#include "utility.h"

#include <stdexcept>


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Valuation::Valuation(Input & inp, Output & out)
{
    opt_ = CreateOption(inp);    
    
	expl_ = BS_Valuation(inp);
    
    is_mth_ = new IS_method(inp, out);
    pl_mth_ = new MC_method(inp, out);
    
    Banner(inp);
}

Valuation::~Valuation()
{
    delete opt_;
    delete is_mth_;
    delete pl_mth_;
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	Banner Generator
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void Valuation::Banner(Input & inp) const
{
	ut::OutputLine(">>>>>>>>>>>>>>>> Option Evaluator running <<<<<<<<<<<<<<<<<<<");
    ut::OutputLine( opt_->Identity() );
    ut::OutputLine( "Exercise Price: ", inp.GetX());
    ut::OutputLine("");
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	getters
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double Valuation::Get_expl() const {return expl_;}

double Valuation::Get_IS_SE() const {return is_mth_->GetOptionSE();}
double Valuation::Get_IS_Value()const {return is_mth_->GetOptionValue();}
double Valuation::Get_IS_Time() const {return is_mth_->GetTime();}

double Valuation::Get_PL_SE() const {return pl_mth_->GetOptionSE();}
double Valuation::Get_PL_Value() const {return pl_mth_->GetOptionValue();}
double Valuation::Get_PL_Time() const {return pl_mth_->GetTime();}

double Valuation::Get_Speed_up() const 
{
	double Pl_se = Get_PL_SE();
	double Pl_time = Get_PL_Time();
	double IS_se = Get_IS_SE();
	double IS_time = Get_IS_Time();
	return(Pl_se*Pl_se*Pl_time)/(IS_se*IS_se*IS_time);
}
        
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	run
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void Valuation::run()
{
    is_mth_->run(*opt_);
    pl_mth_->run(*opt_);
}
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	CreateOption()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

OptionBase * Valuation::CreateOption(const Input & inp)
{
    char o_type = inp.GetOptionType();
    
    switch(o_type)
    {
        case 'c': return new EuroCall(inp);  break;
        case 'p': return new EuroPut(inp);   break;
        default:  throw std::runtime_error("Valuation: CreateOption:  bad option");
    }
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	CreateOption()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double Valuation::BS_Valuation(const Input & inp)
{
    char o_type = inp.GetOptionType();
    double S = inp.GetS0();
	double r = inp.Getr();
	double sig = inp.Getsig();
	double T = inp.GetT();
	double X = inp.GetX();
    
    switch(o_type)
    {
        case 'c': return AnalyticalFormulae::European_call_GBM( S, r, sig, T, X);  break;
        case 'p': return AnalyticalFormulae::European_put_GBM( S, r, sig, T, X);   break;
        default:  throw std::runtime_error("Valuation: BS_Valuation:  bad option");
    }
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
