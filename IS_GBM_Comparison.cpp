//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	IS_GBM_Comparison.cpp
//  Student ID Number: 1721545
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "IS_GBM_Comparison.h"

#include "Input.h"
#include "Output.h"
#include "Valuation.h"
#include "utility.h"

#include "StopWatch.h"
#include <stdexcept>

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

IS_GBM_Comparison::IS_GBM_Comparison(double X, char o_type)
{

	GBM_inp_ = new Input(X, o_type, 'p');      // owns this pointer
    IS_inp_ = new Input(X, o_type, 'i');        // owns this pointer
 
    GBM_out_ = new Output;     // owns this pointer
    IS_out_ = new Output;     // owns this pointer
    
    stw_ = new StopWatch;

    GBM_val_ = new Valuation(*GBM_inp_, *GBM_out_);
    IS_val_ = new Valuation(*IS_inp_, *IS_out_);
}
                                          
IS_GBM_Comparison::~IS_GBM_Comparison()
{
   delete GBM_inp_;
   delete IS_inp_;
   delete GBM_out_;
   delete IS_out_;
   delete stw_;
   delete GBM_val_;
   delete IS_val_;
}
                                           
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	run
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void IS_GBM_Comparison::run()
{
	ut::OutputLine("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
	ut::OutputLine("*****	Comparison Running 	****");
	ut::OutputLine("Exercise price is : ", IS_inp_->GetX());
	
	if (IS_inp_->GetOptionType() == 'c') ut::OutputLine("Option is a Call");
	else if (IS_inp_->GetOptionType() == 'p') ut::OutputLine("Option is a Put");
	else throw std::runtime_error("CreateOption:  bad option");
	
	ut::OutputLine("");
	
    GBM_out_->OutputBanner("Plain Option Evaluator Running");
    
    stw_->StartStopWatch();
    
    GBM_val_->run();
        
    GBM_out_->SetOutput(*GBM_val_, *stw_);
    GBM_out_->DoOutput();
    
    double GBM_time = (GBM_out_->GetTime());
    double GBM_se = GBM_out_->GetSe();
    
	GBM_out_->OutputBanner("");
    
    stw_->Reset();
    
	IS_out_->OutputBanner("IS Option Evaluator Running");
	IS_out_->OutputBanner("");
    
    stw_->StartStopWatch();

    IS_val_->run();
        
    IS_out_->SetOutput(*IS_val_, *stw_);
    IS_out_->DoOutput();
    
    double IS_time = (IS_out_->GetTime());
    double IS_se = IS_out_->GetSe();
    
    double speed_up = (GBM_se*GBM_se*GBM_time)/(IS_se*IS_se*IS_time);
    
    ut::OutputLine("Speed up was:	",speed_up );
}
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
