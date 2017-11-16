//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	Application_Wrapper.cpp
//  By : 1222781
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "Application_Wrapper.h"

#include "GrandInput.h"
#include "Input.h"
#include "Output.h"

#include "Valuation.h"

#include "ErrorHandler.h"

#include <vector>
#include <stdexcept>

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Application_Wrapper::Application_Wrapper()
{

}
                                          
Application_Wrapper::~Application_Wrapper()
{

}
                                           
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	run
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void Application_Wrapper::run()
{
	ErrorHandler the_handler;
	
	GrandInput ginp;	//loading inputs from the depository
	std::vector<Input> inputs =  ginp.GetVec();
	
	double imp_size = inputs.size();
	
    for (double i = 0; i < imp_size; ++i)
	{
	    Input inp = inputs[i];      
	    Output out;
	    
	    Valuation val(inp, out);
	    
		try
	    {
    		val.run();
	    }
	    catch(const std::runtime_error & e)
	    {
	        the_handler.HandleRunTimeError(e);
	    }
	    catch(...)
	    {
	        the_handler.HandleUnknownError();
	    }
	    
        
    	out.SetOutput(val);
    	out.DoOutput();
	}

}
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
