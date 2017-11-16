//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  main.cpp
//  By : 1222781
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "Application_Wrapper.h"
#include "Errorhandler.h"

#include <stdexcept>
#include <iostream>
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  main()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

int main(int argc, char *argv[])
{
	std::cout << "\\\\\\\\\\\\\\\\\\ Importance Sampling Project //////////////////  " <<  std::endl;
	std::cout << "\\\\\\\\\\\\\\\\\\         By: 1222781         //////////////////  " <<  std::endl;
	
    ErrorHandler the_handler;
    
    try
    {
        Application_Wrapper app;    
        app.run();
    }
    catch(const std::runtime_error & e)
    {
        the_handler.HandleRunTimeError(e);
    }
    catch(...)
    {
        the_handler.HandleUnknownError();
    }
                    
    return the_handler.PauseAndReturn();
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  end of file
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

