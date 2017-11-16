//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  IS_process.h
//  By : 1222781
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#ifndef IS_processH
#define IS_processH
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class Input;
class OptionBase;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  class GBM_process
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class IS_process 
{
	
    public:
        explicit IS_process(Input & inp, OptionBase & opt);
        
        double Next_P(double ga_rand) ;
        
        double Get_mu_T() const {return g_alpha_;}
        double Get_sig_T() const {return g_beta_;}
        double Get_g_alpha() const {return g_alpha_;}
        double Get_g_beta() const {return g_beta_;}
        
        
        
    private:         
        double mu_T_;
        double sig_T_; 
        
		double S_0_;
        double sig_;      
        double r_;
        double T_;
        double X_;
        
        OptionBase * opt_;
        
        double g_alpha_;
        double g_beta_;
        
        double f_S(double S) const;
        double Find_mode(); 
        double mode_target(double S);
        double hf(double S) const;
		double g_S(double S) const;

};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
