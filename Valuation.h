//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  Valuation.h
//  By : 1222781
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#ifndef ValuationH
#define ValuationH
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class OptionBase;
class Method_Base;
class IS_process;

class IS_MC_method;
class Input;
class Output;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  class Valuation
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class Valuation
{
    public:
        Valuation(Input & inp, Output & out);
        ~Valuation();
        
        double Get_expl() const; 

		double Get_IS_SE() const; 
		double Get_IS_Value()const; 
		double Get_IS_Time() const; 

		double Get_PL_SE() const; 
		double Get_PL_Value() const; 
		double Get_PL_Time() const;
		
		double Get_Speed_up() const;
        
        void run();
        
    private:
       OptionBase * opt_;

       Method_Base * is_mth_;
       Method_Base * pl_mth_;
       
       OptionBase * CreateOption(const Input & inp);
       
       double BS_Valuation(const Input & inp);
       double expl_;
       
       void Banner(Input & inp) const;
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
