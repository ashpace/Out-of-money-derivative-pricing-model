//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	AnalyticalFormulae.h
//  By : 1222781
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#ifndef AnalyticalFormulaeH
#define AnalyticalFormulaeH
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	name space AnalyticalFormulae
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

namespace AnalyticalFormulae
{
    double Lookback_call_fixed_GBM(double S, double r, double sig, double T, double X, double M);
    double Lookback_put_fixed_GBM(double S, double r, double sig, double T, double X, double M);
    
    double Lookback_call_floating_GBM(double S, double r, double sig, double T, double M);
    double Lookback_put_floating_GBM(double S, double r, double sig, double T, double M);
       
    double European_call_GBM(double S, double r, double sig, double T, double X);    
    double European_put_GBM(double S, double r, double sig, double T, double X);    

    double Barrier_UIC_GBM(double S, double r, double sig, double T, double X, double K);
    double Barrier_UOC_GBM(double S, double r, double sig, double T, double X, double K);
    double Barrier_DIC_GBM(double S, double r, double sig, double T, double X, double K);
    double Barrier_DOC_GBM(double S, double r, double sig, double T, double X, double K);
    double Barrier_DIP_GBM(double S, double r, double sig, double T, double X, double K);
    double Barrier_UIP_GBM(double S, double r, double sig, double T, double X, double K);
    double Barrier_DOP_GBM(double S, double r, double sig, double T, double X, double K);
    double Barrier_UOP_GBM(double S, double r, double sig, double T, double X, double K);

	double AverageGeometricContinuousCall_GBM(double S, double r, double sig, double T, double X);
	double AverageGeometricContinuousPut_GBM(double S, double r, double sig, double T, double X);
	
	double AverageGeometricDiscreteCall_GBM(double S, double r, double sig, double T, double X, long N);
	double AverageGeometricDiscretePut_GBM(double S, double r, double sig, double T, double X, long N);
	double AverageGeometricDiscreteCallPart_GBM(double S, double r, double sig, double X, double tau, double dtau, long N, long M, double A);
	double AverageGeometricDiscretePutPart_GBM(double S, double r, double sig, double X, double tau, double dtau, long N, long M, double A);
}

namespace af = AnalyticalFormulae;    //alias

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

