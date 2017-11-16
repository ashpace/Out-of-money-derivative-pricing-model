//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  AnalyticalFormulae.cpp
//  By : 1222781
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "AnalyticalFormulae.h"

#include "rv_library.h"

#include <cmath>
#include <iostream>

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	global constants
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

const double ROOT_THREE_INVERSE = 0.577350269189626;     // 1/sqr(3)
const double ROOT_PI_INVERSE = 0.564189583547756;        // 1/sqrt(PI)
const double ROOT_TWO_INVERSE = 0.707106781186547;       // 1/sqrt(2)

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	Lookback_call_fixed_GBM()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Lookback_call_fixed_GBM(double S, double r, double sig, double T, double X, double M)
{
    double sgrt = sig * std::sqrt(T);
    double disc	= std::exp(-r * T);
    double PvX = disc * X;
    double PvM = disc * M;
    double b = r / (sig * sig);
    
    double d1 = std::log(S / PvX) / sgrt + 0.5 * sgrt;
    double d2 = d1 - sgrt;
    double d3 = d1 - 2. * sgrt * b;
    
    double a1 = std::log(S / PvM) / sgrt + 0.5 * sgrt;
    double a2 = a1 - sgrt;
    double a3 = a1 - 2. * sgrt * b;
    
    double gamM = disc * pow(M / S, 2. * b);
    double gamX = disc * pow(X / S, 2. * b);
    
    
    if (X > M)
        return S * rv::normal_cdf(d1) - PvX * rv::normal_cdf(d2) 
                     + 0.5 * S * (rv::normal_cdf(d1) - gamX * rv::normal_cdf(d3)) / b;
    else
        return (PvM - PvX) + S * rv::normal_cdf(a1) - PvM * rv::normal_cdf(a2) 
                        + 0.5 * S * (rv::normal_cdf(a1) - gamM * rv::normal_cdf(a3)) / b;
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	Lookback_put_fixed_GBM()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Lookback_put_fixed_GBM(double S, double r, double sig, double T, double X, double M)
{
    double sgrt = sig * std::sqrt(T);
    double disc	= std::exp(-r * T);
    double PvX = disc * X;
    double PvM = disc * M;
    double b = r / (sig * sig);
    
    double d1 = std::log(S / PvX) / sgrt + 0.5 * sgrt;
    double d2 = d1 - sgrt;
    double d3 = d1 - 2. * sgrt * b;
    
    double a1 = std::log(S / PvM) / sgrt + 0.5 * sgrt;
    double a2 = a1 - sgrt;
    double a3 = a1 - 2. * sgrt * b;
    
    double gamM = disc * pow(M / S, 2. * b);
    double gamX = disc * pow(X / S, 2. * b);
    
    if (X < M)
        return PvX * rv::normal_cdf(-d2) - S * rv::normal_cdf(-d1)  
                     + 0.5 * S * (gamX * rv::normal_cdf(-d3) - rv::normal_cdf(-d1)) / b;
    else
        return (PvX - PvM) + PvM * rv::normal_cdf(-a2) - S * rv::normal_cdf(-a1) 
                        + 0.5 * S * (gamM * rv::normal_cdf(-a3) - rv::normal_cdf(-a1)) / b;
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     Lookback_call_floating_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Lookback_call_floating_GBM(double S, double r, double sig, double T, double M)
{
    double sgrt = sig * std::sqrt(T);
    double disc	= std::exp(-r * T);
    double PvM = disc * M;
    double b = r / (sig * sig);
    
    double a1 = std::log(S / PvM) / sgrt + 0.5 * sgrt;
    double a2 = a1 - sgrt;
    double a3 = a1 - 2. * sgrt * b;
    
    double gamM = disc * pow(M / S, 2. * b);
    
    return S * rv::normal_cdf(a1) - PvM * rv::normal_cdf(a2) 
              - 0.5 * S * (rv::normal_cdf(-a1) - gamM * rv::normal_cdf(-a3)) / b;
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     Lookback_put_floating_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Lookback_put_floating_GBM(double S, double r, double sig, double T, double M)
{
    double sgrt = sig * std::sqrt(T);
    double disc	= std::exp(-r * T);
    double PvM = disc * M;
    double b = r / (sig * sig);
    
    double a1 = std::log(S / PvM) / sgrt + 0.5 * sgrt;
    double a2 = a1 - sgrt;
    double a3 = a1 - 2. * sgrt * b;
    
    double gamM = disc * pow(M / S, 2. * b);
    
    return PvM * rv::normal_cdf(-a2) - S * rv::normal_cdf(-a1)  
              + 0.5 * S * (rv::normal_cdf(a1) - gamM * rv::normal_cdf(a3)) / b;
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     European_call_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::European_call_GBM(double S, double r, double sig, double T, double X)
{
    double sgrt = sig * std::sqrt(T);
    double PvX = std::exp(-r * T) * X;
    
    double d1 = std::log(S / PvX) / sgrt + 0.5 * sgrt;
    double d2 = d1 - sgrt;
    
    return S * rv::normal_cdf(d1) - PvX * rv::normal_cdf(d2);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     European_put_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::European_put_GBM(double S, double r, double sig, double T, double X)
{
    double sgrt = sig * std::sqrt(T);
    double PvX = std::exp(-r * T) * X;
    
    double d1 = std::log(S / PvX) / sgrt + 0.5 * sgrt;
    double d2 = d1 - sgrt;
    
    return PvX * rv::normal_cdf(-d2) - S * rv::normal_cdf(-d1);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     Barrier_UIC_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Barrier_UIC_GBM(double S, double r, double sig, double T, double X, double K)
{
    double sgrt = sig*std::sqrt(T);
    double disc = std::exp(-r*T);
    double EofS = S/disc;
    
    double d1 = std::log(EofS/X)/sgrt + 0.5*sgrt;    double d2 = d1 - sgrt;	
	double N_d1 = rv::normal_cdf(-d1);	double N_d2 = rv::normal_cdf(-d2);

    if (X >= K) return disc*(EofS*(1 - N_d1) -  X*(1 - N_d2));

    double arg1 = S*S/K;
    double arg2 = X*S*S/(K*K);
    
    double a1 = std::log(EofS/K)/sgrt + 0.5*sgrt;    double a2 = a1 - sgrt;
    double b1 = std::log(EofS/arg1)/sgrt + 0.5*sgrt; double b2 = b1 - sgrt;    
    double c1 = std::log(EofS/arg2)/sgrt + 0.5*sgrt; double c2 = c1 - sgrt;
		
	double N_a1 = rv::normal_cdf(-a1);	double N_a2 = rv::normal_cdf(-a2);
	double N_b1 = rv::normal_cdf(-b1);	double N_b2 = rv::normal_cdf(-b2);
	double N_c1 = rv::normal_cdf(-c1);	double N_c2 = rv::normal_cdf(-c2);


	double nu = 2.*r/(sig*sig);	  
	double f1 = std::pow(K/S, nu + 1.);
	double f2 = std::pow(K/S, nu - 1.);
		
    double N_1 = EofS * N_d1;
    double N_2 = N_d2;
    
    double DK = EofS * (N_a1 - f1 *  N_b1);
    double DX = EofS * (N_d1 - f1 *  N_c1);

    double FK = N_a2 - f2 * N_b2;
    double FX = N_d2 - f2 * N_c2;

    return disc*((EofS - N_1 - DK + DX) -  X*(1 - N_2 - FK + FX));
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     Barrier_UOC_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Barrier_UOC_GBM(double S, double r, double sig, double T, double X, double K)
{
    if (S >= K) return 0.0;

    if (X >= K) return 0.0;

    double sgrt = sig*std::sqrt(T);
    double disc = std::exp(-r*T);
    double EofS = S/disc;
    
    double arg1 = S*S/K;
    double arg2 = X*S*S/(K*K);
    
    double d1 = std::log(EofS/X)/sgrt + 0.5*sgrt;    double d2 = d1 - sgrt;	
    double a1 = std::log(EofS/K)/sgrt + 0.5*sgrt;    double a2 = a1 - sgrt;
    double b1 = std::log(EofS/arg1)/sgrt + 0.5*sgrt; double b2 = b1 - sgrt;    
    double c1 = std::log(EofS/arg2)/sgrt + 0.5*sgrt; double c2 = c1 - sgrt;
		
	double N_d1 = rv::normal_cdf(-d1);	double N_d2 = rv::normal_cdf(-d2);
	double N_a1 = rv::normal_cdf(-a1);	double N_a2 = rv::normal_cdf(-a2);
	double N_b1 = rv::normal_cdf(-b1);	double N_b2 = rv::normal_cdf(-b2);
	double N_c1 = rv::normal_cdf(-c1);	double N_c2 = rv::normal_cdf(-c2);

	double nu = 2.*r/(sig*sig);	   
	double f1 = std::pow(K/S, nu + 1.);
	double f2 = std::pow(K/S, nu - 1.);
		
    double DK = EofS * (N_a1 - f1 *  N_b1);
    double DX = EofS * (N_d1 - f1 *  N_c1);

    double FK = N_a2 - f2 * N_b2;
    double FX = N_d2 - f2 * N_c2;

    return disc*(DK - DX -  X*(FK - FX));
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     Barrier_DIC_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Barrier_DIC_GBM(double S, double r, double sig, double T, double X, double K)
{
    double sgrt = sig*std::sqrt(T);
    double disc = std::exp(-r*T);
    double EofS = S/disc;
    
    double arg1 = S*S/K;
    double arg2 = X*S*S/(K*K);
    
    double d1 = std::log(EofS/X)/sgrt + 0.5*sgrt;    double d2 = d1 - sgrt;
	double N_d1 = rv::normal_cdf(-d1);	double N_d2 = rv::normal_cdf(-d2);

    if (S <= K) return disc*(EofS*(1. - N_d1) - X*(1. - N_d2));

    double c1 = std::log(EofS/arg2)/sgrt + 0.5*sgrt; double c2 = c1 - sgrt;
	double N_c1 = rv::normal_cdf(-c1);	double N_c2 = rv::normal_cdf(-c2);

	double nu = 2.*r/(sig*sig);	   
	double f1 = std::pow(K/S, nu + 1.);
	double f2 = std::pow(K/S, nu - 1.);

    if (X > K) return disc*(EofS*f1*(1. - N_c1) - X*f2*(1. - N_c2));

    double a1 = std::log(EofS/K)/sgrt + 0.5*sgrt;    double a2 = a1 - sgrt;
    double b1 = std::log(EofS/arg1)/sgrt + 0.5*sgrt; double b2 = b1 - sgrt;    
		
	double N_a1 = rv::normal_cdf(-a1);	double N_a2 = rv::normal_cdf(-a2);
	double N_b1 = rv::normal_cdf(-b1);	double N_b2 = rv::normal_cdf(-b2);

		
    double DX = (X <= K) ?  N_d2 : N_a2 + f2 * (N_c2 - N_b2);
    double EDX = N_a2 + f2 * (1. - N_b2);
    
    double FX = (X <= K) ? EofS * N_d1 : EofS * (N_a1 + f1 * (N_c1 - N_b1));
    double EFX = EofS * (N_a1 + f1 * (1. - N_b1));

    return disc*(EFX - FX -  X*(EDX - DX));
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     Barrier_D0C_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Barrier_DOC_GBM(double S, double r, double sig, double T, double X, double K)
{
    if (S <= K) return 0.0;

    if (X <= K) return 0.0;

    double sgrt = sig*std::sqrt(T);
    double disc = std::exp(-r*T);
    double EofS = S/disc;
    
    double arg2 = X*S*S/(K*K);
    
	double nu = 2.*r/(sig*sig);	   
	double f1 = std::pow(K/S, nu + 1.);
	double f2 = std::pow(K/S, nu - 1.);

    double d1 = std::log(EofS/X)/sgrt + 0.5*sgrt;    double d2 = d1 - sgrt;	
    double c1 = std::log(EofS/arg2)/sgrt + 0.5*sgrt; double c2 = c1 - sgrt;
		
	double N_d1 = rv::normal_cdf(-d1);	double N_d2 = rv::normal_cdf(-d2);
	double N_c1 = rv::normal_cdf(-c1);	double N_c2 = rv::normal_cdf(-c2);
		
    return disc*(EofS*(1. - N_d1) - X*(1. - N_d2))
              - disc*(EofS * f1 * (1. - N_c1) -  X * f2 * (1. - N_c2));
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     Barrier_DIP_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Barrier_DIP_GBM(double S, double r, double sig, double T, double X, double K)
{
    double sgrt = sig*std::sqrt(T);
    double disc = std::exp(-r*T);
    double EofS = S/disc;

    double arg1 = S*S/K;
    double arg2 = X*S*S/(K*K);

	double nu = 2.*r/(sig*sig);
	double f1 = std::pow(K/S, nu + 1.);
	double f2 = std::pow(K/S, nu - 1.);

    double d1 = std::log(EofS/X)/sgrt + 0.5*sgrt;    double d2 = d1 - sgrt;
	double N_d1 = rv::normal_cdf(-d1);	double N_d2 = rv::normal_cdf(-d2);
	
	double p0 = X*N_d2 - EofS*N_d1;

    if (S <= K) return disc*p0;
    
    if (X < K) return disc*p0;

    double a1 = std::log(EofS/K)/sgrt + 0.5*sgrt;    double a2 = a1 - sgrt;
    double b1 = std::log(EofS/arg1)/sgrt + 0.5*sgrt; double b2 = b1 - sgrt;
    double c1 = std::log(EofS/arg2)/sgrt + 0.5*sgrt; double c2 = c1 - sgrt;

	double N_a1 = rv::normal_cdf(-a1);	double N_a2 = rv::normal_cdf(-a2);
	double N_b1 = rv::normal_cdf(-b1);	double N_b2 = rv::normal_cdf(-b2);
	double N_c1 = rv::normal_cdf(-c1);	double N_c2 = rv::normal_cdf(-c2);

    double p1 = X*N_a2 - EofS*N_a1;
    double p2 = X*f2*(1. - N_c2) - EofS*f1*(1. - N_c1);
    double p3 = X*f2*(1. - N_b2) - EofS*f1*(1. - N_b1);

    return disc*(p1 - p2 + p3);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     Barrier_UIP_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Barrier_UIP_GBM(double S, double r, double sig, double T, double X, double K)
{
    double sgrt = sig*std::sqrt(T);
    double disc = std::exp(-r*T);
    double EofS = S/disc;

    double arg1 = S*S/K;
    double arg2 = X*S*S/(K*K);

	double nu = 2.*r/(sig*sig);
	double f1 = std::pow(K/S, nu + 1.);
	double f2 = std::pow(K/S, nu - 1.);

    double d1 = std::log(EofS/X)/sgrt + 0.5*sgrt;    double d2 = d1 - sgrt;
	double N_d1 = rv::normal_cdf(-d1);	double N_d2 = rv::normal_cdf(-d2);
    double p1 = X*N_d2 - EofS*N_d1;

    if (S >= K) return disc*p1;
    
    double c1 = std::log(EofS/arg2)/sgrt + 0.5*sgrt; double c2 = c1 - sgrt;
	double N_c1 = rv::normal_cdf(-c1);	double N_c2 = rv::normal_cdf(-c2);

    if (X < K) return disc*(X*f2*N_c2 - EofS*f1*N_c1);

    double a1 = std::log(EofS/K)/sgrt + 0.5*sgrt;    double a2 = a1 - sgrt;
    double b1 = std::log(EofS/arg1)/sgrt + 0.5*sgrt; double b2 = b1 - sgrt;

	double N_a1 = rv::normal_cdf(-a1);	double N_a2 = rv::normal_cdf(-a2);
	double N_b1 = rv::normal_cdf(-b1);	double N_b2 = rv::normal_cdf(-b2);

    double p2 = X*N_a2 - EofS*N_a1;
    double p3 = X*f2*N_b2 - EofS*f1*N_b1;

    return disc*(p1 - p2 + p3);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     Barrier_DOP_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Barrier_DOP_GBM(double S, double r, double sig, double T, double X, double K)
{
    if (S <= K) return 0.0;

    double sgrt = sig*std::sqrt(T);
    double disc = std::exp(-r*T);
    double EofS = S/disc;

    double arg1 = S*S/K;
    double arg2 = X*S*S/(K*K);

	double nu = 2.*r/(sig*sig);
	double f1 = std::pow(K/S, nu + 1.);
	double f2 = std::pow(K/S, nu - 1.);

    if (X < K) return 0.0;

    double d1 = std::log(EofS/X)/sgrt + 0.5*sgrt;    double d2 = d1 - sgrt;
    double a1 = std::log(EofS/K)/sgrt + 0.5*sgrt;    double a2 = a1 - sgrt;
    double b1 = std::log(EofS/arg1)/sgrt + 0.5*sgrt; double b2 = b1 - sgrt;
    double c1 = std::log(EofS/arg2)/sgrt + 0.5*sgrt; double c2 = c1 - sgrt;

	double N_d1 = rv::normal_cdf(-d1);	double N_d2 = rv::normal_cdf(-d2);
	double N_a1 = rv::normal_cdf(-a1);	double N_a2 = rv::normal_cdf(-a2);
	double N_b1 = rv::normal_cdf(-b1);	double N_b2 = rv::normal_cdf(-b2);
	double N_c1 = rv::normal_cdf(-c1);	double N_c2 = rv::normal_cdf(-c2);

    double p1 = X*N_d2 - EofS*N_d1;
    double p2 = X*N_a2 - EofS*N_a1;
    double p3 = X*f2*(1. - N_c2) - EofS*f1*(1. - N_c1);
    double p4 = X*f2*(1. - N_b2) - EofS*f1*(1. - N_b1);

    return disc*(p1 - p2 + p3 - p4);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     Barrier_UOP_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::Barrier_UOP_GBM(double S, double r, double sig, double T, double X, double K)
{
    if (S >= K) return 0.0;

    double sgrt = sig*std::sqrt(T);
    double disc = std::exp(-r*T);
    double EofS = S/disc;

    double arg1 = S*S/K;
    double arg2 = X*S*S/(K*K);

	double nu = 2.*r/(sig*sig);
	double f1 = std::pow(K/S, nu + 1.);
	double f2 = std::pow(K/S, nu - 1.);

    if (X < K)
    {
        double d1 = std::log(EofS/X)/sgrt + 0.5*sgrt;    double d2 = d1 - sgrt;
        double c1 = std::log(EofS/arg2)/sgrt + 0.5*sgrt; double c2 = c1 - sgrt;

	    double N_d1 = rv::normal_cdf(-d1);	double N_d2 = rv::normal_cdf(-d2);
	    double N_c1 = rv::normal_cdf(-c1);	double N_c2 = rv::normal_cdf(-c2);
	    
        double p1 = X*N_d2 - EofS*N_d1;
        double p3 = X*f2*N_c2 - EofS*f1*N_c1;
        
	    return disc*(p1 - p3);
    }

    double a1 = std::log(EofS/K)/sgrt + 0.5*sgrt;    double a2 = a1 - sgrt;
    double b1 = std::log(EofS/arg1)/sgrt + 0.5*sgrt; double b2 = b1 - sgrt;

	double N_a1 = rv::normal_cdf(-a1);	double N_a2 = rv::normal_cdf(-a2);
	double N_b1 = rv::normal_cdf(-b1);	double N_b2 = rv::normal_cdf(-b2);

    double p2 = X*N_a2 - EofS*N_a1;
    double p4 = X*f2*N_b2 - EofS*f1*N_b1;

    return disc*(p2 - p4);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     AverageGeometricContinuousCall_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::AverageGeometricContinuousCall_GBM(double S, double r, double sig, double T, double X)
{
    double sgrt = sig* ROOT_THREE_INVERSE * std::sqrt(T);
    double div = 0.5 * (r + sig*sig/6.);
    double PvS = std::exp(-r * T) * S;
    double PvX = std::exp(-div * T) * X;
    
    double d1 = std::log(PvS / PvX) / sgrt + 0.5 * sgrt;
    double d2 = d1 - sgrt;
    
    return PvS * rv::normal_cdf(d1) - PvX * rv::normal_cdf(d2);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     AverageGeometricContinuousPut_GBM
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::AverageGeometricContinuousPut_GBM(double S, double r, double sig, double T, double X)
{
    double sgrt = sig* ROOT_THREE_INVERSE * std::sqrt(T);
    double div = 0.5 * (r + sig*sig/6.);
    double PvS = std::exp(-r * T) * S;
    double PvX = std::exp(-div * T) * X;
    
    double d1 = std::log(PvS / PvX) / sgrt + 0.5 * sgrt;
    double d2 = d1 - sgrt;
    
    return PvX * rv::normal_cdf(-d2) - PvS * rv::normal_cdf(-d1);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     AverageGeometricDiscreteCall_GBM
//     No previous average. current time is 0
// 	   N resets at t1,...,tN. t(i+1) - t(i) = dt (t(0) = 0)
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::AverageGeometricDiscreteCall_GBM(double S, double r, double sig, double T, double X, long N)
{
    double dt = T/N;
    double drift = 0.5*(r - 0.5*sig*sig)*(T + dt);
    double sgrt = sig*std::sqrt(dt*(1. + (N - 1.)*(2.*N - 1.)/(6.*N)));
    
    double disc = std::exp(-r * T);
    double S_adj = S*std::exp(drift + 0.5*sgrt*sgrt);
   
    double d1 = (std::log(S/X) + drift)/sgrt + sgrt;
    double d2 = d1 - sgrt;
    
    return disc*(S_adj*rv::normal_cdf(d1) - X*rv::normal_cdf(d2));
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     AverageGeometricDiscretePut_GBM
//     No previous average. current time is 0
// 	   N resets at t1,...,tN. t(i+1) - t(i) = dt (t(0) = 0)
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::AverageGeometricDiscretePut_GBM(double S, double r, double sig, double T, double X, long N)
{
    double dt = T/N;
    double drift = 0.5*(r - 0.5*sig*sig)*(T + dt);
    double sgrt = sig*std::sqrt(dt*(1. + (N - 1.)*(2.*N - 1.)/(6.*N)));
    
    double disc = std::exp(-r * T);
    double S_adj = S*std::exp(drift + 0.5*sgrt*sgrt);
   
    double d1 = (std::log(S/X) + drift)/sgrt + sgrt;
    double d2 = d1 - sgrt;
    
    return disc*(X*rv::normal_cdf(-d2) - S_adj*rv::normal_cdf(-d1));
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     AverageGeometricDiscreteCallPart_GBM
//     Part-way option with an average to date
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::AverageGeometricDiscreteCallPart_GBM
(
    double S,		// current asset value  
	double r, 		// riskless rate
	double sig, 	// sigma
	double X, 		// strike
	double tau, 	// time to maturity
	double dtau, 	// time to next reset date
	long N,			// total number of reset dates
	long M,			// number of previous reset dates
	double A		// average to date,  A = (S1*...*SM)^(1/M)
)
{
	long RN = N - M;									// number of remaining reset dates
	double dt = (tau - dtau)/(RN - 1);					// interval between reset dates
    //double tm = dt - dtau;           					// time of last reset date
    
    double q = double(RN)/N;							// proportion of resets remaining
    double p = (A == 0) ? 1. : pow(A, 1. - q);
    
    double Tbar = dtau + 0.5*(RN - 1.)*dt;
    
    double top = RN*RN*dtau + dt*RN*(RN - 1.)*(2.*RN - 1.)/6.;
    double bot = RN*RN*Tbar;
    
    double sigbar = std::sqrt(sig*sig*top/bot);
    double div = 0.5*(sig*sig - sigbar*sigbar);
    double divq = q*div + (r + 0.5*q*sigbar*sigbar)*(1. - q);
    double disc = std::exp(r*(Tbar - tau));

	double PvS = std::exp(-divq*Tbar)*pow(S, q);
	double PvX = std::exp(-r*Tbar)*X/p;
	double sgrt = sigbar*q*std::sqrt(Tbar);
   
    double d1 = std::log(PvS/PvX)/sgrt + 0.5*sgrt;
    double d2 = d1 - sgrt;
    
    return p*disc*(PvS*rv::normal_cdf(d1) - PvX*rv::normal_cdf(d2));
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//     AverageGeometricDiscretePutPart_GBM
//     Part-way option with an average to date
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double AnalyticalFormulae::AverageGeometricDiscretePutPart_GBM
(
    double S,		// current asset value  
	double r, 		// riskless rate
	double sig, 	// sigma
	double X, 		// strike
	double tau, 	// time to maturity
	double dtau, 	// time to next reset date
	long N,			// total number of reset dates
	long M,			// number of previous reset dates
	double A		// average to date,  A = (S1*...*SM)^(1/M)
)
{
	long RN = N - M;									// number of remaining reset dates
	double dt = (tau - dtau)/(RN - 1);					// interval between reset dates
    //double tm = dt - dtau;           					// time of last reset date
    
    double q = double(RN)/N;							// proportion of resets remaining
    double p = (A == 0) ? 1. : pow(A, 1. - q);
    
    double Tbar = dtau + 0.5*(RN - 1.)*dt;
    
    double top = RN*RN*dtau + dt*RN*(RN - 1.)*(2.*RN - 1.)/6.;
    double bot = RN*RN*Tbar;
    
    double sigbar = std::sqrt(sig*sig*top/bot);
    double div = 0.5*(sig*sig - sigbar*sigbar);
    double divq = q*div + (r + 0.5*q*sigbar*sigbar)*(1. - q);
    double disc = std::exp(r*(Tbar - tau));

	double PvS = std::exp(-divq*Tbar)*pow(S, q);
	double PvX = std::exp(-r*Tbar)*X/p;
	double sgrt = sigbar*q*std::sqrt(Tbar);
   
    double d1 = std::log(PvS/PvX)/sgrt + 0.5*sgrt;
    double d2 = d1 - sgrt;
    
    return p*disc*(PvX*rv::normal_cdf(-d2) - PvS*rv::normal_cdf(-d1));
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

























