/*

Test using gamma function to do low-e pion integral. Compare to 
series expansion. 

All results are quick and return the same value.

Mathcore is fastest. (but all O(1E-5) s)
Large delta, small tee is undefined for spec. function. GSL works best.
For positive lamb1 spec function is undefined. 

-->use series rather than special function. 
*/

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>


#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TObject.h"
#include "TROOT.h"
#include "TFile.h"
#include "TPad.h"
#include "TLegend.h"
#include "Math/SpecFunc.h"
#include "TStopwatch.h"

using namespace std;


const double lamb1 = -0.171;
const double delta = 100.;
const double tee = 0.54;


double seriesint(double x){

	double series = 0;
	double factorialrun = 1.;

	for(int i=1; i<100; i++){	
		series += pow(lamb1, i-1)*pow(x, i)*pow(factorialrun, -1.)/(delta+i);
		factorialrun *= i; 	
	}

	return series;


}

double mathregint(double x){

//	return (-1./lamb1)*pow(-1./(lamb1*x), delta)*TMath::Gamma(delta+1)*TMath::Gamma(delta+1., -1.*lamb1*x);
	return TMath::Gamma(delta+1., -1.*lamb1*x);
}

double mathcoreint(double x){
//	return (-1./lamb1)*pow(-1./(lamb1*x), delta)*ROOT::Math::tgamma(delta+1)*ROOT::Math::inc_gamma(delta+1., -1.*lamb1*x);
	return ROOT::Math::inc_gamma(delta+1., -1.*lamb1*x);
}

double gnufuncint(double x){
	return (-1./lamb1)*pow(-1./(lamb1*x), delta)*gsl_sf_gamma(delta+1.)*gsl_sf_gamma_inc_P(delta+1., -1*lamb1*x);
}

int main(){

 	TStopwatch time;

 	time. Start();
	double seriesexp = seriesint(tee);	
	time.Stop();
	cout << seriesexp <<endl;
 	cout <<time.RealTime() <<" seconds"<<endl;
 	
 	time. Start(1); 	
	double mathreg = mathregint(tee);
	time.Stop();
	cout << mathreg <<endl;
 	cout <<time.RealTime() <<" seconds"<<endl;

 	time. Start(1); 	
	double mathcore = mathcoreint(tee);
	time.Stop();
	cout << mathcore <<endl;
 	cout <<time.RealTime() <<" seconds"<<endl;

 	time. Start(1); 	
	double gnufunc = gnufuncint(tee);
	time.Stop();
	cout << gnufunc <<endl;
	cout <<time.RealTime() <<" seconds"<<endl;
	
	return 0;
}