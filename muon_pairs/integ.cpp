#include <stdlib.h>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>


#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TObject.h"
#include "TROOT.h"
#include "TFile.h"
#include "TPad.h"
#include "TLegend.h"
 
 using namespace std;
 
 struct my_params{double a; double b;};
 
 
 double func(double x, void * p){

//	double pa = *(double *) p; //this works for one number. don't understand. 

	struct my_params *f_params = (struct my_params *)p;
	
	double pa = f_params->a;
	double pb = f_params->b;
	
//	double pb = *(double *) (++p);
	
	return pow(x, 2.)*pa + pb;
}


int main(){


	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

	double result, error;

	double para[2] = {10., 1.};

//	struct my_params para = {10., 1.};
	
	gsl_function F;
	F.function = &func;
	F.params = &para;
	
	gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000, w, &result, &error);

	cout <<result <<endl;
	
	gsl_integration_workspace_free (w);

	return 0;
}