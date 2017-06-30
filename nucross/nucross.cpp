/*
Neutrino CC cross-section.

Antineutrino CC cross-section. 


*/

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <stdlib.h>
#include <fstream>

#include "LHAPDF/LHAPDF.h"

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

int quark = 2;
int aquark = -2;
 
const int SUBSET = 0;
const string NAME = "CTEQ61";

const double mw = 80.4; //GeV
const double mn = 0.938;
const double gf = 1.166E-5; //GeV-2
const double hc = 0.389379*1E-27; //GeV2 cm2

double pi = TMath::Pi();

double nucross(double *x, size_t dim, void* p){

	double ex = x[0];
	double y = x[1];
	double enu = *(double *)p;

	double Q = pow(2.*mn*ex*y*enu, 0.5);	
	
	double factor1 = hc*gf*gf*mn*enu/pi;
	double factor2 = pow(mw, 4.) / pow( 2.*mn*ex*y*enu + mw*mw, 2.);
	
	double down = LHAPDF::xfx(ex, Q, 1);
	double up = LHAPDF::xfx(ex, Q, 2);
	double strange = LHAPDF::xfx(ex, Q, 3);
	
	double adown = LHAPDF::xfx(ex, Q, -1);
	double aup = LHAPDF::xfx(ex, Q, -2);

	double factor3 = up + down + 2.*strange + pow(1.-y, 2.)*(aup + adown);
	
	double integrand = factor1*factor2*factor3;

 	return integrand;
}

double anucross(double *x, size_t dim, void* p){

	double ex = x[0];
	double y = x[1];
	double enu = *(double *)p;

	double Q = pow(2.*mn*ex*y*enu, 0.5);	
	
	double factor1 = hc*gf*gf*mn*enu/pi;
	double factor2 = pow(mw, 4.) / pow( 2.*mn*ex*y*enu + mw*mw, 2.);
	
	double down = LHAPDF::xfx(ex, Q, 1);
	double up = LHAPDF::xfx(ex, Q, 2);
	
	double adown = LHAPDF::xfx(ex, Q, -1);
	double aup = LHAPDF::xfx(ex, Q, -2);
	double astrange = LHAPDF::xfx(ex, Q, -3);

	double factor3 = pow(1.-y, 2.)*(up + down) + (aup + adown + 2.*astrange);
	
	double integrand = factor1*factor2*factor3;

 	return integrand;
}

int main(int argc, char **argv){

//	double enu = 1.E1;
	int n = 30;
	double nueng[n], cross[n];

	LHAPDF::initPDFSet(NAME, LHAPDF::LHPDF, SUBSET);
		
	double xlo = LHAPDF::getXmin(SUBSET);	
	double xhi = LHAPDF::getXmax(SUBSET);
					
	double result, error;	

//------Monte Carlo-------//	

	for(int i=0; i<n; i++){
		nueng[i] = pow(10., 0. + 0.25*i);

		double par1 = nueng[i];

		gsl_monte_function F;
		
		F.f = &anucross;
		F.dim = 2;
		F.params = &par1; 		
		
		double xl[2] = {xlo, 0};
		double xu[2] = {xhi, 1};
	
		
		const gsl_rng_type *T;
		gsl_rng *r;
	 
		size_t calls = 50000;
	   
		gsl_rng_env_setup ();
	 
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
	
		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);  //state alloc and free must be inside loop
	
		gsl_monte_vegas_integrate (&F, xl, xu, 2, calls/10, r, s, &result, &error);
	
		do{
			gsl_monte_vegas_integrate (&F, xl, xu, 2, calls, r, s, &result, &error);		
		}
		while (fabs (s->chisq -1.0) > 0.5);
				
		gsl_monte_vegas_free (s);
	
//		cout <<result <<endl;
		cross[i] = result;	
	}
	
	
	TApplication *app = new TApplication("app", &argc, argv); //needed to run ROOT externally

	TGraph *gr1 = new TGraph(n, nueng, cross);
	gr1->SetLineWidth(2);
	gr1->SetLineColor(2);
	gr1->SetLineStyle(1);

	TCanvas *c1 = new TCanvas("c1","stuff2",500,500 );
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,10);

	c1->cd();

	pad1->SetTitle("");
	pad1->Draw();
	pad1->cd();
	pad1->SetLogx();
	pad1->SetLogy();

	gr1->SetMaximum(1E-31);
	gr1->SetMinimum(1E-39);
	gr1->SetTitle("anu");

	gr1->Draw("Al");

	gr1->GetXaxis()->SetLimits(1., 1E7);
	gr1->GetXaxis()->SetTitle("E (GeV)");
	gr1->GetXaxis()->CenterTitle(1);
	gr1->GetXaxis()->SetTitleOffset(1.2);
	gr1->GetYaxis()->SetTitle("Total Cross-section (cm2)");
	gr1->GetYaxis()->CenterTitle(1);
	gr1->GetYaxis()->SetTitleOffset(1.2);

	TLegend *leg = new TLegend(0.45, 0.7, 0.90, 0.89, "", "brNDC");
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
//		leg->SetShadowColor("1");
	leg->AddEntry(gr1, "Integrand", "l");
//		leg->Draw();

	app->Run();

	
	return 0;
}
