/*
P-function and its derivatives

P(E_nu, E_mu>emin) : pfunction(enu, emu)
dP/dE_mu (E_nu, E_mu): diffpfunction(enu, emu)
d/dE_nu(dP/dE_mu) (E_nu, E_mu): ddpfunction(enu, emu)

implement second derivative d2/dE2_nu(dP/dE_mu) ? maybe not

to test integral equation solver, give form for A_mu(E_mu), get
A_nu(E_nu): neutrinoarea(enu, enu*) (2nd argument is to get right signature for derivative)
dA_nu/dE_nu(E_nu): derivneutrinoarea(enu, enu*) (slow since 4 3D integrals per point)

*/


#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <stdlib.h>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>


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
#include "TStopwatch.h"

using namespace std;
 
const int SUBSET = 0;
const string NAME = "CTEQ61";

const double mw = 80.4; //GeV
const double mn = 0.938;
const double gf = 1.166E-5; //GeV-2
const double hc = 0.389379*1E-27; //GeV2 cm2
const double avogadro = 6.022E23; 
const double alpha = 2.E-3; //GeV g-1 cm2
const double eps = 476.; //GeV

const double pi = TMath::Pi();

const double emin = 1.E2; //detector threshold in GeV

struct my_params2{double a; int b;};

double pfuncint(double *x, size_t dim, void* p){ 

//P(E_nu, E_thresh) integrand, call from pfunction

	double ex = x[0];
	double phi = x[1];
	double emu = x[2];
	
	double enu = *(double *)p;

	double ymax = 1. - emu/enu;
	double ymin = 0.;

	double y = (ymax - ymin)*phi + ymin;

	double Q = pow(2.*mn*ex*y*enu, 0.5);	
	
	double factor1 = hc*gf*gf*mn*enu/pi;
	double factor2 = pow(mw, 4.) / pow( 2.*mn*ex*y*enu + mw*mw, 2.);
	
	double down = LHAPDF::xfx(ex, Q, 1);
	double up = LHAPDF::xfx(ex, Q, 2);
	double strange = LHAPDF::xfx(ex, Q, 3);
	
	double adown = LHAPDF::xfx(ex, Q, -1);
	double aup = LHAPDF::xfx(ex, Q, -2);

	double factor3 = up + down + 2.*strange + pow(1.-y, 2.)*(aup + adown);
	
	double crosssec = factor1*factor2*factor3;

	double range = avogadro/(alpha + alpha*emu/eps);

 	return crosssec*range*(ymax - ymin);
}

double apfuncint(double *x, size_t dim, void* p){

	double ex = x[0];
	double phi = x[1];
	double emu = x[2];
	
	double enu = *(double *)p;

	double ymax = 1. - emu/enu;
	double ymin = 0.;

	double y = (ymax - ymin)*phi + ymin;

	double Q = pow(2.*mn*ex*y*enu, 0.5);	
	
	double factor1 = hc*gf*gf*mn*enu/pi;
	double factor2 = pow(mw, 4.) / pow( 2.*mn*ex*y*enu + mw*mw, 2.);
	
	double down = LHAPDF::xfx(ex, Q, 1);
	double up = LHAPDF::xfx(ex, Q, 2);
	
	double adown = LHAPDF::xfx(ex, Q, -1);
	double aup = LHAPDF::xfx(ex, Q, -2);
	double astrange = LHAPDF::xfx(ex, Q, -3);

	double factor3 = pow(1.-y, 2.)*(up + down) + (aup + adown + 2.*astrange);
	
	double crosssec = factor1*factor2*factor3;

	double range = avogadro/(alpha + alpha*emu/eps);

 	return crosssec*range*(ymax - ymin);
}

double pfunction(double x){
//call from main, P(E_nu, E_mu>Emin), integrate either pfuncint or apfuncint

		double nueng = x;
		
		double xlo = LHAPDF::getXmin(SUBSET);	
		double xhi = LHAPDF::getXmax(SUBSET);

		double result, error;	

		gsl_monte_function F;
		
		F.f = &pfuncint;
		F.dim = 3;
		F.params = &nueng; 		
		
		double xl[3] = {xlo, 0, emin};
		double xu[3] = {xhi, 1, nueng};
			
		const gsl_rng_type *T;
		gsl_rng *r;
	 
		size_t calls = 5000;
	   
		gsl_rng_env_setup ();
	 
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
	
		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
	
		gsl_monte_vegas_integrate (&F, xl, xu, 3, calls/10, r, s, &result, &error);
	
		do{
			gsl_monte_vegas_integrate (&F, xl, xu, 3, calls, r, s, &result, &error);		
		}
		while (fabs (s->chisq -1.0) > 0.5);
				
		gsl_monte_vegas_free (s);

		return result;	
}


double diffpfuncint(double *x, size_t dim, void* p){ 
//call from diffpfunction, dP(E_nu, E_mu)/dE_mu integrand


	double ex = x[0];
	double phi = x[1];
	
//	double enu = *(double *)p;

	struct my_params2 *params = ( struct my_params2 *) p;
	double enu = params->a;
	double emu = params->b;


	double ymax = 1. - emu/enu;
	double ymin = 0.;

	double y = (ymax - ymin)*phi + ymin;

	double Q = pow(2.*mn*ex*y*enu, 0.5);	
	
	double factor1 = hc*gf*gf*mn*enu/pi;
	double factor2 = pow(mw, 4.) / pow( 2.*mn*ex*y*enu + mw*mw, 2.);
	
	double down = LHAPDF::xfx(ex, Q, 1);
	double up = LHAPDF::xfx(ex, Q, 2);
	double strange = LHAPDF::xfx(ex, Q, 3);
	
	double adown = LHAPDF::xfx(ex, Q, -1);
	double aup = LHAPDF::xfx(ex, Q, -2);

	double factor3 = up + down + 2.*strange + pow(1.-y, 2.)*(aup + adown);
	
	double crosssec = factor1*factor2*factor3;

	double range = avogadro/(alpha + alpha*emu/eps);

 	return crosssec*range*(ymax - ymin);
}



double diffpfunction(double x, double z){
//call from main, get dP/dE_mu(E_nu, E_mu)

		double nueng = x;
		double mueng = z;

		if (mueng>nueng) return 0;

		struct my_params2 par1 = {nueng, mueng};		
		
		double xlo = LHAPDF::getXmin(SUBSET);	
		double xhi = LHAPDF::getXmax(SUBSET);

		double result, error;	

		gsl_monte_function F;
		
		F.f = &diffpfuncint;
		F.dim = 2;
		F.params = &par1; 		
		
		double xl[2] = {xlo, 0};
		double xu[2] = {xhi, 1};
			
		const gsl_rng_type *T;
		gsl_rng *r;
	 
		size_t calls = 5000;
	   
		gsl_rng_env_setup ();
	 
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
	
		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
	
		gsl_monte_vegas_integrate (&F, xl, xu, 2, calls/10, r, s, &result, &error);
	
		do{
			gsl_monte_vegas_integrate (&F, xl, xu, 2, calls, r, s, &result, &error);		
		}
		while (fabs (s->chisq -1.0) > 0.5);
				
		gsl_monte_vegas_free (s);

		return result;	
}

double derivdiffpfunction(double x, void *params){
//wrapper for derivative of dP/dE_mu, called from ddpfunction

	double emu = *(double *)params;

	return diffpfunction(x,emu);
}

double ddpfunction(double x, double z){ 
//call from main, E_nu derivative of dP/dE_mu at nueng, mueng

		double nueng = x;
		double mueng = z;

		double result, abserr;		
		gsl_function F;
		F.function = &derivdiffpfunction;
		F.params = &mueng;
		gsl_deriv_forward(&F, nueng, 1e-3, &result, &abserr);//forward b/c E_nu>E_mu

		return result;
}

double neutrinoareaint(double *x, size_t dim, void* p){
//call from neutrinoarea, returns A_mu(E_mu)*dP/dE_mu(E_nu, E_mu)

//	double mueng = x;
//	double nueng = *(double *)p;


	double ex = x[0];
	double phi = x[1];
	double emu = x[2];
	
	double muonarea = 0.25*1.E6*(1+log10(emu)); //in m2, 1 km2 at 1 TeV
	if (emu<emin) muonarea = 0;//can't happen but include anyway
	
	double enu = *(double *)p;

	double ymax = 1. - emu/enu;
	double ymin = 0.;

	double y = (ymax - ymin)*phi + ymin;

	double Q = pow(2.*mn*ex*y*enu, 0.5);	
	
	double factor1 = hc*gf*gf*mn*enu/pi;
	double factor2 = pow(mw, 4.) / pow( 2.*mn*ex*y*enu + mw*mw, 2.);
	
	double down = LHAPDF::xfx(ex, Q, 1);
	double up = LHAPDF::xfx(ex, Q, 2);
	double strange = LHAPDF::xfx(ex, Q, 3);
	
	double adown = LHAPDF::xfx(ex, Q, -1);
	double aup = LHAPDF::xfx(ex, Q, -2);

	double factor3 = up + down + 2.*strange + pow(1.-y, 2.)*(aup + adown);
	
	double crosssec = factor1*factor2*factor3;

	double range = avogadro/(alpha + alpha*emu/eps);

 	return crosssec*range*(ymax - ymin)*muonarea;
}

double neutrinoarea(double nueng, void *p){
//call from main, calls neutrinoareaint, gets A_nu(E_nu)
		
		double xlo = LHAPDF::getXmin(SUBSET);	
		double xhi = LHAPDF::getXmax(SUBSET);

		double result, error;	

		gsl_monte_function F;
		
		F.f = &neutrinoareaint;
		F.dim = 3;
		F.params = &nueng; 		
		
		double xl[3] = {xlo, 0, emin};
		double xu[3] = {xhi, 1, nueng};
			
		const gsl_rng_type *T;
		gsl_rng *r;
	 
		size_t calls = 7500;
	   
		gsl_rng_env_setup ();
	 
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
	
		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
	
		gsl_monte_vegas_integrate (&F, xl, xu, 3, calls/10, r, s, &result, &error);
	
		do{
			gsl_monte_vegas_integrate (&F, xl, xu, 3, calls, r, s, &result, &error);		
		}
		while (fabs (s->chisq -1.0) > 0.5);
				
		gsl_monte_vegas_free (s);

		return result;	
}

double derivneutrinoarea(double nueng, void* p){
//call from main, calls neutrinoarea, gets dA_nu/dE_nu(E_nu)

//		double nueng = x;
//		double mueng = z;

		double result, abserr;		
		gsl_function F;
		F.function = &neutrinoarea;
//		F.params = &mueng;
		gsl_deriv_central(&F, nueng, 1e-5, &result, &abserr);//forward b/c E_nu>E_mu

		return result;
}

int main(int argc, char **argv){

 	TStopwatch time;
 	time. Start();

	int n = 15;
	double de = 0.2;


	double eng[n];
//	double pfunc[n];

	LHAPDF::initPDFSet(NAME, LHAPDF::LHPDF, SUBSET);

	for(int i=0; i<n; i++){
		eng[i] = pow(10., log10(emin)+ 0.001 + de*i); //either nu or mu depending on function

//		pfunc[i] = pfunction(eng[i]);
//		pfunc[i] = diffpfunction(1.E5, eng[i]); //(neutrino energy, muon energy)
//		pfunc[i] = diffpfunction(eng[i], 1.E2); //

//		pfunc[i] = ddpfunction(eng[i], eng[i]); //enu=emu boundary
//		cout <<pfunc[i] <<endl;

//		pfunc[i] = neutrinoarea(eng[i], eng); //2nd argument irrelevant
//		pfunc[i] = derivneutrinoarea(eng[i], eng);


//		cout <<eng[i] <<"	"<<pfunc[i]<<endl;
	}
	

	double realamu[n], amu[n];	
	double anuprim, sum;
	
	realamu[0] = 0.25*1.E6*(1+log10(eng[0]));	
//	amu[0] = realamu[0];
	
	double result, abserr;		
	gsl_function F;
	F.function = &derivneutrinoarea;
	gsl_deriv_central(&F, eng[0], 1e-5, &result, &abserr);	
	
	double doublederivarea = result;
	amu[0] = doublederivarea/ddpfunction(emin, emin);
	
	for(int i=1; i<n; i++){
		realamu[i] = 0.25*1.E6*(1+log10(eng[i]));
		anuprim = derivneutrinoarea(eng[i], eng);	
		sum = 0;
		
		for(int j=1; j<i; j++){ //check
			if (i==1) break;
			sum += amu[j]*ddpfunction(eng[i], eng[j])*de*log(10.)*eng[j];		
		}
	
		amu[i] = (anuprim - 0.5*amu[0]*ddpfunction(eng[i], emin)*de*log(10.)*eng[0] - sum)/(0.5*de*log(10.)*eng[i]*ddpfunction(eng[i], eng[i]));
		cout <<	amu[i] <<"	"<<realamu[i] <<endl;
	}

// 	ifstream infile;
// 	infile.open("p_mu_to_nu.txt");
// 	
// 	double x, y, z;
// 	double ploteng[88], nuprob[88], anuprob[88];
// 	int lines = 0;
// 	
// 	while(!infile.eof()){
// 		infile >> x  >> y >> z;
// //		cout << x << " " << y << endl;	
// //		infile >> x;
// //		cout << x <<endl;
// //		if(!infile.good() ) break;
// 		ploteng[lines] = x;
// //		cout << x <<endl;
// 		nuprob[lines] = pow(10, y);
// 		anuprob[lines] = pow(10, z);
// //		cout << y <<endl;
// 		lines = lines+1;
// 		}
// 	
// 	for(Int_t i=0; i<44; i++){
// //		cout << array2[i]  <<endl;
// 		}
// 
// 	Double_t nu1[40], nu2[46], anu1[40], anu2[46], ploteng1[40], ploteng2[46];
// 
// 	for(Int_t i=0; i<40; i++){
// 		ploteng1[i] = ploteng[i];
// 		nu1[i] = nuprob[i];
// 		anu1[i] = anuprob[i];
// 	}
// 
// 	for(Int_t i=0; i<46; i++){
// 		ploteng2[i] = ploteng[i+40];
// 		nu2[i] = nuprob[i+40];
// 		anu2[i] = anuprob[i+40];
// 	}

 	time.Stop();
 	cout <<time.RealTime() <<" seconds"<<endl;
	cout.precision(8);

	
	TApplication *app = new TApplication("app", &argc, argv); //needed to run ROOT externally

	TGraph *gr1 = new TGraph(n, eng, amu);
	gr1->SetLineWidth(2);
	gr1->SetLineColor(2);
	gr1->SetLineStyle(1);
	
	TGraph *real = new TGraph(n, eng, realamu);
	real->SetLineWidth(2);
	real->SetLineColor(1);
	real->SetLineStyle(2);

// 	TGraph *nup1 = new TGraph(40, ploteng1, nu1);
// 	nup1->SetLineWidth(2);
// 	nup1->SetLineColor(1);
// 	
// 	TGraph *anup1 = new TGraph(40, ploteng1, anu1);
// 	anup1->SetLineWidth(2);
// 	anup1->SetLineColor(1);
// 	anup1->SetLineStyle(2);
// 
// 	TGraph *nup2 = new TGraph(46, ploteng2, nu2);
// 	nup2->SetLineWidth(2);
// 	nup2->SetLineColor(1);
// 	
// 	TGraph *anup2 = new TGraph(46, ploteng2, anu2);
// 	anup2->SetLineWidth(2);
// 	anup2->SetLineColor(1);
// 	anup2->SetLineStyle(2);




	TCanvas *c1 = new TCanvas("c1","stuff2",500,500 );
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,10);

	c1->cd();

	pad1->SetTitle("");
	pad1->Draw();
	pad1->cd();
	pad1->SetLogx();
	pad1->SetLogy();

	gr1->SetMaximum(1E7);
	gr1->SetMinimum(1E5);
	gr1->SetTitle("pfunc");

	gr1->Draw("Al*");
	real->Draw("same");

	gr1->GetXaxis()->SetLimits(1.E2, 1E5);
	gr1->GetXaxis()->SetTitle("E (GeV)");
	gr1->GetXaxis()->CenterTitle(1);
	gr1->GetXaxis()->SetTitleOffset(1.2);
	gr1->GetYaxis()->SetTitle("P(E, >Emin)");
	gr1->GetYaxis()->CenterTitle(1);
	gr1->GetYaxis()->SetTitleOffset(1.2);
	
//	nup1->Draw("same");
//	nup2->Draw("same");
//	anup1->Draw("same");
//	anup2->Draw("same");

	TLegend *leg = new TLegend(0.45, 0.7, 0.90, 0.89, "", "brNDC");
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
//		leg->SetShadowColor("1");
	leg->AddEntry(gr1, "Integrand", "l");
//		leg->Draw();

	app->Run();

	
	return 0;
}
