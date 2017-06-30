/*
Fits and analysis of Markarian 501 data. 

Energies in _MeV_. 

100 MeV to 300 GeV for Fermi
100 GeV to 1.7 TeV for Magic (low state)
*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


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

//const double pi = TMath::Pi();

const double dloge = 0.1;

const double slope = 1.78;
const double cutoff = 1E9;
const double cutslope = 2.;

//const double light = 3.E10;
//const double density = 1.;

//const double kaypi = 0.17;
const double mp = 0.938E-3;
const double mpi = 0.135E-3;
const double mpipl = 0.139E-3;
const double mmu = 0.1056E-3;
const double me = 0.511;

const double r = mmu*mmu/(mpipl*mpipl);


double icenergy[20] = {0.158489,0.251189,0.398107,
							0.630957,1,1.58489,2.51189,3.98107,6.30957,10,15.8489,
							25.1189,39.8107,63.0957,100,158.489,251.189,398.107,
							630.957,1000.}; //energy (TeV) corresponding to areas
							
double ic80[20] = {6.96E-4, 3.69E-3, 1.72E-2, 6.01E-2, 1.57E-1, 
							3.86E-1, 1.01, 2.26, 4.43, 8.14, 1.50E1, 2.50E1, 3.92E1, 
							5.95E1,7.69E1, 1.03E2, 1.33E2, 1.61E2, 1.95E2, 2.15E2}; //Taup average

double markback[80] = {6925.81, 4466.68, 2814.45, 1735.68, 1049.41, 623.042, 363.78, 209.185, 
						118.624, 66.4218, 36.767, 20.1414, 10.9308, 5.8824, 3.14181, 1.66677, 
						0.878952, 0.46104, 0.240691, 0.125131, 0.0648134, 0.0334616, 0.0172254, 
						0.00884444, 0.00453071, 0.00231605, 0.00118165, 0.000601788, 0.000305944,
						0.000155274, 7.86707e-05, 3.97892e-05, 2.00872e-05, 1.01211e-05, 
						5.08891e-06, 2.55291e-06, 1.27753e-06, 6.3757e-07, 3.17248e-07, 
						1.57349e-07, 7.77663e-08, 3.82865e-08, 1.87708e-08, 9.16118e-09, 
						4.44934e-09, 2.14959e-09, 1.03269e-09, 4.93137e-10, 2.33982e-10, 
						1.09716e-10, 5.12337e-11, 2.37495e-11, 1.0922e-11, 4.98162e-12, 
						2.2535e-12, 1.01132e-12, 4.5049e-13, 1.99314e-13, 8.76544e-14, 
						3.83466e-14, 1.67003e-14, 7.2455e-15, 3.13352e-15, 1.35162e-15, 
						5.81758e-16, 2.49957e-16, 1.07243e-16, 4.59593e-17, 1.96777e-17, 
						8.41883e-18, 3.59974e-18, 1.53845e-18, 6.57247e-19, 2.807e-19, 
						1.19854e-19, 5.11655e-20, 2.18392e-20, 9.32055e-21, 3.97745e-21, 
						1.69721e-21};//TeV-1 cm-2 s-1
					
double markeng[80] = {1e-05, 1.25893e-05, 1.58489e-05, 1.99526e-05, 2.51189e-05, 3.16228e-05, 
					3.98107e-05, 5.01187e-05, 6.30957e-05, 7.94328e-05, 0.0001, 0.000125893, 
					0.000158489, 0.000199526, 0.000251189, 0.000316228, 0.000398107, 0.000501187, 
					0.000630957, 0.000794328, 0.001, 0.00125893, 0.00158489, 0.00199526, 
					0.00251189, 0.00316228, 0.00398107, 0.00501187, 0.00630957, 0.00794328, 0.01, 
					0.0125893, 0.0158489, 0.0199526, 0.0251189, 0.0316228, 0.0398107, 0.0501187, 
					0.0630957, 0.0794328, 0.1, 0.125893, 0.158489, 0.199526, 0.251189, 0.316228,
					0.398107, 0.501187, 0.630957, 0.794328, 1, 1.25893, 1.58489, 1.99526, 2.51189, 
					3.16228, 3.98107, 5.01187, 6.30957, 7.94328, 10, 12.5893, 15.8489, 19.9526, 
					25.1189, 31.6228, 39.8107, 50.1187, 63.0957, 79.4328, 100, 125.893, 158.489, 
					199.526, 251.189, 316.228, 398.107, 501.187, 630.957, 794.328}; //TeV corresponding to markback					

double integrator(double *x, int a, int n){

	double event = 0;
	double intfact;
	
// 	for(int j=a; j<n; j++){ //extended simpson's rule 
// 		intfact = 4./3.;
// //		term = j-a; //term is odd 4/3 even 2/3
// 		if ((j-a)%2==0) intfact = 2./3.;
// 
// 		if (j==a) intfact = 1./3.; //1st and last terms
// 		if (j==n-1) intfact = 1./3.;
// 
// 		event += x[j]*intfact;	
// 	}

	for(int j=a; j<n; j++){ //riemann sum
		intfact = 1.;
		if (j==a) intfact = 0.5;
		if (j==n-1) intfact = 0.5;
		event += x[j]*intfact;
	}


	return event;
}


void backflux(double *eng, double *flux, int n){ //return E2dN/dE background

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_akima, 80);
	gsl_spline_init (spline, markeng, markback, 80); 		

	for(int i=0; i<n; i++){
		eng[i] = pow(10, -4 + dloge*i);
		flux[i] = 0;
		if (eng[i]>1E3) continue;
		if (eng[i]>=markeng[0]) flux[i] = gsl_spline_eval(spline, eng[i], acc);	
//		cout <<flux[i] <<endl;
		flux[i]*=eng[i]*eng[i];
	}
	
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
}

void eventrate(double *eng, double *flux, int n){//takes dN/dE flux

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_akima, 20);
	gsl_spline_init (spline, icenergy, ic80, 20); 		

	double area1[n], integrand1[n];
	for(int i=0; i<n; i++){
//		eng1[i] = pow(10, -1 + 0.02*i);
		area1[i] = 0;
		if (eng[i]>=icenergy[0]) area1[i] = gsl_spline_eval(spline, eng[i], acc);
		if (eng[i]>1000) area1[i] = 0;
		integrand1[i] = 3.1557E7*1E4*eng[i]*dloge*log(10)*area1[i]*flux[i];		
	}

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    
    cout <<integrator(integrand1, 0, n) <<" neutrinos/year" <<endl;

}

double pionspec(double eng){ //flux of pi+ 
//	fudge = 6.468891E14;//normalize integral to 1	
//	fudge *= 5.518E-0; //normalize integral to desired number

//	return 1.E-6*pow(eng, -slope);
	return 1.E-6*pow(eng, -slope)*exp(-1.*pow(eng/cutoff, cutslope));
}


double fluxfunc(double eng, void * p){//basic for EdN/dE

	double x = log10(eng*1E6);

//uncorrected for EBL
// 	double p0 = -0.03153; //-0.0274155
// 	double p1 = 0.4489; //0.392597
// 	double p2 = -2.365; //-2.0859
// 	double p3 = 5.725; //5.12713
// 	double p4 = -10.72; //-10.2577

//Kneiske EBL model 
	double p0 =-0.0274155;
	double p1 = 0.392597;
	double p2 =-2.0859;
	double p3 = 5.12713;
	double p4 = -10.2577;

	double flux = p0*x*x*x*x + p1*x*x*x + p2*x*x + (p3-1)*x + p4 ;
	return pow(10, flux);

}


double gaisser1(double t, void *p){ // 1 pion

	double eng = 1./t;
//	double eng = t;
	double nuscale = 29.8E-6;//MeV

	double piinteg = (0.5*mpipl/nuscale)*pionspec(eng)*pow(eng*eng-mpipl*mpipl, -0.5);

	return piinteg/(t*t);
}

double gaisser2(double* eng, size_t dim, void* p){ //1 pion

	//integrals over al(x), bet(e_pi), t(e_mu) 
	
	double al = eng[0];
	double bet = eng[1];
	double t = eng[2];
	
	double enu = *(double *)p;
//	double emu = enu + (1.-t)/t;
	double emu = 1./t;

	double epimax = emu/r - mpipl*mpipl*(1.-r)*0.25/emu ;
	double epimin = emu + mpipl*mpipl*(1.-r)*0.25/emu;

	double epi = (epimax - epimin)*bet + epimin;


	double y = enu/emu;
	double betmu = pow( 1- (1./pow(emu/mmu, 2.)) ,0.5);	
	double xmax = 1.;
	if (2.*y/(1.-betmu)<1.) xmax = 2.*y/(1.-betmu);
//	double xmax = 2.*y/(1.- betmu);	
	double xmin = 2.*y/(1. + betmu);	
	double x = (xmax-xmin)*al + xmin; 

	double pions = pionspec(epi);
	double pidec = pow(epi*epi-mpipl*mpipl, -0.5)*pow(1.-r,-1.);

	double pmu = (1./betmu)*( (2.*epi*r/(emu*(1.-r))) - ((1.+r)/(1.-r)) );
	
	if (y < 0.5*x*(1.-betmu)) return 0;
	if (y > 0.5*x*(1.+betmu)) return 0;
	
	double f0 = 2.*x*x*(3.-2.*x);
	double f1 = 2.*x*x*(1.-2.*x);
	double mudec1 = (1./(betmu*x))*(f0 - (pmu*f1*(2.*y - x))/(betmu*x) ); // dn/dydx, mu+
	double mudec2 = (1./(betmu*x))*(f0 + (pmu*f1*(2.*y - x))/(betmu*x) ); // dn/dydx, mu-
	
	return pions*pidec*(0.5*(mudec1+mudec2)/emu)*(xmax-xmin)*(epimax -epimin)*(1./(t*t));

}


double gaisserwrap(double eng){

	double result, error;	
	double flux;
		
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);		
	gsl_function G;

	double nuscale = 29.8E-6;//MeV
	
	G.function = &gaisser1;

	double lo = 0.5*mpipl*((eng/nuscale) + (nuscale/eng));

	gsl_integration_qags (&G,0., 1./lo, 1e-8, 1e-8, 5000, w, &result, &error);	
//	gsl_integration_qagiu (&G, lo, 1e-8, 1e-8, 5000, w, &result, &error);

	flux = 0.5*2.*result;	//		
	
	gsl_integration_workspace_free (w);

	gsl_monte_function F;
	
	F.f = &gaisser2;
	F.dim = 3;		

	double xl[3] = {0.,0.,0.};
//	double xu[3] = {1.,1.,1.};	
	
	const gsl_rng_type *T;
	gsl_rng *r1;
 
	size_t calls = 2000;
   
	gsl_rng_env_setup ();
 
	T = gsl_rng_default;
	r1 = gsl_rng_alloc (T);

	F.params = &eng; 	

	double xu[3] = {1.,1.,1./(eng + 0.25*mpipl*mpipl/eng)};	
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);  //state alloc and free must be inside loop

	gsl_monte_vegas_integrate (&F, xl, xu, 3, calls/10, r1, s, &result, &error);

	do{
		gsl_monte_vegas_integrate (&F, xl, xu, 3, calls, r1, s, &result, &error);		
	}
	while (fabs (s->chisq -1.0) > 0.5);	
			
	gsl_monte_vegas_free (s);
	
	flux += 0.5*2.*result;	//	multiply by 0.75 to shift to Dermer
	
	return flux*eng*eng;
}

double stecker(double t, void * p){ //stecker's pi-0 decay, takes PIONS as input, 1 pion

	double eng = 1./t;
//	double eng = t;
	double piinteg = 2.*pionspec(eng)*pow(eng*eng-mpi*mpi, -0.5);

	return piinteg/(t*t);
}


double steckerwrap(double eng){

	double result, error;	
		
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
		
	gsl_function F;

	F.function = &stecker;
	double lo;
	
	lo = eng + 0.25*mpi*mpi/eng;
	
//	gsl_integration_qagiu (&F, lo, 1e-5, 1e-5, 5000, w, &result, &error);
	gsl_integration_qags (&F,0., 1./lo, 1e-8, 1e-8, 5000, w, &result, &error);
		
	gsl_integration_workspace_free (w);

	return result*eng*eng;	//	2 pi-0's

}

double gamfluence(){

	double result, error;
		
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
		
	gsl_function F;

	F.function = &fluxfunc;
	
	gsl_integration_qags (&F, 1E-4, 1.7, 1e-8, 1e-8, 5000, w, &result, &error);
	
	gsl_integration_workspace_free (w);

	return result;
}

int main(int argc, char **argv){
	
//	cout << log10(10) <<endl;

// 	double result, error;
// 		
// 	gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
// 		
// 	gsl_function F;
// 
// 	F.function = &flux;
// 	
// 	gsl_integration_qags (&F, 1E2, 1E6, 1e-8, 1e-8, 5000, w, &result, &error);
// 	
// 	gsl_integration_workspace_free (w);
// 
// 	cout << result <<endl;

	const int n = 131;
	int i;
	double eng[n], flux[n], steckflux[n];
	double gaissflux[n];

	for(i=0; i<n; i++){
		eng[i] = pow(10, -4. + dloge*i);
		flux[i] = fluxfunc(eng[i], eng)*eng[i];
		steckflux[i] = steckerwrap(eng[i]);
		gaissflux[i] = gaisserwrap(eng[i]); //E2dN/dE
//		cout <<eng[i]<<"	"<<gaissflux[i] <<"	"<<steckflux[i] <<endl;
	}

//	double normo = flux[15]/steckflux[15];
//	cout <<normo <<endl;

	double diffgams[n], diffneuts[n];

	for(i=0; i<n; i++){
//		steckflux[i] *= normo*1E-6;
//		gaissflux[i] *=normo*1E-6;
//		flux[i]*=1E-6;
	
//		eng[i]*=1E-6;

//		diffgams[i] = steckflux[i]/(eng[i]*eng[i]);
//		diffneuts[i] = gaissflux[i]/(eng[i]*eng[i]);
//		cout <<diffneuts[i] <<endl;
	}
	
//	cout <<"Normalize at "<<eng[5] <<endl;

	cout << "Emax = " <<eng[n-1] <<" TeV" <<endl;
//	eventrate(eng, diffneuts, n);

	cout <<gamfluence() <<endl;

	double gamfluence_uncasc_integrand[n], nufluence_uncasc_integrand[n];
	
	for(i=0; i<n; i++){
	
		gamfluence_uncasc_integrand[i] = steckflux[i]*dloge*log(10);//steckflux is E2dN/dE so divide by eng[i]
		nufluence_uncasc_integrand[i] = gaissflux[i]*dloge*log(10);//gaissflux is E2dN/dE so divide by eng[i]
	
	}

	double nugamrat = integrator(nufluence_uncasc_integrand, 0, n)/integrator(gamfluence_uncasc_integrand, 0, n);
	double gamexprat = gamfluence()/integrator(gamfluence_uncasc_integrand, 0, n);
	cout <<integrator(gamfluence_uncasc_integrand, 0, n) <<endl;
	cout <<integrator(nufluence_uncasc_integrand, 0, n) <<endl;	
	
	for(i=0; i<n; i++){
		steckflux[i]*= gamexprat;
		gaissflux[i]*= gamexprat;

		diffgams[i] = steckflux[i]/(eng[i]*eng[i]);
		diffneuts[i] = gaissflux[i]/(eng[i]*eng[i]);		
	}

	const int m = 90;
	double backgr[m], diffback[m], beng[m];
	backflux(beng, backgr, m);

	for(i=0; i<m; i++){
		diffback[i] = backgr[i]/(beng[i]*beng[i]);	
	
	}


	cout <<"Signal: ";
	eventrate(eng, diffneuts, n);
	cout <<"Background: "; 
	eventrate(beng, diffback, m);

	cout <<"Ratio of neutrino fluence to gamma= "<<nugamrat <<endl;

	TApplication *app = new TApplication("app", &argc, argv); //needed to run ROOT externally


	TGraph *gammas = new TGraph(n, eng, flux);
	gammas->SetLineColor(1);
	gammas->SetLineWidth(2);
	gammas->SetLineStyle(2);

	TGraph *steckgammas = new TGraph(n, eng,steckflux);
	steckgammas->SetLineColor(1);
	steckgammas->SetLineWidth(2);
	steckgammas->SetLineStyle(1);

	TGraph *neutrinos = new TGraph(n, eng,gaissflux);
	neutrinos->SetLineColor(2);
	neutrinos->SetLineWidth(2);
	neutrinos->SetLineStyle(1);

	TGraph *background = new TGraph(m, beng,backgr);
	background->SetLineColor(2);
	background->SetLineWidth(2);
	background->SetLineStyle(2);


	TCanvas *c4 = new TCanvas("events4","stuff4",500,500 );
	TPad *pad4 = new TPad("pad4","pad4",0.0,0.0,1.0,1.0,10);

	c4->cd();

	pad4->SetTitle("");
	pad4->Draw();
	pad4->cd();
	pad4->SetLogx();
	pad4->SetLogy();
	pad4->SetGridx();
	pad4->SetGridy();
	pad4->SetTicks();

	gammas->SetMinimum(1.E-13);
	gammas->SetMaximum(4.E-9);
	gammas->SetTitle("");
	
	gammas->Draw("Al");
	steckgammas->Draw("same");	
	neutrinos->Draw("same");
	background->Draw("same");

	gammas->GetXaxis()->SetLimits(1.E-4,1.E10);
	
	app->Run();

	return 0;

}