/*
Find average energy fraction of higher-energy muon from pair.

i.e. find average splitting. 
*/
#include <stdlib.h>
#include <iostream>
#include <fstream>
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

const double zee = 7.; //7.2 is air avg
const double ay = 14.; //14.7 is air avg
	
const double me = 0.511E-3; //electron mass in GeV
const double mu = 0.105;
const double alpha = 1./(137.036);
const double classrad = 1.859E-30; //classical radius of muon squared in cm2
//const double classrad = 1.859E-30*pow(mu/me, 2.);

const double zeered = pow(zee*alpha, 2.);
const double eff = 1.202*zeered - 1.0369*zeered*zeered + 1.008*pow(zeered, 3.)/(1.+zeered);

const double m = mu;

const double sqrte = 1.6487;

double nu = alpha*zee; 
double coulombcorr = nu*nu*(pow(1. + nu*nu, -1.) + 0.202 - 0.0369*pow(nu, 2.) + 0.0083*pow(nu, 4.) - 0.002*pow(nu, 6.));

// struct my_params1{double a;};
 
 struct my_params2{double a; double b;};

double geanttot(double x, void *p){

	double kay  = *(double *)p;
	

	double xplus = x;
	double xminus = 1.-xplus;
	
	double factor = 4.*alpha*classrad*zee;

	double delta = m*m/(2.*kay*xplus*xminus);
	double deen = 1.54*pow(ay, 0.27);
	double deenprim = pow(deen, 1.-(1./zee));
	double bee = 183.;
	double beeprim = 1429.;
	double sqrte = 1.648;
	
	double formn = log(bee*pow(zee, -1./3.)) + log(m + delta*(deenprim*sqrte - 2.)) - log(deenprim) - log(me + delta*sqrte*bee*pow(zee, -1./3.));
		if(formn<0.) formn = 0;
	
	double forme = log(beeprim*m*pow(zee, -2./3.)) - log(1. + delta*m/(me*me*sqrte)) - log(me + delta*sqrte*beeprim*pow(zee, -2./3.));
		if(forme<0.) forme = 0;
//		if(phot >=  ) forme = 0;


	double integrand = factor*(1. - (4./3.)*xplus*xminus)*(zee*formn + forme);

	return integrand;
}

double geantweigh(double x, void *p){

	double kay  = *(double *)p;
	

	double xplus = x;
	double xminus = 1.-xplus;
	
	double factor = 4.*alpha*classrad*zee;

	double delta = m*m/(2.*kay*xplus*xminus);
	double deen = 1.54*pow(ay, 0.27);
	double deenprim = pow(deen, 1.-(1./zee));
	double bee = 183.;
	double beeprim = 1429.;
	double sqrte = 1.648;
	
	double formn = log(bee*pow(zee, -1./3.)) + log(m + delta*(deenprim*sqrte - 2.)) - log(deenprim) - log(me + delta*sqrte*bee*pow(zee, -1./3.));
		if(formn<0.) formn = 0;
	
	double forme = log(beeprim*m*pow(zee, -2./3.)) - log(1. + delta*m/(me*me*sqrte)) - log(me + delta*sqrte*beeprim*pow(zee, -2./3.));
		if(forme<0.) forme = 0;
//		if(phot >=  ) forme = 0;


	double integrand = xplus*factor*(1. - (4./3.)*xplus*xminus)*(zee*formn + forme);

	return integrand;
}


int main(int argc, char **argv){

	int n = 104;
	double cross[n];
	double weigh[n];
	
	double split[n];
	
	double gamme[n];

	double kay;

	for(int i=0; i<n; i++){
	
		gamme[i] = pow(10, 1. + 0.05*i);
	
		kay = gamme[i];

		double emin = 0.5;
		double emax = 0.5 + pow(0.25 - m/kay, 0.5);

	
		double par1 = kay;
						
		double result, error;	
		
		double xl[1] = {emin};
		double xu[1] = {emax};
		
	//--------Adaptive-------//

		gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
		
		gsl_function F;
		F.function = &geanttot;
		F.params = &par1;
				
		gsl_integration_qags (&F, xl[0], xu[0], 1e-10, 1e-10, 5000, w, &result, &error);
		
		gsl_integration_workspace_free (w);
		cross[i] = result;

		gsl_integration_workspace *z = gsl_integration_workspace_alloc(5000);
		
		F.function = &geantweigh;
		F.params = &par1;
				
		gsl_integration_qags (&F, xl[0], xu[0], 1e-10, 1e-10, 5000, z, &result, &error);
		
		gsl_integration_workspace_free (z);
		weigh[i] = result;

		split[i] = weigh[i]/cross[i];
	
	}
	
	TApplication *app = new TApplication("app", &argc, argv); //needed to run ROOT externally

	TGraph *gr1 = new TGraph(n, gamme, split);
	gr1->SetLineWidth(2);
	gr1->SetLineColor(1);
	gr1->SetLineStyle(1);
	
	TCanvas *c1 = new TCanvas("c1","stuff2",500,500 );
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,10);

	c1->cd();

	pad1->SetTitle("");
	pad1->Draw();
	pad1->cd();
	pad1->SetLogx();
//		pad1->SetLogy();
	pad1->SetTicks();
	pad1->SetGridx();
	pad1->SetGridy();

	gr1->SetMaximum(1.);
	gr1->SetMinimum(0.);
	gr1->SetTitle("");

	gr1->Draw("Al");

	gr1->GetXaxis()->SetLimits(10., 1E6);
	gr1->GetXaxis()->SetTitle("E (GeV)");
	gr1->GetXaxis()->CenterTitle(1);
	gr1->GetXaxis()->SetTitleOffset(1.2);
	gr1->GetYaxis()->SetTitle("<x> of HE muon");
	gr1->GetYaxis()->CenterTitle(1);
	gr1->GetYaxis()->SetTitleOffset(1.2);

	TLegend *leg = new TLegend(0.45, 0.7, 0.90, 0.89, "", "brNDC");
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
//		leg->SetShadowColor("1");
	leg->AddEntry(gr1, "<x>", "l");
//	leg->Draw();

	app->Run();




	return 0;
}
