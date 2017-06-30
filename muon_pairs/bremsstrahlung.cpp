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

const double zee = 26;
const double ay = 55.845;
	
const double me = 0.511E-3; //electron mass in GeV
const double mu = 0.105;
const double alpha = 1./(137.036);
const double classrad = 1.859E-30; //classical radius of muon squared in cm2
//const double classrad = 1.859E-30*pow(mu/me, 2.);

const double zeered = pow(zee*alpha, 2.);
const double eff = 1.202*zeered - 1.0369*zeered*zeered + 1.008*pow(zeered, 3.)/(1.+zeered);

const double m = mu;

const double avogadro = 6.022E23; 

// struct my_params1{double a;};
 
 struct my_params2{double a; double b;};



double cillis(double *x, size_t dim, void *p){

	double vee = x[0];
	double lep = *(double *) p;
	
	double factor = alpha*classrad*4*zee/vee;
	
	double delta = m*m*vee/(2.*lep*(1.-vee));
	
	double form1 = log(189*m*pow(zee, -1./3.)/me) - log(1. + 189*pow(2.718, 0.5)*delta*pow(zee, -1./3.)/me);
		   form1 = form1 + log((2./3.)*pow(zee, -1./3.));
	double forme = log(m/delta) - log( (m*delta/(me*me)) + pow(2.718, 0.5) ) - log( 1. + me/(delta*1429*pow(zee, -2./3.)*pow(2.718, 0.5)) );
	
	double integrand = factor*(zee*form1 + forme)*((4./3.) - (4./3.)*vee + vee*vee);
	
	return integrand*vee*lep*avogadro/ay;

}

double lohmann(double *x, size_t dim, void *p){

	double vee = x[0];
	double lep = *(double *) p;
	
	double factor = alpha*classrad*4*zee*zee/vee;
	
	double delta = m*m*vee/(2.*lep*(1.-vee));
	
	double form1 = log(189*m*pow(zee, -1./3.)/me) - log(1. + 189*pow(2.718, 0.5)*delta*pow(zee, -1./3.)/me);
	double integrand = factor*form1*((4./3.) - (4./3.)*vee + vee*vee);
	
	return integrand*vee*lep*avogadro/ay;

}

double bugaev(double *x, size_t dim, void *p){

	double vee = x[0];
	double lep = *(double *)p;
	
	double delta = m*m*vee/(2.*lep*(1.-vee));
	
	double factor = 4.*classrad*zee*zee*alpha/vee;
	
	double qc = 1.9*mu*pow(zee, -1./3.);
	double beta = pow(1.+ 4.*mu*mu/(qc*qc) , 0.5);
	double bee = (184.15/pow(2.718,0.5))*pow(zee,-1./3.)*(1./me);
	double cee = (1194./pow(2.718,0.5))*pow(zee, -2./3.)*(1./me);
		
	double del1 = log(mu/qc) + (beta/2.)*log((beta+1.)/(beta-1.));
	double del2 = log(mu/qc) + 0.25*(3*beta - pow(beta,3.))*log((beta+1.)/(beta-1.)) + 2*pow(mu/qc, 2.);
	
	double phi01_1 = 0.5*log( pow(mu*bee,2.) / (1+ pow(delta*bee, 2.)) ) + 0.5;
	double phi01_2 = -delta*bee*atan(1./(delta*bee));
	double phi01_3 = (0.5/zee)*log( pow(mu*cee,2.) / (1+ pow(delta*cee, 2.)) ) +(0.5/zee);
	double phi01_4 = -delta*(cee/zee)*atan(1./(delta*cee));
	
	double phi02_1 = 0.5*log( pow(mu*bee,2.) / (1+ pow(delta*bee, 2.)) ) + (1./3.);
	double phi02_2 = 2*pow(delta*bee, 2.)*(1. - delta*bee*atan(1./(delta*bee)) -  0.75*log(1+ pow(delta*bee, -2.)) );
	double phi02_3 = (0.5/zee)*log( pow(mu*cee,2) / (1+ pow(delta*cee, 2.)) ) +(1./(3.*zee));
	double phi02_4 = 2.*pow(delta*cee, 2.)*(1./zee)*(1. - delta*cee*atan(1./(delta*cee)) -  0.75*log(1+ pow(delta*cee, -2.)) );
	
	
	double phi1 = phi01_1 + phi01_2 + phi01_3 + phi01_4 - del1;
	double phi2 = phi02_1 + phi02_2 + phi02_3 + phi02_4 - del2;
	
	double integrand = phi1*(2.-2.*vee + vee*vee) -phi2*(2./3.)*(1.-vee) ;
	integrand = integrand*factor;

	return integrand*vee*lep*avogadro/ay;

}

double geant(double *x, size_t dim, void *p){
//dsigma/dv

	double vee = x[0];
	double lep = *(double *)p;
	
	double phot = vee*lep;

	double factor = (16./3.)*alpha*classrad*zee/vee;

	double delta = m*m*vee/(2.*(lep-phot));
	double deen = 1.54*pow(ay, 0.27);
	double deenprim = pow(deen, 1.-(1./zee));
	double bee = 183.;
	double beeprim = 1429.;
	double sqrte = 1.648;
	
	double formn = log(bee*pow(zee, -1./3.)) + log(m + delta*(deenprim*sqrte - 2.)) - log(deenprim) - log(me + delta*sqrte*bee*pow(zee, -1./3.));
		if(formn<0.) formn = 0;
	
	double forme = log(beeprim*m*pow(zee, -2./3.)) - log(1. + delta*m/(me*me*sqrte)) - log(me + delta*sqrte*beeprim*pow(zee, -2./3.));
		if(forme<0.) forme = 0;
		if(phot >= lep/(1. + mu*mu/(2.*me*lep)) ) forme = 0;


	double integrand = factor*(1.-vee+0.75*vee*vee)*(zee*formn + forme);
	
	return integrand*vee*lep*avogadro/ay;
}

int main(int argc, char **argv){

	//get total energy loss in GeV cm2 g-1 to compare with PDG

	int n = 21;
	double cross[n], ince[n];

//	double zee[3] = {7., 8., 18.}; //specific zees
//	double airfrac[3] = {0.769, 0.218, 0.013};

	for(int i=0; i<n; i++){
	
		ince[i] = pow(10, 1. + 0.2*i);
	
		double lep = ince[i];
	
		double vmin = 1.E-4;
		double vmax = 1.-0.75*pow(2.718, 0.5)*pow(zee, 1./3.)*m/lep;
//		double vmax = lep-m; //geant won't integrate with this limit...use vmax above
	
		double par1 = lep;
						
		double result, error;	

	//------Monte Carlo-------//	
	
		gsl_monte_function F;
		
		F.f = &geant;
		F.dim = 1;
		F.params = &par1; 		
		
		double xl[1] = {vmin};
		double xu[1] = {vmax};
	
		
		const gsl_rng_type *T;
		gsl_rng *r;
	 
		size_t calls = 50000;
	   
		gsl_rng_env_setup ();
	 
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
	
		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);  //state alloc and free must be inside loop
	
		gsl_monte_vegas_integrate (&F, xl, xu, 1, calls/10, r, s, &result, &error);

		do{
			gsl_monte_vegas_integrate (&F, xl, xu, 1, calls, r, s, &result, &error);		
		}
		while (fabs (s->chisq -1.0) > 0.6);

		cross[i] = result*1000; 		
				
		gsl_monte_vegas_free (s);
	
	}

	cout <<cross[n-1] <<endl;
//	cout <<asymptotic()*1.E30 <<endl;
	
	
	TApplication *app = new TApplication("app", &argc, argv); //needed to run ROOT externally

	TGraph *gr1 = new TGraph(n, ince, cross);
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

	gr1->SetMaximum(1.E3);
	gr1->SetMinimum(1.E-1);
	gr1->SetTitle("bug");

	gr1->Draw("Al");

	gr1->GetXaxis()->SetLimits(10., 1E5);
	gr1->GetXaxis()->SetTitle("E (GeV)");
	gr1->GetXaxis()->CenterTitle(1);
	gr1->GetXaxis()->SetTitleOffset(1.2);
	gr1->GetYaxis()->SetTitle("Total Cross-section (ub)");
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