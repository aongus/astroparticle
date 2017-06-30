/*
Examine electron photonuclear interaction, esp to determine the average virtual photon energy
as a function of electron energy.

Energies in GeV. Photon energy >5 GeV 
Cross-section in ub. 

*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
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
#include "TStopwatch.h"

using namespace std;

const int n = 30;
const double emin = 1.E1; //1 geV min
const double emax = 300.E3; //300 TeV max
const double de = log10(emax/emin)/(n-1.);

const double melec = 0.511E-3;//GeV
//const double melec = 0.105;
const double mpi = 0.139;
const double mprot = 0.938;
const double alpha = 1./137.;
const double pi = TMath::Pi();

const double em1sq = 0.54; //GeV2
const double em2sq = 1.8;

const double ay = 14.;//

const double lambda = 0.4;

//handle weighted integral and integral using parameter-->need 2 params


struct my_params2{double a; double b;};


double lohmanncross(double x, void *p){

	double efrac = x;
	
	struct my_params2 *parameters = (my_params2 *)p; 
	double eleceng = parameters->a;
	double fracpower = parameters->b;
	
	double photeng = efrac*eleceng;
	if (photeng < 10.) return 0;
	
	double tee = melec*melec*efrac*efrac/(1.-efrac);
	double kay = 1. - 2./efrac + 2./(efrac*efrac);
	
	double ess = photeng*2.*mprot;
	
//	double photcross = 114.3 + 1.647*log(log(0.0213*photeng));
	double photcross = 67.7*pow(ess, 0.0808) + 129.*pow(ess, -0.4525);

	double ex = 0.00282*photcross*pow(ay, 1./3.);
	double gee = 3.*pow(ex, -3.)*(ex*ex/2. - 1. + exp(-1.*ex)*(1.+ex) );
	
	double cross1 = 0.75*gee*( kay*log(1. + em1sq/tee) - kay*em1sq/(em1sq + tee) - 2.*melec*melec/tee  );
	double cross2 = 0.25*( kay*log(1. + em2sq/tee)- 2.*melec*melec/tee );
	double cross3 = (0.5*melec*melec/tee)*(  0.75*gee*em1sq/(em1sq + tee) + (0.25*em2sq/tee)*log(1. + tee/em2sq)  );
	
	double cross = (cross1 + cross2 + cross3)*alpha*ay*photcross*efrac/(2.*pi);
	
	return cross*pow(efrac, fracpower);
}

double geantcross(double x, void *p){

	double efrac = x;
	
	struct my_params2 *parameters = (my_params2 *)p; 
	double eleceng = parameters->a;
	double fracpower = parameters->b;
	
	double photeng = efrac*eleceng;
	double ess = photeng*2.*mprot;
	double photcross = 67.7*pow(ess, 0.0808) + 129.*pow(ess, -0.4525);
	
	double term1 = eleceng*eleceng*(1.-efrac)*pow(melec, -2.);
	double term2 = 1. + pow(melec*efrac/lambda, 2.)/(1.-efrac);
	double term3 = 1. + eleceng*efrac/lambda;
	double term4 = 1. + 0.5*lambda/mprot + eleceng*efrac/lambda;

	double cross1 = efrac - 1. + (1. - efrac + 0.5*efrac*efrac*(1. + 2*pow(melec/lambda , 2)) )*log(term1*term2/(term3*term4));

	double cross = cross1*alpha*ay*photcross/(pi*efrac);

	return cross*pow(efrac, fracpower);

}


double bugaevnew(double x, void *p){

	double efrac = x;
	
	struct my_params2 *parameters = (my_params2 *)p; 
	double eleceng = parameters->a;
	double fracpower = parameters->b;
	
	double photeng = efrac*eleceng;
//	if (photeng < 1.) return 0;
	
//	double tee = melec*melec*efrac*efrac/(1.-efrac);
	double tee = melec*melec*efrac*efrac/(1.-efrac);
	double haitch = 1. - 2./efrac + 2./(efrac*efrac);
	
	double ess = photeng*2.*mprot;
	
	double photcross = 67.7*pow(ess, 0.0808) + 129.*pow(ess, -0.4525);

	double zee = 0.00282*photcross*pow(ay, 1./3.);
	double gee = (0.5 + ( (1.+ zee)*exp(-1.*zee) -1. )/(zee*zee) )*9./zee;

	double cross1 = haitch*log(1. + em2sq/tee) - 2.*melec*melec/tee;
	double cross2 = gee*( haitch*log(1. + em1sq/tee)  -  (haitch*2.*em1sq/(em1sq + tee)) - 2.*melec*melec/tee );
	double cross3 = (0.5*melec*melec/tee)*(  (2.*gee*em1sq/(em1sq + tee)) + (em2sq/tee)*log(1. + tee/em2sq)  );
	
	double cross = (cross1 + cross2 + cross3)*alpha*ay*photcross*efrac/(8.*pi);
	
	return cross*pow(efrac, fracpower);

}

int main(int argc, char **argv){

	int i;

	double eng[n], weighcross[n], totcross[n], avfrac[n], photoncross[n];


//	double xmin = 0.05;
//	double xmax = 1.-xmin;
	double xmin, xmax;
	
	double result, error;

	for(i=0; i<n; i++){ //define energy and do weighted integral

		eng[i] = pow(10, log10(emin)+ de*i);

//		photoncross[i] = ay*(114.3 + 1.647*log(log(0.0213*eng[i])));
		photoncross[i] = ay*(67.7*pow(eng[i]*2.*mprot, 0.0808) + 129.*pow(eng[i]*2.*mprot, -0.4525));


		struct my_params2 weighpar = {eng[i], 1.};
		
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(50000);
		
		gsl_function F;
		F.function = &geantcross;
		F.params = &weighpar;
		
//		xmin = (mpi + 0.5*mpi*mpi/mprot)/eng[i];
		xmin = 0.4/eng[i];  //somewhat sensitive to this...
//		xmax = 1.- (1. + pow(melec/mprot , 2.))*mprot*0.5/eng[i];
		xmax = 1.-0.5*mprot/eng[i];
				
		gsl_integration_qags (&F, xmin, xmax, 1e-1, 1e-1, 50000, w, &result, &error);
		
		gsl_integration_workspace_free (w);
	
		weighcross[i] = result;	//

	}

	for(i=0; i<n; i++){ //do total integral

		struct my_params2 totpar = {eng[i], 0.};
		
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(50000);

		gsl_function F;		
		F.function = &geantcross;
		F.params = &totpar;
		
//		xmin = (mpi + 0.5*mpi*mpi/mprot)/eng[i];
		xmin = 0.4/eng[i];
//		xmax = 1.- (1. + pow(melec/mprot , 2.))*mprot*0.5/eng[i];
		xmax = 1.-0.5*mprot/eng[i];

			
		gsl_integration_qags (&F, xmin, xmax, 1e-10, 1e-10, 50000, w, &result, &error);
		
		gsl_integration_workspace_free (w);
	
		totcross[i] = result;	

		avfrac[i] = 0;
		if(totcross[i] > 0.) avfrac[i]= weighcross[i]/totcross[i];
	}

//	cout <<avfrac[20] <<endl;
//	cout <<totcross[20] <<endl;

	TApplication *app = new TApplication("app", &argc, argv); //needed to run ROOT externally

//-------Cross-section plot---------// 
	TGraph *enuccross = new TGraph(n, eng, totcross);
	enuccross->SetLineWidth(2);
	enuccross->SetLineColor(2);
	enuccross->SetLineStyle(1);

	TGraph *rephotcross = new TGraph(n, eng, photoncross);
	rephotcross->SetLineWidth(2);
	rephotcross->SetLineColor(1);
	rephotcross->SetLineStyle(1);

	
	TCanvas *c1 = new TCanvas("events1","stuff1",500,500 );
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,10);

	c1->cd();

	pad1->SetTitle("");
	pad1->Draw();
	pad1->cd();
	pad1->SetLogx();
	pad1->SetLogy();
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->SetTicks();

	enuccross->SetMaximum(1E4);
	enuccross->SetMinimum(1E2);
	enuccross->SetTitle("Total Cross-Section");

	enuccross->Draw("Al");
	enuccross->GetXaxis()->SetLimits(emin, emax);
	enuccross->GetXaxis()->SetTitle("Energy (GeV)");
 	enuccross->GetXaxis()->CenterTitle(1);
 	enuccross->GetXaxis()->SetTitleOffset(1.2);
	enuccross->GetYaxis()->SetTitle("cross (ub)");
	enuccross->GetYaxis()->CenterTitle(1);
 	enuccross->GetYaxis()->SetTitleOffset(1.2);

	rephotcross->Draw("same");
	
	TLegend *leg1 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
	leg1->SetFillStyle(0);
	leg1->SetBorderSize(0);
//	leg1->SetShadowColor("1");
 	leg1->AddEntry(enuccross, "e->g*N", "l");
 	leg1->AddEntry(rephotcross, "gN", "l");
 	leg1->Draw();


//------average e transfer plot---------//

	TGraph *averagefrac = new TGraph(n, eng, avfrac);
	averagefrac->SetLineWidth(2);
	averagefrac->SetLineColor(2);
	averagefrac->SetLineStyle(1);
	
	TCanvas *c2 = new TCanvas("events2","stuff2",500,500 );
	TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,1.0,10);

	c2->cd();

	pad2->SetTitle("");
	pad2->Draw();
	pad2->cd();
	pad2->SetLogx();
//	pad2->SetLogy();
	pad2->SetGridx();
	pad2->SetGridy();
	pad2->SetTicks();

	averagefrac->SetMaximum(0.14);
	averagefrac->SetMinimum(0.);
	averagefrac->SetTitle("average photon energy frac of E_e");

	averagefrac->Draw("Al");
	averagefrac->GetXaxis()->SetLimits(emin, emax);
	averagefrac->GetXaxis()->SetTitle("Electron Energy (GeV)");
 	averagefrac->GetXaxis()->CenterTitle(1);
 	averagefrac->GetXaxis()->SetTitleOffset(1.2);
	averagefrac->GetYaxis()->SetTitle("<v>");
	averagefrac->GetYaxis()->CenterTitle(1);
 	averagefrac->GetYaxis()->SetTitleOffset(1.2);



	app->Run();

	return 0;
}