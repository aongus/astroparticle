// muon pair production cross-section using GSL's integrators

/* integrate Tsai 1974 Eq. 2.7 using GSL Vegas
	change variables to get a rectangular region


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
#include "TStopwatch.h"
 
 using namespace std;
 
 
const double pi = TMath::Pi();
const double zee = 7.2;
const double ay = 15.;
const double mnuc = ay*0.938;
	
const double dee = (0.162)*pow(ay, -2./3.); //GeV2
const double mu = 0.105;
const double me = 0.511E-3; //electron mass in GeV
const double alpha = 1./(137.036);

 struct my_params1{double a;};
 
 struct my_params2{double a; double b;};

 struct my_params3{double a; double b; double c;};

double cross_section(double* x, size_t dim, void * p){ //Function of cosplus. Takes three variables as parameters: k, p, costheta.
//	1st integral is over costheta+  = x	
// 	2nd integral is over 2pi*costheta (emission solid angle)
//	3rd integral is over p, muon momentum

// only do coherent scattering, so no change in final nuclear mass
	
	double cosplus = x[0]; //costheta+
	double phi = x[1]; //costheta
	double xi = x[2]; //pee

	double kay = *(double *)p;

	double kthresh = (2.*mu*mu+2.*mu*mnuc)/mnuc;
	if (kay<kthresh) return 0;

//-------------pee

	double pfunc1 = mnuc*kay - mnuc*mu;
	double pfunc2 = mnuc*mnuc + 2.*kay*mnuc;
		
	double pmax = (kay*pfunc1 + (kay+mnuc)*pow(pfunc1*pfunc1 - mu*mu*pfunc2, 0.5))/pfunc2;
	double pmin = (kay*pfunc1 - (kay+mnuc)*pow(pfunc1*pfunc1 - mu*mu*pfunc2, 0.5))/pfunc2;
	
	if(pmin<0.) pmin = 0;
	
	double pee = xi*(pmax - pmin) + pmin;

	double mueng = pow(pee*pee+mu*mu,0.5);

//------------costheta	

	double costhetamax = (1./(kay*pee))*(mu*mnuc - mnuc*(kay - mueng) + kay*mueng) ;

	if(costhetamax<= -1.) costhetamax = -1.;
	if(costhetamax>1.) costhetamax = 1.;

	double costheta = phi*(1.-costhetamax) +costhetamax;
	
	
	double kmin = mnuc*(mu+mueng)/(mnuc - mueng + pee*costheta);
	if (kmin<0.) return 0;
	if (kay<kmin) return 0;

//-----------function

	double dotprod = kay*mueng - kay*pee*costheta;

	double yoo =  mu*mu +mnuc*mnuc + 2*mnuc*(kay-mueng) - 2*dotprod;
			if(yoo<0.) cout <<"yoo is negative" <<endl;
		   yoo =  pow(yoo, 0.5);
	
	
	double kay_s = (kay*mnuc - dotprod)/yoo;
	double ee_plus_s = (yoo*yoo + mu*mu - mnuc*mnuc)/(2.*yoo);
	double pee_plus_s = pow(ee_plus_s*ee_plus_s - mu*mu , 0.5);
			if(ee_plus_s*ee_plus_s - mu*mu<0.) cout <<"p_+s is im" <<endl;
	double pee_init_s = (mnuc/yoo)*pow(kay*kay + pee*pee - 2*pee*kay*costheta , 0.5);
			if(kay*kay + pee*pee - 2*pee*kay*costheta<0.) cout <<"pee_init_s is im" <<endl;
	double ee_s = (dotprod - mu*mu + mueng*mnuc)/yoo;
	double costhetak = (kay_s - ee_s)/(pee_init_s) + (dotprod)/(kay_s*pee_init_s);
	double sinsqthetak = 1. - costhetak*costhetak;
	
	double qsquared = 2.*mu*mu - 2.*dotprod - 2.*ee_plus_s*(kay_s - ee_s) + 2.*pee_plus_s*pee_init_s*cosplus;
	
	double doubleyoo = ee_plus_s - pee_plus_s*cosplus*costhetak;	
	double wai = mu*mu*sinsqthetak + pow(pee_plus_s*cosplus - ee_plus_s*costhetak , 2.);
		if(wai<0.) cout << "wai is negative" <<endl;
		wai = pow(wai, 0.5);
	       
	
	//delta = 0 for coherent scattering
	
	double haitch = -1.*mu*mu*( 0.5*qsquared*(1. - 2.*mueng/mnuc) + 2*mueng*mueng );

	double bee1 = (-2./dotprod)*(  (mu*mu - 0.5*qsquared)*(2.*mueng*(mueng-kay) + 0.5*qsquared*( (kay-2.*mueng)/mnuc + 1.) ) - 0.5*qsquared*kay*kay );
	double bee2 = (qsquared/mnuc)*(mnuc + mueng - kay - 0.5*qsquared/mnuc) + dotprod;
	double bee = bee1 + bee2;
	
	double cee1 = -1.*(mu*mu/(dotprod*dotprod)) *(  2.*(kay- mueng + 0.5*qsquared/mnuc)*(kay-mueng) + 0.5*qsquared );
	double cee2 = qsquared*(1.-mueng/mnuc)/dotprod;
	double cee = cee1 + cee2;

	double bigdee = 1./dotprod;

//-------------formfactors
	
	//-------nuclear
//	only do coherent so no W1
//	double form1 = ;

	double form2;

// 	if(kay< 50.){	
// 	form2 = 2.*mnuc*zee*zee/pow(1. - qsquared/dee, 2.);
// 	}	
// 	else{
	//--------atomic
 	double ay = 184.15*pow(2.718,-0.5)*pow(zee,-1./3.)*(1./me);
 //	double ayprim = 1194.*pow(2.718,-0.5)*pow(zee, -2./3.)*(1./me);
 	
 	double form2el = zee*zee*pow(ay, 4.)*pow(qsquared, 2.)*pow( 1.- ay*ay*qsquared, -2.);
 //	double form2inel = ;
 	form2 = 2*mnuc*form2el; 
//	}
	
	
	double integrand1 = (-1./(2*pi))*pow(alpha, 3.)*pow(0.19732, 2.)*1.E-26*pee*pee*pee_plus_s/(yoo*kay*mueng*qsquared*qsquared);
	double integrand2 = form2*( (haitch*doubleyoo/(wai*wai*wai*kay_s*kay_s)) + (bee/(wai*kay_s)) + cee + bigdee*kay_s*doubleyoo);

	return integrand1*integrand2*(1.-costhetamax)*2.*pi*(pmax - pmin);  
}



 
double integrator(double x){  //do the integral over p

	double kay = x;
	

	double p[1] = {kay};

	double lo[3] = {-1., 0., 0.};
	double hi[3] = {1., 1., 1.};
	
	double result, error;	

	gsl_monte_function F;
	
	F.f = &cross_section;
	F.dim = 3;
	F.params = &p; 		
	
	const gsl_rng_type *T;
	gsl_rng *r;
 
	size_t calls = 100000;
   
	gsl_rng_env_setup ();
 
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);  //state alloc and free must be inside loop

	gsl_monte_vegas_integrate (&F, lo, hi, 3, calls/10, r, s, &result, &error);

	do{
		gsl_monte_vegas_integrate (&F, lo, hi, 3, calls/2, r, s, &result, &error);		
	}
	while (fabs (s->chisq -1.0) > 0.5);
	
	
	gsl_monte_vegas_free (s);
	return result;
}


int main(int argc, char **argv){

// 	double kay = 100.;
// 
// 	double par1 = kay;
// 
// 	TStopwatch time;
// 	time. Start();
// 	
// 	double answer = integrator(par1);
// 
// 	time.Stop();
// 	
// 	cout <<answer <<" in "<<time.RealTime() <<" s"<<endl;
	

//--plot total cross-section versus photon energy---//

	int n = 5;
	double cross[n], gamme[n];
	double kay;

	for(int i=0; i<n; i++){
		gamme[i] = pow(10., 3.5+0.2*i);
		kay = gamme[i];
	
//		cross[i] = 0;
		cross[i] = 1.E30*integrator(kay);
		cout <<cross[i] <<"    "<<kay<<endl;
	}

	int q = 26; 

	double vankoveng2[26] = {0.010, 0.016, 0.025, 0.039, 0.063, 0.1, 0.157, 0.252, 0.396, 0.621, 1., 1.589, 2.525, 4.012, 6.293, 10., 15.890, 24.926, 40.119, 62.934, 100., 158.898, 252.484, 396.064, 629.336, 1000. };
	double vankovcross2[26] = {0.691, 0.829, 0.968, 1.175, 1.382, 1.590, 1.797, 2.074, 2.696, 3.664, 4.839, 6.152, 7.396, 8.502, 9.470, 10.3, 10.783, 11.129, 11.336, 11.475, 11.613, 11.682, 11.751, 11.820, 11.820, 11.889};


	for(int i=0; i<q; i++){
		vankoveng2[i] = vankoveng2[i]*1000;
	}

	int p = 121;
	double berezeng[p], berezcross[p];
	double x, y;
	int lines = 0;
	
  	ifstream infile;
  	infile.open("berez_pairs.data");	  	
 	while(!infile.eof()){
  		infile >> x >> y; 	
  		berezeng[lines] = x*1000.;
  		berezcross[lines] = y;	
  		lines=lines+1;
  	}
	
	TApplication *app = new TApplication("app", &argc, argv); //needed to run ROOT externally

	TGraph *gr1 = new TGraph(n, gamme, cross);
	gr1->SetLineWidth(2);
	gr1->SetLineColor(2);
	gr1->SetLineStyle(1);
	
	TGraph *gr2 = new TGraph(q, vankoveng2, vankovcross2);
	gr2->SetLineWidth(2);
	gr2->SetLineColor(1);
	gr2->SetLineStyle(1);

	TGraph *gr3 = new TGraph(p, berezeng, berezcross);
	gr3->SetLineWidth(2);
	gr3->SetLineColor(3);
	gr3->SetLineStyle(1);

	TCanvas *c1 = new TCanvas("c1","stuff2",500,500 );
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,10);

	c1->cd();

	pad1->SetTitle("");
	pad1->Draw();
	pad1->cd();
	pad1->SetLogx();
//		pad1->SetLogy();

	gr1->SetMaximum(30.);
	gr1->SetMinimum(0.);
	gr1->SetTitle("");

	gr1->Draw("Al");

	gr1->GetXaxis()->SetLimits(10., 1E6);
	gr1->GetXaxis()->SetTitle("E (GeV)");
	gr1->GetXaxis()->CenterTitle(1);
	gr1->GetXaxis()->SetTitleOffset(1.2);
	gr1->GetYaxis()->SetTitle("Total Cross-section (ub)");
	gr1->GetYaxis()->CenterTitle(1);
	gr1->GetYaxis()->SetTitleOffset(1.2);

	gr2->Draw("same");
	gr3->Draw("same");

	TLegend *leg = new TLegend(0.45, 0.7, 0.90, 0.89, "", "brNDC");
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
//		leg->SetShadowColor("1");
	leg->AddEntry(gr1, "Integrand", "l");
//		leg->Draw();

	app->Run();
 
 
 	return 0;

}