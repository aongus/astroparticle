/* muon pair production cross-section using GSL's integrators

 integrate Tsai 1974 Eq. 3.5
 doesn't seem to work well, some flaw in the approximation causes the 
 cross-section to collapse above 10 TeV; also much too big for low energies?
 
 know it doesn't work for large lepton masses. know approximation is for small production
 angles so don't know how much integrating to all angles hurts. 
 
 However, do see "knee" in Vankov's cross-section...
 also peak is at Berezinsky's asymptotic value...

*/

#include <stdlib.h>
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
 
 
	const double pi = TMath::Pi();
	const double zee = 7.2;
	const double ay = 14.5;
	const double mnuc = ay*0.938;
		
	const double dee = (0.164)*pow(ay, -2./3.); //GeV2
//	const double mu = 0.1056;
	const double mu = 0.105;
	const double me = 0.511E-3; //electron mass in GeV
	const double alpha = 1./(137.036);
	const double zeered = pow(zee*alpha, 2.);


	const double eff = 1.202*zeered - 1.0369*zeered*zeered + 1.008*pow(zeered, 3.)/(1.+zeered);


 struct my_params1{double a;};
 
 struct my_params2{double a; double b;};



double crossdiff(double x, void * p){

	double theta = x;
//	double ell = x;

	struct my_params2 *G_params = (struct my_params2 *)p;	
	
	double kay = G_params->a;
	double pee = G_params->b;
	
	double mueng = pow(pee*pee+mu*mu,0.5);
//	double theta = mu*pow(ell,0.5)/mueng;

	double kthresh = (2.*mu*mu+2.*mu*mnuc)/mnuc;
	if (kay<kthresh) return 0;

	double kmin = mnuc*(mu+mueng)/(mnuc - mueng + pee*cos(theta));
	if (kmin<0.) return 0;
	if (kay<kmin) return 0;

	double ex, cee, bee, ecks, integrand1, integrand2, factor, integrand;

	ex = mueng/kay;

	double ell = pow(ex*kay*theta/mu, 2.);
	
	cee = pow(mu*(1.+ell), 2.)/dee;
	bee = pow( mu*mu*(1.+ell)/(2.*kay*ex*(1.-ex)) , 2.)/dee; //kay is photon energy
	ecks = pow(zee, 2.)*( (1.+2.*bee)*log( (1.+(1./bee)) / (1.+(1./cee)) ) - (1.-(bee/cee))*(1.+2.*cee)/(1.+cee)  );
			
	integrand1 = (2.*ex*ex - 2.*ex +1.)/pow(1.+ell, 2.);
	integrand2 = 4.*ell*ex*(1.-ex)/pow(1.+ell, 4.);
		
	factor = 0.38937E-27*(2.*pow(alpha, 3.)*pow(ex*kay, 2.)*(1./pi)*(1./kay)*pow(mu, -4.)); //1st factor converts GeV2 to cm2
	integrand = (integrand1 + integrand2)*(ecks - 2.*zee*zee*eff)*factor*2.*pi*sin(theta);

//	factor = 1.859E-30*2.*(alpha/kay)*sin(theta)/theta; //not trustful of this since the sintheta/theta factor is missing
//	integrand = (integrand1 + integrand2)*(ecks - 2.*zee*zee*eff)*factor;

	return integrand;  
}

double crossdiffint(double x, void * p){

	double pee = x;

	struct my_params1 *F_params = (struct my_params1 *)p;	
	double kay = F_params->a;	

	double mueng = pow(pee*pee+mu*mu,0.5);

	double thetamax;
	double thetamin = 0.;	
	double costhetamax = (1./(kay*pee))*(mu*mnuc - mnuc*(kay - mueng) + kay*mueng) ;

	if(costhetamax>= -1.) {
		if(costhetamax<= 1.){
			thetamax = acos(costhetamax);
		}
		else{	
			thetamax = 0.;
		}
	}	
	else{
		thetamax = pi;
	}

//	double lmin = thetamin;
//	double lmax = pow(mueng*thetamax/mu, 2.);

	struct my_params2 par2 = {kay, pee};
	double result2, error2;	

	//---------Gaussian Quadrature----------//

	int n_divisions2 = 50000;
	gsl_integration_workspace *v = gsl_integration_workspace_alloc(n_divisions2);

	gsl_function G;
	G.function = &crossdiff;

	G.params = &par2;
			
	gsl_integration_qags (&G, thetamin, thetamax, 1e-10, 1e-10, n_divisions2, v, &result2, &error2);
//	gsl_integration_qags (&G, lmin, lmax, 1e-5, 1e-7, n_divisions2, v, &result2, &error2);

	gsl_integration_workspace_free (v);

	//----------Monte Carlo-----------//
	
// 	gsl_monte_function G;
// 
// 	G.f = &crossdiff;
// 	G.dim = 1;
// 	G.params = &par2; 
//      
//     double zl[1] = { thetamin };
//     double zu[1] = { thetamax };
//      
//     const gsl_rng_type *T;
//     gsl_rng *r;
//      
//     size_t calls = 500;
// 	   
// 	gsl_rng_env_setup ();
//      
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);
// 
//     gsl_monte_vegas_state *k = gsl_monte_vegas_alloc (1);
//      
//     gsl_monte_vegas_integrate (&G, zl, zu, 1, calls, r, k, &result2, &error2);
//    
//     gsl_monte_vegas_free (k);
 	

	return result2;
}



int main(int argc, char **argv){

	int n = 50;
	double cross[n], gamme[n];

	for(int i=0; i<n; i++){
	
		gamme[i] = pow(10, 1. + 0.1*i);
		double kay = gamme[i];
		
//		double kay = 200.; //parameter to pass to integral (k)
	
		double pfunc1 = mnuc*kay - mnuc*mu;
		double pfunc2 = mnuc*mnuc + 2.*kay*mnuc;
	
	//	cout <<acos(1.) <<endl;
		
		double pmax = (kay*pfunc1 + (kay+mnuc)*pow(pfunc1*pfunc1 - mu*mu*pfunc2, 0.5))/pfunc2;
		double pmin = (kay*pfunc1 - (kay+mnuc)*pow(pfunc1*pfunc1 - mu*mu*pfunc2, 0.5))/pfunc2;
	
		if(pmin<0.) pmin = 0;
	
	//	cout <<pmax <<"	"<<pmin <<endl;

	//--------Gaussian Quadrature----------//
	
		struct my_params1 par1 = {kay};
		double result, error;	
	
		int n_divisions1 = 50000;
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(n_divisions1);
	
		gsl_function F;
		F.function = &crossdiffint;
	
	
		F.params = &par1;
				
		gsl_integration_qags (&F, pmin, pmax, 1e-10, 1e-10, n_divisions1, w, &result, &error);
	
//	 	cout <<result <<endl;
		cross[i] = result*1E30;
		
		gsl_integration_workspace_free (w);
	
	
	//------Monte Carlo-------//	
	
// 		gsl_monte_function F;
// 		
// 		F.f = &crossdiffint;
// 		F.dim = 1;
// 		F.params = &par1; 		
// 		
// 		double xl[1] = {pmin};
// 		double xu[1] = {pmax};
// 	
// 		
// 		const gsl_rng_type *T;
// 		gsl_rng *r;
// 	 
// 		size_t calls = 500;
// 	   
// 		gsl_rng_env_setup ();
// 	 
// 		T = gsl_rng_default;
// 		r = gsl_rng_alloc (T);
// 	
// 		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);  //state alloc and free must be inside loop
// 	
// 		gsl_monte_vegas_integrate (&F, xl, xu, 1, calls, r, s, &result, &error);
// 		cross[i] = result*1E30;		
// 				
// 		gsl_monte_vegas_free (s);
	
	}

	int q = 26; 

	double vankoveng2[26] = {0.010, 0.016, 0.025, 0.039, 0.063, 0.1, 0.157, 0.252, 0.396, 0.621, 1., 1.589, 2.525, 4.012, 6.293, 10., 15.890, 24.926, 40.119, 62.934, 100., 158.898, 252.484, 396.064, 629.336, 1000. };
	double vankovcross2[26] = {0.691, 0.829, 0.968, 1.175, 1.382, 1.590, 1.797, 2.074, 2.696, 3.664, 4.839, 6.152, 7.396, 8.502, 9.470, 10.3, 10.783, 11.129, 11.336, 11.475, 11.613, 11.682, 11.751, 11.820, 11.820, 11.889};

//	double newcross[26];
//	double lambda = 0.010; //30 MeV?
//	double newfactor = (28./9.)*zee*zee*alpha*1.859E-30;//in cm2

	for(int i=0; i<q; i++){
		vankoveng2[i] = vankoveng2[i]*1000;
//		newcross[i] = log(2*lambda*vankoveng2[i]/(mu*mu)) - (171./42.) -(3./35.)*pow(lambda/mu, 2.)*(1. + 8.*0.116*zeered);
//		newcross[i] = newfactor*newcross[i]*1.E30;
	}

	ifstream infile;
	infile.open("berez_pairs.data");
	
	int p = 121; //number of lines
	
	double x, y;
	double berez_eng[p], berez_cross[p];
	int lines = 0;
	
	while(!infile.eof()){ //read in until the end of the file is reached
		infile >> x  >> y;
		berez_eng[lines] = x*1000.;
		berez_cross[lines] = y;
		lines = lines+1;
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

// 	TGraph *gr3 = new TGraph(q, vankoveng2, newcross);
// 	gr3->SetLineWidth(2);
// 	gr3->SetLineColor(4);
// 	gr3->SetLineStyle(1);

	TGraph *gr3 = new TGraph(p, berez_eng, berez_cross);
	gr3->SetLineWidth(2);
	gr3->SetLineColor(4);
	gr3->SetLineStyle(1);

	TCanvas *c1 = new TCanvas("c1","stuff2",500,500 );
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,10);

	c1->cd();

	pad1->SetTitle("");
	pad1->Draw();
	pad1->cd();
	pad1->SetLogx();
//		pad1->SetLogy();

	gr1->SetMaximum(50.);
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