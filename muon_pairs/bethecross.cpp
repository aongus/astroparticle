/*
2 approximations for muon pair production:
formula bhcross is the Bethe-Heitler HE cross-section assuming exponential screening
formula wlcross is the Wheeler-Lamb cross-section

bhcross is in agreement with Vankov above 400 GeV if we set Z=7.9 (really z_eff = 7.2)

wlcross does not agree since there is no screening to modulate. formula is probably overcooked 
	assuming light elements (inelastic corrections) and electrons. 


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

const double zee = 7; //7.2 is air avg
const double ay = 14; //14.7 is air avg
	
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


double bhcross( double x, void *p){

//these formulas work in units of lepton mass. Therefore convert _afterwards_ to get units right.

	double eplus = x;
	double kay = *(double *) p;
//	struct my_params2 *G_params = (struct my_params2 *)p;	
	
//	double kay = G_params->a;
//	double zee = G_params->b;

	double eminus = kay-eplus;



	double bee = 2*eplus*eminus*pow(zee, 1./3.)/(111.*kay);
//	 bee = bee*pow(me/(m*m*m*m), -2.);
	double emm = pow(kay/(2*eplus*eminus), 2.) + pow( pow(zee, 1./3.)/111., 2.) ;
	emm = 1./emm;
	
	double factor = 2.*alpha*zee*zee*classrad/(pow(kay, 3.));

	double integrand1 = (eplus*eplus + eminus*eminus + (2./3.)*eplus*eminus)*( log(emm) + 1. - (2./bee)*atan(bee) );
	double integrand2 = eplus*( (2./(bee*bee))*log(1. + bee*bee) + 4.*(2.-bee*bee)*atan(bee)/(3.*bee*bee*bee) - 8./(3.*bee*bee) + 2./9.  );

	double integrand = factor*(integrand1 - integrand2);

	return integrand*1.25;
}

double wlcross( double x, void *p){

	//asymptotes to 2 ub by 1 GeV because can't put me/mu into screening. 
	//basically can use total formula since they agree completely at GeV energies. 
	
//	double wltotal = alpha*zee*classrad*(22.6 - 2.08*log(zee));
	
	double eplus = x;
	double kay = *(double *) p;
//	struct my_params2 *G_params = (struct my_params2 *)p;	
	
//	double kay = G_params->a;
//	double zee = G_params->b;

	double eminus = kay-eplus;
	
	double factor = alpha*zee*classrad*pow(kay, -3.);

	double integrand1 = (eplus*eplus + eminus*eminus)*( 29.1 - (8./3.)*log(zee));
	double integrand2 = ((2./3.)*eplus*eminus)*(28.4 - (8./3.)*log(zee));

	double integrand = factor*(integrand1 + integrand2);

	return integrand;
}

double bhscreened(double* x, size_t dim, void *p){

	//Bethe-Heitler for complete screening

	double eplus = x[0];
	double kay  = *(double *) p;
	
	double eminus = kay - eplus;

	double factor = 4.*alpha*zee*zee*classrad*pow(kay, -3.);
	
	
	double integrand1 = (eplus*eplus + eminus*eminus + (2./3.)*eplus*eminus)*log(111.7*pow(zee, -1./3.)*(mu/me));
	double integrand2 = (1./9.)*eplus*eminus;
	
	double integrand = factor*(integrand1 + integrand2);

	return integrand;
}

double bhunscreened(double x){
	//bethe-heitler bare point nucleus, no screening
	
	double kay  = x;
	
	double factor = (28./9.)*alpha*zee*zee*classrad;

	double formula = factor*(log(2.*kay/m) - 2.595 - coulombcorr); //coulombcorr is negligible

	return formula;
}

double asymptotic(){

	double factor = (7./9.)*alpha*classrad;
	
 	double ay = 184.15*pow(2.718,-0.5)*pow(zee,-1./3.)*(1./me);
 	double ayprim = 1194.*pow(2.718,-0.5)*pow(zee, -2./3.)*(1./me);
	
	double psi0 = 2.*(1. + log( pow(ay*m, 2.)*pow(zee, 2./3.) ) );
	double phi0 = 2.*(1. + log( pow(ayprim*m, 2.)*pow(zee, 4./3.) ) );
	
	double int1 = zee*zee*(psi0 - (4./3.)*log(zee) - 4.*eff);
	double int2 = zee*(phi0 - (8./3.)*log(zee));
	double int3 = (-2./21.)*(zee*zee + zee);

	double answer =  factor*(int1 + int2 + int3);
	
	return answer;
}

double bugaev(double x, void *p){

	double kay  = *(double *)p;
	
	double vee = x;
	
	double delta = m*m/(2.*kay*vee*(1.-vee));
	
	double factor = 4.*classrad*zee*zee*alpha;
	
	double qc = 1.9*mu*pow(zee, -1./3.);
	double beta = pow(1.+ 4.*mu*mu/(qc*qc) , 0.5);

	double bee = (184.15/pow(2.718,0.5))*pow(zee,-1./3.)*(1./me);
	double cee = (1194./pow(2.718,0.5))*pow(zee, -2./3.)*(1./me);

//	double bee = (182.5/pow(2.718,0.5))*pow(zee,-1./3.)*(1./me);
//	double cee = (1429./pow(2.718,0.5))*pow(zee, -2./3.)*(1./me);

// 	double deen = 1.54*pow(ay, 0.27);
// 	double bee = (182.5/deen)*pow(zee,-1./3.)*(1./me);
// 	double cee = (1429./deen)*pow(zee, -2./3.)*(1./me);
	
		
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
	
	double integrand = phi1*(2.*vee*vee-2.*vee + 1) + phi2*(2./3.)*vee*(1.-vee) ;

	integrand = integrand*factor;

	return integrand;
}

double geantel(double* x, size_t dim, void *p){

	double kay  = *(double *)p;
	

	double xplus = x[0];
	double xminus = 1.-xplus;
	

	double factor = 4.*alpha*zee*zee*classrad;

	double bee = 183.;
	double deen = 1.54*pow(ay, 0.27);
	double delta = m*m/(2.*kay*xplus*xminus);

	double w_asymptotic = bee*m*pow(zee, -1./3.)/(deen*me);

	double w = w_asymptotic*(1. + (delta/m)*(deen*sqrte - 2.))/(1. + bee*pow(zee, -1./3.)*sqrte*delta/me);

	double integrand = factor*(1. - (4./3.)*xplus*xminus)*log(w);

	return integrand;
}

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
//	forme = 0;

	double integrand = factor*(1. - (4./3.)*xplus*xminus)*(zee*formn + forme);

	return integrand;
}

double mufrome(double x){ //e+e-->mu+mu-

	double ethresh = 43.69;//GeV

	double positeng = x;
	
	if (positeng < ethresh) return 0;
	
	double xi = ethresh/positeng;
	
	double crosssec = (3.14195*classrad/3.)*xi*(1. + 0.5*xi)*pow(1.-xi, 0.5);
	
	return crosssec;
}

int main(int argc, char **argv){

	int n = 104;
	double cross[n];
	double gamme[n];

	double unscreened[n];
	double elecmuons[n];
	double kay;

//	double zee[3] = {7., 8., 18.}; //specific zees
//	double ay[3] = {14., 15.9994, 39.948};
//	double airfrac[3] = {0.769, 0.218, 0.013};

	for(int i=0; i<n; i++){
	
		//all energies in units of muon mass
		gamme[i] = pow(10, 1. + 0.05*i);
	
//		double kay = gamme[i]/m;
		kay = gamme[i];
		
//		unscreened[i] = bhunscreened(kay)*1.E30;
		elecmuons[i] = mufrome(kay)*1.E30*zee;
//		cout << elecmuons[i] <<"	"<<gamme[i] <<endl;
		
//		double kay = 200.; //parameter to pass to integral (k)
	
//		double emin = 1.0;
//		double emax = kay-1.0;

		double emin = 0.5 - pow(0.25 - m/kay, 0.5);
		double emax = 0.5 + pow(0.25 - m/kay, 0.5);

	
		double par1 = kay;
						
		double result, error;	
		
		double xl[1] = {emin};
		double xu[1] = {emax};


	//------Monte Carlo-------//	
	
// 		gsl_monte_function F;
// 		
// 		F.f = &bugaev;
// 		F.dim = 1;
// 		F.params = &par1; 			
// 		
// 		const gsl_rng_type *T;
// 		gsl_rng *r;
// 	 
// 		size_t calls = 50000;
// 	   
// 		gsl_rng_env_setup ();
// 	 
// 		T = gsl_rng_default;
// 		r = gsl_rng_alloc (T);
// 	
// 		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);  //state alloc and free must be inside loop
// 	
// 		gsl_monte_vegas_integrate (&F, xl, xu, 1, calls/10, r, s, &result, &error);
// 
// 		do{
// 			gsl_monte_vegas_integrate (&F, xl, xu, 1, calls, r, s, &result, &error);		
// 		}
// 		while (fabs (s->chisq -1.0) > 0.5);	
// 				
// 		gsl_monte_vegas_free (s);
		
	//--------Adaptive-------//

		gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
		
		gsl_function F;
		F.function = &geanttot;
		F.params = &par1;
				
		gsl_integration_qags (&F, xl[0], xu[0], 1e-10, 1e-10, 5000, w, &result, &error);
		
		gsl_integration_workspace_free (w);
		cross[i] = result*1E30;	//result is in cm2, convert to ub

/*
		par1 = kay/m;

		gsl_integration_workspace *z = gsl_integration_workspace_alloc(5000);		
		F.function = &bhcross;
		F.params = &par1;
				
		gsl_integration_qags (&F, 1.0, kay/m - 1.0, 1e-10, 1e-10, 5000, w, &result, &error);	
	
		gsl_integration_workspace_free (z);

		unscreened[i] = result*1E30;
*/

		unscreened[i] = bhunscreened(kay)*1E30;	
	}
/*
	for(int i=0; i<n; i++){
		cout <<elecmuons[i] <<", ";
	}
	cout <<endl;
*/
	
// 	cout <<cross[n-1] <<"	"<< gamme[n-1]<<endl;
// 	cout <<asymptotic()*1.E30 <<endl;

//	double kaytest = 1.E6;//in GeV
//	cout <<geanttot(0.5, &kaytest)*1E30*1.E3/kaytest <<endl; //convert to per TeV

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
	gr1->SetLineColor(1);
	gr1->SetLineStyle(1);
	
	TGraph *gr2 = new TGraph(q, vankoveng2, vankovcross2);
	gr2->SetLineWidth(2);
	gr2->SetLineColor(1);
	gr2->SetLineStyle(2);

	TGraph *gr3 = new TGraph(p, berezeng, berezcross);
	gr3->SetLineWidth(2);
	gr3->SetLineColor(3);
	gr3->SetLineStyle(1);

 	TGraph *gr4 = new TGraph(n, gamme, unscreened);
 	gr4->SetLineWidth(2);
 	gr4->SetLineColor(1);
 	gr4->SetLineStyle(3);
	
	TGraph *etomu = new TGraph(n, gamme, elecmuons);
	etomu->SetLineWidth(2);
	etomu->SetLineColor(1);
	etomu->SetLineStyle(4);

	TCanvas *c1 = new TCanvas("c1","stuff2",500,500 );
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,10);

	c1->cd();

	pad1->SetTitle("");
	pad1->Draw();
	pad1->cd();
	pad1->SetLogx();
//		pad1->SetLogy();
	pad1->SetTicks();


	gr1->SetMaximum(25.);
	gr1->SetMinimum(0.);
//	gr1->SetTitle("");
	gr1->SetTitle("Muon pair production cross-section on N (Z=7)");

	gr1->Draw("Al");

	gr1->GetXaxis()->SetLimits(10., 1E6);
	gr1->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
	gr1->GetXaxis()->CenterTitle(1);
	gr1->GetXaxis()->SetTitleOffset(1.2);
	gr1->GetYaxis()->SetTitle("Total Cross-section (#mub)");
	gr1->GetYaxis()->CenterTitle(1);
	gr1->GetYaxis()->SetTitleOffset(1.2);

	gr2->Draw("same");
//	gr3->Draw("same");
	gr4->Draw("same");
	etomu->Draw("same");

	TLegend *leg = new TLegend(0.45, 0.7, 0.90, 0.89, "", "brNDC");
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
//		leg->SetShadowColor("1");
	leg->AddEntry(gr1, "GEANT+inelastic", "l");
	leg->AddEntry(gr4, "No screening", "l");
	leg->AddEntry(gr2, "Stanev & Vankov (1985)", "l");
	leg->AddEntry(etomu, "Electronic muon pairs" , "l");
	leg->Draw();

	TLegend *leg1 = new TLegend(0.45, 0.7, 0.60, 0.85, "", "brNDC");
	leg1->SetFillStyle(0);
	leg1->SetBorderSize(0);
	leg1->AddEntry(gr1, "GEANT", "");
//	leg1->Draw();	

	TLegend *leg2 = new TLegend(0.45, 0.4, 0.60, 0.55, "", "brNDC");
	leg2->SetFillStyle(0);
	leg2->SetBorderSize(0);
	leg2->AddEntry(gr2, "EFF", "");
//	leg2->Draw();	

	TLegend *leg3 = new TLegend(0.25, 0.35, 0.40, 0.5, "", "brNDC");
	leg3->SetFillStyle(0);
	leg3->SetBorderSize(0);
	leg3->AddEntry(gr4, "EFF+ESC", "");
//	leg3->Draw();	

	TLegend *leg4 = new TLegend(0.45, 0.3, 0.60, 0.45, "", "brNDC");
	leg4->SetFillStyle(0);
	leg4->SetBorderSize(0);
	leg4->AddEntry(etomu, "e^{+}e^{-}#rightarrow #mu^{+}#mu^{-}", "");
//	leg4->Draw();	

	app->Run();




	return 0;
}