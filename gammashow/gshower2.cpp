/*

Rewrite of gshower.cpp:

(1) Delete extraneous routines (HE, E-2) (done)
(2) Begin in-ice (done)

Binning, number of points is fixed by effective area data. 30 points from ~100 GeV.

Probably want to use this code for event rates & significance. Use gshower for fluxes (more points)

Note: IC80 effective areas derived from integrated event rates. Therefore sensitive to technique used
to solve the integral equation. Currently technique is Riemann sum. Possibly redo using Trapeziod?

Trapezoid/Simpson's rule give very ugly effective areas (highly oscillatory). Therefore stick
with riemann sum for now. 

(3) Split muon pairs for energy loss, rebin. (done)

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

const double pi = TMath::Pi();

//-----Basic constants

	const double specindex = 2.2; //from 1.5 to 2.6 in steps of 0.1 (2.2 for VJ)
	const double emax = 300.;  //cutoff for spectrum (300 or 50 for VJ)

	const int n = 31; //481, 961, 1921 (31?)
	const double de = 0.1; //fixed by Kappes' background binning

	const double zenith =0.707;//cosine of zenith; 0.342(70);  .866(30), 0.707(45), 0.5(60)

	const double photonorm = 1.89E-11;//Gamma normalization (1.89E-11 for VJ)
	const double mcfact = 1.; //HHS pion MC is factor ~1.25 bigger than analytical at HE where decay, E-loss unimportant

	const double equilfact = 0.567; //DHH have 0.5, cascadeconstants give 0.567

//	const double rcross = 505.; //mb
	const double rcross = 650.; 	
	const double radleng = 2.4E4/rcross; //rad length in g/cm2


//------Pion constants
	const double pithresheng = 0.115; //TeV; above this, pi interaction>decay
	const double rat = pow(105.658/139.57,2.);//ratio^2 of muon and pion masses

	const double avmass = 14.5; //I think the numbers are cooked so that 14.5 is the right factor--want ~1.5mb total
//	const double avmass = pow(14.5, 0.91);
	const double zpipi = 0.3;
	const double piacross = 198.; //mb
	const double gammncross = 0.1; //mb
	const double zgammpi = 2./3.;

	const double piinelast = 0.25; //paper says 1/4. Know pi->mu =0.8!
	
	double ratfact = (1.-pow(rat,3.))/(3.-3.*rat);
	double loratfact = (1.-pow(rat,2.))/(2.-2.*rat);

	double pilengabs = 2.4E4/((1-zpipi)*piacross); //Lambda_pi in g/cm2
	double gammpileng = 2.4E4/(avmass*gammncross);

	double interactfact = pilengabs*zgammpi/gammpileng;
	double evofact = radleng/pilengabs;


// Meson constants

	const double kmass = 493.677; //K+-: 493.677, pi+-:139.57, masses easier in MeV
	const double kthresheng = 0.850; //TeV; K:850 , pi: 0.115
	const double krat = pow(105.658/kmass,2.);//ratio^2 of muon and meson masses (both in MeV)


//	const double kngammcross = 0.01; //mb; gp->k+x 10ub? 5ub?
	const double kagammcross = gammncross*avmass*pow(139.57/kmass, 2.);
	const double zgammk = 2./3.;

	const double kinelast = 0.5*(1.+krat)*0.32; //muon <x>

	double klengabs = 200.; //attenuation length in g/cm2; K:180, pi:160
	double gammkleng = 2.4E4/kagammcross;
	
	double kratfact = (1.-pow(krat,3.))/(3.-3.*krat);
	double lokratfact = (1.-pow(krat,2.))/(2.-2.*krat);

//	double kinteractfact = (1./(1.-zkk))*kgammncross*zgammk/kncross;
//	double interactfact = (1./(1.-zpipi))*gammncross*zgammpi/piacross;
//	double evofact = (1-zpipi)*(avmass*pincross/rcross);
//	double kevofact = (1-zkk)*(kavmass*pincross/rcross);

	double kevofact = radleng/klengabs;
	double kinteractfact = zgammk*klengabs/gammkleng;


//-----General spectra constants

	const double lambda1[17] = {0.813, 0.576, 0.389, 0.235, 0.108, 0.0, -0.092, -0.171, -0.239, -0.298, -0.350, -0.395, -0.435, -0.470, -0.5, -0.526, 0.550};
	const double lambda2[17] = {-2.201, -2.055, -1.953, -1.878, -1.824, -1.7868, -1.760, -1.744, -1.734, -1.732, -1.734, -1.741, -1.751, -1.762, -1.780, -1.797, -1.816};
	
	const int lambindex = 17 - (32 - 10*specindex); //call lamb1,2[lambindex] for right parameters

	const double radchange = 9./9.; // 1 for radiation=brem, 9/7 for radiation = pairprod

	const double lamb1 = lambda1[lambindex]*radchange;
	const double lamb2 = lambda2[lambindex]*radchange;

	const double paircs = 0.7733*radchange; //pair prod prob per rad length


	const double alph = paircs + lamb1;
	const double bet = paircs + lamb2;

	const double genconstant = (paircs + lamb1)*(paircs + lamb2)/(lamb2-lamb1); //per (rad length)^2

	
	const double icedepth = 680./zenith; //g/cm2  
	const double icedepthnum = icedepth/radleng;//number of radiation lengths from atmo to ice

	const double gammnpileng = gammpileng/radleng; //in radiation lengths
	const double gammnkleng = gammkleng/radleng;

	const double pileng = pilengabs/radleng ;//in radiation lengths, E-2.7 zpipi
	const double kleng = klengabs/radleng;

	const double genpiratfact = (1.-pow(rat, specindex+1.))/((specindex+1.)*(1.-rat));
	const double genkratfact = (1.-pow(krat, specindex+1.))/((specindex+1.)*(1.-krat));

	const double genmcfact = 1.2;
	const double radcorrect = 8./7.; //correct length in tmax to avg length not rad or pair

	const int terms = 100;

//------Muon pair constants
	const double zee = 7.2; //7.2 is air avg
	const double ay = 14.7; //14.7 is air avg
		
	const double me = 0.511E-6; //electron mass in TeV
	const double mu = 0.105E-3; //muon mass in TeV
	const double alpha = 1./(137.036);
	const double classrad = 1.859E-3; //classical radius of muon squared in mb
	//const double classrad = 1.859E-30*pow(mu/me, 2.);
	
	const double zeered = pow(zee*alpha, 2.);
	const double eff = 1.202*zeered - 1.0369*zeered*zeered + 1.008*pow(zeered, 3.)/(1.+zeered);
	
	const double m = mu;
	
	const double sqrte = 1.6487;
	
	//pair misreconstruction
	
	const int q = 240;
	const double pairde = de/8.;
	
//-------Ice constants

	const double range = (1400.E2/zenith)*0.92; //range is in g/cm2 = depth*density
	const double alphar = 2.59E-6;
	const double beta = 3.63E-6;
	const double drat = exp(-1.*beta*range); //dE_ice/dE_surf
	
	
//---------Event Rate & Significance

	const double angres = 1.; //angular resolution in degrees, diameter!
	const double sourcesize = 2.; //source diameter (0 point, 2 velajr)
	const double bindi = 1.6*pow(angres*angres + sourcesize*sourcesize,0.5);//degrees

	const double solid = 2.*pi*(1. - cos( 0.5*bindi*TMath::DegToRad() ));
	const double binarea = solid*129600./(4.*pi*pi);

//	const double binarea = binrad*binrad*pi; //bin in square deg
//	const double binarea = 1.;	
	
//	const double solid = binarea*4.*3.14159*(1./41252.96); //1st # is # of sq. deg

	const double years = 10.;
	const double icetime = 3.1557E7; // s, 1 year
//	const double icetime = 10.; 

					//muon areas in km2

	const double ic80_0[30] = {0.00921524,0.0541854,0.358172,0.829176,0.91783,0.859948,0.835427,0.825464,
								0.814068,0.806992,0.811272,0.810732,0.82457,0.823537,0.869003,
								0.899253,0.951343,0.757751,0.773057,0.958224,1.2731,0.776643,
								0.634353,0.877282,1.622,0,1.85936,0,0,0};

	const double ic80_30[30] = {0.0491908,0.189943,0.446932,0.782435,0.939407,0.898962,0.854311,0.8539,
								0.849247,0.867357,0.859868,0.858154,0.854355,0.897646,0.924585,
								0.926626,0.874329,0.87613,1.05927,0.85799,0.703524,1.10115,0.50501,
								1.23983,0.572418,0,1.9653,0,0,0};

	const double ic80_45[30] = {0.110642,0.275404,0.52986,0.846921,1.03893,1.01591,0.957614,0.949949,
								0.933375,0.93715,0.947172,0.928199,0.962862,0.943725,0.960085,
								0.975073,1.09565,0.875846,1.18016,0.99432,1.04397,0.844461,
								0.965642,2.12907,1.30839,2.41873,2.24054,0,0,0};

//below is test effective area to see what effect of correcting HE effective area is
// 	const double ic80_45[30] = {0.110642,0.275404,0.52986,0.846921,1.03893,1.01591,0.957614,0.949949,
// 								0.933375,0.93715,0.947172,0.928199,0.962862,0.943725,0.960085,
// 								0.975073,1.09565,0.875846,1.18016,0.99432,1.04397,0.844461,
// 								0.965642,1.12907,1.30839,1.41873,1.24054,0,0,0};



struct my_params2{double a; int b;};	


double genpionint(double* x, size_t dim, void* p){

	double depth = x[0];
	double peng = x[1];

	double z = *(double *)p;

//	vector<double> *factorial = (vector<double> *)p;
//	double z = factorial->front();

	double thresheng, massrat, mesonleng;

	if(z==1.0){
		thresheng = pithresheng;
		massrat = rat;
		mesonleng = pileng;	
	}
	else{
		thresheng = kthresheng;
		massrat = krat;
		mesonleng = kleng;	
	}
	
	double decay = depth*peng*zenith/thresheng;
	double edist = 1./(peng-peng*massrat);
	double photons = photonorm*pow(peng, -1.*specindex);

	double delta = depth/decay;
//	delta = 0; //HE

	double dec1 = 0;
	double dec2 = 0;

	double factorialrun = 1.;

	if(specindex==2.){
		dec1 = depth/(delta+1.);
	
 		for(int i=1; i<terms; i++){
			dec2 += pow(depth, i)*pow(lamb2, i-1.)/((factorialrun)*(delta + i));
 			factorialrun *= i;
 		}	
	}

	else{
		for(int i=1; i<terms; i++){ //
			dec1 += pow(depth, i)*pow(lamb1, i-1.)/((factorialrun)*(delta + i)); 
			dec2 += pow(depth, i)*pow(lamb2, i-1.)/((factorialrun)*(delta + i)); 	
			factorialrun *= i;
		}
	}

	double decayfact = dec1/alph - dec2/bet;
	
	double int1 = (1./alph)*(exp(lamb1*depth) - exp(-1.*depth/mesonleng))/(lamb1 + (1./mesonleng));
	double int2 = (1./bet)*(exp(lamb2*depth) - exp(-1.*depth/mesonleng))/(lamb2 + (1./mesonleng));
	
	double interact = int1-int2;

	
	if(decayfact>interact) return interact*edist*photons/decay;
	return decayfact*edist*photons/decay;

}

double genpion(double x, double z){


	double emu = x;

	double zgamm, gammnleng, inelast, massrat;

	if(z==1.0){
		zgamm = zgammpi;
		gammnleng = gammnpileng;
		inelast = piinelast;
		massrat = rat;	
	}

	else{
		zgamm = zgammk;
		gammnleng = gammnkleng;
		inelast = kinelast;
		massrat = krat;				
	}

	double cascadeconstants = genconstant*(zgamm/gammnleng); 

	double tmin = 0;
	//tmax is a number of radiation lengths here. 
	double tmax1 = icedepthnum; //
	double tmax2 = log(emax*inelast/emu)*radcorrect;
	
	double tmax = tmax2;
	if(tmax2>tmax1) tmax = tmax1;

	if (tmax <0.) return 0.;
							
	double result, error;	
		
	double xl[2] = {tmin, emu};
	double xu[2] = {tmax, emu/massrat};	
	
	gsl_monte_function F;
	
	F.f = &genpionint;
	F.dim = 2;
	F.params = &z; 			
	
	const gsl_rng_type *T;
	gsl_rng *r;
 
	size_t calls = 1000;
   
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
	
	double engdep_int = result;	//	

	return cascadeconstants*engdep_int*genmcfact;
}

double cross_sec(double x, void*p){ //GEANT muon pair cross-section*x^gamma

//	double emu  = *(double *)p;

	struct my_params2 *params = ( struct my_params2 *) p;
	double emu = params->a;
	double tot = params->b; //misreconstruction parameter

	double kay = emu/x;
		if (tot==1) kay = emu; //right? i think so...

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
			if (integrand<0.) integrand = 0.;


//	if (specindex==2.) return integrand*pow(xplus, specindex-1)*log(emax*x/emu);
	if (tot==1) return integrand;
	return integrand*pow(xplus, specindex-1); //mb
}

double stanev_cross(double x, void *p){
 //	double emu  = *(double *)p;

	struct my_params2 *params = ( struct my_params2 *) p;
	double emu = params->a;
	double tot = params->b; //misreconstruction parameter

	double kay = emu/x;
		if (tot==1) kay = emu; //right? i think so...

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
			if (integrand<0.) integrand = 0.;


//	if (specindex==2.) return integrand*pow(xplus, specindex-1)*log(emax*x/emu);
	if (tot==1) return integrand;
	return integrand*pow(xplus, specindex-1); //mb


}

double genpair(double x, double y, int a){ //a is a marker for total cross or weighted cross (re. misreconstruction)

	double emu = x;
	double photons = y;

	double num_mu = 2.;
	if(a==1) num_mu = 1.; //a=1 is misreconstruct

	//tmax is a number of radiation lengths here. 
	double tmax1 = icedepthnum;
	double tmax2 = log(emax/(emu*num_mu))*radcorrect;
	
	double tmax = tmax2;
	if(tmax2>tmax1) tmax = tmax1;
	if(tmax<0.) return 0.;

	double depthint1 = 0;
	double depthint2 = 0;

	if(specindex==2){ 		
		depthint1 = tmax;
		depthint2 = exp(lamb2*tmax)/lamb2 - 1./lamb2;
	}
	else{
		depthint1 = exp(lamb1*tmax)/lamb1 - 1./lamb1;
		depthint2 = exp(lamb2*tmax)/lamb2 - 1./lamb2;
	}

	double depthint = genconstant*(depthint1/alph - depthint2/bet);

	double genpairconst = num_mu*photons*radleng/2.4E4;
	
	double xmin = 1E-4;
	double xmax = 1.-xmin;
	
//	double par1 = emu;
	struct my_params2 par1 = {emu, a};
						
	double result, error;	
		
	double xl[1] = {xmin};
	double xu[1] = {xmax};

	gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
		
	gsl_function F;
	F.function = &cross_sec;
	F.params = &par1;
				
	gsl_integration_qags (&F, xl[0], xu[0], 1e-10, 1e-10, 5000, w, &result, &error);
		
	gsl_integration_workspace_free (w);
	
	double cross_int = result;	//
	
	return genpairconst*depthint*cross_int;

}


double propagateeng(double x, int i){ //detector<->surface

	double eng = x;
	
	double seng;
	
	if(i==1){ //detector->surface
		seng = (1./beta)*( ((alphar+beta*eng)*exp(beta*range)) - alphar);
	}
	else{ //surface->detector
		seng = (1./beta)*( ((alphar+beta*eng)*exp(-beta*range)) - alphar);
	}	
	
	return seng;
}

double eventrate(double* x, double* y, double* dz ,int a){ // int (xy dz) from a to n
	
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
// 		event += x[j]*y[j]*dz[j]*intfact;	
// 	}

	for(int j=a; j<n; j++){ //riemann sum
		intfact = 1.;
		event += x[j]*y[j]*dz[j]*intfact;
	}

	return event;
}

void printout(double* x, int a){
	cout <<"Printing array..." <<endl;
	for(int i=0; i<a; i++){
		cout <<x[i] <<", ";
	}
	cout << endl;
}

int main(int argc, char **argv){

	int i;

 	TStopwatch time;
 	time. Start();

	double detenergy[n], energy[n], photons[n];
	double pions[n], pairs[n], signal[n];
	double gaisser[n];
	double kaons[n], ksignal[n];

	double icarea[n], dice_energy[n];
	double signal_events[n], back_events[n];
	double ksignal_events[n];

	double sigma[n], sigmarker5[n], sigmarker3[n];
	double ksigma[n];
	
	for(i=0; i<n; i++){
	
		detenergy[i] = pow(10, -0.95 + de*i);  //energy at detector, binning matches Kappes' event rates
		energy[i] = propagateeng(detenergy[i], 1);  //corresponding energy at the surface
//		energy[i] = pow(10, -1 + de*i);//to make flux plot

		photons[i] = photonorm*pow(energy[i], -1.*specindex);
		
		pions[i] = genpion(energy[i], 1.0);
		kaons[i] = genpion(energy[i], 2.0)*0.6354; //1.0 is pions, anything else is kaons
		pairs[i] = genpair(energy[i], photons[i], 0); //1 is misreconstruct (correct)

//		pions[i] = 1.;
//		kaons[i] = 1.;
//		pairs[i] = 1.;
		
		signal[i] = pions[i]+pairs[i];
		ksignal[i] = pions[i]+pairs[i]+kaons[i];

//		signal[i] = pairs[i];
//		ksignal[i] = pairs[i];
		
		gaisser[i] =solid*1000*0.14*pow(1000*energy[i],-2.715)* ( (1./(1.+(1.06*energy[i]*zenith/0.115))) + (0.054/(1.+(1.13*energy[i]*zenith/0.85))) ) ;

//		icarea[i] = ic80_0[i];
//		icarea[i] = ic80_30[i];
		icarea[i] = ic80_45[i];
		
//		ratio[i] = detenergy[i]/energy[i];
//		if (ratio[i]==0) dice_energy[i] = 0;
//		else dice_energy[i] = log(10.)*detenergy[i]*de;
		dice_energy[i] = log(10.)*detenergy[i]*de;
	}

	double pairdetenergy[q];
	double mui1, mui2;
	double muf1, muf2;
	double pairenergy[q];
	double pairs2[q], photons2[q];
	int j;
	double pairs2_rebin[n];
	
	for(i=0; i<q; i++){ //muon pairs
		pairenergy[i] = pow(10, -0.3 + i*pairde ); //start at 500 GeV?
		photons2[i] = photonorm*pow(pairenergy[i], -1.*specindex);				
		pairs2[i] = genpair(pairenergy[i], photons2[i], 1); //1 is misreconstruct (correct)

		mui1 = pairenergy[i]*0.75;
		mui2 = pairenergy[i]*0.25;
		
		muf1 = propagateeng(mui1, 2); //energy at detector
		muf2 = propagateeng(mui2, 2);
		
		if(muf1<0.01) muf1 = 0; //range out at 10 GeV?
		if(muf2<0.01) muf2 = 0;
		
		pairdetenergy[i] = muf1+ muf2;
	
	}

	int k=0;
	int p=0;
	int l;
	double detenergybound;
	double ksig_rebin[n], sig_rebin[n];

	for(i=0; i<n; i++){
		
		detenergybound = detenergy[n-1];	
		if(i!=(n-1)) detenergybound = 0.5*detenergy[i+1]+0.5*detenergy[i];
		
		pairs2_rebin[i] = 0;
		p = k;
		l = k;
		
		for(j=k; j<q; j++){
			if( pairdetenergy[j] > detenergybound) break;		
			pairs2_rebin[i] += pairs2[j];
			p++;
		}	
		k = p;
		pairs2_rebin[i] = pairs2_rebin[i]/(p-l);
		
		ksig_rebin[i] = pions[i]+pairs2_rebin[i]+kaons[i];
		sig_rebin[i] = pions[i]+pairs2_rebin[i];
	}
	
	for(i=0; i<n; i++){
	
		signal_events[i] = 0.72*eventrate( sig_rebin, icarea, dice_energy, i)*1.E10*icetime/drat;
		ksignal_events[i] = 0.72*eventrate( ksig_rebin, icarea, dice_energy, i)*1.E10*icetime/drat;
		back_events[i] = eventrate(gaisser, icarea, dice_energy, i)*1.E10*icetime/drat;

		sigma[i] = years*signal_events[i]/pow(years*back_events[i], 0.5);
		ksigma[i] = years*ksignal_events[i]/pow(years*back_events[i], 0.5);
		
		sigmarker5[i] = 5.;
		sigmarker3[i] = 3.;
	}

	double pifrac[n], pairfrac[n], kfrac[n];
	
	for(i=0; i<n; i++){ //fractional event rates due to different components

		pifrac[i] = 0.72*eventrate( pions, icarea, dice_energy, i)*1.E10*icetime/drat;	
		pairfrac[i] = 0.72*eventrate( pairs2_rebin, icarea, dice_energy, i)*1.E10*icetime/drat;	
		kfrac[i] = 0.72*eventrate( kaons, icarea, dice_energy, i)*1.E10*icetime/drat;
		
		pifrac[i] = pifrac[i]/ksignal_events[i];
		pairfrac[i] = pairfrac[i]/ksignal_events[i];
		kfrac[i] = kfrac[i]/ksignal_events[i];	
	}


	for(i=0; i<n; i++){
	
		pions[i] = pions[i]*pow(energy[i], 2.);
		kaons[i] = kaons[i]*pow(energy[i], 2.);
		pairs[i] = pairs[i]*pow(energy[i], 2.);		
		gaisser[i] = gaisser[i]*pow(energy[i], 2.);
		ksignal[i] = ksignal[i]*pow(energy[i], 2.);		
	}

//	printout(detenergy,n);
	printout(sigma,n);
	printout(ksigma,n);
//	printout(back_events,n);
//	printout(gaisser,n);

 	time.Stop();
 	
 	cout <<time.RealTime() <<" seconds"<<endl;

	cout.precision(8);

//	cout <<back_events[0] <<endl; 
	cout <<ksignal_events[0]*(1./0.72) <<endl;
//	cout <<signal_events[0] <<endl;

//	cout <<sigma[0] <<endl;
//	cout <<ksigma[0] <<endl;	
	
	TApplication *app = new TApplication("app", &argc, argv); //needed to run ROOT externally

//--------Flux Plot------//
	
	char Fluxtitle[50];
	sprintf(Fluxtitle,"dN/dE Muons, zenith=%.0f deg, norm=%.1E, index=-%.1f, Emax=%.0E TeV",acos(zenith)/0.01745329, photonorm, specindex, emax);
//	sprintf(Fluxtitle, "Vela Jr");
//	sprintf(Fluxtitle, "");
	
	TGraph *pidiff = new TGraph(n, energy, pions);
	pidiff->SetLineWidth(2);
	pidiff->SetLineColor(1);
	pidiff->SetLineStyle(3);

	TGraph *pairdiff = new TGraph(n, energy, pairs);
	pairdiff->SetLineWidth(2);
	pairdiff->SetLineColor(1);
	pairdiff->SetLineStyle(4);

	TGraph *kdiff = new TGraph(n, energy, kaons);
	kdiff->SetLineWidth(2);
	kdiff->SetLineColor(1);
	kdiff->SetLineStyle(6);
	
	TGraph *backdiff = new TGraph(n, energy, gaisser);
	backdiff->SetLineWidth(2);
	backdiff->SetLineColor(1);
	backdiff->SetLineStyle(2);

// 	TGraph *pairdiff2 = new TGraph(q, pairenergy, pairs2);
// 	pairdiff2->SetLineWidth(2);
// 	pairdiff2->SetLineColor(4);
// 	pairdiff2->SetLineStyle(1);

	TGraph *sigdiff = new TGraph(n, energy, ksignal);
	sigdiff->SetLineWidth(2);
	sigdiff->SetLineColor(1);
	sigdiff->SetLineStyle(1);

	TCanvas *c1 = new TCanvas("events1","stuff1",500,500 );
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,10, 10);

	c1->cd();

	pad1->SetTitle("");
	pad1->Draw();
	pad1->cd();
	pad1->SetLogx();
	pad1->SetLogy();
//	pad1->SetGridx();
//	pad1->SetGridy();
	pad1->SetTicks();

	pidiff->SetMaximum(1E-8);
	pidiff->SetMinimum(1E-17);

//	pidiff->SetTitle(Fluxtitle);
	pidiff->SetTitle("");

	pidiff->Draw("Al");
	pidiff->GetXaxis()->SetLimits(1.E-1, 1.E2);
	pidiff->GetXaxis()->SetTitle("E_{#mu} (TeV)");
 	pidiff->GetXaxis()->CenterTitle(1);
 	pidiff->GetXaxis()->SetTitleOffset(1.2);
	pidiff->GetYaxis()->SetTitle("E^{2} dN/dE (TeV cm^{-2} s^{-1})");
	pidiff->GetYaxis()->CenterTitle(1);
 	pidiff->GetYaxis()->SetTitleOffset(1.6);

	pairdiff->Draw("same");
	kdiff->Draw("same");
	backdiff->Draw("same");
	sigdiff->Draw("same");
//	pairdiff2->Draw("same");

	TLegend *leg1 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
	leg1->SetFillStyle(0);
	leg1->SetBorderSize(0);
//	leg1->SetShadowColor("1");
 	leg1->AddEntry(sigdiff, "Total signal", "l");	
 	leg1->AddEntry(pidiff, "Pions", "l");
 	leg1->AddEntry(pairdiff, "Pairs", "l");
 	leg1->AddEntry(kdiff, "Kaons", "l");
 	leg1->AddEntry(backdiff, "Background", "l");
 	leg1->Draw();
	
// 	TLegend *leg11 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
// 	leg11->SetFillStyle(0);
// 	leg11->SetBorderSize(0);
// //	leg11->SetShadowColor("1");
//  	leg11->AddEntry(sigdiff, "Total signal", "");	
//  	leg11->Draw();
// 
// 	TLegend *leg12 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
// 	leg12->SetFillStyle(0);
// 	leg12->SetBorderSize(0);
// //	leg12->SetShadowColor("1");
//  	leg12->AddEntry(pidiff, "Pions", "");
//  	leg12->Draw();
// 
// 	TLegend *leg13 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
// 	leg13->SetFillStyle(0);
// 	leg13->SetBorderSize(0);
// //	leg13->SetShadowColor("1");
//  	leg13->AddEntry(pairdiff, "Pairs", "");
//  	leg13->Draw();
//  	
//  	TLegend *leg14 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
// 	leg14->SetFillStyle(0);
// 	leg14->SetBorderSize(0);
// //	leg14->SetShadowColor("1");
//  	leg14->AddEntry(kdiff, "Kaons", "");
//  	leg14->Draw();
//  	
//  		TLegend *leg15 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
// 	leg15->SetFillStyle(0);
// 	leg15->SetBorderSize(0);
// //	leg15->SetShadowColor("1");
//  	leg15->AddEntry(backdiff, "Background", "");
//  	leg15->Draw();
 	
//----------Significance Plot--------//

	char Sigtitle[50];
//	sprintf(Sigtitle,"Sigma, %.0f years, %.1f deg angres, zenith=%.0f deg, norm=%.1E, index=-%.1f, Emax=%.0E TeV",years, angres,acos(zenith)/0.01745329, photonorm, specindex, emax);
//	sprintf(Sigtitle,"Pt. Pevatron, Emax= %.0f, %.0f years, bin = %.1f deg2", emax, years, binarea);
//	sprintf(Sigtitle,"Pt. Pevatron, %.1E norm, Emax= %.0f, %.0f years, %.1f deg angres, E-%.1f", emax, years, angres, specindex);
	sprintf(Sigtitle,"%.1E norm, Emax= %.0f, %.0f years, %.1f deg angres, %.1f deg source diam, E-%.1f",photonorm, emax, years, angres, sourcesize, specindex);
	

	TGraph *signif = new TGraph(n, detenergy, sigma);
	signif->SetLineColor(1);
	signif->SetLineWidth(2);
	signif->SetLineStyle(1);	

	TGraph *ksignif = new TGraph(n, detenergy, ksigma);
	ksignif->SetLineColor(1);
	ksignif->SetLineWidth(2);
	ksignif->SetLineStyle(2);

	TGraph *discovery = new TGraph(n, detenergy, sigmarker5);
	discovery->SetLineColor(2);
	discovery->SetLineWidth(2);
	discovery->SetLineStyle(1);	

	TGraph *disco3 = new TGraph(n, detenergy, sigmarker3);
	disco3->SetLineColor(2);
	disco3->SetLineWidth(2);
	disco3->SetLineStyle(1);	

	TCanvas *c3 = new TCanvas("events3","stuff3",500,500 );
	TPad *pad3 = new TPad("pad3","pad3",0.0,0.0,1.0,1.0,10);

	c3->cd();

	pad3->SetTitle("");
	pad3->Draw();
	pad3->cd();
	pad3->SetLogx();
	pad3->SetLogy();
	pad3->SetGridx();
	pad3->SetGridy();
	pad3->SetTicks();

	signif->SetMinimum(1E-1);
	signif->SetMaximum(10);
	signif->SetTitle(Sigtitle);
	
	signif->Draw("Al");
	ksignif->Draw("same");
	
	discovery->Draw("same");
	disco3->Draw("same");
	
	signif->GetXaxis()->SetLimits(0.1,100);
	signif->GetXaxis()->SetTitle("In-ice Energy(TeV)");	
	signif->GetXaxis()->CenterTitle(1);
 	signif->GetXaxis()->SetTitleOffset(1.2);
	signif->GetYaxis()->SetTitle("Sigma");	
 	signif->GetYaxis()->CenterTitle(1);
 	signif->GetYaxis()->SetTitleOffset(1.2);
	
	
	TLegend *leg3 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
	leg3->SetFillStyle(0);
	leg3->SetBorderSize(0);
//	leg3->SetShadowColor("1");
 	leg3->AddEntry(signif, "No kaons", "l");
 	leg3->AddEntry(ksignif, "Kaons", "l");
	leg3->Draw();


	
//-------Plot 5--Event rates---------//

	char Evtitle[50];
//	sprintf(Evtitle,"Sigma, %.0f years, %.1f deg angres, zenith=%.0f deg, norm=%.1E, index=-%.1f, Emax=%.0E TeV",years, angres,acos(zenith)/0.01745329, photonorm, specindex, emax);
	sprintf(Evtitle,"Event rate per year, %.1f deg2 bin, zenith=%.0f deg", binarea ,acos(zenith)/0.01745329);
	

	TGraph *signalev = new TGraph(n, detenergy, signal_events);
	signalev->SetLineColor(1);
	signalev->SetLineWidth(2);
	signalev->SetLineStyle(1);	

	TGraph *testsignal = new TGraph(n, detenergy, ksignal_events);
	testsignal->SetLineColor(1);
	testsignal->SetLineWidth(2);
	testsignal->SetLineStyle(2);	
	testsignal->SetMarkerStyle(27);
	testsignal->SetMarkerSize(2);

	TGraph *backev = new TGraph(n, detenergy, back_events);
	backev->SetLineColor(1);
	backev->SetLineWidth(2);
	backev->SetLineStyle(1);

	double corsika_eng[n], diffae[n];
	
	double corsika_0deg[n] = {9299584,9273409,9123834,8203224,6307365,4527258,3179188,2170061,1436419,926749.3,584848.3,360542,218892.5,130375.6,77354.68,44481.73,24836.1,13004.07,7701.981,4688.162,2623.138,1116.229,613.926,390.6802,223.2458,55.81146,55.81146,0,0,0};
	double corsika_45deg[n] = {6130330,5999118,5667765,5045635,4115593,3095973,2246244,1595817,1096416,733864.9,476462.4,299819.1,186521.9,112069.4,67141.18,39681.95,23273.38,12613.39,7757.793,4074.236,2344.081,1339.475,892.9833,613.926,279.0573,167.4344,55.81146,0,0,0};
	double corsika_30deg[n] = {7949114,7850161,7471703,6625267,5277978,3874543,2764062,1932528,
									1307774,861226.6,546896.5,340170.8,207897.7,125966.5,73782.75,41914.4,
									23329.19,13283.13,7590.358,3739.368,2009.212,1227.852,558.1146,390.6802,
									167.4344,111.6229,111.6229,55.81146,55.81146,0};


	for(i=0; i<n; i++){ 
		corsika_eng[i] = pow(10, -0.95 + i*0.1);	
		diffae[i] = -corsika_30deg[i+1]+corsika_30deg[i];			
	}


	TGraph *mcback_0 = new TGraph(n, corsika_eng, corsika_0deg);
//	TGraph *mcback_0 = new TGraph(n, corsika_eng, diffae);
	mcback_0->SetLineColor(4);
	mcback_0->SetLineWidth(2);	
	mcback_0->SetLineStyle(2);
	
	TGraph *mcback_30 = new TGraph(n, corsika_eng, corsika_30deg);
	mcback_30->SetLineColor(4);
	mcback_30->SetLineWidth(2);	
	mcback_30->SetLineStyle(2);
	
	TGraph *mcback_45 = new TGraph(n, corsika_eng, corsika_45deg);
	mcback_45->SetLineColor(4);
	mcback_45->SetLineWidth(2);	
	mcback_45->SetLineStyle(2);

	//--Kappes' points --taken from histogram so only plot points
	//--pions only?, 30deg, 10^-11 E-2, 30deg or 0deg?
	double kappeseng[17] = {0.115, 0.139, 0.177, 0.222, 0.279, 0.350, 0.440, 0.559, 0.702, 0.882, 1.121, 1.408, 1.768, 2.221, 2.854, 3.504, 6.267};
	double kappessig[17] = {550.740, 532.194, 472.166, 404.802, 347.048, 282.653, 255.085, 230.206, 187.491, 145.064, 112.238, 82.496, 51.985, 38.210, 28.085, 13.693, 4.661};
	
	TGraph *kappes = new TGraph(17, kappeseng, kappessig);
	kappes->SetLineColor(1);
	kappes->SetLineWidth(2);
	kappes->SetLineStyle(2);	
	kappes->SetMarkerStyle(27);
	kappes->SetMarkerSize(2);


	TCanvas *c5 = new TCanvas("events5","stuff5",500,500 );
	TPad *pad5 = new TPad("pad5","pad5",0.0,0.0,1.0,1.0,10);

	c5->cd();

	pad5->SetTitle("");
	pad5->Draw();
	pad5->cd();
	pad5->SetLogx();
	pad5->SetLogy();
	pad5->SetGridx();
	pad5->SetGridy();
	pad5->SetTicks();

//	signalev->SetMinimum(1E1);
//	signalev->SetMaximum(5E7);

	signalev->SetMinimum(1E0);
	signalev->SetMaximum(1E7);

	signalev->SetTitle(Evtitle);
	
	signalev->Draw("Al");
	backev->Draw("same");
	testsignal->Draw("same");
//	kappes->Draw("*");
	
//	mcback_0->Draw("same");
//	mcback_30->Draw("same");
//	mcback_45->Draw("same");
	
	signalev->GetXaxis()->SetLimits(0.1,40);
	signalev->GetXaxis()->SetTitle("E_{#mu} in-ice (TeV)");	
	signalev->GetXaxis()->CenterTitle(1);
 	signalev->GetXaxis()->SetTitleOffset(1.2);
	signalev->GetYaxis()->SetTitle("Events/year");	
 	signalev->GetYaxis()->CenterTitle(1);
 	signalev->GetYaxis()->SetTitleOffset(1.2);
	
	
	TLegend *leg5 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
	leg5->SetFillStyle(0);
	leg5->SetBorderSize(0);
//	leg5->SetShadowColor("1");
 	leg5->AddEntry(signalev, "No kaons", "l");
 	leg5->AddEntry(testsignal, "kaons", "l");
// 	leg5->AddEntry(mcback_45, "IC MC", "l");
	leg5->Draw();
	
	
	//----------Fractional Event rates plot--------//

	char Fractitle[50];
//	sprintf(Sigtitle,"Sigma, %.0f years, %.1f deg angres, zenith=%.0f deg, norm=%.1E, index=-%.1f, Emax=%.0E TeV",years, angres,acos(zenith)/0.01745329, photonorm, specindex, emax);
//	sprintf(Sigtitle,"Pt. Pevatron, Emax= %.0f, %.0f years, bin = %.1f deg2", emax, years, binarea);
//	sprintf(Sigtitle,"Pt. Pevatron, %.1E norm, Emax= %.0f, %.0f years, %.1f deg angres, E-%.1f", emax, years, angres, specindex);
	sprintf(Fractitle,"%.1E norm, Emax= %.0f, %.0f years, %.1f deg angres, %.1f deg source diam, E-%.1f",photonorm, emax, years, angres, sourcesize, specindex);
	

	TGraph *piev = new TGraph(n, detenergy, pifrac);
	piev->SetLineColor(1);
	piev->SetLineWidth(2);
	piev->SetLineStyle(1);	

	TGraph *pairev = new TGraph(n, detenergy, pairfrac);
	pairev->SetLineColor(1);
	pairev->SetLineWidth(2);
	pairev->SetLineStyle(2);

	TGraph *kev = new TGraph(n, detenergy, kfrac);
	kev->SetLineColor(1);
	kev->SetLineWidth(2);
	kev->SetLineStyle(3);	

	TCanvas *c4 = new TCanvas("events4","stuff4",500,500 );
	TPad *pad4 = new TPad("pad4","pad4",0.0,0.0,1.0,1.0,10);

	c4->cd();

	pad4->SetTitle("");
	pad4->Draw();
	pad4->cd();
	pad4->SetLogx();
//	pad4->SetLogy();
//	pad4->SetGridx();
	pad4->SetGridy();
	pad4->SetTicks();

	piev->SetMinimum(0.);
	piev->SetMaximum(1.);
	piev->SetTitle("");
	
	piev->Draw("Al");
	pairev->Draw("same");
	
	kev->Draw("same");
	
	piev->GetXaxis()->SetLimits(0.1,40);
	piev->GetXaxis()->SetTitle("E_{#mu} in-ice (TeV)");	
	piev->GetXaxis()->CenterTitle(1);
 	piev->GetXaxis()->SetTitleOffset(1.2);
	piev->GetYaxis()->SetTitle("Event rate fraction");	
 	piev->GetYaxis()->CenterTitle(1);
 	piev->GetYaxis()->SetTitleOffset(1.2);
	
	
	TLegend *leg4 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
	leg4->SetFillStyle(0);
	leg4->SetBorderSize(0);
//	leg4->SetShadowColor("1");
 	leg4->AddEntry(piev, "Pions", "l");
 	leg4->AddEntry(pairev, "Pairs", "l");
 	leg4->AddEntry(kev, "Kaons", "l");
	leg4->Draw();

	
	
	
		
	app->Run();





	return 0;

}
