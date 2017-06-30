/*
Calculate IC80 muon area from Kappes event rate and Gaisser flux

Issue: what quadrature method? 

Kappes event rate is not smooth->Fluctuations at HE. 
Therefore it may be problematic to solve this with an integral equation b/c have to
work backwards from HE. 
Simpson's rule gives oscillating area because odd<->even number of bins. 

Actually, b/c have number of events per bin, can we use standard integral eqn methods? does
this even make sense?

extrapolate event rates to more bins? see if this helps? hard to do b/c log bins
Fuck it, just use riemann sum. 
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

	const int n = 28;
	
//	const double angres = 1.128; //angular resolution
	const double solid = 1.*4.*3.14159*(1./41252.96); //1st # is # of sq. deg
		
	const double years = 1.;
	const double icetime = 3.1557E7; // 1 year
	
	const double deng = 0.1*log(10.)*icetime;
	const double exdeng = 0.05*log(10.)*icetime;
	

	const double alphar = 2.59E-6; //2, 2.59
	const double beta = 3.63E-6; //4.2, 3.63

	const double normconst = 1.;

	const double E1=115;  // all Dima stuff in GeV, mwe
	const double E2=850;
	const double pK=0.054;
	const double a=0.260;
	const double b=0.360e-3;
	const double as=-0.522;
	const double bs=8.893e-3;
	const double A=0.14;
//	const double x0 = 1./();//atmo is 10.8 mwe?
	
	const double corsika_0deg[29] = {9299584.,9273409.,9123834.,8203224.,6307365.,4527258.,3179188.,2170061.,
									1436419,926749.3,584848.3,360542.,218892.5,130375.6,77354.68,44481.73,
									24836.1,13004.07,7701.981,4688.162,2623.138,1116.229,613.926, 
									390.6802,223.2458,55.81146,55.81146,0,0};

	const double corsika_45deg[29] = {6130330,5999118,5667765,5045635,4115593,3095973,2246244,1595817,
									1096416,733864.9,476462.4,299819.1,186521.9,112069.4,67141.18,
									39681.95,23273.38,12613.39,7757.793,4074.236,2344.081,1339.475,
									892.9833,613.926,279.0573,167.4344,55.81146,0.,0.};

	const double corsika_30deg[29] = {7949114,7850161,7471703,6625267,5277978,3874543,2764062,1932528,
									1307774,861226.6,546896.5,340170.8,207897.7,125966.5,73782.75,41914.4,
									23329.19,13283.13,7590.358,3739.368,2009.212,1227.852,558.1146,390.6802,
									167.4344,111.6229,111.6229,55.81146,55.81146};



double gaisser(double x, double y){ //energy and zenith, energy in TeV 

	double eng = x;
	double zenith = y;

	double flux = solid*1000.*0.14*pow(1000.*eng,-2.715)* ( (1./(1.+(1.06*eng*zenith/0.115))) + (0.054/(1.+(1.13*eng*zenith/0.85))) ) ;

	return flux;
}

double surface_eng(double x, double y){ //energy and zenith

	double ef = x;
	double range = (1400.E2/y)*0.92;

	double ei = (ef+ alphar/beta)*exp(beta*range) - alphar/beta;


	return ei;
}

double fluct_eng(double x, double y){ //energy and zenith--return flux not energy--Chirkin

	double ef = x*1000.;

	double xs = 1./(-0.017326+0.114236*pow(y, 1.15043) + 0.0200854*pow(1.-y*y, 1.16714));
//	double xs = x0/y;

	double ei = ((a+b*ef)*exp(b*xs)-a)/b;

	double zfi, zff, aux, aux1, aux2, auf1, auf2, f1, f2, f3, f4, de;

//	0.701, 2.715

	aux1=(bs*a-as*b);
	aux1*=2*aux1;
	aux2=2*b*b*b;
	zfi=(bs*ei*b*(-2*bs*a+4*as*b+bs*ei*b)+aux1*log(a+ei*b))/aux2;
	zff=(bs*ef*b*(-2*bs*a+4*as*b+bs*ef*b)+aux1*log(a+ef*b))/aux2;

	aux=1.1*xs;
	auf1=aux/E1;
	auf2=aux/E2;
	aux1=1/(1+auf1*ei);
	aux2=1/(1+auf2*ei);
	aux=aux1+pK*aux2;

	f1=pow(ei, -2.715-1)*aux;
	f2=pow(ei, -2.715-1)*(-(2.715+1)*aux/ei-auf1*aux1*aux1-pK*auf2*aux2*aux2);
	f3=(-2.715-1)*f2/ei+pow(ei, -2.715-1)*((2.715+1)*aux/(ei*ei)
		      +((2.715+1)/ei)*(auf1*aux1*aux1+pK*auf2*aux2*aux2)
		      +2*auf1*auf1*aux1*aux1*aux1+2*pK*auf2*auf2*aux2*aux2*aux2);
	aux=f2/f1;
	f4=f3/f1-aux*aux;

	aux1=as+bs*ef;
	aux2=as+bs*ei;
	de=exp(b*xs)-(f2/(2*f1))*(aux1*aux1-aux2*aux2)/(a+b*ef)-(f4/2)*(zff-zfi);

	ei=ei-(f2/(2*f1))*(zff-zfi);
	ei = ei/1000.; //convert to TeV

//	cout <<ei/x <<"   "<<de <<endl;

	double gaiss = gaisser(ei, y)*de;

	return gaiss*normconst;
}

double areadet_riem(double eng, double flux, int i, double zenith){ //zenith=1, 0.866, 0.707

	double events = corsika_0deg[i] - corsika_0deg[i+1];
	
	if (zenith < 0.89){ events = corsika_30deg[i] - corsika_30deg[i+1]; }	
	if (zenith < 0.720){ events = corsika_45deg[i] - corsika_45deg[i+1]; }	
	
	
	double icarea = 1E-10*events/(flux*eng*0.1*log(10.)*icetime);

	
	return icarea; //area in km2
}

int main(int argc, char **argv){
	
	double zenith = 1;
	double range = (1400.E2/zenith)*0.92; //range is in g/cm2 = depth*density
	double drat = exp(beta*range); //dE_surf/dE_ice
	int i;

	double eng[n];
	double area[n];

	double surfeng;
	double iceflux[n];
	

	for(i=0; i<n; i++){

		eng[i] = pow(10, -0.95 + i*0.1);//TeV
		surfeng = surface_eng(eng[i], zenith);	
		iceflux[i] = gaisser(surfeng, zenith)*drat;	

		area[i] = 0;
		area[i] = areadet_riem(eng[i], iceflux[i], i, zenith);
//		cout <<area[i] <<endl;		
	}
	
//	cout <<iceflux[31] <<endl;

	double extrapev[2*29-1], eng2[2*n], iceflux2[2*n], areaex[2*n];
	double events;

	extrapev[0] = corsika_0deg[0];

	for(i=1; i<2*29-1; i++){ //pretty crude, divides log bin by 2
		
		if((i-1)%2==0){ 
			extrapev[i] = corsika_0deg[(i-1)/2]; //original event points
		} 	

		else{
			extrapev[i] = 0.5*corsika_0deg[(i/2)-1] + 0.5*corsika_0deg[i/2];
		}	
	}

	for(i=0; i<2*n; i++){
		eng2[i] = pow(10, -0.95 + i*0.05);//TeV
		surfeng = surface_eng(eng2[i], zenith);	
		iceflux2[i] = gaisser(surfeng, zenith)*drat;	

		events = extrapev[i] - extrapev[i+1];
		areaex[i] = 1E-10*events/(iceflux2[i]*eng2[i]*exdeng);
//		cout <<areaex[i] <<endl;	
	}


	int j, k;
 	double sum, intfact;
 
 	for(j=n-1; j>-1; j--){
 	
 //		j = n-1-i;
 		
 		if(j==n-1){ //one bin, take area.
 			area[j] = corsika_0deg[j]/(iceflux[j]*eng[j]*deng);
 			area[j] = 0.;
 		}
 
 		else{
 			if(j==n-2){ //two bins, trapezoidal rule.
 				area[j] = (corsika_0deg[j]-0.5*deng*area[n-1]*iceflux[n-1]*eng[n-1])/(0.5*iceflux[j]*eng[j]*deng);
 			}	
 
 			else{	
// 				if((n-1-j-0)%2!=0){ //odd # of bins, Simpson's rule
// 					sum = (1./3.)*area[n-1]*iceflux[n-1]*eng[n-1];
// 				
// 					for(k=n-2; k>j; k--){	
// 						intfact = 4./3.;
// 	//					term = n-1-k;
// 						if ((n-1-k)%2==0) intfact = 2./3.;
// 	//					intfact = 1.;
// 						sum += intfact*area[k]*iceflux[k]*eng[k];				
// 					}
// 					area[j] = (corsika_0deg[j] - sum*deng)/((1./3.)*deng*iceflux[j]*eng[j]);	
// 				}
//				else{ //even # of bins, trap?
 					sum = (1./2.)*area[n-1]*iceflux[n-1]*eng[n-1];
 				
 					for(k=n-2; k>j; k--){	
 						intfact = 1.;
 						sum += intfact*area[k]*iceflux[k]*eng[k];				
 					}
 					area[j] = (corsika_0deg[j] - sum*deng)/((1./2.)*deng*iceflux[j]*eng[j]);									
 				
// 				}
 			}
 		}
 		area[j] = 0.25*1E-10*area[j];
 //		cout <<area[j] <<endl;
 	}
// 
// //	cout <<area[-1] <<endl;
// 
// 	double intflux[n];
// 	
// 	for(i=0; i<n; i++){
// 	
// 		sum = 0;
// 		for(j=i; j<n; j++){
// 		
// 			intfact = 1.;
// 			
// 			if((n-1-i)%2!=0){		
// 				intfact = 4./3.;
// 				if((i-j)%2==0) {intfact = 2./3.;}
// 				if(j==i) {intfact = 1./3.;}
// 				if(j==n-1) {intfact = 1./3.;}
// 			}	
// 			else{
// 				intfact = 1.;
// 				if(j==i){intfact = 0.5;}
// 				if(j==n-1){intfact = 0.5;}		
// 			}
// 
// 			sum += intfact*area[j]*iceflux[j]*eng[j]*deng;
// 		
// 		}
// 	
// 		intflux[i] = sum;
// //		cout <<intflux[i] <<endl;
// 	}
// 
// 
// //	int j, k;
// //	double sum, intfact;
// 
// 	for(j=2*n-1; j>-1; j--){
// 	
// //		j = n-1-i;
// 		
// 		if(j==2*n-1){ //one bin, take area.
// 			areaex[j] = extrapev[j]/(iceflux2[j]*eng2[j]*exdeng);
// 			areaex[j] = 0.;
// 		}
// 
// 		else{
// 			if(j==2*n-2){ //two bins, trapezoidal rule.
// 				areaex[j] = (extrapev[j]-0.5*exdeng*areaex[2*n-1]*iceflux2[2*n-1]*eng2[2*n-1])/(0.5*iceflux2[j]*eng2[j]*exdeng);
// 			}	
// 
// 			else{	
// 				if((2*n-1-j-0)%2!=0){ //odd # of bins, Simpson's rule
// 					sum = (1./3.)*areaex[2*n-1]*iceflux2[2*n-1]*eng2[2*n-1];
// 				
// 					for(k=2*n-2; k>j; k--){	
// 						intfact = 4./3.;
// 	//					term = n-1-k;
// 						if ((2*n-1-k)%2==0) intfact = 2./3.;
// 	//					intfact = 1.;
// 						sum += intfact*areaex[k]*iceflux2[k]*eng2[k];				
// 					}
// 					areaex[j] = (extrapev[j] - sum*exdeng)/((1./3.)*exdeng*iceflux2[j]*eng2[j]);	
// 				}
// 				else{ //even # of bins, trap?
// 					sum = (1./2.)*areaex[2*n-1]*iceflux2[2*n-1]*eng2[2*n-1];
// 				
// 					for(k=2*n-2; k>j; k--){	
// 						intfact = 1.;
// 						sum += intfact*areaex[k]*iceflux2[k]*eng2[k];				
// 					}
// 					areaex[j] = (extrapev[j] - sum*exdeng)/((1./2.)*exdeng*iceflux2[j]*eng2[j]);									
// 				
// 				}
// 			}
// 		}
// //		area[j] = area[j];
// //		cout <<area[j] <<endl;
// 	}
// 
// //	cout <<area[-1] <<endl;
// 
// //	double intflux[n];
// 	
// 	for(i=0; i<2*n; i++){
// 	
// 		sum = 0;
// 		for(j=i; j<2*n; j++){
// 		
// 			intfact = 1.;
// 			
// 			if((2*n-1-i)%2!=0){		
// 				intfact = 4./3.;
// 				if((i-j)%2==0) {intfact = 2./3.;}
// 				if(j==i) {intfact = 1./3.;}
// 				if(j==n-1) {intfact = 1./3.;}
// 			}	
// 			else{
// 				intfact = 1.;
// 				if(j==i){intfact = 0.5;}
// 				if(j==n-1){intfact = 0.5;}		
// 			}
// 
// 			sum += intfact*areaex[j]*iceflux2[j]*eng2[j]*exdeng;
// 		
// 		}
// 	
// 		intflux[i] = sum;
// //		cout <<intflux[i] <<endl;
// 	}


	TGraph *effarea = new TGraph(n, eng, area);
	effarea->SetLineWidth(2);
	effarea->SetLineColor(2);

	
 	TH1 *area1 = new TH1D("area1","",30,-1,2);
	
 	for(i=0; i<n; i++){ 
		area1->SetBinContent(i+1, area[n-1-i]);	
	}

 	TH1 *area1ex = new TH1D("area1ex","",60,-1,2);
	
 	for(i=0; i<2*n-1; i++){ 
		area1ex->SetBinContent(i+1, areaex[i]*1E-10);	
	}
	

	zenith = 0.707;
	range = (1400.E2/zenith)*0.92; //range is in g/cm2 = depth*density
	drat = exp(beta*range); //dE_surf/dE_ice


	for( i=0; i<n; i++){
		surfeng = surface_eng(eng[i], zenith);	
		iceflux[i] = gaisser(surfeng, zenith)*drat;	
		area[i] = areadet_riem(eng[i], iceflux[i], i, zenith);

	}


 	TH1 *area2 = new TH1D("area2","",30,-1,2);

 	for(i=0; i<n; i++){ 
		area2->SetBinContent(i+1, area[i]);	
	}

	zenith = 0.866;
	range = (1400.E2/zenith)*0.92; //range is in g/cm2 = depth*density
	drat = exp(beta*range); //dE_surf/dE_ice


	for( i=0; i<n; i++){
		surfeng = surface_eng(eng[i], zenith);	
		iceflux[i] = gaisser(surfeng, zenith)*drat;	
		area[i] = areadet_riem(eng[i], iceflux[i], i, zenith);

	}


 	TH1 *area3 = new TH1D("area3","",30,-1,2);

 	for(i=0; i<n; i++){ 
		area3->SetBinContent(i+1, area[i]);	
	}

	TApplication *app = new TApplication("app", &argc, argv); //needed to run ROOT externally

	TCanvas *c3 = new TCanvas("events3","stuff3",500,500 );
	TPad *pad3 = new TPad("pad3","pad3",0.0,0.0,1.0,1.0,10);

	c3->cd();

	pad3->SetTitle("");
	pad3->Draw();
	pad3->cd();
//	pad3->SetLogx();
//	pad3->SetLogy();
	pad3->SetGridx();
	pad3->SetGridy();
	pad3->SetTicks();

	area1->SetMinimum(0);
   	area1->SetMaximum(2.5);
   	area1->SetStats(0);
    area1->SetLineColor(1);
	area1->SetLineWidth(2);  
	area1->SetMarkerStyle(29);
	area1->SetMarkerColor(1);	
	area1->SetMarkerSize(1.5);	
//    area1->SetTitle("IC80 effective area");
    area1->SetTitle("");
    area1->Draw("ph");

	area2->SetLineWidth(2);    
	area2->SetLineColor(1);
	area2->SetMarkerStyle(26);
	area2->SetMarkerColor(1);	
	area2->SetMarkerSize(1.5);
//	area2->Draw("phsame");
	
	area1ex->SetLineWidth(2);    
	area1ex->SetLineColor(3);
//	area1ex->Draw("same");

	area3->SetLineWidth(2);	
	area3->SetLineColor(1);
	area3->SetMarkerStyle(27);
	area3->SetMarkerColor(1);
	area3->SetMarkerSize(1.5);
//	area3->Draw("phsame");	

	area1->GetXaxis()->SetTitle(" log_{10}(E_{#mu}/TeV)");	
	area1->GetXaxis()->CenterTitle(1);
 	area1->GetXaxis()->SetTitleOffset(1.1);
	area1->GetYaxis()->SetTitle("A_{eff}(km^{2})");	
 	area1->GetYaxis()->CenterTitle(1);
 	area1->GetYaxis()->SetTitleOffset(1.1);

	TLegend *leg5 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
	leg5->SetFillStyle(0);
	leg5->SetBorderSize(0);
//	leg5->SetShadowColor("1");
 	leg5->AddEntry(area1, "0^{o}", "p");
 	leg5->AddEntry(area3, "30^{o}", "p");
 	leg5->AddEntry(area2, "45^{o}", "p");
	leg5->Draw();


	app->Run();

	return 0;
}