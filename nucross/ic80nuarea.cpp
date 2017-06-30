/*
IC80 area interpolation and derivative
*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
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

int main (int argc, char **argv){
 	TStopwatch time;
 	time. Start();

	const int n = 128;

    double xi, yi;

//  	 const double energy[n] = {0.158489,0.251189,0.398107,
//  							0.630957,1,1.58489,2.51189,3.98107,6.30957,10,15.8489,
//  							25.1189,39.8107,63.0957,100,158.489,251.189,398.107,
//  							630.957,1000.}; //energy (TeV) corresponding to areas
// 							
// 	 const double ic80[n] = {6.96E-4, 3.69E-3, 1.72E-2, 6.01E-2, 1.57E-1, 
// 							3.86E-1, 1.01, 2.26, 4.43, 8.14, 1.50E1, 2.50E1, 3.92E1, 
// 							5.95E1,7.69E1, 1.03E2, 1.33E2, 1.61E2, 1.95E2, 2.15E2}; //Taup average, m2




	const double energy[n] = {0.13342, 0.13695, 0.14180, 0.14683, 0.15204, 0.15881, 0.16301, 0.17027, 0.18097, 0.18739, 0.19744, 0.20632, 0.21541, 0.22305, 0.23501, 0.24761, 0.25863, 0.27013, 0.28710, 0.30514, 0.31871, 0.33579, 0.35379, 0.37602, 0.40313, 0.42475, 0.44751, 0.47150, 0.49678, 0.52798, 0.56115, 0.59640, 0.63940, 0.67367, 0.71599, 0.76095, 0.83011, 0.88995, 0.95408, 1.04080, 1.12556, 1.21722, 1.31635, 1.41122, 1.52615, 1.63614, 1.75405, 1.88046, 1.98122, 2.12399, 2.25735, 2.42002, 2.59441, 2.85489, 3.06059, 3.25272, 3.45691, 3.83718, 4.14961, 4.48742, 4.98097, 5.38645, 5.77450, 6.57905, 7.11461, 7.82875, 8.68977, 9.39708, 10.43063, 11.57772, 12.96327, 15.15930, 16.39306, 18.83955, 20.91118, 24.03203, 27.14243, 30.38995, 34.62293, 38.42996, 43.40369, 49.88034, 58.83804, 68.20815, 80.45621, 94.90279, 113.90621, 132.04379, 158.48534, 180.55542, 212.97652, 242.63443, 278.83148, 320.42776, 368.23001, 415.87106, 473.77826, 539.74731, 642.21436, 770.80145, 933.20081, 1149.62305, 1321.11804, 1571.89160, 1854.07971, 2168.00269, 2448.48755, 2913.25000, 3557.81787, 4160.22559, 4780.83154, 5541.94971, 6313.54541, 7192.58008, 8337.62012, 9920.16895, 11803.1484, 13563.8056, 15860.3701, 18385.4199, 21875.2343, 25802.4121, 31511.2168, 36211.9023, 42712.8203, 49512.5859, 59943.4570, 71321.5156};

	const double ic80[n] = {0.00035, 0.00039, 0.00045, 0.00052, 0.00060, 0.00071, 0.00081, 0.00098, 0.00125, 0.00145, 0.00175, 0.00216, 0.00253, 0.00297, 0.00360, 0.00416, 0.00473, 0.00538, 0.00632, 0.00753, 0.00843, 0.00989, 0.01143, 0.01363, 0.01626, 0.01909, 0.02170, 0.02507, 0.02897, 0.03400, 0.04121, 0.04838, 0.05771, 0.06667, 0.07953, 0.09187, 0.10959, 0.12865, 0.14392, 0.17168, 0.20154, 0.23282, 0.27331, 0.31573, 0.37063, 0.42816, 0.48675, 0.56231, 0.62909, 0.71518, 0.82619, 0.93925, 1.06778, 1.27372, 1.42500, 1.59424, 1.78359, 2.12759, 2.45783, 2.74975, 3.17660, 3.55389, 3.91276, 4.74287, 5.30618, 6.12984, 7.08137, 7.79646, 9.15221, 10.40480, 12.21417, 14.57015, 15.78637, 18.23700, 20.40325, 23.95153, 27.22963, 30.46412, 35.19307, 39.37338, 44.76214, 50.07951, 57.85414, 65.77285, 73.58710, 79.73154, 87.78601, 95.11545, 108.13538, 117.16359, 131.08244, 142.02650, 149.03011, 153.89182, 158.91248, 164.09572, 172.18651, 180.67661, 195.76331, 212.11066, 219.03432, 229.84084, 233.56436, 237.35217, 237.36505, 237.37743, 241.22255, 241.23663, 253.13591, 257.23950, 265.63129, 269.93576, 274.30853, 278.75333, 278.76700, 269.98785, 265.71140, 265.72247, 265.73605, 274.40585, 274.42160, 283.37598, 297.35443, 312.01593, 322.19696, 322.21176, 312.06769, 312.08527};

	const int m = 200;
	
	const double dej = log10(energy[n-1]/energy[0])/(m-1.);
    double xia[m], yia[m], yiaprim[m];
     
     
   	{
  		gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
     
        gsl_spline_init (spline, energy, ic80, n);
 
        for(int j=0; j<m; j++){
        	xi = pow(10, -0.9 + dej*j);
        	yi = gsl_spline_eval (spline, xi, acc);
			xia[j] = xi;
			yia[j] = yi;
			yiaprim[j] = gsl_spline_eval_deriv (spline, xi, acc);
//			cout<<yiaprim[j]<<endl;
        	}
         
      	gsl_spline_free (spline);
        gsl_interp_accel_free (acc);
	}

 	time.Stop();
 	
 	cout <<time.RealTime() <<" seconds"<<endl;

	cout.precision(8);

 	TApplication *app = new TApplication("app", &argc, argv); //needed to run ROOT externally
 
 	TGraph *data = new TGraph(n, energy, ic80);
	data->SetLineColor(1);
	data->SetLineWidth(2);
	data->SetLineStyle(1);	

	TGraph *smooth = new TGraph(m, xia, yia);
	smooth->SetLineColor(2);
	smooth->SetLineWidth(2);
	smooth->SetLineStyle(1);	
 
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

	data->SetMinimum(1E-4);
	data->SetMaximum(1E3);

	data->SetTitle("");
	
	data->Draw("A*");
	smooth->Draw("same");
	
	data->GetXaxis()->SetLimits(1E-1,1E3);
//	data->GetXaxis()->SetTitle("E_{#mu} in-ice (TeV)");	
	data->GetXaxis()->CenterTitle(1);
 	data->GetXaxis()->SetTitleOffset(1.2);
//	data->GetYaxis()->SetTitle("Events/year");	
 	data->GetYaxis()->CenterTitle(1);
 	data->GetYaxis()->SetTitleOffset(1.2);
 








//  	TGraph *data2 = new TGraph(n, energy, lamb1prim);
// 	data2->SetLineColor(1);
// 	data2->SetLineWidth(2);
// 	data2->SetLineStyle(1);	

	TGraph *smoothp = new TGraph(m, xia, yiaprim);
	smoothp->SetLineColor(2);
	smoothp->SetLineWidth(2);
	smoothp->SetLineStyle(1);	

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

//	data2->SetMinimum(1E0);
//	data2->SetMaximum(1E7);

	smoothp->SetTitle("");
	
	smoothp->Draw("Al");
//	smoothp->Draw("same");
	
// 	data2->GetXaxis()->SetLimits(0,10);
// //	data2->GetXaxis()->SetTitle("E_{#mu} in-ice (TeV)");	
// 	data2->GetXaxis()->CenterTitle(1);
//  	data2->GetXaxis()->SetTitleOffset(1.2);
// //	data2->GetYaxis()->SetTitle("Events/year");	
//  	data2->GetYaxis()->CenterTitle(1);
//  	data2->GetYaxis()->SetTitleOffset(1.2);

 
 
 	app->Run();




  	return 0;
    }
