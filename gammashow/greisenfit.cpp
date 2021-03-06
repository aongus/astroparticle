{
//find analytical fits to Greisen & Rossi tabulated parameters

	gROOT->Reset();

 double ess[38] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 
 					1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 
 					2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 
 					9.0, 10.0};
 
 double lamb1[38] = {1.E5, 3.789, 2.270, 1.569, 1.127, 0.813, 0.576, 0.389, 0.235, 0.108, 0.0, -0.092, -0.171, -0.239, -0.298, -0.350, -0.395, -0.435, -0.470, -0.5, -0.526, -0.550, -0.570, -0.589, -0.605, -0.619, -0.632, -0.643, -0.654, -0.663, -0.671, -0.720, -0.742, -0.752, -0.759, -0.763, -0.765, -0.766};

 double lamb1prim[38] = {-1.E5, -25.005, -9.488, -5.415, -3.654, -2.693, -2.093, -1.685, -1.389, -1.1660, -0.9908, -0.8501, -0.7333, -0.6362, -0.5531, -0.4825, -0.4214, -0.3691, -0.3238, -0.2841, -0.2498, -0.2202, -0.1943, -0.1719, -0.1523, -0.1354, -0.1205, -0.1077, -0.0964, -0.0863, -0.0777, -0.0307, -0.0146, -0.0080, -0.0048, -0.0031, -0.0021, -0.0015};
 
 double lamb2[38] = {-1.E5, -4.715, -3.330, -2.749, -2.415, -2.201, -2.055, -1.953, -1.878, -1.824, -1.7868, -1.760, -1.744, -1.734, -1.732, -1.734, -1.741, -1.751, -1.762, -1.780, -1.797, -1.816, -1.837, -1.859, -1.882, -1.902, -1.928, -1.951, -1.977, -2.003, -2.026, -2.264, -2.480, -2.669, -2.837, -2.988, -3.123, -3.246};

 double ess2[31] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 
 					1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 
 					2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0};
 double lamb1doubprim[31] = {1.E5, 1.E3, 75, 26, 12.5, 7.6, 4.95, 3.50, 2.55, 1.97, 1.5634, 1.275, 1.060, 0.893, 0.764, 0.655, 0.565, 0.487, 0.423, 0.370, 0.320, 0.277, 0.241, 0.210, 0.182, 0.159, 0.138, 0.120, 0.107, 0.093, 0.080};


  	delete (TF1*) gROOT->FindObject("curve");
   	
   	TF1 *f1 = new TF1("curve", "[0] + [1]*pow(x,[2])+[3]*x+[4]*pow(x,2)", 0.05,10);

	f1->SetParLimits(0, -1., 0.);
	f1->SetParLimits(1, -10., 10.);
	f1->SetParLimits(2, -1., -0.1.);
	f1->SetParLimits(3, -10., 10.);
	f1->SetParLimits(4, -10., -10.);


//   	delete (TF1*) gROOT->FindObject("hifit");	
//   	TF1 *hi = new TF1("hifit", "[0]+ [1]*pow(x, [2]*log(x) + [3])  ", 7., 1E3);
//   	double hip[4] = {-3.32101, 6.48719E-1, -8.98896E-2, 1.46431};
//   	hifit->SetParameters(hip[0],hip[1],hip[2],hip[3]);


	TCanvas *c2 = new TCanvas("events2","stuff2",500,500 );
	pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,1.0,10);


	pad2->SetTitle("");
	pad2->Draw();
	pad2->cd();
//	pad2->SetLogx();
//	pad2->SetLogy();

	TGraph *gr1 = new TGraph(38, ess, lamb1);
	gr1->SetMaximum(4);
	gr1->SetMinimum(-1);
	gr1->Draw("A*");

	double a = gr1->Fit("pol6","q","l",0.05,10);
	cout <<a <<endl;

//	->Draw("same");

}