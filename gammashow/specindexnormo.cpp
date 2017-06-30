{
/*
Normalization for 5 sigma after 10 years

*/

gROOT->Reset();

const int n=12;

const double binfact = 2.222; //1.6 *radius gets 0.72 of signal
//const double binfact = 1.;


double indices[12] = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6};

double knormo[12] = {10.28, 4.82, 2.98, 2.15, 1.71, 1.45, 1.26, 1.13, 1.03, 0.96, 0.89, 0.84};
double normo[12] = {9.66, 4.44, 2.71, 1.94, 1.54, 1.30, 1.13, 1.01, 0.93, 0.86, 0.80, 0.76 };



double knormo5[n], normo5[n],knormo3[n],normo3[n];
for(int i=0; i<12; i++){

	knormo5[i] = (5./knormo[i])*1.E-11*binfact;
	normo5[i] = (5./normo[i])*1.E-11*binfact;
	
	knormo3[i] = (3./knormo[i])*1.E-11*binfact;
	normo3[i] = (3./normo[i])*1.E-11*binfact;
	
}

	

	TGraph *signif5 = new TGraph(n, indices, normo5);
	signif5->SetLineColor(1);
	signif5->SetLineWidth(2);
	signif5->SetLineStyle(1);	

	TGraph *ksignif5 = new TGraph(n, indices, knormo5);
	ksignif5->SetLineColor(1);
	ksignif5->SetLineWidth(2);
	ksignif5->SetLineStyle(2);

	TGraph *signif3 = new TGraph(n, indices, normo3);
	signif3->SetLineColor(1);
	signif3->SetLineWidth(2);
	signif3->SetLineStyle(1);	

	TGraph *ksignif3 = new TGraph(n, indices, knormo3);
	ksignif3->SetLineColor(1);
	ksignif3->SetLineWidth(2);
	ksignif3->SetLineStyle(2);


	TCanvas *c3 = new TCanvas("events3","stuff3",500,500 );
	TPad *pad3 = new TPad("pad3","pad3",0.0,0.0,1.0,1.0,10);

	c3->cd();

	pad3->SetTitle("");
	pad3->Draw();
	pad3->cd();
//	pad3->SetLogx();
	pad3->SetLogy();
	pad3->SetGridx();
	pad3->SetGridy();
	pad3->SetTicks();

	signif5->SetMinimum(1E-12);
	signif5->SetMaximum(1E-9);
	signif5->SetTitle("");
	
	signif5->Draw("Al");
	ksignif5->Draw("same");
	
	signif3->Draw("same");	
	ksignif3->Draw("same");
	
	signif5->GetXaxis()->SetLimits(1.5,2.6);
	signif5->GetXaxis()->SetTitle("Spectral Index");	
	signif5->GetXaxis()->CenterTitle(1);
 	signif5->GetXaxis()->SetTitleOffset(1.2);
	signif5->GetYaxis()->SetTitle("#phi_{#gamma} at 1 TeV (TeV^{-1} cm^{-2} s^{-1})");	
 	signif5->GetYaxis()->CenterTitle(1);
 	signif5->GetYaxis()->SetTitleOffset(1.5);
	
	
	TLegend *leg3 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
	leg3->SetFillStyle(0);
	leg3->SetBorderSize(0);
//	leg3->SetShadowColor("1");
 	leg3->AddEntry(signif5, "No kaons", "l");
 	leg3->AddEntry(ksignif5, "Kaons", "l");
	leg3->Draw();

	TLegend *leg1 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
	leg1->SetFillStyle(0);
	leg1->SetBorderSize(0);
//	leg1->SetShadowColor("1");
 	leg1->AddEntry(signif5, "5#sigma", "");
	leg1->Draw();

	TLegend *leg2 = new TLegend(0.25, 0.7, 0.90, 0.89, "", "brNDC");
	leg2->SetFillStyle(0);
	leg2->SetBorderSize(0);
//	leg2->SetShadowColor("1");
 	leg2->AddEntry(signif5, "3#sigma", "");
	leg2->Draw();


}