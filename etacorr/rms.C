Double_t powerlaw(Double_t *x, Double_t *par){
  Float_t xx =x[0];
  Double_t f = par[0]*TMath::Power(xx,par[1]);
  return f;
}

Double_t Get_a(float val){
	double loga = TMath::Log10(val) + 0.5*TMath::Log10(0.1);
	double a = TMath::Power(10.,loga);
	return a;
}

void rms(){

	float x04[5] = {0.1,0.32,1.0,3.2,10.0};
	float r04[5] = {6.406,3.599,2.006,1.122,0.626};
	float x09[4] = {0.1,0.32,1.0,3.2};
	float r09[4] = {3.060,1.720,0.953,0.5394};
	float x12[1] = {0.1};
	float r12[1] = {2.334};
	float x16[1] = {0.1};
	float r16[1] = {1.7787};
	float x21[1] = {0.1};
	float r21[1] = {1.3706};
	float x25[1] = {0.1};
	float r25[1] = {1.1984};
	float x31[1] = {0.1};
	float r31[1] = {0.9343};
	float x36[1] = {0.1};
	float r36[1] = {0.80495};
	float x49[1] = {0.1};
	float r49[1] = {0.5952};
	float x81[1] = {0.1};
	float r81[1] = {0.3605};

// 	r36[0] = 0.85;
// 	r49[0] = 0.64;
// 	r81[0] = 0.41;
	
	TGraph *g04	= new TGraph(5,x04,r04);
		g04->SetMarkerStyle(20);
	TGraph *g09	= new TGraph(4,x09,r09);
		g09->SetMarkerStyle(20);
	TGraph *g12	= new TGraph(1,x12,r12);
		g12->SetMarkerStyle(20);
	TGraph *g16	= new TGraph(1,x16,r16);
		g16->SetMarkerStyle(20);
	TGraph *g21	= new TGraph(1,x21,r21);
		g21->SetMarkerStyle(20);
	TGraph *g25	= new TGraph(1,x25,r25);
		g25->SetMarkerStyle(20);
	TGraph *g31	= new TGraph(1,x31,r31);
		g31->SetMarkerStyle(20);
	TGraph *g36	= new TGraph(1,x36,r36);
		g36->SetMarkerStyle(20);
	TGraph *g49	= new TGraph(1,x49,r49);
		g49->SetMarkerStyle(20);
	TGraph *g81	= new TGraph(1,x81,r81);
		g81->SetMarkerStyle(20);
		
	TGraph *ga	= new TGraph();
		ga->SetMarkerStyle(21);

	double pars[2]	= {0};
	TF1 *f04	= new TF1("f04",powerlaw,0.07,12,2);
		f04->SetLineColor(2);
		f04->SetParName(0,"a");
		f04->SetParName(1,"k");
		g04->Fit("f04","0NQR");
		f04->GetParameters(pars);		
		parse = f04->GetParErrors();		
		cout<<"NTRK = 04 \t"
			<<"a = "<<pars[0]<<" +- "<<parse[0]<<" \t "
			<<"k = "<<pars[1]<<" +- "<<parse[1]<<endl;
		ga->SetPoint(0,4,pars[0]);
		//
	TF1 *f09	= new TF1("f09",powerlaw,0.07,12,2);
		f09->SetLineColor(2);
		g09->Fit("f09","0NQR");
		f09->GetParameters(pars);		
		parse = f09->GetParErrors();		
		cout<<"NTRK = 09 \t"
			<<"a = "<<pars[0]<<" +- "<<parse[0]<<" \t "
			<<"k = "<<pars[1]<<" +- "<<parse[1]<<endl;
		ga->SetPoint(1,9,pars[0]);
		//
		ga->SetPoint(2,12,Get_a(r12[0]));
		ga->SetPoint(3,16,Get_a(r16[0]));
		ga->SetPoint(4,21,Get_a(r21[0]));
		ga->SetPoint(5,25,Get_a(r25[0]));
		ga->SetPoint(6,31,Get_a(r31[0]));
		ga->SetPoint(7,36,Get_a(r36[0]));
		ga->SetPoint(8,49,Get_a(r49[0]));
		ga->SetPoint(9,81,Get_a(r81[0]));
	TF1 *fa	= new TF1("fa",powerlaw,0.3,100,2);
		fa->SetParName(0,"a0");
		fa->SetParName(1,"k0");
		fa->SetLineColor(6);
		ga->Fit("fa","0NQR");
		fa->SetRange(0.3,100.);
		fa->GetParameters(pars);		
		parse = fa->GetParErrors();		
		cout<<"NTRK DEP: \t"
			<<"a = "<<pars[0]<<" +- "<<parse[0]<<" \t "
			<<"k = "<<pars[1]<<" +- "<<parse[1]<<endl;
	

	gStyle->SetPadRightMargin(0.01);
	gStyle->SetPadTopMargin(0.02);
	gStyle->SetPadBottomMargin(0.085);
	gStyle->SetPadLeftMargin(0.14);
	gStyle->SetStatW(0.3);
	gStyle->SetStatH(0.12);
	gStyle->SetStatY(0.98);
	gStyle->SetStatX(0.97);
	gStyle->SetLabelSize(0.05,"xyzt");
	gStyle->SetLabelOffset(-0.02,"X");
	gStyle->SetOptFit(11);
	
	TLatex *eq1 = new TLatex();
	TLatex *eq2 = new TLatex();
	
	TCanvas *canvas = new TCanvas("canvas","canvas",30,30,800,(8.5/11.)*800);
	canvas->cd(); canvas->Divide(2,1,0.0001,0.0001);
		canvas->cd(1);
			gPad->SetLogx(1);
			gPad->SetLogy(1);
			TH1F *frame = (TH1F*)gPad->DrawFrame(0.07,0.2,12.0,10.0);	
				frame->GetXaxis()->SetTitle("N_{evt} (Millions)");
				frame->GetYaxis()->SetTitle("#sigma(#LTR_{2}#GT)");
				frame->GetXaxis()->SetTitleSize(0.06);
				frame->GetYaxis()->SetTitleSize(0.06);
				frame->GetXaxis()->SetTitleOffset(0.65);
				frame->GetYaxis()->SetTitleOffset(1.0);
			g04->Draw("LP");
			g04->Fit("f04","QR");
			f04->Draw("same");
			g09->Draw("LP");
			f09->Draw("same");
			g12->Draw("LP");
			g16->Draw("LP");
			g21->Draw("LP");
			g25->Draw("LP");
			g36->Draw("LP");
			g49->Draw("LP");
			g81->Draw("LP");
			eq1->DrawLatex(1.0,3.5,"#sigma = [a]*N_{evt}^{[k]}");
		canvas->cd(2);
			gPad->SetLogx(1);
			gPad->SetLogy(1);
			TH1F *frame2 = (TH1F*)gPad->DrawFrame(3,0.04,110,4.0);	
				frame2->GetXaxis()->SetTitle("N_{trk}");
				frame2->GetYaxis()->SetTitle("a");
				frame2->GetXaxis()->SetTitleSize(0.06);
				frame2->GetYaxis()->SetTitleSize(0.06);
				frame2->GetXaxis()->SetTitleOffset(0.65);
				frame2->GetYaxis()->SetTitleOffset(1.0);
			ga->Draw("LP");
			fa->Draw("same");
			ga->Fit("fa","QR");
			eq2->DrawLatex(25.0,1.1,"a = [a0]*N_{trk}^{[k0]}");
	canvas->cd(); canvas->Update();
	canvas->Print("rms.pdf");



}