const int NTRIALS	=   1000;
const int NEV		= 1000;
const int NSG		= 1000;
void sgtrial();
void sgerrcalc(float*);

void sgerr(){

	gRandom->SetSeed(0);

	gStyle->SetPadRightMargin(0.01);
	gStyle->SetPadTopMargin(0.002);
	gStyle->SetPadBottomMargin(0.085);
	gStyle->SetPadLeftMargin(0.07);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.12);
	gStyle->SetStatY(0.92);
	gStyle->SetStatX(0.98);
	gStyle->SetLabelSize(0.05,"xyzt");
	//gStyle->SetLabelOffset(-0.02,"X");
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(11);

	sgtrial();

	int kdraw=0;
	if (kdraw){
		TGraphErrors* gresrat1	= new TGraphErrors();
			gresrat1->SetMarkerStyle(20);
			gresrat1->SetTitle("SGerror/TrueError vs NSG");
			int ipnt=-1;
			gresrat1->SetPoint(++ipnt,0.99*2.,0.813);  gresrat1->SetPointError(ipnt,0,0.580);
			gresrat1->SetPoint(++ipnt,1.00*2.,0.765);  gresrat1->SetPointError(ipnt,0,0.570);
			gresrat1->SetPoint(++ipnt,1.01*2.,0.828);  gresrat1->SetPointError(ipnt,0,0.594);
			gresrat1->SetPoint(++ipnt,0.99*3.,0.861);  gresrat1->SetPointError(ipnt,0,0.440);
			gresrat1->SetPoint(++ipnt,1.00*3.,0.879);  gresrat1->SetPointError(ipnt,0,0.456);
			gresrat1->SetPoint(++ipnt,1.01*3.,0.880);  gresrat1->SetPointError(ipnt,0,0.430);
			gresrat1->SetPoint(++ipnt,0.99*4.,0.930);  gresrat1->SetPointError(ipnt,0,0.409);
			gresrat1->SetPoint(++ipnt,1.00*4.,0.930);  gresrat1->SetPointError(ipnt,0,0.389);
			gresrat1->SetPoint(++ipnt,1.01*4.,0.927);  gresrat1->SetPointError(ipnt,0,0.352);
			gresrat1->SetPoint(++ipnt,0.99*5.,0.926);  gresrat1->SetPointError(ipnt,0,0.33);
			gresrat1->SetPoint(++ipnt,1.00*5.,0.930);  gresrat1->SetPointError(ipnt,0,0.34);
			gresrat1->SetPoint(++ipnt,1.01*5.,0.932);  gresrat1->SetPointError(ipnt,0,0.33);
			gresrat1->SetPoint(++ipnt,1.00*10.,0.976); gresrat1->SetPointError(ipnt,0,0.24);
			gresrat1->SetPoint(++ipnt,1.00*25.,0.990); gresrat1->SetPointError(ipnt,0,0.149);
			gresrat1->SetPoint(++ipnt,1.00*50.,0.997); gresrat1->SetPointError(ipnt,0,0.097);
//			TF1 *fresrat1	= new TF1("fresrat1","1.-exp([0]+[1]*x)",1.5,55.);
			TF1 *fresrat1	= new TF1("fresrat1","1.-([0]/x)",1.5,55.);
				fresrat1->SetLineWidth(2);
				fresrat1->SetLineColor(16);
//			TLatex *eq = new TLatex(16,1.45,"1.-exp([p0]+[p1]*x)");
			TLatex *eq = new TLatex(16,1.45,"1. - [p0]/#surdNSG");
			
		TGraphErrors* gresrat2	= new TGraphErrors();
			gresrat2->SetMarkerStyle(20);
			gresrat2->SetMarkerColor(6);
			gresrat2->SetLineColor(6);
			gresrat2->SetTitle("SGerror/TrueError vs NSG");
			ipnt=-1;
			gresrat2->SetPoint(++ipnt,0.99*2. ,0.789);  gresrat2->SetPointError(ipnt,0,0.602);
			gresrat2->SetPoint(++ipnt,1.00*2. ,0.818);  gresrat2->SetPointError(ipnt,0,0.63);
			gresrat2->SetPoint(++ipnt,1.01*2. ,0.807);  gresrat2->SetPointError(ipnt,0,0.62);
			gresrat2->SetPoint(++ipnt,1.02*2. ,0.826);  gresrat2->SetPointError(ipnt,0,0.64);
			gresrat2->SetPoint(++ipnt,0.99*5. ,0.946);  gresrat2->SetPointError(ipnt,0,0.34);
			gresrat2->SetPoint(++ipnt,1.01*5. ,0.950);  gresrat2->SetPointError(ipnt,0,0.36);
			gresrat2->SetPoint(++ipnt,0.99*10.,0.971);  gresrat2->SetPointError(ipnt,0,0.235);
			gresrat2->SetPoint(++ipnt,1.00*10.,0.977);  gresrat2->SetPointError(ipnt,0,0.223);
			gresrat2->SetPoint(++ipnt,1.01*10.,0.960);  gresrat2->SetPointError(ipnt,0,0.24);
			gresrat2->SetPoint(++ipnt,0.99*25.,0.985);  gresrat2->SetPointError(ipnt,0,0.16);
			gresrat2->SetPoint(++ipnt,1.01*25.,0.996);  gresrat2->SetPointError(ipnt,0,0.15);
			gresrat2->SetPoint(++ipnt,0.99*50.,0.995);  gresrat2->SetPointError(ipnt,0,0.107);
			gresrat2->SetPoint(++ipnt,1.01*50.,0.989);  gresrat2->SetPointError(ipnt,0,0.0963);

		TGraphErrors* gresrat3	= new TGraphErrors();
			gresrat3->SetMarkerStyle(20);
			gresrat3->SetMarkerColor(2);
			gresrat3->SetLineColor(2);
			gresrat3->SetTitle("SGerror/TrueError vs NSG");
			ipnt=-1;
			gresrat3->SetPoint(++ipnt,0.99*3. ,0.869);  gresrat3->SetPointError(ipnt,0,0.44);
			gresrat3->SetPoint(++ipnt,1.00*3. ,0.920);  gresrat3->SetPointError(ipnt,0,0.45);
			gresrat3->SetPoint(++ipnt,1.01*3. ,0.870);  gresrat3->SetPointError(ipnt,0,0.44);
			gresrat3->SetPoint(++ipnt,0.99*4. ,0.937);  gresrat3->SetPointError(ipnt,0,0.39);
			gresrat3->SetPoint(++ipnt,1.00*4. ,0.901);  gresrat3->SetPointError(ipnt,0,0.38);
			gresrat3->SetPoint(++ipnt,1.01*4. ,0.936);  gresrat3->SetPointError(ipnt,0,0.38);
			gresrat3->SetPoint(++ipnt,0.99*25.,0.995);  gresrat3->SetPointError(ipnt,0,0.135);
			gresrat3->SetPoint(++ipnt,1.01*25.,0.972);  gresrat3->SetPointError(ipnt,0,0.135);

		TLegend* leg	= new TLegend(0.7,0.24,0.98,0.38);
			leg->AddEntry(gresrat1,"N_{ev}=2,500,000","P");

		TCanvas *canvasb = new TCanvas("canvasb","canvasb",680,30,600,(8.5/11.)*600);
		canvasb->cd(); canvasb->Divide(1,1,0.0001,0.0001);
			canvasb->cd(1);
				TH1F *frame2 = (TH1F*)gPad->DrawFrame(1.4,0.05,55,1.55);	
					//frame2->GetXaxis()->SetTitle("NSG");
					//frame2->GetYaxis()->SetTitle("");
					frame2->GetXaxis()->SetTitleSize(0.06);
					frame2->GetYaxis()->SetTitleSize(0.06);
					frame2->GetXaxis()->SetTitleOffset(0.65);
					frame2->GetYaxis()->SetTitleOffset(1.0);
				gPad->SetLogx(1);
				gresrat1->Fit("fresrat1","WQR");
				TLine* line = new TLine(1.4,1.0,55,1.0);
						line->SetLineStyle(2);
						line->SetLineWidth(2);
						line->SetLineColor(17);
						line->Draw("same");
				eq->Draw("same");
				//fresrat1->Draw("same");
				leg->Draw("same");
				gresrat1->Draw("P");
		canvasb->cd(); canvasb->Update();
		canvasb->Print("sgerr1.pdf");
		//delete canvasb;

		leg->AddEntry(gresrat2,"N_{ev}=50,000","P");
		leg->AddEntry(gresrat3,"N_{ev}=500","P");

		TCanvas *canvasc = new TCanvas("canvasc","canvasc",710,60,600,(8.5/11.)*600);
		canvasc->cd(); canvasc->Divide(1,1,0.0001,0.0001);
			canvasc->cd(1);
				TH1F *frame3 = (TH1F*)gPad->DrawFrame(1.4,0.05,55,1.55);	
					//frame3->GetXaxis()->SetTitle("NSG");
					//frame3->GetYaxis()->SetTitle("");
					frame3->GetXaxis()->SetTitleSize(0.06);
					frame3->GetYaxis()->SetTitleSize(0.06);
					frame3->GetXaxis()->SetTitleOffset(0.65);
					frame3->GetYaxis()->SetTitleOffset(1.0);
				gPad->SetLogx(1);
				line->Draw("same");
				fresrat1->Draw("same");
				gresrat1->Draw("P");
				gresrat2->Draw("P");
				gresrat3->Draw("P");
				eq->Draw("same");
				leg->Draw("same");
		canvasc->cd(); canvasc->Update();
		canvasc->Print("sgerr2.pdf");
		//delete canvasc;
	}

}


void sgtrial(){
	
	TH1D* herrrat	= new TH1D("herrrat",Form("U_{SG}/U_{true}, %d trials",NTRIALS),500,0.0,5.0);
	float res[2]	= {0};
	for (int i=0;i<NTRIALS;i++){
		sgerrcalc(res);
		cout<<i
			<<"  true uncert = "<<res[0]
			<<"   SG estimate = "<<res[1]
			<<"   ratio = "<<res[1]/res[0]
			<<endl;
		float rat	= res[1]/res[0];
		herrrat->Fill(rat);
	}
		
	gStyle->SetOptStat(111111);
	TCanvas *canvas = new TCanvas("canvas","canvas",30,30,600,(8.5/11.)*600);
	canvas->cd(); canvas->Divide(1,1,0.0001,0.0001);
		canvas->cd(1);
			herrrat->Draw();
	canvas->cd(); canvas->Update();
	//canvas->Print("sgerr.pdf");
	//delete canvas;
				
}

//-----------------------------------
void sgerrcalc(float* result){

	double pars[2]	= {0};
	TF1 *fparent	= new TF1("fparent","gaus",-10,10);
		fparent->SetParameter(0,1);
		fparent->SetParameter(1,2);
		fparent->SetParameter(2,1);
	TH1D* hparent	= new TH1D("hparent","hparent",200,-10,10);
		
	double true_S2		=  0;	// used to calculate the true uncertainty
	double true_S		=  0;
	double true_N		=  0;
	float sg_N[NSG]		= {0};	// used to calculate the SG std dev
	float sg_S[NSG]		= {0};
	float sg_S2[NSG]	= {0};
	for (int iev=0;iev<NEV;iev++){
		float val	= fparent->GetRandom();
		hparent->Fill(val);
		//
		true_N		+= 1.;
		true_S		+= val;
		true_S2		+= val*val;
		//
		int isg		 = iev%NSG;
		sg_N[isg]	+=	1.;	
		sg_S[isg]	+=	val;	
		sg_S2[isg]	+=	val*val;	
	}

	//---- calculate population mean value, and true uncertainty of this mean...
	float avg	= true_S/true_N;
	float var	= true_S2/true_N - avg*avg;
	float rms	= sqrt(var);
	float avge	= rms/sqrt(NEV);
	result[0]	= avge;
	//cout<<"MEAN = "<<hparent->GetMean()<<"   avg = "<<avg<<endl;
	//cout<<"RMS  = "<<hparent->GetRMS() <<"   rms = "<<rms<<endl;
	//cout<<"avge	= "<<avge<<endl;

	//---- calculate mean value in each SG...
	//---- increment variables needed to calculate the std dev of these SG mean values...
	float sgn=0;
	float sgs=0;
	float sgs2=0;
	for (int isg=0;isg<NSG;isg++){
		float aavg		= sg_S[isg]/sg_N[isg];						// mean value in this SG
		//float avar		= sg_S2[isg]/sg_N[isg] - aavg*aavg;		// not needed...
		//float arms		= sqrt(avar);							// not needed...
		//float aavge		= arms/sqrt(sg_N[isg]);					// not needed...
		//
		sgn		+=	1.;			// increment variables needed to calculate the std dev over the SGs
		sgs		+=	aavg;	
		sgs2	+=	aavg*aavg;	
		//
	}
	if (sgn!=NSG){ cout<<sgn<<" "<<NSG<<endl; exit(0); }
	float favg	= sgs/sgn;
	float fvar	= sgs2/sgn - favg*favg;
	float frms	= sqrt(fvar);
	float favge	= frms/sqrt(sgn-1.);		// SG estimate of the uncertainty of the population mean...
	//cout<<"SG avge	= "<<favge<<endl;
	result[1]	= favge;

	delete hparent; hparent=0;
	delete fparent; fparent=0;
}