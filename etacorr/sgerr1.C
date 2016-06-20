//
const int NTRIALS	=   1000;
const int NEV		= 100000;
const int NSG		=      5;
//
void sgtrial();
void sgerrcalc(float*);
TF1 *fparent;

//----------------------------------------------------------------
void sgerr(){
	//
	gRandom->SetSeed(0);
	//
	gStyle->SetPadRightMargin(0.01);
	gStyle->SetPadTopMargin(0.002);
	gStyle->SetPadBottomMargin(0.085);
	gStyle->SetPadLeftMargin(0.07);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.12);
	gStyle->SetStatY(0.98);
	gStyle->SetStatX(0.98);
	gStyle->SetLabelSize(0.05,"xyzt");
	//gStyle->SetLabelOffset(-0.02,"X");
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(11);
	//
	fparent	= new TF1("fparent","gaus",-10,10);
		fparent->SetParameter(0,1);		// parent norm
		fparent->SetParameter(1,0);		// parent mean
		fparent->SetParameter(2,1);		// parent std dev
	//
	sgtrial();
	//
}


//----------------------------------------------------------------
void sgtrial(){
	//	
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
	//	
	TLatex *txt1 = new TLatex(); txt1->SetNDC(); txt1->SetTextSize(0.04);
	gStyle->SetOptStat(111111);
	TCanvas *canvas = new TCanvas("canvas","canvas",30,30,600,(8.5/11.)*600);
	canvas->cd(); canvas->Divide(1,1,0.0001,0.0001);
		canvas->cd(1);
			herrrat->Draw();
			txt1->DrawLatex(0.4,0.88,Form("NEV=%d, NSG=%d",NEV,NSG));
	canvas->cd(); canvas->Update();
	//canvas->Print("sgerr.pdf");
	//delete canvas;
	//			
}

//----------------------------------------------------------------
void sgerrcalc(float* result){
	//
	double pars[2]	= {0};
	//TH1D* hparent	= new TH1D("hparent","hparent",200,-10,10);
	//		
	double true_S2		=  0;	// used to calculate the true uncertainty
	double true_S		=  0;
	double true_N		=  0;
	float sg_N[NSG]		= {0};	// used to calculate the SG std dev
	float sg_S[NSG]		= {0};
	float sg_S2[NSG]	= {0};
	for (int iev=0;iev<NEV;iev++){
		float val	= fparent->GetRandom();
		//hparent->Fill(val);
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
	//
	//---- calculate population mean value, and true uncertainty of this mean...
	float avg	= true_S/true_N;
	float var	= true_S2/true_N - avg*avg;
	float rms	= sqrt(var);
	float avge	= rms/sqrt(NEV);
	result[0]	= avge;
	//cout<<"MEAN = "<<hparent->GetMean()<<"   avg = "<<avg<<endl;
	//cout<<"RMS  = "<<hparent->GetRMS() <<"   rms = "<<rms<<endl;
	//cout<<"avge	= "<<avge<<endl;
	//
	//---- calculate mean value in each SG...
	//---- increment variables needed to calculate the std dev of these SG mean values...
	float sgn	= 0;
	float sgs	= 0;
	float sgs2	= 0;
	for (int isg=0;isg<NSG;isg++){
		float aavg		= sg_S[isg]/sg_N[isg];						// mean value in this SG
		//float avar		= sg_S2[isg]/sg_N[isg] - aavg*aavg;		// not needed...
		//float arms		= sqrt(avar);							// not needed...
		//float aavge		= arms/sqrt(sg_N[isg]);					// not needed...
		//
		sgn		+=	1.;			// increment variables needed for the std dev of the SG means
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
	//
	//delete hparent; hparent=0;
}