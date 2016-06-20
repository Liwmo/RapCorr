// #include <iostream>
// #include <TH1.h>
// #include <TH2.h>
// #include <TStyle.h>
// #include <TCanvas.h>
// #include <TProfile.h>
// #include <TRandom.h>
// #include <TFile.h>
// #include "TNtuple.h"
// #include "TTree.h"
// #include "TMath.h"
// #include "TF1.h"
// #include "TLatex.h"

//
const int NTRIALS	=   2500;
const int NEV		= 100000;
const int NSG		=    100;
//
void sgtrial();
void sgerrcalc(double*);
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
	TH1D* herrrat	= new TH1D("herrrat",Form("U_{SG}^{2}/U_{true}^{2}, %d trials",NTRIALS),4500,0.0,5.0);
	double res[2]	= {0};
	for (int i=0;i<NTRIALS;i++){
		sgerrcalc(res);
		if(i%100==0){cout<<i
			<<"  true uncert = "<<res[0]
			<<"   SG estimate = "<<res[1]
			<<"   ratio = "<<res[1]/res[0]
				 <<endl;
		}
		double rat	= res[1]/res[0];
		herrrat->Fill(rat);
	}
	cout<<"std dev of variance ratio = "<<sqrt(herrrat->GetRMS())<<endl;
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
void sgerrcalc(double* result){
	//
	double pars[2]	= {0};
	//TH1D* hparent	= new TH1D("hparent","hparent",200,-10,10);
	//		
	double true_S2		=  0;	// used to calculate the true uncertainty
	double true_S		=  0;
	double true_N		=  0;
	double sg_N[NSG]		= {0};	// used to calculate the SG std dev
	double sg_S[NSG]		= {0};
	double sg_S2[NSG]	= {0};
	for (int iev=0;iev<NEV;iev++){
		double val	= fparent->GetRandom();
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
	double avg	= true_S/true_N;
	double var	= true_S2/true_N - avg*avg;
	double rms	= sqrt(var);
	//	double avge	= rms/sqrt(NEV-1.);  !============
		double avge	= rms/sqrt(true_N-1.);  //============
		avge=1./sqrt(true_N); // comment this out if use an estimate (Might lead to some bias in the ratio!)
	
		//	result[0]	= avge;
		result[0]	= avge*avge;
	//cout<<"MEAN = "<<hparent->GetMean()<<"   avg = "<<avg<<endl;
	//cout<<"RMS  = "<<hparent->GetRMS() <<"   rms = "<<rms<<endl;
	//cout<<"avge	= "<<avge<<endl;
	//
	//---- calculate mean value in each SG...
	//---- increment variables needed to calculate the std dev of these SG mean values...
	double sgn	= 0;
	double sgs	= 0;
	double sgs2	= 0;
	for (int isg=0;isg<NSG;isg++){
		double aavg		= sg_S[isg]/sg_N[isg];						// mean value in this SG
		//double avar		= sg_S2[isg]/sg_N[isg] - aavg*aavg;		// not needed...
		//double arms		= sqrt(avar);							// not needed...
		//double aavge		= arms/sqrt(sg_N[isg]);					// not needed...
		//
		sgn		+=	1.;			// increment variables needed for the std dev of the SG means
		sgs		+=	aavg;	
		sgs2	+=	aavg*aavg;	
		//
	}
	if (sgn!=NSG){ cout<<sgn<<" "<<NSG<<endl; exit(0); }
	double favg	= sgs/sgn;
	double fvar	= sgs2/sgn - favg*favg;
	double frms	= sqrt(fvar);
	double favge	= frms/sqrt(sgn-1.);		// SG estimate of the uncertainty of the population mean...
	//cout<<"SG avge	= "<<favge<<endl;
	//result[1]	= favge;
	result[1]	= favge*favge;
	
	//
	//delete hparent; hparent=0;
}
