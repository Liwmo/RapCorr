int PLOTON	= 0;

#include <string>
#include "stdio.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TH3D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TLatex.h"
#include "TPad.h"
#include "TCutG.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TPaletteAxis.h"
#include "TColor.h"
#include "TObjectTable.h"
#include "TKey.h"
#include "TRandom.h"
#include "TApplication.h"
using namespace std;

float 	ec();
void	ecIncrement(int,float*,TH1D*,TH2D*); 
void	ecCalculate(TH1D*,TH2D*,TH2D*,TH2D*,TH1D*); 
float	GetC2Baseline(TH1D*);

void	AddToRegisterG(TGraph *g);
void	AddToRegisterH(TH1D *h);
void	AddToRegisterH2(TH2D *h);
int		nRegisteredF,nRegisteredG,nRegisteredH,nRegisteredH2;
TGraph	*RegisterG[10000];
TH1D	*RegisterH[1000];
TH2D	*RegisterH2[150];

#if !defined(ARRAY_SIZE)
    #define ARRAY_SIZE(x) (sizeof((x)) / sizeof((x)[0]))
#endif

//-----------------------------------------------------------------------------------
//
int main(int argc, char **argv){

	//---- set stack size to 63MB (default=8.192kB)
	const	rlim_t kStackSize = 63L * 1024L * 1024L;
	struct	rlimit rl;
	int		result = getrlimit(RLIMIT_STACK, &rl);
	//cout<<"Initial RLIMIT_STACK = "<<rl.rlim_cur/1024L/1024L<<" MB"<<endl;
	if (result!=0){
		fprintf(stderr, "Problem ...getrlimit returned %d\n", result);
	} else {
		if (rl.rlim_cur < kStackSize){
			rl.rlim_cur = kStackSize;
			result 		= setrlimit(RLIMIT_STACK, &rl);
			if (result != 0){
				fprintf(stderr, "Problem ...setrlimit returned %d\n", result);
			}
		}
	}
	result = getrlimit(RLIMIT_STACK, &rl);
	//cout<<"Present RLIMIT_STACK = "<<rl.rlim_cur/1024L/1024L<<" MB"<<endl;

//	TApplication theApp("App", &argc, argv);
	TCanvas *mydummycanvas=new TCanvas();

//	gRandom->SetSeed(123456);
	gRandom->SetSeed(0);

	const int NRUN	= 100;
	float avg		=   0;
	for (int irun=0;irun<NRUN;irun++){
		float rms	= ec();
		avg	+= rms/((float)NRUN);
	}	
	cout<<endl;
	cout<<"NRUN = "<<NRUN<<"    <rms> = "<<avg<<endl;
	cout<<endl;

	return 0;
}


//---------------------------------------------------------------------------------
//
float ec(){

	TH1::AddDirectory(kFALSE);
	char 		buf[200],bufb[200],bufc[200];
	int 		ican=-1,iframe=-1,ileg=-1;
	TCanvas 	*ccan[100];
	TH1F		*frame[100];
	TLegend		*legend[100];
	int			icolstart_energy,icolstart_sets;
	TString 	plotfilebase,rootoutfile;
	TString 	plotfileO,plotfile,plotfileC,plotfilePDF;

	//---- set color table...
//  	const Int_t NRGBs = 5;
//  	      Int_t NCont = NE;
// 	//Double_t stops[NRGBs]	= { 0.00, 0.34, 0.61, 0.84, 1.00 };
//  	Double_t stops[NRGBs]	= { 0.00, 0.10, 0.61, 0.99, 1.00 };
//  	Double_t red[NRGBs]		= { 0.00, 0.00, 0.87, 1.00, 0.51 };
//  	Double_t green[NRGBs]	= { 0.00, 0.81, 1.00, 0.20, 0.00 };
//  	Double_t blue[NRGBs]	= { 0.51, 1.00, 0.12, 0.00, 0.00 };
//  	icolstart_energy	= TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NE);
//  	icolstart_sets 		= TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NSET/2);

	//---- log-variable binning ...
	const int anb	= 200;
	float ax1		= 1.0E-10;
	float ax2		= 1.0;
	float logax1	= TMath::Log10(ax1);
	float logax2	= TMath::Log10(ax2);
	float abw		= (logax2-logax1)/((float)anb);
	double xb[anb+1];
	xb[0]			= ax1;
	for (int i=1;i<=anb;i++){
		xb[i]		= TMath::Power(10,logax1+i*abw);
	}

	if (PLOTON){
		rootoutfile		= TString(Form("./root/etacorr.root"));
		plotfilebase	= TString(Form("./ps/etacorr"));
		cout<<"root file    = "<<rootoutfile.Data()<<endl;
		cout<<"plot file    = "<<plotfilebase.Data()<<endl;
		plotfileO		= plotfilebase + TString(".ps(");
		plotfile		= plotfilebase + TString(".ps");
		plotfileC		= plotfilebase + TString(".ps]");
		plotfilePDF		= plotfilebase + TString(".pdf");
	}

	//---- booking
	const int   NBETA	= 25;
	const float ETAMAX	= 1.0;
	TH1D *hmultGen		= new TH1D("hmultGen","hmultGen",80,-0.5,79.5);
	TH1D *hetaGen		= new TH1D("hetaGen","Generated #eta",NBETA,-ETAMAX,ETAMAX);
	TH1D *hdeta			= new TH1D("hdeta","#Delta#eta",2*NBETA-1,-2.*ETAMAX,2.*ETAMAX);
	TH1D *heta1D		= new TH1D("heta1D","heta1D",NBETA,-ETAMAX,ETAMAX);
	TH2D *heta2D		= new TH2D("heta2D","#eta_{2} vs #eta_{1}",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH2D *hetaT			= new TH2D("hetaT","hetaT",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH2D *hR2			= new TH2D("hR2","R_{2} vs (#eta_{1},#eta_{2})",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH2D *hC2			= new TH2D("hC2","C_{2} vs (#eta_{1},#eta_{2})",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH1D *hcorr1D		= new TH1D("hcorr1D","hcorr1D",anb,xb);
	TH2D *hm1			= new TH2D("hm1","hm1",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH1D *hR2deta		= new TH1D("hR2deta","#LTR_{2}#GT vs #eta_{1}-#eta_{2}",2*NBETA-1,-2.*ETAMAX,2.*ETAMAX);
	TH1D *hR2detabase	= new TH1D("hR2detabase","#LTR_{2}#GT-Baseline vs #eta_{1}-#eta_{2}",2*NBETA-1,-2.*ETAMAX,2.*ETAMAX);
	hmultGen	->Reset();	
	hetaGen		->Reset();
	hdeta		->Reset();	
	heta1D		->Reset();
	heta2D		->Reset();
	hetaT		->Reset();	
	hR2			->Reset();
	hC2			->Reset();
	hcorr1D		->Reset();
	hm1			->Reset();
	hR2deta		->Reset();
	hR2detabase	->Reset();

	int kfreqnb		= 500;
	float freqbw	= 0.0001; 
	TH1D *hR2detabasevaldist	= new TH1D("hR2detabasevaldist","Frequency #LTR_{2}#GT-Baseline",
			kfreqnb,-kfreqnb*freqbw/2., kfreqnb*freqbw/2.);	// bw=0.0002

//	const int NEVT	= 10000000;
//	const int NEVT	=  3200000;
//	const int NEVT	=  1000000;
//	const int NEVT	=   320000;
	const int NEVT	=   100000;
	      int NTRK	=       12;
	      int NTRKP	=        0;		// number of smeared deta pairs added
	float *Aeta		= new float[NTRK];

	//---- Generate & Calculate
	for (int ievt=0;ievt<NEVT;ievt++){
		//
//		NTRK = (int)(61 + 8.1*gRandom->Gaus());
		if (hmultGen){
			hmultGen->Fill(NTRK);
		} else {
			cout<<hmultGen<<endl;
		}
		//
		for (int itrk=0;itrk<NTRK-2*NTRKP;itrk++){		// random sampled...
			float eta = -1.0 + 2.0*gRandom->Rndm();			
			hetaGen->Fill(eta);
			Aeta[itrk]	= eta;
		}
		for (int itrk=0;itrk<NTRKP;itrk++){
			float eta1	=  0.2 + 0.6*gRandom->Rndm();	
			float eta2;
			if (eta1<0){ eta2 = eta1 + 0.7 + gRandom->Gaus(0.0,0.2); } else
			if (eta1>0){ eta2 = eta1 - 0.7 + gRandom->Gaus(0.0,0.2); }
			hetaGen->Fill(eta1);
			hetaGen->Fill(eta2);
			Aeta[NTRK-2*NTRKP+itrk]		= eta1;
			Aeta[NTRK-2*NTRKP+itrk+1]	= eta2;
		}
		//
		//for (int i=0;i<NTRK;i++){cout<<Aeta[i]<<" "; } cout<<endl;
		std::random_shuffle(Aeta,Aeta+NTRK);
		//for (int i=0;i<NTRK;i++){cout<<Aeta[i]<<" "; } cout<<endl;cout<<endl;
		//
		for (int i=0;i<NTRK;i++){
			heta1D->Fill(Aeta[i]);
			for (int j=0;j<NTRK;j++){
				if (i!=j){
					hdeta->Fill(Aeta[j]-Aeta[i]);
					heta2D->Fill(Aeta[i],Aeta[j]);
					hR2->Fill(Aeta[i],Aeta[j]);
				}
			}
		}
		//
	}
	//
	heta1D->Scale(1./NEVT);
	heta2D->Scale(1./NEVT);
	hR2->Scale(1./NEVT);
	//
	int nbin	= heta1D->GetNbinsX();
	for (int ibin=1;ibin<=nbin;ibin++){
		float valx1	= heta1D->GetBinCenter(ibin);
		float valn1	= heta1D->GetBinContent(ibin);
		for (int jbin=1;jbin<=nbin;jbin++){
			float valx2	= heta1D->GetBinCenter(jbin);
			float valn2	= heta1D->GetBinContent(jbin);
			hetaT->Fill(valx1,valx2,valn1*valn2);
			hm1->Fill(valx1,valx2,-1);
		}
	}
	//
	//X hC2 = (TH2D*)heta2D->Clone("hC2");
	//X hC2->SetTitle("hC2");
	//X hC2->Add(hetaT,-1.0);
	//TH2D *hR2	= (TH2D*)heta2D->Clone("hR2");
	//hR2->SetTitle("R_{2} vs (#eta_{1},#eta_{2})");
	hR2->Divide(hetaT);
	hR2->Add(hm1);
	//
//cout<<" c "<<endl;
	//
	//hR2->Print("range");
	//
	//---- get and apply Baseline
	float base = GetC2Baseline(hmultGen);
	for (int ix=0;ix<NBETA;ix++){
		for (int iy=0;iy<NBETA;iy++){
			float oval		 = hR2->GetBinContent(ix+1,iy+1);
			hR2->SetBinContent(ix+1,iy+1,oval-base);
		}
	}
	//
//cout<<" d "<<endl;
	//
	float dAvg_deta[2*NBETA-1]	= {0};
	float dAvg_N[2*NBETA-1]		= {0};
	float dAvg_S[2*NBETA-1]		= {0};
	for (int ix=0;ix<NBETA;ix++){
		for (int iy=0;iy<NBETA;iy++){
			int k=ix-iy+NBETA-1;
			dAvg_N[k]		+= 1.0;
			dAvg_S[k]		+= hR2->GetBinContent(ix+1,iy+1);
			dAvg_deta[k]	+= hR2->GetXaxis()->GetBinCenter(ix+1) - hR2->GetYaxis()->GetBinCenter(iy+1);
			//
			//cout<<ix<<" "<<iy<<" "<<k<<" "<<dAvg_N[k]
			//		<<" \t "<<hR2->GetBinContent(ix+1,iy+1)
			//		<<" \t "<<hR2->GetXaxis()->GetBinCenter(ix+1)
			//		<<" "   <<hR2->GetYaxis()->GetBinCenter(iy+1)
			//		<<" \t "<<dAvg_S[k]<<" "<<dAvg_deta[k]<<endl;
		}	
	}
	for (int k=0;k<2*NBETA-1;k++){
		if (dAvg_N[k]){
			hR2deta->Fill(dAvg_deta[k]/dAvg_N[k],dAvg_S[k]/dAvg_N[k]);
		}
	}
	//
	int nbinx	= heta2D->GetNbinsX();
	int nbiny	= heta2D->GetNbinsY();
	for (int ibin=1;ibin<=nbinx;ibin++){
		for (int jbin=1;jbin<=nbiny;jbin++){
			float val	= hC2->GetBinContent(ibin,jbin);
			hcorr1D->Fill(val);			
		}
	}	
	//
	//---- calculate integral...	
	float integral	= 0.0;
	for (int ibin=1;ibin<=hR2deta->GetNbinsX();ibin++){
		float deta	= hR2deta->GetBinCenter(ibin);
		float val	= hR2deta->GetBinContent(ibin);
		float vale	= hR2deta->GetBinError(ibin);		// INCORRECT ERRORS......... 
		float bw	= hR2deta->GetBinWidth(ibin);
		hR2detabase->SetBinContent(ibin,val);			// straight copy now...
		hR2detabase->SetBinError(ibin,vale);			// INCORRECT ERRORS......... 
		if (deta>=0.0){
			hR2detabasevaldist->Fill(val);
		}
		//hR2detabase->SetBinError(ibin,0.0);				// INCORRECT ERRORS......... 
		//hR2deta->SetBinError(ibin,0.0);					// INCORRECT ERRORS......... 
		integral	+= (val-base)*bw;
	}
	//
	float rms = 1000.*hR2detabasevaldist->GetRMS();
	//
	cout<<"Baseline = "<<base
		<<"\t <R2>-baseline Integral = "<<integral
		<<"\t RMS = "<<rms
		<<"\t OVER/UNDER = "
		<<hR2detabasevaldist->GetBinContent(0)<<" "
		<<hR2detabasevaldist->GetBinContent(hR2detabasevaldist->GetNbinsX()+1)
		<<endl;


	//---- plotting setup
	//
	gStyle->SetPaperSize(TStyle::kUSLetter);
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetTitleXSize(0.055);
	gStyle->SetTitleYSize(0.055);
	gStyle->SetTitleOffset(0.85,"X");
	gStyle->SetTitleOffset(1.2,"Y");
	gStyle->SetOptStat(111110);
	gStyle->SetStatStyle(0); 
	gStyle->SetTitleStyle(0); 
	gStyle->SetStatX(0.94);
	gStyle->SetStatY(0.92);
	gStyle->SetStatH(0.26);
	gStyle->SetStatW(0.4);
	//gStyle->SetTitleSize(0.6);
	gStyle->SetErrorX(0.0001);
	gStyle->SetPadRightMargin(0.06);
	gStyle->SetPadTopMargin(0.08);
	gStyle->SetPadBottomMargin(0.05);
	gStyle->SetPadLeftMargin(0.08);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleY(1.0);
	gStyle->SetTitleW(0.75);
	gStyle->SetTitleH(0.075);
	gStyle->SetTitleTextColor(1);
	gStyle->SetTitleSize(0.1,"T");
	gStyle->SetPalette(1);
	gStyle->SetHistMinimumZero(kFALSE);
	//
	gStyle->SetHatchesSpacing(2);
	gStyle->SetHatchesLineWidth(2);

	double ylimits[2]	= {0};

	//cout<<"plot start"<<endl;

	if (PLOTON){

		//---- 
		//
		++ican;
		sprintf(buf,"ccan%d",ican);
		ccan[ican] = new TCanvas(buf,buf,30*ican,30*ican,800,(8.5/11.)*800);
		ccan[ican]->cd(); ccan[ican]->Divide(3,2,0.0001,0.0001);
			ccan[ican]->cd(1);
	//cout<<"1"<<endl;
				hetaGen->SetMinimum(0.5);
				hetaGen->SetStats(0);
				hetaGen->Draw();
			ccan[ican]->cd(2);
	//cout<<"2"<<endl;
				heta2D->SetStats(0);
				heta2D->Draw("colz");
			ccan[ican]->cd(3);
	//cout<<"3"<<endl;
	//			heta2D->SetStats(0);
	//			heta2D->Draw("lego2");
				hmultGen->Draw();
			ccan[ican]->cd(4);
	//cout<<"4"<<endl;
	//			hC2->SetStats(0);
	//			hC2->Draw("colz");
	//			hdeta->Draw();
			hR2detabasevaldist->Draw();
			ccan[ican]->cd(5);
	//cout<<"5"<<endl;
	//			hR2->SetStats(0);
	//			hR2->SetMinimum(-0.1);
	//			hR2->SetMaximum( 0.1);
				hR2->SetStats(0);
				hR2->Draw("colz");
			ccan[ican]->cd(6);
	//cout<<"6"<<endl;
				hR2detabase->SetStats(0);
				hR2detabase->SetMinimum(-0.005);
				hR2detabase->SetMaximum( 0.005);
				hR2detabase->SetMarkerStyle(20);
				hR2detabase->SetMarkerSize(1);
				hR2detabase->SetMarkerColor(4);
				hR2detabase->SetLineColor(4);
				hR2detabase->Draw("hist");
				//hR2deta->SetMarkerStyle(20);
				//hR2deta->SetMarkerSize(1);
				//hR2deta->SetMarkerColor(16);
				//hR2deta->SetLineColor(16);
				//hR2deta->Draw("P same");
	//cout<<"7"<<endl;
		ccan[ican]->cd(); 
	//cout<<"8"<<endl;
		ccan[ican]->Update();
	//cout<<"9"<<endl;
	//cout<<plotfileO<<endl;
	//cout<<plotfileO.Data()<<endl;
		ccan[ican]->Print(plotfileO.Data());
		//

	}	// PLOTON
	
//	cout<<"closeout"<<endl;
	
	if (PLOTON){
	 	cout<<" You plotted "<<ican+1<<" canvasses......."<<endl;
	 	ccan[ican]->Print(plotfileC.Data());
	 	sprintf(buf,"/usr/bin/pstopdf %s -o %s",plotfile.Data(),plotfilePDF.Data());
	 	//cout<<" "<<buf<<endl;
	 	gSystem->Exec(buf);
	 	//sprintf(buf,"/bin/cp %s /Library/WebServer/WebPages/files/",plotfilePDF.Data());
	 	//cout<<" "<<buf<<endl;
	 	//gSystem->Exec(buf);
	}
	//
	//cout<<"rootfile"<<endl;
	//
// 	TFile *fout = new TFile(rootoutfile.Data(),"RECREATE");
// 	cout<<"Writing "<<nRegisteredG<<" tgraph's to "<<rootoutfile.Data()<<endl;
// 	for (int i=0;i<nRegisteredG;i++){
// 		RegisterG[i]->Write();
// 	}
// 	cout<<"Writing "<<nRegisteredH<<" th1d's to "<<rootoutfile.Data()<<endl;
// 	for (int i=0;i<nRegisteredH;i++){
// 		RegisterH[i]->Write();
// 	}
// 	cout<<"Writing "<<nRegisteredH2<<" th2d's to "<<rootoutfile.Data()<<endl;
// 	for (int i=0;i<nRegisteredH2;i++){
// 		RegisterH2[i]->Write();
// 	}
// 	fout->Close();

	//cout<<"done"<<endl;
	
	return rms;
	
}

//----------------------------------------------------------------------------
void ecIncrement(int ntrk, float* Aeta, TH1D* heta1D, TH2D* heta2D){
	//
}
		
//----------------------------------------------------------------------------
void ecCalculate(TH1D* heta1D, TH2D* heta2D, TH2D* hetaT, TH2D* hR2, TH2D *hC2, TH1D* hcorr1D){
	//

}

//----------------------------------------------------------------------------
void AddToRegisterG(TGraph *g){
	RegisterG[nRegisteredG]	= g;
	++nRegisteredG;
}
//----------------------------------------------------------------------------
void AddToRegisterH(TH1D *h){
	RegisterH[nRegisteredH]	= h;
	++nRegisteredH;
}
//----------------------------------------------------------------------------
void AddToRegisterH2(TH2D *h2){
	RegisterH2[nRegisteredH2]	= h2;
	++nRegisteredH2;
}

//----------------------------------------------------------------------------
float GetC2Baseline(TH1D *hmult){
	if (!hmult){ exit(0); }
	float nent 		= (float)hmult->GetEntries();
	if (!nent){ exit(0); }
	TH1D *hmultwork	= (TH1D*)hmult->Clone("hmultwork");
	hmultwork->Scale(1./nent);
	int nbx			= hmultwork->GetNbinsX();
	float sumnum	= 0;
	float sumden	= 0;
	for (int ibx=1;ibx<=nbx;ibx++){
		float ni	 = hmultwork->GetBinCenter(ibx);
		float yi	 = hmultwork->GetBinContent(ibx);		
		sumnum		+= yi*ni*(ni-1);
		sumden		+= yi*ni;
	}
	float result 	= 0;
	if (sumden!=0.0){
		result		= sumnum/sumden/sumden - 1.;
	}
	return result;
}


