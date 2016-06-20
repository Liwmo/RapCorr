int PLOT_ON	= 1;

#include <string>
#include "stdio.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
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
float	GetC2Baseline(TH1D*);
float	GetC3Baseline(TH1D*);
float 	WhiteOut(TH2D*,TH2D*);

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

int main(int argc, char **argv) {

	//---- set stack size to 63MB (default=8.192kB)
	const	rlim_t kStackSize = 63L * 1024L * 1024L;
	struct	rlimit rl;
	int		result = getrlimit(RLIMIT_STACK, &rl);
	//cout<<"Initial RLIMIT_STACK = "<<rl.rlim_cur/1024L/1024L<<" MB"<<endl;
	if(result!=0) { 
		fprintf(stderr, "Problem ...getrlimit returned %d\n", result);
	} else {
		if(rl.rlim_cur < kStackSize) {
			rl.rlim_cur = kStackSize;
			result 		= setrlimit(RLIMIT_STACK, &rl);
			if(result != 0) {
				fprintf(stderr, "Problem ...setrlimit returned %d\n", result);
			}
		}
	}
	result = getrlimit(RLIMIT_STACK, &rl);
	//cout<<"Present RLIMIT_STACK = "<<rl.rlim_cur/1024L/1024L<<" MB"<<endl;

//	TApplication theApp("App", &argc, argv);
	TCanvas *mydummycanvas = new TCanvas();

	gRandom->SetSeed(123456);
//	gRandom->SetSeed(0);

//	const int NRUN	= 100;
//	float avg		=   0;
//	for (int irun=0;irun<NRUN;irun++){
		float rms	= ec();
//		avg	+= rms/((float)NRUN);
//	}	
//	cout<<endl;
//	cout<<"NRUN = "<<NRUN<<"    <rms> = "<<avg<<endl;
//	cout<<endl;

	return 0;
}


//---------------------------------------------------------------------------------

float ec(){

	TH1::AddDirectory(kFALSE);
	char 		buf[200],bufb[200],bufc[200];
	int 		iCanvas=-1,iFrame=-1,iLeg=-1;
	TCanvas 	*canvases[100];
	TH1F		*frame[100];
	TLegend		*legend[100];
	int			icolstart_energy,icolstart_sets;
	TString 	plotFileBase,rootOutFile;
	TString 	plotFile0,plotFile,plotFileC,plotFilePDF;

	const int numBins = 200;
	float axis1		  = 1.0E-10;
	float axis2		  = 1.0;
	float logAxis1	  = TMath::Log10(axis1);
	float logAxis2	  = TMath::Log10(axis2);
	float binWidth	  = (logAxis2 - logAxis1) / ((float)numBins);
	double xb[numBins+1];
	xb[0] = axis1;
	for (int i = 1; i <= numBins; i++){
		xb[i] = TMath::Power(10, logAxis1 + i * binWidth);
	}

	if(PLOT_ON) {
		rootOutFile		= TString(Form("./root/etacorr.root"));
		plotFileBase	= TString(Form("./ps/etacorr"));
		cout << "root file = " << rootOutFile.Data() << endl;
		cout << "plot file = " << plotFileBase.Data() << endl;
		plotFile0		= plotFileBase + TString(".ps(");
		plotFile		= plotFileBase + TString(".ps");
		plotFileC		= plotFileBase + TString(".ps]");
		plotFilePDF		= plotFileBase + TString(".pdf");
	}

	const int   NBETA	= 35;
	const float ETAMAX	= 0.7;
	TH1D *hMultGen		= new TH1D("hMultGen","hMultGen",80,-0.5,79.5);
	TH1D *hEtaGen		= new TH1D("hEtaGen","Generated #eta",NBETA,-ETAMAX,ETAMAX);
	TH1D *hdEta			= new TH1D("hdEta","#Delta#eta",2*NBETA-1,-2.*ETAMAX,2.*ETAMAX);
	TH1D *hEta1D		= new TH1D("hEta1D","hEta1D",NBETA,-ETAMAX,ETAMAX);
	TH2D *hEta2D		= new TH2D("hEta2D","#eta_{2} vs #eta_{1}",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH2D *hEtaT2D		= new TH2D("hEtaT2D","hEtaT2D",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH2D *hR2			= new TH2D("hR2","R_{2} vs (#eta_{1},#eta_{2})",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH2D *hC2			= new TH2D("hC2","C_{2} vs (#eta_{1},#eta_{2})",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH1D *hCorr1D		= new TH1D("hCorr1D","hCorr1D",numBins,xb);
	TH2D *hm2D			= new TH2D("hm2D","hm2D",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH1D *hR2dEta		= new TH1D("hR2dEta","#LTR_{2}#GT vs #eta_{1}-#eta_{2}",2*NBETA-1,-2.*ETAMAX,2.*ETAMAX);
	TH1D *hR2dEta_N		= new TH1D("hR2dEta_N","NBINS vs #eta_{1}-#eta_{2}",2*NBETA-1,-2.*ETAMAX,2.*ETAMAX);
	TH1D *hR2dEtaBase	= new TH1D("hR2dEtaBase","#LTR_{2}#GT-Baseline vs #eta_{1}-#eta_{2}",2*NBETA-1,-2.*ETAMAX,2.*ETAMAX);

	hMultGen	->Reset();	
	hEtaGen		->Reset();
	hdEta		->Reset();	
	hEta1D		->Reset();
	hEta2D		->Reset();
	hEtaT2D		->Reset();	
	hR2			->Reset();
	hC2			->Reset();
	hCorr1D		->Reset();
	hm2D		->Reset();
	hR2dEta		->Reset();
	hR2dEta_N	->Reset();
	hR2dEtaBase	->Reset();

	TH3D *hEta3D		= new TH3D("hEta3D","hEta3D",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH3D *hR3			= new TH3D("hR3","hR3",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH3D *hEtaT3D		= new TH3D("hEtaT3D","hEtaT3D",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH3D *hC2rho1		= new TH3D("hC2rho1","hC2rho1",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH3D *hm3D			= new TH3D("hm3D","hm3D",NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX,NBETA,-ETAMAX,ETAMAX);
	TH2D *hR3dEta		= new TH2D("hR3dEta","#LTR_{3}#GT vs (#Delta#eta_{12},#Delta#eta_{13})",2*NBETA-1,-2.*ETAMAX,2.*ETAMAX,2*NBETA-1,-2.*ETAMAX,2.*ETAMAX);
	TH2D *hR3dEta_N		= new TH2D("hR3dEta_N","NBINS vs (#Delta#eta_{12},#Delta#eta_{13})",2*NBETA-1,-2.*ETAMAX,2.*ETAMAX,2*NBETA-1,-2.*ETAMAX,2.*ETAMAX);
	
	//==================================================
	
	int kfreqnb		= 500;
	float freqbw	= 0.0001; 
	TH1D *hR2dEtaBaseValDist	= new TH1D("hR2dEtaBaseValDist","Frequency #LTR_{2}#GT-Baseline",
			kfreqnb,-kfreqnb*freqbw/2., kfreqnb*freqbw/2.);	// bw=0.0002

	const int N_EVENTS   	= 10000000;
	      int N_TRACK    	=        4;
	      int N_TRACK_PAIRS	=        0;		// number of smeared dEta pairs added

	//==================================================
	
	float *etaArr = new float[N_TRACK];
	
	for(int iEvent = 0; iEvent < N_EVENTS;iEvent++) {
		if(iEvent % (N_EVENTS/10) == 0){ 
			cout << "processing " << iEvent << " of " << N_EVENTS << endl; 
		}
		
		//N_TRACK = (int)(61 + 8.1*gRandom->Gaus());
		hMultGen->Fill(N_TRACK);
		
		for(int iTrack = 0; iTrack < N_TRACK - 2*N_TRACK_PAIRS; iTrack++) {		// random sampled...
			float eta = -1.0 + 2.0 * gRandom->Rndm();			
			hEtaGen->Fill(eta);
			etaArr[iTrack] = eta;
		}

		// for (int iTrack=0;iTrack<N_TRACK_PAIRS;iTrack++){
		// 	float eta1	=  0.2 + 0.6*gRandom->Rndm();	
		// 	float eta2;
		// 	if (eta1<0){ eta2 = eta1 + 0.7 + gRandom->Gaus(0.0,0.2); } else
		// 	if (eta1>0){ eta2 = eta1 - 0.7 + gRandom->Gaus(0.0,0.2); }
		// 	hEtaGen->Fill(eta1);
		// 	hEtaGen->Fill(eta2);
		// 	etaArr[N_TRACK-2*N_TRACK_PAIRS+iTrack]		= eta1;
		// 	etaArr[N_TRACK-2*N_TRACK_PAIRS+iTrack+1]	= eta2;
		// }
		//
		//for (int i=0;i<N_TRACK;i++){cout<<etaArr[i]<<" "; } cout<<endl;
		std::random_shuffle(etaArr, etaArr + N_TRACK);
		//for (int i=0;i<N_TRACK;i++){cout<<etaArr[i]<<" "; } cout<<endl;cout<<endl;
		//
		//---- calculate rho2 -> hEtaT2D, hR2
		//---- calculate rho3 -> hEtaT3D, hR3

		for(int i = 0; i < N_TRACK; i++) {
			hEta1D->Fill(etaArr[i]);
			
			for(int j = 0; j < N_TRACK; j++) {
				if(i != j) {
					hdEta->Fill(etaArr[j] - etaArr[i]);
					hEta2D->Fill(etaArr[i], etaArr[j]);
					hR2->Fill(etaArr[i], etaArr[j]);
					
					for (int k = 0; k < N_TRACK; k++) {
						if (i != k && j != k) {
							hEta3D->Fill(etaArr[i], etaArr[j], etaArr[k]);
							hR3->Fill(etaArr[i], etaArr[j], etaArr[k]);
						}
					}
				}
			}
		}
		
	}

	cout << "Normalizing..." << endl;
	hEta1D->Scale(1. / N_EVENTS);
	hEta2D->Scale(1. / N_EVENTS);
	hEta3D->Scale(1. / N_EVENTS);
	hR2->Scale(1. / N_EVENTS);
	hR3->Scale(1. / N_EVENTS);
	//---- calculate rho1*rho1      -> hEtaT2D
	//---- calculate rho1*rho1*rho1 -> hEtaT3D
	//---- fill TH2D all with -1    -> hm_2D
	//---- fill TH3D all with -1    -> hm_3D
	cout << "Calculating rho1*rho1 & rho1*rho1*rho1..." << endl;
	int nbin = hEta1D->GetNbinsX();
	for(int ibin = 1; ibin <= nbin; ibin++) {
		float valx1	= hEta1D->GetBinCenter(ibin);
		float valn1	= hEta1D->GetBinContent(ibin);

		for(int jbin = 1; jbin <= nbin; jbin++) {
			float valx2	= hEta1D->GetBinCenter(jbin);
			float valn2	= hEta1D->GetBinContent(jbin);
			hEtaT2D->Fill(valx1, valx2, valn1 * valn2);
			hm2D->Fill(valx1, valx2, -1.0);

			for(int kbin = 1; kbin <= nbin; kbin++) {
				float valx3	= hEta1D->GetBinCenter(kbin);
				float valn3	= hEta1D->GetBinContent(kbin);
				hEtaT3D->Fill(valx1, valx2, valx3, valn1 * valn2 * valn3);
				hm3D->Fill(valx1, valx2, valx3, +2.);
			}
		}
	}

	cout << "Calculating C2rho1..." << endl;
	float valxi, valxj, valxk, valnij, valnk;
	for(int ibin = 1; ibin <= nbin; ibin++) {
		for(int jbin = 1; jbin <= nbin; jbin++) {
			valxi	= hEta2D->GetXaxis()->GetBinCenter(ibin);
			valxj	= hEta2D->GetYaxis()->GetBinCenter(jbin);
			valnij	= hEta2D->GetBinContent(ibin, jbin);

			for(int kbin = 1; kbin <= nbin; kbin++) {
				valxk = hEta1D->GetBinCenter(kbin);
				valnk = hEta1D->GetBinContent(kbin);
				hC2rho1->Fill(valxi, valxj, valxk, valnij * valnk);
			}
		}
	}
	
	cout<<"Calculating R2 and R3..."<<endl;
	hR2->Divide(hEtaT2D);
	hR2->Add(hm2D);
	
	hR3->Divide(hEtaT3D);
	hC2rho1->Divide(hEtaT3D);
	hC2rho1->Scale(3.0);
	hR3->Add(hC2rho1,-1.);
	hR3->Add(hm3D);
	
	//hR2->Print("range");
	
	cout<<"Calculating Baselines..."<<endl;
	float baseC2 = GetC2Baseline(hMultGen);
	float baseC3 = GetC3Baseline(hMultGen);
	for(int ibx = 1; ibx <= NBETA; ibx++) {
		for (int iby = 1; iby <= NBETA; iby++) {
			float oval2	= hR2->GetBinContent(ibx, iby);
			hR2->SetBinContent(ibx, iby, oval2 - baseC2);

			for(int ibz = 1; ibz <= NBETA; ibz++) {
				float oval3	= hR3->GetBinContent(ibx, iby, ibz);
				hR3->SetBinContent(ibx, iby, ibz, oval3 - baseC3);
			}
		}
	}

	cout<<"Averaging in dEta slices..."<<endl;
	for(int ibx = 1; ibx <= NBETA; ibx++) {								// x-bin
		for(int iby = 1; iby <= NBETA; iby++) {							// y-bin
			float dEtaXY, dEtaXZ;									
			
			dEtaXY = hR2->GetXaxis()->GetBinCenter(ibx) 		
			       - hR2->GetYaxis()->GetBinCenter(iby);
			hR2dEta->Fill(dEtaXY,hR2->GetBinContent(ibx, iby));
			hR2dEta_N->Fill(dEtaXY, 1.0);
			
			for(int ibz = 1; ibz <= NBETA; ibz++) {					// z-bin
				dEtaXY	= hR3->GetXaxis()->GetBinCenter(ibx) 		//
			            - hR3->GetYaxis()->GetBinCenter(iby);		// this is the dy for this xy diagonal
				dEtaXZ	= hR3->GetXaxis()->GetBinCenter(ibx) 		//
			        	- hR3->GetZaxis()->GetBinCenter(ibz);		// this is the dy for this xz diagonal
			    //cout<<ibx<<" "<<iby<<" "<<ibz<<" \t"<<dEtaXY<<" "<<dEtaXZ<<" \t"<<hR3->GetBinContent(ibx,iby,ibz)<<endl;
			    hR3dEta->Fill(dEtaXY, dEtaXZ, hR3->GetBinContent(ibx, iby, ibz));
			    hR3dEta_N->Fill(dEtaXY, dEtaXZ, 1.0);
			}
		}	
	}
	hR2dEta->Divide(hR2dEta_N);
	hR3dEta->Divide(hR3dEta_N);
	
	int nbinx = hEta2D->GetNbinsX();
	int nbiny = hEta2D->GetNbinsY();
	for(int ibin = 1; ibin <= nbinx; ibin++) {
		for(int jbin = 1; jbin <= nbiny; jbin++) {
			float val = hC2->GetBinContent(ibin, jbin);
			hCorr1D->Fill(val);			
		}
	}	
	
	//---- calculate integral...	
	float integralC2 = 0.0;
	for(int ibin = 1; ibin <= hR2dEta->GetNbinsX(); ibin++) {
		float dEta	 = hR2dEta->GetBinCenter(ibin);
		float val	 = hR2dEta->GetBinContent(ibin);
		float valErr = hR2dEta->GetBinError(ibin);		// INCORRECT ERRORS VALS......... 
		float bWidth = hR2dEta->GetBinWidth(ibin);
		hR2dEtaBase->SetBinContent(ibin, val);			// straight copy now...
		hR2dEtaBase->SetBinError(ibin, valErr);			// INCORRECT ERRORS VALS......... 
		if(dEta >= 0.0) {
			hR2dEtaBaseValDist->Fill(val);
		}
		//hR2dEtaBase->SetBinError(ibin,0.0);				// INCORRECT ERROR VALS......... 
		//hR2dEta->SetBinError(ibin,0.0);					// INCORRECT ERRORS VALS......... 
		integralC2	+= val * bWidth;
	}
	
	float rms = 1000. * hR2dEtaBaseValDist->GetRMS();
	
	cout << "C2 Baseline = " << baseC2
		 << "\t <R2>-baseline Integral = " << integralC2
		 << "\t 1000. * RMS = " << rms
		 << "\t OVER/UNDER = "
		 << hR2dEtaBaseValDist->GetBinContent(0) << " "
		 << hR2dEtaBaseValDist->GetBinContent(hR2dEtaBaseValDist->GetNbinsX() + 1)
		 << endl;
	cout << "C3 Baseline = " << baseC3
		//<<"\t <R2>-baseline Integral = "<<integralC3
		//<<"\t RMS = "<<rms
		//<<"\t OVER/UNDER = "
		//<<hR2dEtaBaseValDist->GetBinContent(0)<<" "
		//<<hR2dEtaBaseValDist->GetBinContent(hR2dEtaBaseValDist->GetNbinsX()+1)
		 << endl;

	//---- plotting setup
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
	
	gStyle->SetHatchesSpacing(2);
	gStyle->SetHatchesLineWidth(2);

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

	double ylimits[2]	= {0};

	if(PLOT_ON) {

		//---- R2
		++iCanvas;
		sprintf(buf, "canvases%d", iCanvas);
		canvases[iCanvas] = new TCanvas(buf, buf, 30 * iCanvas, 30 * iCanvas, 800, (8.5 / 11.) * 800);
		canvases[iCanvas]->cd(); canvases[iCanvas]->Divide(3,2,0.0001,0.0001);
			canvases[iCanvas]->cd(1);
				hEtaGen->SetMinimum(0.5);
				hEtaGen->SetStats(0);
				hEtaGen->Draw();
			canvases[iCanvas]->cd(2);
				hEta2D->SetStats(0);
				hEta2D->Draw("colz");
			canvases[iCanvas]->cd(3);
				hMultGen->Draw();
			canvases[iCanvas]->cd(4);
				hR2dEtaBaseValDist->Draw();
			canvases[iCanvas]->cd(5);
				hR2->SetStats(0);
				hR2->Draw("colz");
			canvases[iCanvas]->cd(6);
				hR2dEtaBase->SetStats(0);
				hR2dEtaBase->SetMinimum(-0.005);
				hR2dEtaBase->SetMaximum( 0.005);
				hR2dEtaBase->SetMarkerStyle(20);
				hR2dEtaBase->SetMarkerSize(1);
				hR2dEtaBase->SetMarkerColor(4);
				hR2dEtaBase->SetLineColor(4);
				hR2dEtaBase->Draw("hist");
				//hR2dEta->SetMarkerStyle(20);
				//hR2dEta->SetMarkerSize(1);
				//hR2dEta->SetMarkerColor(16);
				//hR2dEta->SetLineColor(16);
				//hR2dEta->Draw("P same");
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Update();
		canvases[iCanvas]->Print(plotFile0.Data());
		
		//---- R3
		gStyle->SetOptStat(0);
		gStyle->SetPadRightMargin(0.15);
		++iCanvas;
		sprintf(buf,"canvases%d",iCanvas);
		canvases[iCanvas] = new TCanvas(buf,buf,30*iCanvas,30*iCanvas,800,(8.5/11.)*800);
		canvases[iCanvas]->cd(); canvases[iCanvas]->Divide(3,2,0.0001,0.0001);
			canvases[iCanvas]->cd(1);
				hR2dEtaBase->Draw("HIST");
			canvases[iCanvas]->cd(2);
				//float min = WhiteOut(hR3dEta,hR3dEta_N);
				//hR3dEta->SetMinimum(min);
				hR3dEta->Draw("COLZ1");
			canvases[iCanvas]->cd(3);
				hR3dEta_N->Draw("COLZ1");
			canvases[iCanvas]->cd(4);
				//hR3->Project3D("xy")->Draw("colz");
			canvases[iCanvas]->cd(5);
				//hR3->Project3D("yz")->Draw("colz");
			canvases[iCanvas]->cd(6);
				//hR3->Project3D("zx")->Draw("colz");
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Update();
		canvases[iCanvas]->Print(plotFile.Data());
		gStyle->SetPadRightMargin(0.06);
		//
	}	// PLOT_ON
	
	
	if (PLOT_ON){
	 	cout << " You plotted " << iCanvas + 1 << " canvasses......." << endl;
	 	canvases[iCanvas]->Print(plotFileC.Data());
	 	sprintf(buf, "/usr/bin/pstopdf %s -o %s", plotFile.Data(), plotFilePDF.Data());
	 	cout << " " << buf << endl;
	 	gSystem->Exec(buf);
	 	//sprintf(buf,"/bin/cp %s /Library/WebServer/WebPages/files/",plotFilePDF.Data());
	 	//cout<<" "<<buf<<endl;
	 	//gSystem->Exec(buf);
	}
	
// 	TFile *fout = new TFile(rootOutFile.Data(),"RECREATE");
// 	cout<<"Writing "<<nRegisteredG<<" tgraph's to "<<rootOutFile.Data()<<endl;
// 	for (int i=0;i<nRegisteredG;i++){
// 		RegisterG[i]->Write();
// 	}
// 	cout<<"Writing "<<nRegisteredH<<" th1d's to "<<rootOutFile.Data()<<endl;
// 	for (int i=0;i<nRegisteredH;i++){
// 		RegisterH[i]->Write();
// 	}
// 	cout<<"Writing "<<nRegisteredH2<<" th2d's to "<<rootOutFile.Data()<<endl;
// 	for (int i=0;i<nRegisteredH2;i++){
// 		RegisterH2[i]->Write();
// 	}
// 	fout->Close();
	
	return rms;
	
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
	if(!hmult) { exit(0); }
	float nent = (float)hmult->GetEntries();
	if(!nent) { exit(0); }
	TH1D *hmultwork	= (TH1D*)hmult->Clone("hmultwork");
	hmultwork->Scale(1. / nent);
	int nBinsX		= hmultwork->GetNbinsX();
	float sumnum	= 0;
	float sumden	= 0;
	for(int ibx = 1; ibx <= nBinsX; ibx++) {
		float ni	 = hmultwork->GetBinCenter(ibx);
		//if (ni<1.0) continue;
		float yi	 = hmultwork->GetBinContent(ibx);		
		sumnum		+= yi * ni * (ni - 1.);
		sumden		+= yi * ni;
	}
	float result = 0;
	if(sumden != 0.0) {
		result = sumnum / sumden / sumden - 1.;
	}
	delete hmultwork; hmultwork = 0;
	return result;
}

//----------------------------------------------------------------------------
float GetC3Baseline(TH1D *hmult){
	if (!hmult){ exit(0); }
	float nent 		= (float)hmult->GetEntries();
	if (!nent){ exit(0); }
	TH1D *hmultwork	= (TH1D*)hmult->Clone("hmultwork");
	hmultwork->Scale(1./nent);
	int nBinsX		= hmultwork->GetNbinsX();
	float sumnum3	= 0;
	float sumnum2	= 0;
	float sumden	= 0;
	for (int ibx=1;ibx<=nBinsX;ibx++){
		float ni	 = hmultwork->GetBinCenter(ibx);		
		//if (ni<2.0) continue;
		float yi	 = hmultwork->GetBinContent(ibx);		
		sumnum3		+= yi*ni*(ni-1.)*(ni-2.);
		sumnum2		+= yi*ni*(ni-1.);
		sumden		+= yi*ni;
	}
	float result 	= 0;
	float term3		= 0;
	float term2		= 0;
	if (sumden!=0.0){
		term3		= sumnum3/sumden/sumden/sumden;
		term2		= sumnum2/sumden/sumden;
		result		= term3 - 3*term2 + 2.;
	}
	delete hmultwork; hmultwork=0;
	return result;
}

//----------------------------------------------------------------------------
float WhiteOut(TH2D* hist,TH2D* histN) {
   Int_t bin 	= 0;
   float min	= 9999999;
   for(Int_t i = 1; i <= hist->GetNbinsX(); ++i) {
      for(Int_t j = 1; j <= hist->GetNbinsY(); ++j) {
         bin = hist->GetBin(i, j);
         float n = histN->GetBinContent(bin);
         float val = hist->GetBinContent(bin);
         if (val<min) min=val;
		 if (n==0){
		 	hist->SetBinContent(bin,-999);
		 }
      }
   }
   if (min<0.0){ min *= 1.05; } else
   if (min>0.0){ min /= 1.05; }
   return min;
}