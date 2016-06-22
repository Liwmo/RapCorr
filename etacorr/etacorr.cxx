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

float 	etaCorrelations();
double * getXValsForLogAxis(int, float, float);
float	GetC2Baseline(TH1D*);
float	GetC3Baseline(TH1D*);

#if !defined(ARRAY_SIZE)
    #define ARRAY_SIZE(x) (sizeof((x)) / sizeof((x)[0]))
#endif

int main(int argc, char **argv) {
	gRandom->SetSeed(123456);
	float rms = etaCorrelations();
	return 0;
}

float etaCorrelations() {

	TH1::AddDirectory(kFALSE);
	char buf[200];
	int iCanvas = -1;
	TCanvas *canvases[100];
	TString plotFileBase, rootOutFile;
	TString plotFile0, plotFile, plotFileC, plotFilePDF;

	const int numBins = 200;
	float axis1 = 1.0E-10;
	float axis2 = 1.0;
	double * binXVals = new double[numBins + 1];
	binXVals = getXValsForLogAxis(numBins, axis1, axis2);

	if(PLOT_ON) {
		rootOutFile	= TString(Form("./root/etacorr.root"));
		plotFileBase = TString(Form("./ps/etacorr"));
		cout << "root file = " << rootOutFile.Data() << endl;
		cout << "plot file = " << plotFileBase.Data() << endl;
		plotFile0 = plotFileBase + TString(".ps(");
		plotFile = plotFileBase + TString(".ps");
		plotFileC = plotFileBase + TString(".ps]");
		plotFilePDF	= plotFileBase + TString(".pdf");
	}

	const int NBETA	= 35;
	const float ETAMAX = 0.7;
	TH1D *hMultGen = new TH1D("hMultGen", "hMultGen", 80, -0.5, 79.5);
	TH1D *hEtaGen = new TH1D("hEtaGen", "Generated #eta", NBETA, -ETAMAX, ETAMAX);
	TH1D *hdEta	= new TH1D("hdEta", "#Delta#eta", 2 * NBETA - 1, -2. * ETAMAX, 2. * ETAMAX);
	TH1D *hEta1D = new TH1D("hEta1D", "hEta1D", NBETA, -ETAMAX, ETAMAX);
	TH2D *hEta2D = new TH2D("hEta2D", "#eta_{2} vs #eta_{1}", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH2D *hEtaT2D = new TH2D("hEtaT2D", "hEtaT2D", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH2D *hR2 = new TH2D("hR2", "R_{2} vs (#eta_{1},#eta_{2})", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH2D *hC2 = new TH2D("hC2", "C_{2} vs (#eta_{1},#eta_{2})", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH1D *hCorr1D = new TH1D("hCorr1D", "hCorr1D", numBins, binXVals);
	TH2D *hm2D = new TH2D("hm2D", "hm2D", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH1D *hR2dEta = new TH1D("hR2dEta", "#LTR_{2}#GT vs #eta_{1}-#eta_{2}", 2 * NBETA - 1, -2. * ETAMAX, 2. * ETAMAX);
	TH1D *hR2dEta_N = new TH1D("hR2dEta_N", "NBINS vs #eta_{1}-#eta_{2}", 2 * NBETA - 1, -2. * ETAMAX, 2. * ETAMAX);
	TH1D *hR2dEtaBase = new TH1D("hR2dEtaBase", "#LTR_{2}#GT-Baseline vs #eta_{1}-#eta_{2}", 2 * NBETA - 1, -2. * ETAMAX, 2. * ETAMAX);

	hMultGen->Reset();	
	hEtaGen->Reset();
	hdEta->Reset();	
	hEta1D->Reset();
	hEta2D->Reset();
	hEtaT2D->Reset();	
	hR2->Reset();
	hC2->Reset();
	hCorr1D->Reset();
	hm2D->Reset();
	hR2dEta->Reset();
	hR2dEta_N->Reset();
	hR2dEtaBase->Reset();

	TH3D *hEta3D = new TH3D("hEta3D", "hEta3D", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH3D *hR3 = new TH3D("hR3", "hR3", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH3D *hEtaT3D = new TH3D("hEtaT3D", "hEtaT3D", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH3D *hC2rho1 = new TH3D("hC2rho1","hC2rho1", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH3D *hm3D = new TH3D("hm3D", "hm3D", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH2D *hR3dEta = new TH2D("hR3dEta", "#LTR_{3}#GT vs (#Delta#eta_{12}, #Delta#eta_{13})", 2 * NBETA - 1, -2. * ETAMAX, 2.* ETAMAX, 2* NBETA - 1, -2. * ETAMAX, 2. * ETAMAX);
	TH2D *hR3dEta_N	= new TH2D("hR3dEta_N", "NBINS vs (#Delta#eta_{12}, #Delta#eta_{13})", 2 * NBETA - 1, -2. * ETAMAX, 2. * ETAMAX, 2 * NBETA - 1, -2. * ETAMAX, 2. * ETAMAX);
		
	int freqNumBins	= 500;
	float freqBinWidth = 0.0001; 
	TH1D *hR2dEtaBaseValDist = new TH1D("hR2dEtaBaseValDist", "Frequency #LTR_{2}#GT-Baseline",
			freqNumBins, -freqNumBins * freqBinWidth / 2., freqNumBins * freqBinWidth / 2.);	// bw=0.0002

	const int N_EVENTS = 10000000;
	      int N_TRACK = 4;
	      int N_TRACK_PAIRS	= 0; // number of smeared dEta pairs added
	
	float *etaArr = new float[N_TRACK];
	
	for(int iEvent = 0; iEvent < N_EVENTS; iEvent++) {
		if(iEvent % (N_EVENTS / 10) == 0){ 
			cout << "processing " << iEvent << " of " << N_EVENTS << endl; 
		}
		hMultGen->Fill(N_TRACK);
		
		for(int iTrack = 0; iTrack < N_TRACK - 2 * N_TRACK_PAIRS; iTrack++) {		// random sampled...
			float eta = -1.0 + 2.0 * gRandom->Rndm();			
			hEtaGen->Fill(eta);
			etaArr[iTrack] = eta;
		}

		std::random_shuffle(etaArr, etaArr + N_TRACK);

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
			valxi = hEta2D->GetXaxis()->GetBinCenter(ibin);
			valxj = hEta2D->GetYaxis()->GetBinCenter(jbin);
			valnij = hEta2D->GetBinContent(ibin, jbin);

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
	
	float integralC2 = 0.0;
	for(int ibin = 1; ibin <= hR2dEta->GetNbinsX(); ibin++) {
		float dEta	 = hR2dEta->GetBinCenter(ibin);
		float val	 = hR2dEta->GetBinContent(ibin);
		float valErr = hR2dEta->GetBinError(ibin);		// INCORRECT ERRORS VALS
		float bWidth = hR2dEta->GetBinWidth(ibin);
		hR2dEtaBase->SetBinContent(ibin, val);			// straight copy now...
		hR2dEtaBase->SetBinError(ibin, valErr);			// INCORRECT ERRORS VALS.
		if(dEta >= 0.0) {
			hR2dEtaBaseValDist->Fill(val);
		}
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
				hR3dEta->Draw("COLZ1");
			canvases[iCanvas]->cd(3);
				hR3dEta_N->Draw("COLZ1");
			canvases[iCanvas]->cd(4);
			canvases[iCanvas]->cd(5);
			canvases[iCanvas]->cd(6);
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Update();
		canvases[iCanvas]->Print(plotFile.Data());
		gStyle->SetPadRightMargin(0.06);
	}	
	
	if (PLOT_ON){
	 	cout << " You plotted " << iCanvas + 1 << " canvases......." << endl;
	 	canvases[iCanvas]->Print(plotFileC.Data());
	 	sprintf(buf, "/usr/bin/pstopdf %s -o %s", plotFile.Data(), plotFilePDF.Data());
	 	cout << " " << buf << endl;
	 	gSystem->Exec(buf);
	}
	
	delete[] binXVals;
	return rms;
	
}

double * getXValsForLogAxis(int numBins, float axis1, float axis2) {
	float logAxis1 = TMath::Log10(axis1);
	float logAxis2 = TMath::Log10(axis2);
	float binWidth = (logAxis2 - logAxis1) / ((float)numBins);
	double * binXVals = new double[numBins + 1];
	binXVals[0] = axis1;

	for(int i = 1; i <= numBins; i++) {
		binXVals[i] = TMath::Power(10, logAxis1 + i * binWidth);
	}

	return binXVals;
}


float GetC2Baseline(TH1D *hmult){
	if(!hmult) { 
		exit(0); 
	}
	float nent = (float)hmult->GetEntries();
	if(!nent) {
		exit(0); 
	}
	TH1D *hMultWork	= (TH1D*)hmult->Clone("hMultWork");
	hMultWork->Scale(1. / nent);
	int nBinsX = hMultWork->GetNbinsX();
	float sumNum = 0;
	float sumDen = 0;
	for(int ibx = 1; ibx <= nBinsX; ibx++) {
		float ni = hMultWork->GetBinCenter(ibx);
		float yi = hMultWork->GetBinContent(ibx);		
		sumNum += yi * ni * (ni - 1.);
		sumDen += yi * ni;
	}
	float result = 0;
	if(sumDen != 0.0) {
		result = sumNum / sumDen / sumDen - 1.;
	}
	delete hMultWork; hMultWork = 0;
	return result;
}

float GetC3Baseline(TH1D *hmult){
	if(!hmult) { 
		exit(0); 
	}
	float nent = (float)hmult->GetEntries();
	if(!nent) { 
		exit(0); 
	}
	TH1D *hMultWork	= (TH1D*)hmult->Clone("hMultWork");
	hMultWork->Scale(1./nent);
	int nBinsX = hMultWork->GetNbinsX();
	float sumNum3 = 0;
	float sumNum2 = 0;
	float sumDen = 0;
	for(int ibx = 1; ibx <= nBinsX; ibx++) {
		float ni = hMultWork->GetBinCenter(ibx);		
		float yi = hMultWork->GetBinContent(ibx);		
		sumNum3 += yi * ni * (ni - 1.) * (ni - 2.);
		sumNum2	+= yi * ni * (ni - 1.);
		sumDen += yi * ni;
	}
	float result = 0;
	float term3 = 0;
	float term2	= 0;
	if(sumDen != 0.0) {
		term3 = sumNum3 / sumDen / sumDen / sumDen;
		term2 = sumNum2 / sumDen / sumDen;
		result = term3 - 3 * term2 + 2.;
	}
	delete hMultWork; hMultWork=0;
	return result;
}
