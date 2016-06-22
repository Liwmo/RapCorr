const int PLOT_ON = 1;

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

float etaCorrelations();
double * getXValsForLogAxis(int, float, float);
void resetHistograms(TH1**, int);
void setupOutputFilePaths(TString&, TString&, TString&, TString&);
float * fillRandomRapidities(TH1D*, int, int); 
void fill1DRapidityDist(TH1D*, TH1D*, float*, int);
void fill2DRapidityDist(TH2D*, TH2D*, float*, int);
void fill3DRapidityDist(TH3D*, TH3D*, float*, int);
void normalizeHistograms(TH1**, int, int);
void fill2DTensorProduct(TH2D*, TH1D*);
void fill3DTensorProduct(TH3D*, TH1D*);
void fillConstant2DHistogram(TH2D*, float, TH1D*);
void fillConstant3DHistogram(TH3D*, float, TH1D*);
void fillC2rho1Histogram(TH3D*, TH2D*, TH1D*);
void calculateR2Histogram(TH2D*, TH2D*, TH2D*);
void calculateR3Histogram(TH3D*, TH3D*, TH3D*, TH3D*);
void applyC2BaselineAdjustment(TH2D*, float, int);
void applyC3BaselineAdjustment(TH3D*, float, int);
void setStyle();
float getC2Baseline(TH1D*);
float getC3Baseline(TH1D*);

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
	TString plotFile0, plotFile, plotFileC, plotFilePDF;

	const int numBins = 200;
	float axis1 = 1.0E-10;
	float axis2 = 1.0;
	double * binXVals = new double[numBins + 1];
	binXVals = getXValsForLogAxis(numBins, axis1, axis2);
	setupOutputFilePaths(plotFile0, plotFile, plotFileC, plotFilePDF);

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
	TH2D *hConst2D = new TH2D("hConst2D", "hConst2D", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH1D *hR2dEta = new TH1D("hR2dEta", "#LTR_{2}#GT vs #eta_{1}-#eta_{2}", 2 * NBETA - 1, -2. * ETAMAX, 2. * ETAMAX);
	TH1D *hR2dEta_N = new TH1D("hR2dEta_N", "NBINS vs #eta_{1}-#eta_{2}", 2 * NBETA - 1, -2. * ETAMAX, 2. * ETAMAX);
	TH1D *hR2dEtaBase = new TH1D("hR2dEtaBase", "#LTR_{2}#GT-Baseline vs #eta_{1}-#eta_{2}", 2 * NBETA - 1, -2. * ETAMAX, 2. * ETAMAX);

	int freqNumBins	= 500;
	float freqBinWidth = 0.0001; 
	TH1D *hR2dEtaBaseValDist = new TH1D("hR2dEtaBaseValDist", "Frequency #LTR_{2}#GT-Baseline",
			freqNumBins, -freqNumBins * freqBinWidth / 2., freqNumBins * freqBinWidth / 2.);
	TH1 * histoResets[13] = {hMultGen, hEtaGen, hdEta, hEta1D, hEta2D, hEtaT2D, 
		hR2, hC2, hCorr1D, hConst2D, hR2dEta, hR2dEta_N, hR2dEtaBase};
	resetHistograms(histoResets, 13);

	TH3D *hEta3D = new TH3D("hEta3D", "hEta3D", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH3D *hR3 = new TH3D("hR3", "hR3", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH3D *hEtaT3D = new TH3D("hEtaT3D", "hEtaT3D", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH3D *hC2rho1 = new TH3D("hC2rho1","hC2rho1", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH3D *hConst3D = new TH3D("hConst3D", "hConst3D", NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX, NBETA, -ETAMAX, ETAMAX);
	TH2D *hR3dEta = new TH2D("hR3dEta", "#LTR_{3}#GT vs (#Delta#eta_{12}, #Delta#eta_{13})", 2 * NBETA - 1, -2. * ETAMAX, 2.* ETAMAX, 2* NBETA - 1, -2. * ETAMAX, 2. * ETAMAX);
	TH2D *hR3dEta_N	= new TH2D("hR3dEta_N", "NBINS vs (#Delta#eta_{12}, #Delta#eta_{13})", 2 * NBETA - 1, -2. * ETAMAX, 2. * ETAMAX, 2 * NBETA - 1, -2. * ETAMAX, 2. * ETAMAX);

	const int N_EVENTS = 10000000;
	      int N_TRACKS = 4;
	      int N_TRACK_PAIRS	= 0;
	
	float *etaArr = new float[N_TRACKS];
	for(int iEvent = 0; iEvent < N_EVENTS; iEvent++) {
		if(iEvent % (N_EVENTS / 10) == 0) { 
			cout << "processing " << iEvent << " of " << N_EVENTS << endl; 
		}
		hMultGen->Fill(N_TRACKS);
		etaArr = fillRandomRapidities(hEtaGen, N_TRACKS, N_TRACK_PAIRS);
		std::random_shuffle(etaArr, etaArr + N_TRACKS);
		fill1DRapidityDist(hEta1D, hdEta, etaArr, N_TRACKS);
		fill2DRapidityDist(hEta2D, hR2, etaArr, N_TRACKS);
		fill3DRapidityDist(hEta3D, hR3, etaArr, N_TRACKS);		
	}

	cout << "Normalizing..." << endl;
	TH1 * histosNormalize[5] = {hEta1D, hEta2D, hEta3D, hR2, hR3};
	normalizeHistograms(histosNormalize, N_EVENTS, 5);

	cout << "Calculating rho1^2 & rho1^3 Tensor Products..." << endl;
	fill2DTensorProduct(hEtaT2D, hEta1D);
	fill3DTensorProduct(hEtaT3D, hEta1D);
	fillConstant2DHistogram(hConst2D, -1.0, hEta1D);
	fillConstant3DHistogram(hConst3D, +2.0, hEta1D);

	cout << "Calculating C2rho1..." << endl;
	fillC2rho1Histogram(hC2rho1, hEta2D, hEta1D);
	
	cout<<"Calculating R2 and R3..."<<endl;
	calculateR2Histogram(hR2, hEtaT2D, hConst2D);
	calculateR3Histogram(hR3, hEtaT3D, hC2rho1, hConst3D);
		
	cout<<"Calculating Baselines..."<<endl;
	float baseC2 = getC2Baseline(hMultGen);
	float baseC3 = getC3Baseline(hMultGen);
	applyC2BaselineAdjustment(hR2, baseC2, NBETA);
	applyC3BaselineAdjustment(hR3, baseC3, NBETA);

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

	setStyle();

	if(PLOT_ON) {

		//---- R2
		iCanvas++;
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
	delete[] etaArr;
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

void resetHistograms(TH1 ** histograms, int size) {
	for(int i = 0; i < size; i++) {
		histograms[i]->Reset();
	}
}

void setupOutputFilePaths(TString &plotFile0, TString &plotFile, 
	TString &plotFileC, TString &plotFilePDF) {

	TString plotFileBase, rootOutFile;

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
}

float * fillRandomRapidities(TH1D* hEtaGen, int N_TRACKS, int N_TRACK_PAIRS) {
	float *etaArr = new float[N_TRACKS];
	for(int iTrack = 0; iTrack < N_TRACKS - 2 * N_TRACK_PAIRS; iTrack++) {
		float eta = -1.0 + 2.0 * gRandom->Rndm();			
		hEtaGen->Fill(eta);
		etaArr[iTrack] = eta;
	}
	return etaArr;
}

void fill1DRapidityDist(TH1D *hEta1D, TH1D *hdEta, float * etaArr, int N_TRACKS) {
	for(int i = 0; i < N_TRACKS; i++) {
		hEta1D->Fill(etaArr[i]);
		
		for(int j = 0; j < N_TRACKS; j++) {
			if(i != j) {
				hdEta->Fill(etaArr[j] - etaArr[i]);
			}
		}
	}
}

void fill2DRapidityDist(TH2D *hEta2D, TH2D *hR2, float * etaArr, int N_TRACKS) {
	for(int i = 0; i < N_TRACKS; i++) {			
		for(int j = 0; j < N_TRACKS; j++) {
			if(i != j) {
				hEta2D->Fill(etaArr[i], etaArr[j]);
				hR2->Fill(etaArr[i], etaArr[j]);
			}
		}
	}
}

void fill3DRapidityDist(TH3D *hEta3D, TH3D *hR3, float * etaArr, int N_TRACKS) {
	for(int i = 0; i < N_TRACKS; i++) {
		for(int j = 0; j < N_TRACKS; j++) {
			for (int k = 0; k < N_TRACKS; k++) {
				if (i != j && i != k && j != k) {
					hEta3D->Fill(etaArr[i], etaArr[j], etaArr[k]);
					hR3->Fill(etaArr[i], etaArr[j], etaArr[k]);
				}
			}
		}
	}
}

void normalizeHistograms(TH1 **histos, int normConst, int size) {
	for(int i = 0; i < size; i++) {
		histos[i]->Scale(1. / normConst);
	}
}

void fill2DTensorProduct(TH2D *hEtaT2D, TH1D *hEta1D) {
	int nbin = hEta1D->GetNbinsX();
	for(int ibin = 1; ibin <= nbin; ibin++) {
		float valx1	= hEta1D->GetBinCenter(ibin);
		float valn1	= hEta1D->GetBinContent(ibin);
		for(int jbin = 1; jbin <= nbin; jbin++) {
			float valx2	= hEta1D->GetBinCenter(jbin);
			float valn2	= hEta1D->GetBinContent(jbin);
			hEtaT2D->Fill(valx1, valx2, valn1 * valn2);
		}
	}
}

void fill3DTensorProduct(TH3D *hEtaT3D, TH1D *hEta1D) {
	int nbin = hEta1D->GetNbinsX();
	for(int ibin = 1; ibin <= nbin; ibin++) {
		float valx1	= hEta1D->GetBinCenter(ibin);
		float valn1	= hEta1D->GetBinContent(ibin);
		for(int jbin = 1; jbin <= nbin; jbin++) {
			float valx2	= hEta1D->GetBinCenter(jbin);
			float valn2	= hEta1D->GetBinContent(jbin);
			for(int kbin = 1; kbin <= nbin; kbin++) {
				float valx3	= hEta1D->GetBinCenter(kbin);
				float valn3	= hEta1D->GetBinContent(kbin);
				hEtaT3D->Fill(valx1, valx2, valx3, valn1 * valn2 * valn3);
			}
		}
	}
}

void fillConstant2DHistogram(TH2D *hConst2D, float val, TH1D *hEta1D) {
	int nbin = hEta1D->GetNbinsX();
	for(int ibin = 1; ibin <= nbin; ibin++) {
		float valx1	= hEta1D->GetBinCenter(ibin);
		for(int jbin = 1; jbin <= nbin; jbin++) {
			float valx2	= hEta1D->GetBinCenter(jbin);
			hConst2D->Fill(valx1, valx2, val);
		}
	}
}

void fillConstant3DHistogram(TH3D* hConst3D, float val, TH1D *hEta1D) {
	int nbin = hEta1D->GetNbinsX();
	for(int ibin = 1; ibin <= nbin; ibin++) {
		float valx1	= hEta1D->GetBinCenter(ibin);
		for(int jbin = 1; jbin <= nbin; jbin++) {
			float valx2	= hEta1D->GetBinCenter(jbin);
			for(int kbin = 1; kbin <= nbin; kbin++) {
				float valx3	= hEta1D->GetBinCenter(kbin);
				hConst3D->Fill(valx1, valx2, valx3, val);
			}
		}
	}
}

void fillC2rho1Histogram(TH3D *hC2rho1, TH2D *hEta2D, TH1D *hEta1D) {
	int nbin = hEta1D->GetNbinsX();
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
}

void calculateR2Histogram(TH2D *hR2, TH2D *hEtaT2D, TH2D *hConst2D) {
	hR2->Divide(hEtaT2D);
	hR2->Add(hConst2D);
}

void calculateR3Histogram(TH3D *hR3, TH3D *hEtaT3D, TH3D *hC2rho1, TH3D *hConst3D) {
	hR3->Divide(hEtaT3D);
	hC2rho1->Divide(hEtaT3D);
	hC2rho1->Scale(3.0);
	hR3->Add(hC2rho1, -1.0);
	hR3->Add(hConst3D);
}

float getC2Baseline(TH1D *hmult){
	if(!hmult) { exit(0); }
	float nent = (float)hmult->GetEntries();
	if(!nent) {exit(0); }
	TH1D *hMultWork	= (TH1D*)hmult->Clone("hMultWork");
	hMultWork->Scale(1. / nent);
	int nBinsX = hMultWork->GetNbinsX();

	float sumNum = 0;
	float sumDen = 0;
	float result = 0;

	for(int ibx = 1; ibx <= nBinsX; ibx++) {
		float ni = hMultWork->GetBinCenter(ibx);
		float yi = hMultWork->GetBinContent(ibx);		
		sumNum += yi * ni * (ni - 1.);
		sumDen += yi * ni;
	}
	if(sumDen != 0.0) {
		result = sumNum / sumDen / sumDen - 1.;
	}
	delete hMultWork; hMultWork = 0;
	return result;
}

float getC3Baseline(TH1D *hmult){
	if(!hmult) { exit(0); }
	float nent = (float)hmult->GetEntries();
	if(!nent) { exit(0); }
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

void applyC2BaselineAdjustment(TH2D *hR2, float baseC2, int NBETA) {
	for(int ibx = 1; ibx <= NBETA; ibx++) {
		for (int iby = 1; iby <= NBETA; iby++) {
			float oval2	= hR2->GetBinContent(ibx, iby);
			hR2->SetBinContent(ibx, iby, oval2 - baseC2);
		}
	}
}
void applyC3BaselineAdjustment(TH3D *hR3, float baseC3, int NBETA) {
	for(int ibx = 1; ibx <= NBETA; ibx++) {
		for (int iby = 1; iby <= NBETA; iby++) {
			for(int ibz = 1; ibz <= NBETA; ibz++) {
				float oval3	= hR3->GetBinContent(ibx, iby, ibz);
				hR3->SetBinContent(ibx, iby, ibz, oval3 - baseC3);
			}
		}
	}
}

void setStyle() {
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
}
