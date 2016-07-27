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
#include "TChain.h"
#include "TBranch.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
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

void setStyle();

int main(int argc, char **argv) {
	TCanvas *c1 = new TCanvas("c1","graphs",200,10,600,400);
	c1->SetGrid();
	TMultiGraph *mg = new TMultiGraph();
   	mg->SetTitle("y vs #eta");

	TGraph ** g = new TGraph*[3];
	for(int k = 0; k < 3; k++) {
		g[k] = new TGraph();
		g[k]->SetMarkerStyle(20);
	}
	//const int n = 11;

	double pt[3] = {0.4, 1.0, 2.0};
	double eta[11] = {-1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5};
	double E, pz, p, y, theta;
	double m = .9383E+00;

	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 11; j++) {
			theta = 2 * TMath::ATan(TMath::Exp(-eta[j])); 
			cout << theta << endl;
			if(eta[j] < 0) {
				double modifiedEta = std::fabs(eta[j]);
				theta = 2 * TMath::ATan(TMath::Exp(-modifiedEta)); 
				theta = TMath::Pi() - theta;
			}
			int n = g[i]->GetN();
			pz = pt[i] / TMath::Tan(theta);
			cout << pz << endl;
			p = sqrt(pz * pz + pt[i] * pt[i]);
			cout << p << endl;
			E = sqrt(p * p + m * m);
			cout << E << endl;
			y = 0.5 * TMath::Log( (E + pz) / (E - pz) );
			cout << y << endl;
			g[i]->SetPoint(n, eta[j], y);
		}
		cout << endl;
		if(i == 0) g[i]->SetLineColor(kGreen);
		if(i == 1) g[i]->SetLineColor(kBlue);
		if(i == 2) g[i]->SetLineColor(kRed);
		mg->Add(g[i]);
	}

	g[0]->SetLineColor(kGreen);
	g[1]->SetLineColor(kBlue);
	g[2]->SetLineColor(kRed);
    TLegend * legend = new TLegend(0.1, 0.9, 0.28, 0.75);
	legend->SetTextSize(0.03);
    legend->AddEntry(g[0], "p_{T} = 0.4", "l");
    legend->AddEntry(g[1], "p_{T} = 1.0", "l");
    legend->AddEntry(g[2], "p_{T} = 2.0", "l");

	setStyle();
	TString plotFileBase = TString(Form("./ps/etacorr"));
	TString plotFile = plotFileBase + TString(".pdf");
	c1->cd(); 
	mg->Draw("ALP");
	legend->Draw();
	c1->cd(); 
	c1->Update();
	c1->Print(plotFile.Data());

	return 0;
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
