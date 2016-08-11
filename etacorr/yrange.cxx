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
	TMultiGraph *mgProton = new TMultiGraph();
   	mgProton->SetTitle("y_{p} vs #eta_{p} (Proton)");

   	TMultiGraph *mgKaon = new TMultiGraph();
   	mgKaon->SetTitle("y_{p} vs #eta_{p} (Kaon)");

	TMultiGraph *mgPion = new TMultiGraph();
   	mgPion->SetTitle("y_{p} vs #eta_{p} (Pion)");

	TGraph ** g = new TGraph*[9];
	for(int k = 0; k < 9; k++) {
		g[k] = new TGraph();
		g[k]->SetMarkerStyle(20);
	}
	const int num = 16;

	double pt[3] = {0.4, 1.0, 2.0};
	double eta[num] = {-1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5};
	double E, pz, p, y, theta;
	double mProton = .9383E+00;
	double mKaon = .4977E+00;
	double mPion = .1350E+00;
	double mass[3] = {mProton, mKaon, mPion};

	int count = 0;
	for(int m = 0; m < 3; m++) {
		for(int i = 0; i < 3; i++) {
			for(int j = 0; j < num; j++) {
				theta = 2 * TMath::ATan(TMath::Exp(-eta[j])); 
				if(eta[j] < 0) {
					double modifiedEta = std::fabs(eta[j]);
					theta = 2 * TMath::ATan(TMath::Exp(-modifiedEta)); 
					theta = TMath::Pi() - theta;
				}
				int n = g[count]->GetN();
				pz = pt[i] / TMath::Tan(theta);
				p = sqrt(pz * pz + pt[i] * pt[i]);
				E = sqrt(p * p + mass[m] * mass[m]);
				y = 0.5 * TMath::Log( (E + pz) / (E - pz) );
				g[count]->SetPoint(n, eta[j], y);
			}
			if(i == 0) g[count]->SetLineColor(kGreen);
			if(i == 1) g[count]->SetLineColor(kBlue);
			if(i == 2) g[count]->SetLineColor(kRed);
			if(m == 0) mgProton->Add(g[count]);
			if(m == 1) mgKaon->Add(g[count]);
			if(m == 2) mgPion->Add(g[count]);
			count++;
		}
	}

	double yLo_proton_pt1 = g[0]->Eval(-1.094);
	double yHi_proton_pt1 = g[0]->Eval(1);
	double yLo_proton_pt2 = g[1]->Eval(-1.094);
	double yHi_proton_pt2 = g[1]->Eval(1);
	double yLo_proton_pt3 = g[2]->Eval(-1.094);
	double yHi_proton_pt3 = g[2]->Eval(1);

	double yLo_kaon_pt1 = g[3]->Eval(-1.094);
	double yHi_kaon_pt1 = g[3]->Eval(1);
	double yLo_kaon_pt2 = g[4]->Eval(-1.094);
	double yHi_kaon_pt2 = g[4]->Eval(1);
	double yLo_kaon_pt3 = g[5]->Eval(-1.094);
	double yHi_kaon_pt3 = g[5]->Eval(1);

	double yLo_pion_pt1 = g[6]->Eval(-1.094);
	double yHi_pion_pt1 = g[6]->Eval(1);
	double yLo_pion_pt2 = g[7]->Eval(-1.094);
	double yHi_pion_pt2 = g[7]->Eval(1);
	double yLo_pion_pt3 = g[8]->Eval(-1.094);
	double yHi_pion_pt3 = g[8]->Eval(1);
        /*-1.094 to 0.899*/
	cout << "Proton" << endl;
	cout << "pt: 0.4 || " << yLo_proton_pt1 << " < y < " << yHi_proton_pt1 << endl;
	cout << "pt: 1.0 || " << yLo_proton_pt2 << " < y < " << yHi_proton_pt2 << endl;
	cout << "pt: 2.0 || " << yLo_proton_pt3 << " < y < " << yHi_proton_pt3 << endl;

	cout << "Kaon" << endl;
	cout << "pt: 0.4 || " << yLo_kaon_pt1 << " < y < " << yHi_kaon_pt1 << endl;
	cout << "pt: 1.0 || " << yLo_kaon_pt2 << " < y < " << yHi_kaon_pt2 << endl;
	cout << "pt: 2.0 || " << yLo_kaon_pt3 << " < y < " << yHi_kaon_pt3 << endl;

	cout << "Pion" << endl;
	cout << "pt: 0.4 || " << yLo_pion_pt1 << " < y < " << yHi_pion_pt1 << endl;
	cout << "pt: 1.0 || " << yLo_pion_pt2 << " < y < " << yHi_pion_pt2 << endl;
	cout << "pt: 2.0 || " << yLo_pion_pt3 << " < y < " << yHi_pion_pt3 << endl;

	setStyle();
	TString plotFileBase = TString(Form("./ps/etacorr"));
	TString plotFile0 = plotFileBase + TString(".pdf(");
	TString plotFile = plotFileBase + TString(".pdf");
	TString plotFileC = plotFileBase + TString(".pdf)");
	
	TCanvas *c1 = new TCanvas("c1","graphs",200,10,600,400);
    TLegend * legend1 = new TLegend(0.1, 0.9, 0.28, 0.75);
	legend1->SetTextSize(0.03);
    legend1->AddEntry(g[0], "p_{T} = 0.4", "l");
    legend1->AddEntry(g[1], "p_{T} = 1.0", "l");
    legend1->AddEntry(g[2], "p_{T} = 2.0", "l");
	c1->SetGrid();
	c1->cd(); 
	mgProton->Draw("ALP");
	legend1->Draw();
	c1->cd(); 
	c1->Update();
	c1->Print(plotFile0.Data());

	TCanvas *c2 = new TCanvas("c1","graphs",200,10,600,400);
    TLegend * legend2 = new TLegend(0.1, 0.9, 0.28, 0.75);
	legend2->SetTextSize(0.03);
    legend2->AddEntry(g[3], "p_{T} = 0.4", "l");
    legend2->AddEntry(g[4], "p_{T} = 1.0", "l");
    legend2->AddEntry(g[5], "p_{T} = 2.0", "l");
	c2->SetGrid();
	c2->cd(); 
	mgKaon->Draw("ALP");
	legend2->Draw();
	c2->cd(); 
	c2->Update();
	c2->Print(plotFile.Data());

	TCanvas *c3 = new TCanvas("c1","graphs",200,10,600,400);
    TLegend * legend3 = new TLegend(0.1, 0.9, 0.28, 0.75);
	legend3->SetTextSize(0.03);
    legend3->AddEntry(g[6], "p_{T} = 0.4", "l");
    legend3->AddEntry(g[7], "p_{T} = 1.0", "l");
    legend3->AddEntry(g[8], "p_{T} = 2.0", "l");
	c3->SetGrid();
	c3->cd(); 
	mgPion->Draw("ALP");
	legend3->Draw();
	c3->cd(); 
	c3->Update();
	c3->Print(plotFileC.Data());

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
