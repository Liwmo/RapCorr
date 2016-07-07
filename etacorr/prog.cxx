#include "rapcorr.h"

RapCorr::RapCorr() {
	nEvents	= 0;
	numBins = 18;
	yLower = -0.72;
	yUpper = 0.72;
	binWidth = (yUpper - yLower) / numBins;
	
	numBinsDY = 2 * numBins - 1;
	yLowerDY = yLower - yUpper + binWidth / 2.;
	yUpperDY = yUpper - yLower - binWidth / 2.;

	maxMult	= 200;
	hMultiplicity = new TH1D("hmult", "hmult", maxMult, -0.5, ((float)maxMult)-0.5);
	
	runR2 = false;
	runR3 = false; 
}

RapCorr::~RapCorr() {
	delete hMultiplicity; hMultiplicity = 0;
	delete hRapidity1D; hRapidity1D = 0;
	delete hRapidity2D; hRapidity2D = 0;
	delete hTensorProduct2D; hTensorProduct2D = 0;
	delete hR2; hR2 = 0;
	delete hConstant2D; hConstant2D = 0;
	delete hR2_dRapidity; hR2_dRapidity = 0;
	delete hR2_dRapidity_N; hR2_dRapidity_N = 0;
	delete hR2_dRapidityBaseline; hR2_dRapidityBaseline = 0;

	TH3D *hRapidity3D;
	TH3D *hR3;
	TH3D *hTensorProduct3D;
	TH3D *hC2rho1;
	TH3D *hConstant3D;
	TH2D *hR3_dRapidity;
	TH2D *hR3_dRapidity_N;
}

void RapCorr::book() {
	TH1::AddDirectory(kFALSE);

	TH1D *hRapidity1D = new TH1D("hRapidity1D", "Rapidity 1D", numBins, -yLower, yUpper);
	TH2D *hRapidity2D = new TH2D("hRapidity2D", "#y_{2} vs #y_{1}", numBins, -yLower, yUpper, numBins, yLower, yUpper);
	TH2D *hTensorProduct2D = new TH2D("hTensorProduct2D", "Tensor Product 2D", numBins, yLower, yUpper, numBins, yLower, yUpper);
	TH2D *hR2 = new TH2D("hR2", "R_{2} vs (#y_{1},#y_{2})", numBins, yLower, yUpper, numBins, yLower, yUpper);
	TH2D *hConstant2D = new TH2D("hConstant2D", "hConstant2D", numBins, yLower, yUpper, numBins, yLower, yUpper);
	TH1D *hR2_dRapidity = new TH1D("hR2_dRapidity", "#LTR_{2}#GT vs #y_{1}-#y_{2}", numBinsDY, yLowerDY, yUpperDY);
	TH1D *hR2_dRapidity_N = new TH1D("hR2_dRapidity_N", "NBINS vs #y_{1}-#y_{2}", numBinsDY, yLowerDY, yUpperDY);
	TH1D *hR2_dRapidityBaseline = new TH1D("hR2_dRapidityBaseline", "#LTR_{2}#GT-Baseline vs #y_{1}-#y_{2}", numBinsDY, yLowerDY, yUpperDY);

	TH3D *hRapidity3D = new TH3D("hRapidity3D", "Rapidity 3D", numBins, yLower, yUpper, numBins, yLower, yUpper, numBins, yLower, yUpper);
	TH3D *hR3 = new TH3D("hR3", "hR3", numBins, yLower, yUpper, numBins, yLower, yUpper, numBins, yLower, yUpper);
	TH3D *hTensorProduct3D = new TH3D("hTensorProduct3D", "Tensor Product 3D", numBins, yLower, yUpper, numBins, yLower, yUpper, numBins, yLower, yUpper);
	TH3D *hC2rho1 = new TH3D("hC2rho1","hC2rho1", numBins, yLower, yUpper, numBins, yLower, yUpper, numBins, yLower, yUpper);
	TH3D *hConstant3D = new TH3D("hConstant3D", "Constant 3D", numBins, yLower, yUpper, numBins, yLower, yUpper, numBins, yLower, yUpper);
	TH2D *hR3_dRapidity = new TH2D("hR3_dRapidity", "#LTR_{3}#GT vs (#Delta#y_{12}, #Delta#y_{13})", numBinsDY, yLowerDY, yUpperDY, numBinsDY, yLowerDY, yUpperDY);
	TH2D *hR3_dRapidity_N	= new TH2D("hR3_dRapidity_N", "NBINS vs (#Delta#y_{12}, #Delta#y_{13})", numBinsDY, yLowerDY, yUpperDY, numBinsDY, yLowerDY, yUpperDY);
}

void RapCorr::increment(double *rapidities, int nTracks) {
	nEvents++;
	if(nTracks >= maxMult) {
		cout << "Warning: Increase maxMult! Mult:" << nTracks << ">" << maxMult << endl; 
	}
	hMultiplicity->Fill(nTracks);
	if(runR2) {
		fill1DRapidityDist(rapidities, nTracks);
		fill2DRapidityDist(rapidities, nTracks);	
	}
	if(runR3) {
		fill3DRapidityDist(rapidities, nTracks);		
	}
}

void RapCorr::calculate() {
	int iCanvas = -1;
	TCanvas *canvases[100];
	TString plotFile0, plotFile, plotFileC, plotFilePDF;
	setupOutputFilePaths(plotFile0, plotFile, plotFileC, plotFilePDF);

	TH1 * histoResets[9] = {hMultiplicity, hRapidity1D, hRapidity2D, hTensorProduct2D, 
		hR2, hConstant2D, hR2_dRapidity, hR2_dRapidity_N, hR2_dRapidityBaseline};
	resetHistograms(histoResets, 9);

	TH1 * histosNormalize[5] = {hRapidity1D, hRapidity2D, hRapidity3D, hR2, hR3};
	normalizeHistograms(histosNormalize, nEvents, 5);

	if(runR2) {
		fill2DTensorProduct();
		fillConstant2DHistogram(-1.0);
		calculateR2Histogram();
		float baseC2 = getC2Baseline(hMultiplicity);
		applyC2BaselineAdjustment(baseC2);
		fillR2dRapidityHistogram();
		fillR2dRapidityBaseHistogram();
		double integral = calculateIntegral(baseC2);
		drawR2HistogramsToFile(canvases, iCanvas, plotFile0, integral);
	}

	if(runR3) {
		fill3DTensorProduct();
		fillConstant3DHistogram(+2.0);
		fillC2rho1Histogram();
		calculateR3Histogram();
		float baseC3 = getC3Baseline(hMultiplicity);
		applyC3BaselineAdjustment(baseC3);
		fillR3dRapidityHistogram();
		drawR3HistogramsToFile(canvases, iCanvas, plotFile);
	}
		
	setStyle();
	executeFilePlots(canvases, iCanvas, plotFileC, plotFile, plotFilePDF);
}

void RapCorr::resetHistograms(TH1 ** histograms, int size) {
	for(int i = 0; i < size; i++) {
		histograms[i]->Reset();
	}
}

void RapCorr::setupOutputFilePaths(TString &plotFile0, TString &plotFile, 
	TString &plotFileC, TString &plotFilePDF) {

	TString plotFileBase, rootOutFile;

	rootOutFile	= TString(Form("./root/etacorr.root"));
	plotFileBase = TString(Form("./ps/etacorr"));
	cout << "root file = " << rootOutFile.Data() << endl;
	cout << "plot file = " << plotFileBase.Data() << endl;
	plotFile0 = plotFileBase + TString(".ps(");
	plotFile = plotFileBase + TString(".ps");
	plotFileC = plotFileBase + TString(".ps]");
	plotFilePDF	= plotFileBase + TString(".pdf");
}

void RapCorr::fill1DRapidityDist(double *rapidityArr, int nTracks) {
	for(int i = 0; i < nTracks; i++) {
		hRapidity1D->Fill(rapidityArr[i]);
	}
}

void RapCorr::fill2DRapidityDist(double *rapidityArr, int nTracks) {
	for(int i = 0; i < nTracks; i++) {			
		for(int j = 0; j < nTracks; j++) {
			if(i != j) {
				hRapidity2D->Fill(rapidityArr[i], rapidityArr[j]);
				hR2->Fill(rapidityArr[i], rapidityArr[j]);
			}
		}
	}
}

void RapCorr::fill3DRapidityDist(double *rapidityArr, int nTracks) {
	for(int i = 0; i < nTracks; i++) {
		for(int j = 0; j < nTracks; j++) {
			for (int k = 0; k < nTracks; k++) {
				if (i != j && i != k && j != k) {
					hRapidity3D->Fill(rapidityArr[i], rapidityArr[j], rapidityArr[k]);
					hR3->Fill(rapidityArr[i], rapidityArr[j], rapidityArr[k]);
				}
			}
		}
	}
}


void RapCorr::normalizeHistograms(TH1 **histos, int normConst, int size) {
	for(int i = 0; i < size; i++) {
		histos[i]->Scale(1. / normConst);
	}
}

void RapCorr::fill2DTensorProduct() {
	int nbin = hRapidity1D->GetNbinsX();
	for(int ibin = 1; ibin <= nbin; ibin++) {
		float valx1	= hRapidity1D->GetBinCenter(ibin);
		float valn1	= hRapidity1D->GetBinContent(ibin);
		for(int jbin = 1; jbin <= nbin; jbin++) {
			float valx2	= hRapidity1D->GetBinCenter(jbin);
			float valn2	= hRapidity1D->GetBinContent(jbin);
			hTensorProduct2D->Fill(valx1, valx2, valn1 * valn2);
		}
	}
}

void RapCorr::fill3DTensorProduct() {
	int nbin = hRapidity1D->GetNbinsX();
	for(int ibin = 1; ibin <= nbin; ibin++) {
		float valx1	= hRapidity1D->GetBinCenter(ibin);
		float valn1	= hRapidity1D->GetBinContent(ibin);
		for(int jbin = 1; jbin <= nbin; jbin++) {
			float valx2	= hRapidity1D->GetBinCenter(jbin);
			float valn2	= hRapidity1D->GetBinContent(jbin);
			for(int kbin = 1; kbin <= nbin; kbin++) {
				float valx3	= hRapidity1D->GetBinCenter(kbin);
				float valn3	= hRapidity1D->GetBinContent(kbin);
				hTensorProduct3D->Fill(valx1, valx2, valx3, valn1 * valn2 * valn3);
			}
		}
	}
}


void RapCorr::fillConstant2DHistogram(float val) {
	int nbin = hRapidity1D->GetNbinsX();
	for(int ibin = 1; ibin <= nbin; ibin++) {
		float valx1	= hRapidity1D->GetBinCenter(ibin);
		for(int jbin = 1; jbin <= nbin; jbin++) {
			float valx2	= hRapidity1D->GetBinCenter(jbin);
			hConstant2D->Fill(valx1, valx2, val);
		}
	}
}

void RapCorr::fillConstant3DHistogram(float val) {
	int nbin = hRapidity1D->GetNbinsX();
	for(int ibin = 1; ibin <= nbin; ibin++) {
		float valx1	= hRapidity1D->GetBinCenter(ibin);
		for(int jbin = 1; jbin <= nbin; jbin++) {
			float valx2 = hRapidity1D->GetBinCenter(jbin);
			for(int kbin = 1; kbin <= nbin; kbin++) {
				float valx3	= hRapidity1D->GetBinCenter(kbin);
				hConstant3D->Fill(valx1, valx2, valx3, val);
			}
		}
	}
}

void RapCorr::fillC2rho1Histogram() {
	int nbin = hRapidity1D->GetNbinsX();
	float valxi, valxj, valxk, valnij, valnk;
	for(int ibin = 1; ibin <= nbin; ibin++) {
		for(int jbin = 1; jbin <= nbin; jbin++) {
			valxi = hRapidity2D->GetXaxis()->GetBinCenter(ibin);
			valxj = hRapidity2D->GetYaxis()->GetBinCenter(jbin);
			valnij = hRapidity2D->GetBinContent(ibin, jbin);

			for(int kbin = 1; kbin <= nbin; kbin++) {
				valxk = hRapidity1D->GetBinCenter(kbin);
				valnk = hRapidity1D->GetBinContent(kbin);
				hC2rho1->Fill(valxi, valxj, valxk, valnij * valnk);
			}
		}
	}
}



void RapCorr::calculateR2Histogram() {
	hR2->Divide(hTensorProduct2D);
	hR2->Add(hConstant2D);
}

void RapCorr::calculateR3Histogram() {
	hR3->Divide(hTensorProduct3D);
	hC2rho1->Divide(hTensorProduct3D);
	hC2rho1->Scale(3.0);
	hR3->Add(hC2rho1, -1.0);
	hR3->Add(hConstant3D);
}

void RapCorr::fillR2dRapidityHistogram() {
	for(int ibx = 1; ibx <= numBins; ibx++) {								
		for(int iby = 1; iby <= numBins; iby++) {						
			float dEtaXY;									
			dEtaXY = hR2->GetXaxis()->GetBinCenter(ibx) 		
			       - hR2->GetYaxis()->GetBinCenter(iby);
			hR2_dRapidity->Fill(dEtaXY, hR2->GetBinContent(ibx, iby));
			hR2_dRapidity_N->Fill(dEtaXY, 1.0);
		}	
	} 
	hR2_dRapidity->Divide(hR2_dRapidity_N);
}

void RapCorr::fillR3dRapidityHistogram() {
	for(int ibx = 1; ibx <= numBins; ibx++) {						
		for(int iby = 1; iby <= numBins; iby++) {						
			for(int ibz = 1; ibz <= numBins; ibz++) {					
				float dEtaXY, dEtaXZ;									
				dEtaXY	= hR3->GetXaxis()->GetBinCenter(ibx) 		
			            - hR3->GetYaxis()->GetBinCenter(iby);
				dEtaXZ	= hR3->GetXaxis()->GetBinCenter(ibx) 		
			        	- hR3->GetZaxis()->GetBinCenter(ibz);
			    hR3_dRapidity->Fill(dEtaXY, dEtaXZ, hR3->GetBinContent(ibx, iby, ibz));
			    hR3_dRapidity_N->Fill(dEtaXY, dEtaXZ, 1.0);
			}
		}	
	}
	hR3_dRapidity->Divide(hR3_dRapidity_N);
}


float RapCorr::getC2Baseline(TH1D *hmult){
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

float RapCorr::getC3Baseline(TH1D *hmult){
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


void RapCorr::applyC2BaselineAdjustment(float baseC2) {
	for(int ibx = 1; ibx <= numBins; ibx++) {
		for (int iby = 1; iby <= numBins; iby++) {
			float oval2	= hR2->GetBinContent(ibx, iby);
			hR2->SetBinContent(ibx, iby, oval2 - baseC2);
		}
	}
}

void RapCorr::applyC3BaselineAdjustment(float baseC3) {
	for(int ibx = 1; ibx <= numBins; ibx++) {
		for (int iby = 1; iby <= numBins; iby++) {
			for(int ibz = 1; ibz <= numBins; ibz++) {
				float oval3	= hR3->GetBinContent(ibx, iby, ibz);
				hR3->SetBinContent(ibx, iby, ibz, oval3 - baseC3);
			}
		}
	}
}


void RapCorr::fillR2dRapidityBaseHistogram() { 
	for(int ibin = 1; ibin <= hR2_dRapidity->GetNbinsX(); ibin++) {
		float val = hR2_dRapidity->GetBinContent(ibin);
		float valErr = hR2_dRapidity->GetBinError(ibin);	// INCORRECT ERRORS VALS
		hR2_dRapidityBaseline->SetBinContent(ibin, val);	// straight copy now...
		hR2_dRapidityBaseline->SetBinError(ibin, valErr);	// INCORRECT ERRORS VALS.
	}
}

double RapCorr::calculateIntegral(float baseline) {
	double integral = 0;
	int numBins = hR2_dRapidityBaseline->GetNbinsX();
	for(int i = 1; i < numBins; i++) {
		double value = hR2_dRapidityBaseline->GetBinContent(i);
		integral += (value - baseline);
	}
	return integral;
}


void RapCorr::setStyle() {
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

void RapCorr::drawR2HistogramsToFile(TCanvas **canvases, int &iCanvas, TString plotFile0, double integral) {
		TLatex * text = new TLatex();
		text->SetTextSize(0.05);
		text->SetNDC();
		char buf[200];
		iCanvas++;
		sprintf(buf, "canvases%d", iCanvas);
		canvases[iCanvas] = new TCanvas(buf, buf, 30 * iCanvas, 30 * iCanvas, 800, (8.5 / 11.) * 800);
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Divide(3,2,0.0001,0.0001);
			canvases[iCanvas]->cd(1);
				hMultiplicity->Draw();
			canvases[iCanvas]->cd(2);
				hRapidity1D->SetMinimum(0.5);
				hRapidity1D->Draw();
			canvases[iCanvas]->cd(3);
				hRapidity2D->SetStats(0);
				hRapidity2D->Draw("colz");
			canvases[iCanvas]->cd(4);
				hR2->SetStats(0);
				hR2->Draw("colz");
			canvases[iCanvas]->cd(5);
				hR2_dRapidityBaseline->SetStats(0);
				hR2_dRapidityBaseline->SetMinimum(-0.005);
				hR2_dRapidityBaseline->SetMaximum(0.005);
				hR2_dRapidityBaseline->SetMarkerStyle(20);
				hR2_dRapidityBaseline->SetMarkerSize(1);
				hR2_dRapidityBaseline->SetMarkerColor(4);
				hR2_dRapidityBaseline->SetLineColor(4);
				hR2_dRapidityBaseline->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral));
			canvases[iCanvas]->cd(6);
				hR2_dRapidityBaseline->Draw();
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Update();
		canvases[iCanvas]->Print(plotFile0.Data());
}

void RapCorr::drawR3HistogramsToFile(TCanvas **canvases, int &iCanvas, TString plotFile) {
		char buf[200];
		iCanvas++;	
		gStyle->SetOptStat(0);
		gStyle->SetPadRightMargin(0.15);
		iCanvas++;
		sprintf(buf, "canvases%d", iCanvas);
		canvases[iCanvas] = new TCanvas(buf, buf, 30 * iCanvas, 30 * iCanvas, 800, (8.5 / 11.) * 800);
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Divide(3, 2, 0.0001, 0.0001);
			canvases[iCanvas]->cd(1);
				hR2_dRapidityBaseline->Draw("HIST");
			canvases[iCanvas]->cd(2);
				hR3_dRapidity->Draw("COLZ1");
			canvases[iCanvas]->cd(3);
				hR3_dRapidity_N->Draw("COLZ1");
			canvases[iCanvas]->cd(4);
			canvases[iCanvas]->cd(5);
			canvases[iCanvas]->cd(6);
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Update();
		canvases[iCanvas]->Print(plotFile.Data());
		gStyle->SetPadRightMargin(0.06);
}

void RapCorr::executeFilePlots(TCanvas **canvases, int iCanvas, TString plotFileC, 
	TString plotFile, TString plotFilePDF) {
		char buf[200];
	 	canvases[iCanvas]->Print(plotFileC.Data());
	 	sprintf(buf, "/usr/bin/pstopdf %s -o %s", plotFile.Data(), plotFilePDF.Data());
	 	cout << " " << buf << endl;
	 	gSystem->Exec(buf);
}
