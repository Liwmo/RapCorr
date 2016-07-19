#include "rapcorr.h"

RapCorr::RapCorr(int nb, float low, float high) {
	nEvents	= 0;
	numBins = nb;
	yLower = low;
	yUpper = high;
	binWidth = (yUpper - yLower) / numBins;
	
	numBinsDY = 2 * numBins - 1;
	yLowerDY = yLower - yUpper + binWidth / 2.;
	yUpperDY = yUpper - yLower - binWidth / 2.;

	maxMult	= 200;
	hMultiplicity = new TH1D("hmult", "Multiplicity", maxMult, -0.5, ((float)maxMult)-0.5);
	
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
	delete hR2_dRapidity_Error; hR2_dRapidity_Error = 0;

	delete hRapidity3D; hRapidity3D = 0;
	delete hR3; hR3 = 0;
	delete hTensorProduct3D; hTensorProduct3D = 0;
	delete hC2rho1; hC2rho1 = 0;
	delete hConstant3D; hConstant3D = 0;
	delete hR3_dRapidity; hR3_dRapidity = 0;
	delete hR3_dRapidity_N; hR3_dRapidity_N = 0;
}

void RapCorr::book() {
	TH1::AddDirectory(kFALSE);
	hRapidity1D = new TH1D("hRapidity1D", "Rapidity 1D", numBins, yLower, yUpper);
	hRapidity2D = new TH2D("hRapidity2D", "y_{2} vs y_{1}", numBins, yLower, yUpper, numBins, yLower, yUpper);
	hTensorProduct2D = new TH2D("hTensorProduct2D", "Tensor Product 2D", numBins, yLower, yUpper, numBins, yLower, yUpper);
	hR2 = new TH2D("hR2", "R_{2} vs (y_{1},y_{2})", numBins, yLower, yUpper, numBins, yLower, yUpper);
	hConstant2D = new TH2D("hConstant2D", "hConstant2D", numBins, yLower, yUpper, numBins, yLower, yUpper);
	hR2_dRapidity = new TH1D("hR2_dRapidity", "#LTR_{2}#GT vs y_{1}-y_{2}", numBinsDY, yLowerDY, yUpperDY);
	hR2_dRapidity_N = new TH1D("hR2_dRapidity_N", "NBINS vs y_{1}-y_{2}", numBinsDY, yLowerDY, yUpperDY);
	hR2_dRapidity_Error = new TH1D("hR2_dRapidity_Error", "#LTR_{2}#GT Error", numBinsDY, yLowerDY, yUpperDY);

	hRapidity3D = new TH3D("hRapidity3D", "Rapidity 3D", numBins, yLower, yUpper, numBins, yLower, yUpper, numBins, yLower, yUpper);
	hR3 = new TH3D("hR3", "hR3", numBins, yLower, yUpper, numBins, yLower, yUpper, numBins, yLower, yUpper);
	hTensorProduct3D = new TH3D("hTensorProduct3D", "Tensor Product 3D", numBins, yLower, yUpper, numBins, yLower, yUpper, numBins, yLower, yUpper);
	hC2rho1 = new TH3D("hC2rho1","hC2rho1", numBins, yLower, yUpper, numBins, yLower, yUpper, numBins, yLower, yUpper);
	hConstant3D = new TH3D("hConstant3D", "Constant 3D", numBins, yLower, yUpper, numBins, yLower, yUpper, numBins, yLower, yUpper);
	hR3_dRapidity = new TH2D("hR3_dRapidity", "#LTR_{3}#GT vs (#Delta#y_{12}, #Delta#y_{13})", numBinsDY, yLowerDY, yUpperDY, numBinsDY, yLowerDY, yUpperDY);
	hR3_dRapidity_N	= new TH2D("hR3_dRapidity_N", "NBINS vs (#Delta#y_{12}, #Delta#y_{13})", numBinsDY, yLowerDY, yUpperDY, numBinsDY, yLowerDY, yUpperDY);
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

	TH1 * histosNormalize[5] = {hRapidity1D, hRapidity2D, hRapidity3D, hR2, hR3};
	normalizeHistograms(histosNormalize, nEvents, 5);

	if(runR2) {
		fill2DTensorProduct();
		fillConstant2DHistogram(-1.0);
		calculateR2Histogram();
		float baseC2 = getC2Baseline(hMultiplicity);
		applyC2BaselineAdjustment(baseC2);
		fillR2dRapidityHistogram();
		integral = calculateIntegral(baseC2);
	}

	if(runR3) {
		fill3DTensorProduct();
		fillConstant3DHistogram(+2.0);
		fillC2rho1Histogram();
		calculateR3Histogram();
		float baseC3 = getC3Baseline(hMultiplicity);
		applyC3BaselineAdjustment(baseC3);
		fillR3dRapidityHistogram();
	}
}

void RapCorr::fill1DRapidityDist(double *rapidities, int nTracks) {
	for(int i = 0; i < nTracks; i++) {
		hRapidity1D->Fill(rapidities[i]);
	}
}

void RapCorr::fill2DRapidityDist(double *rapidities, int nTracks) {
	for(int i = 0; i < nTracks; i++) {			
		for(int j = 0; j < nTracks; j++) {
			if(i != j) {
				hRapidity2D->Fill(rapidities[i], rapidities[j]);
				hR2->Fill(rapidities[i], rapidities[j]);
			}
		}
	}
}

void RapCorr::fill3DRapidityDist(double *rapidities, int nTracks) {
	for(int i = 0; i < nTracks; i++) {
		for(int j = 0; j < nTracks; j++) {
			for (int k = 0; k < nTracks; k++) {
				if (i != j && i != k && j != k) {
					hRapidity3D->Fill(rapidities[i], rapidities[j], rapidities[k]);
					hR3->Fill(rapidities[i], rapidities[j], rapidities[k]);
				}
			}
		}
	}
}


void RapCorr::normalizeHistograms(TH1 **histos, int normConst, int size) {
	for(int i = 0; i < size; i++) {
		histos[i]->Sumw2();
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
			double val = hR2->GetBinContent(ibx, iby);
			dEtaXY = hR2->GetXaxis()->GetBinCenter(ibx) 		
			       - hR2->GetYaxis()->GetBinCenter(iby);
			hR2_dRapidity->Fill(dEtaXY, val);
			hR2_dRapidity_Error->Fill(dEtaXY, val * val);
			hR2_dRapidity_N->Fill(dEtaXY, 1.0);
		}	
	} 
	hR2_dRapidity->Divide(hR2_dRapidity_N);
	hR2_dRapidity_Error->Divide(hR2_dRapidity_N);
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
	if(!hmult) { cout << "RapCorr::getC2Baseline - no hmult - exiting" << endl; exit(0); }
	float nent = (float)hmult->GetEntries();
	if(!nent) { cout << "RapCorr::getC2Baseline - no entries - exiting" << endl; exit(0); }
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

double RapCorr::calculateIntegral(float baseline) {
	double integral = 0;
	int numBins = hR2_dRapidity->GetNbinsX();
	for(int i = 1; i < numBins; i++) {
		double value = hR2_dRapidity->GetBinContent(i);
		integral += value;
	}
	return integral;
}
