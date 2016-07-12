#ifndef __RapCorr_h__
#define __RapCorr_h__

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

struct TrackInfo {
	int *geantID;
	float *px, *py, *pz;
	TString name;
	float mass, charge, lifetime;
    float eta, rapidity, phi;
    float pTotal, pt, baryonNo;
};

class RapCorr{
	public:
		RapCorr();		
		RapCorr(int, float, float);		
		virtual ~RapCorr();	
		virtual void book();	
		virtual void setupOutputFilePaths(TString&, TString&, TString&, TString&);
		virtual void calculate();
		virtual void increment(double*, int);
	
		virtual void setRunR2(bool);
		virtual void setRunR3(bool);
		virtual bool getRunR2();
		virtual bool getRunR3();
		virtual void setNumBins(int i);
		virtual void setYLower(float y);
		virtual void setYUpper(float y);
		virtual void setMaxMult(int i);
		virtual TH1D* getMultiplicity();
		virtual TH1D* getRapidity1D();
		virtual TH2D* getRapidity2D();
		virtual TH2D* getTensorProduct2D();
		virtual TH2D* getR2();
		virtual TH1D* getR2dRapidity();
		virtual TH3D* getRapidity3D();
		virtual TH3D* getTensorProduct3D();
		virtual TH3D* getR3();
		virtual TH2D* getR3dRapidity();
		virtual TH2D* getR3dRapidityN();
		virtual double getIntegral();

	private:
		int	maxMult;
		int numBins;
		float yLower;
		float yUpper;
		float binWidth;
		int numBinsDY;
		float yLowerDY;
		float yUpperDY;
		bool runR2;
		bool runR3;
		int nEvents;
		double integral;

		TH1D *hMultiplicity;
		TH1D *hRapidity1D;
		TH2D *hRapidity2D;
		TH2D *hTensorProduct2D;
		TH2D *hR2;
		TH2D *hConstant2D;
		TH1D *hR2_dRapidity;
		TH1D *hR2_dRapidity_N;

		TH3D *hRapidity3D;
		TH3D *hR3;
		TH3D *hTensorProduct3D;
		TH3D *hC2rho1;
		TH3D *hConstant3D;
		TH2D *hR3_dRapidity;
		TH2D *hR3_dRapidity_N;

		virtual void fill1DRapidityDist(double*, int);
		virtual void fill2DRapidityDist(double*, int);
		virtual void fill3DRapidityDist(double*, int);
		virtual void normalizeHistograms(TH1**, int, int);
		virtual void fill2DTensorProduct();
		virtual void fill3DTensorProduct();
		virtual void fillConstant2DHistogram(float);
		virtual void fillConstant3DHistogram(float);
		virtual void fillC2rho1Histogram();
		virtual void calculateR2Histogram();
		virtual void calculateR3Histogram();
		virtual float getC2Baseline(TH1D*);
		virtual float getC3Baseline(TH1D*);
		virtual void applyC2BaselineAdjustment(float);
		virtual void applyC3BaselineAdjustment(float);
		virtual void fillR2dRapidityHistogram();
		virtual void fillR3dRapidityHistogram();
		virtual double calculateIntegral(float);
		virtual void setStyle();
		virtual void drawR2HistogramsToFile(TCanvas**, int&, TString, double);
		virtual void drawR3HistogramsToFile(TCanvas**, int&, TString);
		virtual void executeFilePlots(TCanvas**, int, TString, TString, TString);
};

inline void RapCorr::setRunR2(bool b){ runR2=b; }
inline void RapCorr::setRunR3(bool b){ runR3=b; }
inline bool RapCorr::getRunR2(){ return runR2; }
inline bool RapCorr::getRunR3(){ return runR3; }
inline void RapCorr::setNumBins(int i){ numBins = i; cout<<"Changed Nbins to  "<<numBins<<endl; }
inline void RapCorr::setYLower(float y){ yLower = y; cout<<"Changed Ylower to "<<yLower<<endl; }
inline void RapCorr::setYUpper(float y){ yUpper = y; cout<<"Changed Yupper to "<<yUpper<<endl; }
inline void RapCorr::setMaxMult(int i){ maxMult=i; cout<<"Changed MaxMult to "<<maxMult<<endl; }
inline TH1D* RapCorr::getMultiplicity(){ return hMultiplicity; }
inline TH1D* RapCorr::getRapidity1D(){ return hRapidity1D; }
inline TH2D* RapCorr::getRapidity2D(){ return hRapidity2D; }
inline TH2D* RapCorr::getTensorProduct2D(){ return hTensorProduct2D; }
inline TH2D* RapCorr::getR2(){ return hR2; }
inline TH1D* RapCorr::getR2dRapidity(){ return hR2_dRapidity; }
inline TH3D* RapCorr::getRapidity3D(){ return hRapidity3D; }
inline TH3D* RapCorr::getTensorProduct3D(){ return hTensorProduct3D; }
inline TH3D* RapCorr::getR3(){ return hR3; }
inline TH2D* RapCorr::getR3dRapidity(){ return hR3_dRapidity; }
inline TH2D* RapCorr::getR3dRapidityN(){ return hR3_dRapidity_N; }
inline double RapCorr::getIntegral(){ return integral; }


#endif
