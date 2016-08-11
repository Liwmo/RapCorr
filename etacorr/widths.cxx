#include "rapcorr.h"

void plotRapidityCorrelations();
int addUrQMDEventsFromPath(TChain*);
double * fillRapidities(TrackInfo&, int, int&, double); 
void getInfo(int, float, float, float, TString&, float&, float&, float&,
		    	float&, float&, float&, float&, float&, float&);
void setupOutputFilePaths(TString&, TString&, TString&);
void setStyle();
void drawR2HistogramsToFile(TCanvas**, int&, TString, TH1**, double);
void drawR2WindowWidthsToFile(TCanvas**, int&, TString, TString, TH1D**, double*);
void executeFilePlots(TCanvas**, int, TString);

int main(int argc, char **argv) {
	gRandom->SetSeed(123456);
	plotRapidityCorrelations();
	return 0;
}

void plotRapidityCorrelations() {
	float mass, charge, lifetime, eta, rapidity, phi, pTotal, pt, baryonNo; 
	TString name;

	Int_t ievt;
	Float_t parimp;
	Float_t gecm;
	Float_t evinfo[5];
	Int_t mpart;
	Int_t igid[11000];   
	Float_t gpx[11000];    
	Float_t gpy[11000];    
	Float_t gpz[11000];    
	
	TBranch *b_ievt;     
	TBranch *b_parimp;   
	TBranch *b_gecm;     
	TBranch *b_evinfo;   
	TBranch *b_mpart;    
	TBranch *b_igid;     
	TBranch *b_gpx;      
	TBranch *b_gpy;      
	TBranch *b_gpz;      

	TChain *chain = new TChain("h1", "events chain");
	int nb, nentries = addUrQMDEventsFromPath(chain);

	chain->SetMakeClass(1);
	chain->SetBranchAddress("ievt", &ievt, &b_ievt);
	chain->SetBranchAddress("parimp", &parimp, &b_parimp);
	chain->SetBranchAddress("gecm", &gecm, &b_gecm);
	chain->SetBranchAddress("evinfo", evinfo, &b_evinfo);
	chain->SetBranchAddress("mpart", &mpart, &b_mpart);
	chain->SetBranchAddress("igid", igid, &b_igid);
	chain->SetBranchAddress("gpx", gpx, &b_gpx);
	chain->SetBranchAddress("gpy", gpy, &b_gpy);
	chain->SetBranchAddress("gpz", gpz, &b_gpz);

	const int N_COMPARE = 9;
	RapCorr ** rapCorrWidths = new RapCorr*[N_COMPARE];	
	const int bins[N_COMPARE] = {8, 20, 28, 56, 88, 116, 144, 200, 320};
	const float widths[N_COMPARE] = {0.10, 0.25, 0.35, 0.70, 1.10, 1.45, 1.80, 2.50, 4.00};
	for(int i = 0; i < N_COMPARE; i++) {
		rapCorrWidths[i] = new RapCorr(bins[i], -widths[i], widths[i]);
		rapCorrWidths[i]->setMaxMult(210);
		rapCorrWidths[i]->setRunR2(true);
		rapCorrWidths[i]->book();
	}
	int N_EVENTS = 200000;
	if(N_EVENTS > nentries) {
		N_EVENTS = nentries;
	}
	TrackInfo trackInfo = {igid, gpx, gpy, gpz, name, mass, charge, 
	lifetime, eta, rapidity, phi, pTotal, pt, baryonNo};

	const int MAX_MULT = 210;
	const int PROTON = 14;
	double ** rapidityWidthsArr = new double*[N_COMPARE];
    for(int i = 0; i < N_COMPARE; i++) {
    	rapidityWidthsArr[i] = new double[MAX_MULT];
   		memset(rapidityWidthsArr[i], 0, MAX_MULT * sizeof(double)); 
    }

	for(Long64_t iEvent = 0; iEvent < N_EVENTS; iEvent++) {
		if(iEvent % (N_EVENTS / 10) == 0) { 
			cout << "processing " << iEvent << " of " << N_EVENTS << endl; 
		}
               
		nb = chain->GetEntry(iEvent);
		int N_TRACKS = mpart;
	    int N_PROTON_TRACKS = 0;

		if(parimp > 3.2) {
			continue; // allow only 0-5% central collisions...  
		} 

		for(int i = 0; i < N_COMPARE; i++) {
			int protonCount = 0; 

			for(int iTrack = 0; iTrack < N_TRACKS; iTrack++) {
				getInfo(trackInfo.geantID[iTrack], trackInfo.px[iTrack], trackInfo.py[iTrack], trackInfo.pz[iTrack], 
					trackInfo.name, trackInfo.mass, trackInfo.charge, trackInfo.lifetime, trackInfo.eta, trackInfo.rapidity, trackInfo.phi, trackInfo.pTotal, trackInfo.pt, trackInfo.baryonNo);
				if(trackInfo.geantID[iTrack] == PROTON && fabs(trackInfo.rapidity) <= widths[i]) { 
					rapidityWidthsArr[i][protonCount] = trackInfo.rapidity;
					protonCount++;
				}
			}
				N_PROTON_TRACKS = protonCount;
				rapCorrWidths[i]->increment(rapidityWidthsArr[i], N_PROTON_TRACKS);
		}

	}

	double * widthIntegrals = new double[N_COMPARE];
	for(int i = 0; i < N_COMPARE; i++) {
		rapCorrWidths[i]->calculate();
		widthIntegrals[i] = rapCorrWidths[i]->getIntegral();
	}

	TH1D ** widthHistograms = new TH1D*[N_COMPARE];
	for(int i = 0; i < N_COMPARE; i++) {
	    widthHistograms[i] = rapCorrWidths[i]->getR2dRapidity();
	}

	int iCanvas = -1;
	TCanvas *canvases[100];
	TString plotFile0, plotFile, plotFileC;
	setupOutputFilePaths(plotFile0, plotFile, plotFileC);
	drawR2WindowWidthsToFile(canvases, iCanvas, plotFile0, plotFile, widthHistograms, widthIntegrals);
	for(int i = 0; i < N_COMPARE; i++) {
		TH1 **rapidityHistograms = new TH1*[6];
		rapidityHistograms[0] = (TH1D*) rapCorrWidths[i]->getMultiplicity();
		rapidityHistograms[1] = (TH1D*) rapCorrWidths[i]->getRapidity1D();
		rapidityHistograms[2] = (TH2D*) rapCorrWidths[i]->getRapidity2D();
		rapidityHistograms[3] = (TH2D*) rapCorrWidths[i]->getTensorProduct2D();
		rapidityHistograms[4] = (TH2D*) rapCorrWidths[i]->getR2(); 
		rapidityHistograms[5] = (TH1D*) rapCorrWidths[i]->getR2dRapidity();	
		drawR2HistogramsToFile(canvases, iCanvas, plotFile, rapidityHistograms, widthIntegrals[i]);
	}
	executeFilePlots(canvases, iCanvas, plotFileC);
}

int addUrQMDEventsFromPath(TChain *chain) {
    TString path = TString("/nfs/rhi/UrQMD/events_2016/007/");
	TString	filenames = TString("urqmd_19_00*.root");
	TString input = path + filenames;
	cout << input.Data() << endl;
	chain->Add(input.Data());
	int neventtree = chain->GetEntries();
	cout << "N_events = " << neventtree << endl;
	return neventtree;
}

double * fillRapidities(TrackInfo &info, int N_TRACKS, int &N_PROTON_TRACKS, double window) {
	double *rapidityArr = new double[N_TRACKS];
	int protonCount = 0; 
	const int PROTON = 14;

	for(int iTrack = 0; iTrack < N_TRACKS; iTrack++) {
		getInfo(info.geantID[iTrack], info.px[iTrack], info.py[iTrack], info.pz[iTrack], 
			info.name, info.mass, info.charge, info.lifetime, info.eta, info.rapidity, info.phi, info.pTotal, info.pt, info.baryonNo);
		if(info.geantID[iTrack] == PROTON && fabs(info.rapidity) <= window) { 
			rapidityArr[protonCount] = info.rapidity;
			protonCount++;
		}
	}
	N_PROTON_TRACKS = protonCount;
	return rapidityArr;
}

void getInfo(int geantID, float px, float py, float pz, 
        TString &name, float &mass, float &charge, float &lifetime,
        float &eta, float &rapidity, float &phi, 
        float &pTotal, float &pt, float &baryonNo) {

        mass = charge = lifetime = -9; 
        if(geantID == 1)      {name = TString("GAMMA");         mass = .0000E+00; charge =  0.; lifetime = .10000E+16; baryonNo =  0;}
        else if(geantID == 2) {name = TString("POSITRON");      mass = .5110E-03; charge =  1.; lifetime = .10000E+16; baryonNo =  0;}
        else if(geantID == 3) {name = TString("ELECTRON");      mass = .5110E-03; charge = -1.; lifetime = .10000E+16; baryonNo =  0;}
        else if(geantID == 4) {name = TString("NEUTRINO");      mass = .0000E+00; charge =  0.; lifetime = .10000E+16; baryonNo =  0;}
        else if(geantID == 5) {name = TString("MUON+");         mass = .1057E+00; charge =  1.; lifetime = .21970E-05; baryonNo =  0;}
        else if(geantID == 6) {name = TString("MUON-");         mass = .1057E+00; charge = -1.; lifetime = .21970E-05; baryonNo =  0;}
        else if(geantID == 7) {name = TString("PION0");         mass = .1350E+00; charge =  0.; lifetime = .84000E-16; baryonNo =  0;}
        else if(geantID == 8) {name = TString("PION+");         mass = .1396E+00; charge =  1.; lifetime = .26030E-07; baryonNo =  0;}
        else if(geantID == 9) {name = TString("PION-");         mass = .1396E+00; charge = -1.; lifetime = .26030E-07; baryonNo =  0;}
        else if(geantID == 10) {name = TString("KAON0LONG");    mass = .4977E+00; charge =  0.; lifetime = .51700E-07; baryonNo =  0;}
        else if(geantID == 11) {name = TString("KAON+");        mass = .4937E+00; charge =  1.; lifetime = .12370E-07; baryonNo =  0;}
        else if(geantID == 12) {name = TString("KAON-");        mass = .4937E+00; charge = -1.; lifetime = .12370E-07; baryonNo =  0;}
        else if(geantID == 13) {name = TString("NEUTRON");      mass = .9396E+00; charge =  0.; lifetime = .88700E+03; baryonNo =  1;}
        else if(geantID == 14) {name = TString("PROTON");       mass = .9383E+00; charge =  1.; lifetime = .10000E+16; baryonNo =  1;}
        else if(geantID == 15) {name = TString("ANTIPROTON");   mass = .9383E+00; charge = -1.; lifetime = .10000E+16; baryonNo = -1;}
        else if(geantID == 16) {name = TString("KAON 0 SHORT"); mass = .4977E+00; charge =  0.; lifetime = .89260E-10; baryonNo =  0;}
        else if(geantID == 17) {name = TString("ETA");          mass = .5475E+00; charge =  0.; lifetime = .54850E-18; baryonNo =  0;}
        else if(geantID == 18) {name = TString("LAMBDA");       mass = .1116E+01; charge =  0.; lifetime = .26320E-09; baryonNo =  1;}
        else if(geantID == 19) {name = TString("SIGMA+");       mass = .1189E+01; charge =  1.; lifetime = .79900E-10; baryonNo =  1;}
        else if(geantID == 20) {name = TString("SIGMA0");       mass = .1193E+01; charge =  0.; lifetime = .74000E-19; baryonNo =  1;}
        else if(geantID == 21) {name = TString("SIGMA-");       mass = .1197E+01; charge = -1.; lifetime = .14790E-09; baryonNo =  1;}
        else if(geantID == 22) {name = TString("XI0");          mass = .1315E+01; charge =  0.; lifetime = .29000E-09; baryonNo =  1;}
        else if(geantID == 23) {name = TString("XI-");          mass = .1321E+01; charge = -1.; lifetime = .16390E-09; baryonNo =  1;}
        else if(geantID == 24) {name = TString("OMEGA-");       mass = .1672E+01; charge = -1.; lifetime = .82200E-10; baryonNo =  1;}
        else if(geantID == 25) {name = TString("ANTINEUTRON");  mass = .9396E+00; charge =  0.; lifetime = .88700E+03; baryonNo = -1;}
        else if(geantID == 26) {name = TString("ANTILAMBDA");   mass = .1116E+01; charge =  0.; lifetime = .26320E-09; baryonNo = -1;}
        else if(geantID == 27) {name = TString("ANTISIGMA-");   mass = .1189E+01; charge = -1.; lifetime = .79900E-10; baryonNo = -1;}
        else if(geantID == 28) {name = TString("ANTISIGMA0");   mass = .1193E+01; charge =  0.; lifetime = .74000E-19; baryonNo = -1;}
        else if(geantID == 29) {name = TString("ANTISIGMA+");   mass = .1197E+01; charge =  1.; lifetime = .14790E-09; baryonNo = -1;}
        else if(geantID == 30) {name = TString("ANTIXI0");      mass = .1315E+01; charge =  0.; lifetime = .29000E-09; baryonNo = -1;}
        else if(geantID == 31) {name = TString("ANTIXI+");      mass = .1321E+01; charge =  1.; lifetime = .16390E-09; baryonNo = -1;}
        else if(geantID == 32) {name = TString("ANTIOMEGA+");   mass = .1672E+01; charge =  1.; lifetime = .82200E-10; baryonNo = -1;}
        else if(geantID == 45) {name = TString("DEUTERON");     mass = .1876E+01; charge =  1.; lifetime = .10000E+16; baryonNo =  2;}
        else if(geantID == 46) {name = TString("TRITON");       mass = .2809E+01; charge =  1.; lifetime = .10000E+16; baryonNo =  3;}
        else if(geantID == 47) {name = TString("ALPHA");        mass = .3727E+01; charge =  2.; lifetime = .10000E+16; baryonNo =  4;}
        else if(geantID == 48) {name = TString("GEANTINO");     mass = .0000E+00; charge =  0.; lifetime = .10000E+16; baryonNo =  0;}
        else if(geantID == 49) {name = TString("HE3");          mass = .2809E+01; charge =  2.; lifetime = .10000E+16; baryonNo =  3;}
        else if(geantID == 50) {name = TString("Cerenkov");     mass = .0000E+00; charge =  0.; lifetime = .10000E+16; baryonNo =  0;}
         
        if(mass < 0) { cout << "Unknown particle " << geantID << endl; }
        
        pt = sqrt(px * px + py * py);
        pTotal = sqrt(pt * pt + pz * pz);
        float E_Total = sqrt(mass * mass + pTotal * pTotal);
        rapidity = 0.5 * TMath::Log( (E_Total + pz) / (E_Total - pz) );

        if(pz != 0.0 && pTotal > 0.001) { 
             TVector3 pVect(px, py, pz);
             eta = pVect.PseudoRapidity(); 
        } 
        else {
             eta = 0;
        }
}

void setupOutputFilePaths(TString &plotFile0, TString &plotFile, 
	TString &plotFileC) {

	TString plotFileBase, rootOutFile;

	rootOutFile	= TString(Form("./root/etacorr.root"));
	plotFileBase = TString(Form("./ps/etacorr"));
	cout << "root file = " << rootOutFile.Data() << endl;
	cout << "plot file = " << plotFileBase.Data() << endl;
	plotFile0 = plotFileBase + TString(".pdf(");
	plotFile = plotFileBase + TString(".pdf");
	plotFileC = plotFileBase + TString(".pdf]");
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

void drawR2HistogramsToFile(TCanvas **canvases, int &iCanvas, TString plotFile, TH1 **histograms, double integral) {
		TLatex * text = new TLatex();
		text->SetTextSize(0.05);
		text->SetNDC();
		char buf[200];
		iCanvas++;
		sprintf(buf, "canvases%d", iCanvas);
		canvases[iCanvas] = new TCanvas(buf, buf, 30 * iCanvas, 30 * iCanvas, 800, (8.5 / 11.) * 800);
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Divide(3, 2, 0.0001, 0.0001);
			canvases[iCanvas]->cd(1);
				histograms[0]->Draw();
			canvases[iCanvas]->cd(2);
				histograms[1]->SetMinimum(0.5);
				histograms[1]->Draw();
			canvases[iCanvas]->cd(3);
				histograms[2]->SetStats(0);
				histograms[2]->Draw("colz");
			canvases[iCanvas]->cd(4);
				histograms[3]->SetStats(0);
				histograms[3]->Draw("colz");
			canvases[iCanvas]->cd(5);
				histograms[4]->SetStats(0);
				histograms[4]->Draw("colz");
			canvases[iCanvas]->cd(6);
				histograms[5]->SetStats(0);
				histograms[5]->SetMinimum(-0.005);
				histograms[5]->SetMaximum(0.005);
				histograms[5]->SetMarkerStyle(20);
				histograms[5]->SetMarkerSize(1);
				histograms[5]->SetMarkerColor(4);
				histograms[5]->SetLineColor(4);
				histograms[5]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral));
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Update();
		canvases[iCanvas]->Print(plotFile.Data());
}

void drawR2WindowWidthsToFile(TCanvas **canvases, int &iCanvas, TString plotFile0, TString plotFile, TH1D** histograms, double *integral) {
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
				histograms[0]->SetStats(0);
				histograms[0]->SetMinimum(-0.01);
				histograms[0]->SetMaximum(0.01);
				histograms[0]->SetLineColor(kBlue);
				histograms[0]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[0]));
			canvases[iCanvas]->cd(2);
				histograms[1]->SetStats(0);
				histograms[1]->SetMinimum(-0.01);
				histograms[1]->SetMaximum(0.01);
				histograms[1]->SetLineColor(kRed);
				histograms[1]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[1]));
			canvases[iCanvas]->cd(3);
				histograms[2]->SetStats(0);
				histograms[2]->SetMinimum(-0.01);
				histograms[2]->SetMaximum(0.01);
				histograms[2]->SetLineColor(kGreen);
				histograms[2]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[2]));
			canvases[iCanvas]->cd(4);
				histograms[3]->SetStats(0);
				histograms[3]->SetMinimum(-0.01);
				histograms[3]->SetMaximum(0.01);
				histograms[3]->SetLineColor(kOrange);
				histograms[3]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[3]));
			canvases[iCanvas]->cd(5);
				histograms[4]->SetStats(0);
				histograms[4]->SetMinimum(-0.01);
				histograms[4]->SetMaximum(0.01);
				histograms[4]->SetLineColor(kGray);
				histograms[4]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[4]));
			canvases[iCanvas]->cd(6);
				histograms[5]->SetStats(0);
				histograms[5]->SetMinimum(-0.01);
				histograms[5]->SetMaximum(0.01);
				histograms[5]->SetLineColor(kViolet);
				histograms[5]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[5]));
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Update();
		canvases[iCanvas]->Print(plotFile0.Data());
		iCanvas++;
		canvases[iCanvas] = new TCanvas(buf, buf, 30 * iCanvas, 30 * iCanvas, 800, (8.5 / 11.) * 800);
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Divide(3,2,0.0001,0.0001);
			canvases[iCanvas]->cd(1);
				histograms[6]->SetStats(0);
				histograms[6]->SetMinimum(-0.01);
				histograms[6]->SetMaximum(0.01);
				histograms[6]->SetLineColor(kYellow);
				histograms[6]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[6]));
			canvases[iCanvas]->cd(2);
				histograms[7]->SetStats(0);
				histograms[7]->SetMinimum(-0.01);
				histograms[7]->SetMaximum(0.01);
				histograms[7]->SetLineColor(kMagenta);
				histograms[7]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[7]));
			canvases[iCanvas]->cd(3);
				histograms[8]->SetStats(0);
				histograms[8]->SetMinimum(-0.01);
				histograms[8]->SetMaximum(0.01);
				histograms[8]->SetLineColor(kBlack);
				histograms[8]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[8]));
			canvases[iCanvas]->cd(4);
				histograms[8]->SetLineColor(kBlack);
				histograms[8]->Draw("hist");
				histograms[7]->SetLineColor(kMagenta);
				histograms[7]->Draw("same");
				histograms[6]->SetLineColor(kYellow);
				histograms[6]->Draw("same");
				histograms[5]->SetLineColor(kViolet);
				histograms[5]->Draw("same");
				histograms[4]->SetLineColor(kGray);
				histograms[4]->Draw("same");
				histograms[3]->SetLineColor(kOrange);
				histograms[3]->Draw("same");
				histograms[2]->SetLineColor(kGreen);
				histograms[2]->Draw("same");
				histograms[1]->SetLineColor(kRed);
				histograms[1]->Draw("same");
				histograms[0]->SetLineColor(kBlue);
				histograms[0]->Draw("same");
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Update();
		canvases[iCanvas]->Print(plotFile.Data());
}

void executeFilePlots(TCanvas **canvases, int iCanvas, TString plotFileC) {
		char buf[200];
	 	canvases[iCanvas]->Print(plotFileC.Data());
	 	gSystem->Exec(buf);
}
