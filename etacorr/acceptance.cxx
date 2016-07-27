#include "rapcorr.h"

void plotRapidityCorrelations();
int addUrQMDEventsFromPath(TChain*);
float Efficiency(int i, float f, int j, float g, float h);
float GetEffEta(float e);
float getDetectorEta(float, float);
void getInfo(int, float, float, float, TString&, float&, float&, float&,
		    	float&, float&, float&, float&, float&, float&);
void setupOutputFilePaths(TString&, TString&, TString&);
void setStyle();
void drawR2HistogramsToFile(TCanvas**, int&, TString, TH1**, double);
void drawR2AcceptancesToFile(TCanvas**, int&, TString, TH1D**, double*);
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

	float roots		= 7.7;
	int   icentbin 	= 15; 		
	float eff		= 0.0;

	TH2D *hproton_etaeff = new TH2D("hproton_etaeff","hproton_etaeff",100,-2.0,2.0,100,0.0,10.0);

	const int N_COMPARE = 3;
	RapCorr ** rapCorr = new RapCorr*[N_COMPARE];	
	for(int i = 0; i < N_COMPARE; i++) {
		rapCorr[i] = new RapCorr(40, -1.00, 1.00);
		rapCorr[i]->setRunR2(true);
		rapCorr[i]->book();
	}
	int N_EVENTS = 10000000;
	if(N_EVENTS > nentries) {
		N_EVENTS = nentries;
	}
	TrackInfo trackInfo = {igid, gpx, gpy, gpz, name, mass, charge, 
	lifetime, eta, rapidity, phi, pTotal, pt, baryonNo};

	const int MAX_MULT = 200;
	const int PROTON = 14;
	double ** rapidityArr = new double*[N_COMPARE];
    for(int i = 0; i < N_COMPARE; i++) {
    	rapidityArr[i] = new double[MAX_MULT];
   		memset(rapidityArr[i], 0, MAX_MULT * sizeof(double)); 
    }

	for(Long64_t iEvent = 0; iEvent < N_EVENTS; iEvent++) {
		if(iEvent % (N_EVENTS / 10) == 0) { 
			cout << "processing " << iEvent << " of " << N_EVENTS << endl; 
		}
               
		nb = chain->GetEntry(iEvent);
		int N_TRACKS = mpart;

		if(parimp > 3.2) {
			continue; // allow only 0-5% central collisions...  
		} 


		for(int i = 0; i < N_COMPARE; i++) {
			int protonCount = 0; 
			for(int iTrack = 0; iTrack < N_TRACKS; iTrack++) {
				getInfo(trackInfo.geantID[iTrack], trackInfo.px[iTrack], trackInfo.py[iTrack], trackInfo.pz[iTrack], 
					trackInfo.name, trackInfo.mass, trackInfo.charge, trackInfo.lifetime, trackInfo.eta, trackInfo.rapidity, trackInfo.phi, trackInfo.pTotal, trackInfo.pt, trackInfo.baryonNo);
				
				if(trackInfo.geantID[iTrack] == PROTON && fabs(trackInfo.rapidity) <= 1.00) {
					
					if(i == 0) { //efficiency
						float zVertex = -30.0 + 60.0 * gRandom->Rndm();
						float detectorEta = getDetectorEta(trackInfo.eta, zVertex);
						//cout << "EtaPart: " << trackInfo.eta << " Vertex: " << zVertex << " EtaDet: " << detectorEta << endl;

						float eff = Efficiency(trackInfo.geantID[iTrack], roots, icentbin, trackInfo.pt, detectorEta);
						if(gRandom->Rndm() <= eff) {
							rapidityArr[i][protonCount] = trackInfo.rapidity;
							protonCount++;
						}	
					}

					if(i == 1) { //pt > 0.2
						if(trackInfo.pt > 0.2) {
							rapidityArr[i][protonCount] = trackInfo.rapidity;
							protonCount++;
						}
					}

					if(i == 2) { //perfect
						rapidityArr[i][protonCount] = trackInfo.rapidity;
						protonCount++;
					}

				}
			}
			rapCorr[i]->increment(rapidityArr[i], protonCount);
		}
	}

	double * integrals = new double[N_COMPARE];
	for(int i = 0; i < N_COMPARE; i++) {
		rapCorr[i]->calculate();
		integrals[i] = rapCorr[i]->getIntegral();
	}

	TH1D ** histograms = new TH1D*[N_COMPARE];
	for(int i = 0; i < N_COMPARE; i++) {
	    histograms[i] = rapCorr[i]->getR2dRapidity();
	}

	int iCanvas = -1;
	TCanvas *canvases[100];
	TString plotFile0, plotFile, plotFileC;
	setupOutputFilePaths(plotFile0, plotFile, plotFileC);
	drawR2AcceptancesToFile(canvases, iCanvas, plotFile0, histograms, integrals);
	for(int i = 0; i < N_COMPARE; i++) {
		TH1 **rapidityHistograms = new TH1*[6];
		rapidityHistograms[0] = (TH1D*) rapCorr[i]->getMultiplicity();
		rapidityHistograms[1] = (TH1D*) rapCorr[i]->getRapidity1D();
		rapidityHistograms[2] = (TH2D*) rapCorr[i]->getRapidity2D();
		rapidityHistograms[3] = (TH2D*) rapCorr[i]->getTensorProduct2D();
		rapidityHistograms[4] = (TH2D*) rapCorr[i]->getR2(); 
		rapidityHistograms[5] = (TH1D*) rapCorr[i]->getR2dRapidity();	
		drawR2HistogramsToFile(canvases, iCanvas, plotFile, rapidityHistograms, integrals[i]);
	}
	executeFilePlots(canvases, iCanvas, plotFileC);
}

int addUrQMDEventsFromPath(TChain *chain) {
    TString path = TString("/nfs/rhi/UrQMD/events_2016/007/");
	TString	filenames = TString("urqmd_19_*.root");
	TString input = path + filenames;
	cout << input.Data() << endl;
	chain->Add(input.Data());
	int neventtree = chain->GetEntries();
	cout << "N_events = " << neventtree << endl;
	return neventtree;
}

float Efficiency(int igid, float ecm, int icentbin, float ptloc, float etaloc){
	float	eff_pt,eff_eta,eff;
	float	rmax 		= 14.6;
	float	roots[8]	= {7.7,11.5,14.5,19.6,27.0,39.0,62.4,200};
	float rmeval;
	float a0,a1,b0,b1,c0,c1,a,b,c;
	float par_a0_bypart[6]={ 0.8991210, 0.9118330, 0.8933340, 0.8894710, 0.9023340, 0.9033230};
	float par_a1_bypart[6]={-0.0004252,-0.0004348,-0.0004258,-0.0004039,-0.0004311,-0.0004405};
	float par_b0_bypart[6]={ 0.2982580, 0.2919930, 0.2805300, 0.2838100, 0.1432090, 0.1447280};
	float par_b1_bypart[6]={ 0.0000310, 0.0000497, 0.0001448, 0.0001197, 0.0000274, 0.0000352};
	float par_c0_bypart[6]={ 5.6147900, 6.7754998, 1.2439700, 1.2842700, 4.1700602, 4.5686498};
	float par_c1_bypart[6]={-0.0074872,-0.0073951, 0.0002585, 0.0000616,-0.0029782,-0.0037452};
	float avgrm_ds19_[16]={4.45983,6.47677,9.41367,13.4289,18.3822,24.4147,31.8356,40.8872,51.3251,64.2972,78.804,95.7218,116.15,140.037,168.328,206.456};
	float avgrm_ds20_[16]={5.68122,8.14528,11.6093,16.0707,22.0018,29.4685,38.4022,49.3212,62.2412,77.1492,94.5203,114.854,138.652,166.871,200.417,245.22};
	float avgrm_ds23_[16]={6.48961,9.43501,13.9132,19.4061,25.9052,34.3486,45.3418,58.3398,73.8073,91.8219,112.777,137.746,166.693,200.552,240.366,294.058};
	float avgrm_ds25_[16]={6.94922,10.424,15.3716,21.8801,29.8557,39.3737,50.842,65.3122,82.319,102.282,125.774,152.732,183.684,220.041,263.317,319.92};
	float avgrm_ds18_[16]={8.01674,11.9431,17.4254,23.907,32.3186,43.2726,56.2182,71.6365,90.5218,112.457,137.795,167.654,201.932,241.655,288.174,348.394};
	float avgrm_ds17_[16]={8.40675,12.8959,18.3962,25.3682,34.8266,46.8351,60.8408,77.8056,97.8285,121.29,149.276,181.273,218.222,261.155,310.917,374.983};
	float avgrm_ds16_[16]={11.9048,17.3617,24.759,34.6954,47.1183,63.0139,82.4465,105.334,132.719,165.08,201.967,244.275,293.054,347.794,409.85,479.722};
	float avgrmuse[16];
	int iroots 		= -1;
	for (int irs=0;irs<8;irs++){
		if (std::fabs(roots[irs]-ecm)<0.2) {
			iroots	= irs;
		}
	}
	for (int icb=0;icb<16;icb++){
		if (iroots==0) {						// 7.7 
			avgrmuse[icb]	= avgrm_ds19_[icb];
		} else if (iroots==1) {					// 11.5
			avgrmuse[icb]	= avgrm_ds20_[icb];
		} else if (iroots==2||iroots==3) {		// 14.5 & 19.6
			avgrmuse[icb]	= avgrm_ds23_[icb];
		} else if (iroots==4) {					// 27.
			avgrmuse[icb]	= avgrm_ds25_[icb];
		} else if (iroots==5) {					// 39.
			avgrmuse[icb]	= avgrm_ds18_[icb];
		} else if (iroots==6) {					// 62.4
			avgrmuse[icb]	= avgrm_ds17_[icb];
		} else if (iroots==7) {					// 200.
			avgrmuse[icb]	= avgrm_ds16_[icb];
		}
	}
	eff = 0;
//---- basic cuts
	if (ptloc<0.2){ return eff; }
	if (std::fabs(etaloc)>1.3){ return eff; }
	//
//---- identify the particle
	int kp		= -1;
	if        (igid==8) {		// pi+
		kp    = 1;
	} else if (igid==9) {		// pi-
		kp    = 0;
	} else if (igid==11) {		// K+
		kp    = 3;
	} else if (igid==12) {		// K-
		kp    = 2;
	} else if (igid==14) {		// p
		kp    = 5;
	} else if (igid==15) {		// pbar
		kp    = 4;
	}
	if (kp<0){ return false; }
	//
//---- get the equivalent refmult
	if (icentbin>=0&&icentbin<16) {
		rmeval	= avgrmuse[icentbin];
	} else { 
		cout<<"confused about centbin.... "<<icentbin<<endl;
		exit(0);
	}
	a0			= par_a0_bypart[kp];
	a1			= par_a1_bypart[kp];
	b0			= par_b0_bypart[kp];
	b1			= par_b1_bypart[kp];
	c0			= par_c0_bypart[kp];
	c1			= par_c1_bypart[kp];
	a			= a0 + a1*rmeval;
	b			= b0 + b1*rmeval;
	c			= c0 + c1*rmeval;
	eff_pt	= a*TMath::Exp( -std::pow(b/ptloc,c) );
	eff_eta	= GetEffEta(etaloc);
	eff 	= eff_pt*eff_eta/0.91;
	if (eff_eta<=0.0) eff_eta = 0;
	//cout<<ptloc<<" "<<etaloc<<" "<<eff_pt<<" "<<eff_eta<<" "<<eff<<endl;
	return eff;
}

float GetEffEta(float etaloc) {
	const int proton_etaeff_N	= 15;
	float proton_etaeff_eta[proton_etaeff_N] = {-1.25761,-1.15457,-1.05152,-0.957845,-0.765808,-0.339578,-0.0538642,
					0.348946,0.653396,0.798595,0.957845,1.05152	,1.15457,1.25761,1.36066};
	float proton_etaeff_eff[proton_etaeff_N] = {0.0512,0.5152,0.7872,0.8976,0.91,0.91,0.91,0.91,
					0.91,0.8768,0.8496,0.728,0.5056,0.056,0.016};
	float effeta = 0.0;
	if (std::fabs(etaloc)>1.3){ return effeta; }
	for (int ie=0;ie<proton_etaeff_N-1;ie++){
		if (etaloc<proton_etaeff_eta[ie+1]){
			float slope = (proton_etaeff_eff[ie+1]-proton_etaeff_eff[ie])/(proton_etaeff_eta[ie+1]-proton_etaeff_eta[ie]);
			effeta  = proton_etaeff_eff[ie] + (etaloc-proton_etaeff_eta[ie])*slope;
			effeta	/= 0.91;
			break;
		}
	}
	return effeta; 
}

float getDetectorEta(float particleEta, float zVertex) {
	float radius = 30;
	float particleTheta, detectorTheta = 0;
	particleTheta = 2 * TMath::ATan(TMath::Exp(-particleEta)); 
	if(particleEta < 0) {
		float modifiedEta = std::fabs(particleEta);
		particleTheta = 2 * TMath::ATan(TMath::Exp(-modifiedEta)); 
		particleTheta = TMath::Pi() - particleTheta;
	}
	float z = radius / TMath::Tan(particleTheta);
	detectorTheta = TMath::ATan(radius / (z - zVertex));
	if(radius / (z - zVertex) < 0) {
		float modifiedDen = std::fabs(z - zVertex);
		detectorTheta = TMath::ATan(radius / modifiedDen);
		detectorTheta = TMath::Pi() - detectorTheta;
	}
	float detectorEta = -TMath::Log( TMath::Tan(detectorTheta / 2) );
	return detectorEta;
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

void drawR2AcceptancesToFile(TCanvas **canvases, int &iCanvas, TString plotFile0, TH1D** histograms, double *integral) {
		TLatex * text = new TLatex();
        TLegend * legend = new TLegend(0.1, 0.9, 0.38, 0.7);
		text->SetTextSize(0.05);
		text->SetNDC();
		char buf[200];
		iCanvas++;
		sprintf(buf, "canvases%d", iCanvas);
		canvases[iCanvas] = new TCanvas(buf, buf, 30 * iCanvas, 30 * iCanvas, 800, (8.5 / 11.) * 800);
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Divide(2, 2, 0.0001, 0.0001);
			canvases[iCanvas]->cd(1);
                legend->AddEntry(histograms[0], "Z_{vtx} Acceptance", "l");
				histograms[0]->SetStats(0);
				histograms[0]->SetMinimum(-0.01);
				histograms[0]->SetMaximum(0.01);
				histograms[0]->SetLineColor(kBlue);
				histograms[0]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[0]));
			canvases[iCanvas]->cd(2);
                legend->AddEntry(histograms[1], "P_{t} > 0", "l");
				histograms[1]->SetStats(0);
				histograms[1]->SetMinimum(-0.01);
				histograms[1]->SetMaximum(0.01);
				histograms[1]->SetLineColor(kRed);
				histograms[1]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[1]));
			canvases[iCanvas]->cd(3);
                legend->AddEntry(histograms[2], "Perfect", "l");
				histograms[2]->SetStats(0);
				histograms[2]->SetMinimum(-0.01);
				histograms[2]->SetMaximum(0.01);
				histograms[2]->SetLineColor(kGreen);
				histograms[2]->Draw("hist");
				text->DrawLatex(0.2, 0.8, Form("integral=%5.3f", integral[2]));
			canvases[iCanvas]->cd(4);
				histograms[2]->SetLineColor(kGreen);
				histograms[2]->Draw("hist");
				histograms[1]->SetLineColor(kRed);
				histograms[1]->Draw("same");
				histograms[0]->SetLineColor(kBlue);
				histograms[0]->Draw("same");
                legend->Draw();
		canvases[iCanvas]->cd(); 
		canvases[iCanvas]->Update();
		canvases[iCanvas]->Print(plotFile0.Data());
}

void executeFilePlots(TCanvas **canvases, int iCanvas, TString plotFileC) {
		char buf[200];
	 	canvases[iCanvas]->Print(plotFileC.Data());
	 	gSystem->Exec(buf);
}
