//-----------------------------------------------------------
//
//	just a little program to test my urqmd files at PDSF
//	w.j. llope, 04/21/2016
//	
//-----------------------------------------------------------

#include <string>
#include "stdio.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>

void GetInfo(int geid, float px, float py, float pz, TString &name, 
	            float &amass, float &achg, float &alife,
                float &aeta, float &arap, float &aphi, 
                float &aptot, float &apt, float &abn);

void UrQMD_simpleread(){
	float amass, achg, alife, aeta, arap, aphi, aptot, apt, abn; 
	TString aname;

	Int_t		ievt;
	Float_t     parimp;
	Float_t     gecm;
	Float_t     evinfo[5];
	Int_t       mpart;
	Int_t       igid[11000];   //[mpart]
	Float_t     gpx[11000];    //[mpart]
	Float_t     gpy[11000];    //[mpart]
	Float_t     gpz[11000];    //[mpart]
	
	TBranch     *b_ievt;     //!
	TBranch     *b_parimp;   //!
	TBranch     *b_gecm;     //!
	TBranch     *b_evinfo;   //!
	TBranch     *b_mpart;    //!
	TBranch     *b_igid;     //!
	TBranch     *b_gpx;      //!
	TBranch     *b_gpy;      //!
	TBranch     *b_gpz;      //!

    TString path		= TString("./events/");
	TString	filenames	= TString("urqmd_23_0099_*.root");
	TString input		= path + filenames;
	cout << input.Data() << endl;
	TChain	*chain	= new TChain("h1", "events chain");
	chain->Add(input.Data());
	int neventtree	= chain->GetEntries();
	cout << "Nevt = " << neventtree << endl;
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

	TH1D *hpid	= new TH1D("hpid", "PID ditribution", 50, 0.5, 50.5);
	TH1D *hpt	= new TH1D("hpt", "pt ditribution", 50, 0., 5.);

	//---- start event loop...
	int nb,nentries	= neventtree;
	for(Long64_t jentry = 0; jentry < nentries; jentry++) {
		nb = chain->GetEntry(jentry);
		if(jentry < 10) {
			cout << "entry=" << jentry
				 << "  ievt=" << ievt
				 << "  Ecm=" << gecm
				 << "  mpart=" << mpart
				 << "  b=" << parimp
				 << "  refmult= " << evinfo[0]
				 << "  refmult2=" << evinfo[1]
				 << "  refmult3=" << evinfo[2]
				 << "  itotcoll=" << evinfo[3]
				 << endl;
		}

		//---- start track loop...
		for(int itrk = 0; itrk < mpart; itrk++) {
			GetInfo(igid[itrk], gpx[itrk], gpy[itrk], gpz[itrk], 
    		aname, amass, achg, alife, aeta, arap, aphi, aptot, apt, abn);
			hpt->Fill(std::sqrt(gpx[itrk] * gpx[itrk] + gpy[itrk] * gpy[itrk]));
			hpid->Fill(igid[itrk]);
		}
	}
	cout << "finished..." << endl;


	//---- plotting setup
	gStyle->SetPaperSize(TStyle::kUSLetter);
	gStyle->SetLabelSize(0.05, "X");
	gStyle->SetLabelSize(0.05, "Y");
	gStyle->SetTitleXSize(0.055);
	gStyle->SetTitleYSize(0.055);
	gStyle->SetTitleOffset(0.85, "X");
	gStyle->SetTitleOffset(1.2, "Y");
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
	gStyle->SetTitleSize(0.1, "T");
	gStyle->SetPalette(1);
	gStyle->SetHistMinimumZero(kFALSE);
	
	gStyle->SetHatchesSpacing(2);
	gStyle->SetHatchesLineWidth(2);
	
	TCanvas * canvas = new TCanvas();
	canvas->Divide(2, 1, 0.01, 0.01);
	canvas->cd(1);
	hpt->Draw();
	canvas->cd(2);
	hpid->Draw();
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
}

void GetInfo(int geid, float px, float py, float pz, 
        TString &name, float &amass, float &achg, float &alife,
        float &aeta, float &arap, float &aphi, 
        float &aptot, float &apt, float &abn) {

        amass = achg = alife = -9; 
        // GEANT ID            PARTICLE                            MASS          CHARGE        LIFETIME       BARYON #?
        if(geid == 1)      {name=TString("GAMMA");         amass = .0000E+00; achg =  0.; alife = .10000E+16; abn =  0;}
        else if(geid == 2) {name=TString("POSITRON");      amass = .5110E-03; achg =  1.; alife = .10000E+16; abn =  0;}
        else if(geid == 3) {name=TString("ELECTRON");      amass = .5110E-03; achg = -1.; alife = .10000E+16; abn =  0;}
        else if(geid == 4) {name=TString("NEUTRINO");      amass = .0000E+00; achg =  0.; alife = .10000E+16; abn =  0;}
        else if(geid == 5) {name=TString("MUON+");         amass = .1057E+00; achg =  1.; alife = .21970E-05; abn =  0;}
        else if(geid == 6) {name=TString("MUON-");         amass = .1057E+00; achg = -1.; alife = .21970E-05; abn =  0;}
        else if(geid == 7) {name=TString("PION0");         amass = .1350E+00; achg =  0.; alife = .84000E-16; abn =  0;}
        else if(geid == 8) {name=TString("PION+");         amass = .1396E+00; achg =  1.; alife = .26030E-07; abn =  0;}
        else if(geid == 9) {name=TString("PION-");         amass = .1396E+00; achg = -1.; alife = .26030E-07; abn =  0;}
        else if(geid == 10) {name=TString("KAON0LONG");    amass = .4977E+00; achg =  0.; alife = .51700E-07; abn =  0;}
        else if(geid == 11) {name=TString("KAON+");        amass = .4937E+00; achg =  1.; alife = .12370E-07; abn =  0;}
        else if(geid == 12) {name=TString("KAON-");        amass = .4937E+00; achg = -1.; alife = .12370E-07; abn =  0;}
        else if(geid == 13) {name=TString("NEUTRON");      amass = .9396E+00; achg =  0.; alife = .88700E+03; abn =  1;}
        else if(geid == 14) {name=TString("PROTON");       amass = .9383E+00; achg =  1.; alife = .10000E+16; abn =  1;}
        else if(geid == 15) {name=TString("ANTIPROTON");   amass = .9383E+00; achg = -1.; alife = .10000E+16; abn = -1;}
        else if(geid == 16) {name=TString("KAON 0 SHORT"); amass = .4977E+00; achg =  0.; alife = .89260E-10; abn =  0;}
        else if(geid == 17) {name=TString("ETA");          amass = .5475E+00; achg =  0.; alife = .54850E-18; abn =  0;}
        else if(geid == 18) {name=TString("LAMBDA");       amass = .1116E+01; achg =  0.; alife = .26320E-09; abn =  1;}
        else if(geid == 19) {name=TString("SIGMA+");       amass = .1189E+01; achg =  1.; alife = .79900E-10; abn =  1;}
        else if(geid == 20) {name=TString("SIGMA0");       amass = .1193E+01; achg =  0.; alife = .74000E-19; abn =  1;}
        else if(geid == 21) {name=TString("SIGMA-");       amass = .1197E+01; achg = -1.; alife = .14790E-09; abn =  1;}
        else if(geid == 22) {name=TString("XI0");          amass = .1315E+01; achg =  0.; alife = .29000E-09; abn =  1;}
        else if(geid == 23) {name=TString("XI-");          amass = .1321E+01; achg = -1.; alife = .16390E-09; abn =  1;}
        else if(geid == 24) {name=TString("OMEGA-");       amass = .1672E+01; achg = -1.; alife = .82200E-10; abn =  1;}
        else if(geid == 25) {name=TString("ANTINEUTRON");  amass = .9396E+00; achg =  0.; alife = .88700E+03; abn = -1;}
        else if(geid == 26) {name=TString("ANTILAMBDA");   amass = .1116E+01; achg =  0.; alife = .26320E-09; abn = -1;}
        else if(geid == 27) {name=TString("ANTISIGMA-");   amass = .1189E+01; achg = -1.; alife = .79900E-10; abn = -1;}
        else if(geid == 28) {name=TString("ANTISIGMA0");   amass = .1193E+01; achg =  0.; alife = .74000E-19; abn = -1;}
        else if(geid == 29) {name=TString("ANTISIGMA+");   amass = .1197E+01; achg =  1.; alife = .14790E-09; abn = -1;}
        else if(geid == 30) {name=TString("ANTIXI0");      amass = .1315E+01; achg =  0.; alife = .29000E-09; abn = -1;}
        else if(geid == 31) {name=TString("ANTIXI+");      amass = .1321E+01; achg =  1.; alife = .16390E-09; abn = -1;}
        else if(geid == 32) {name=TString("ANTIOMEGA+");   amass = .1672E+01; achg =  1.; alife = .82200E-10; abn = -1;}
        else if(geid == 45) {name=TString("DEUTERON");     amass = .1876E+01; achg =  1.; alife = .10000E+16; abn =  2;}
        else if(geid == 46) {name=TString("TRITON");       amass = .2809E+01; achg =  1.; alife = .10000E+16; abn =  3;}
        else if(geid == 47) {name=TString("ALPHA");        amass = .3727E+01; achg =  2.; alife = .10000E+16; abn =  4;}
        else if(geid == 48) {name=TString("GEANTINO");     amass = .0000E+00; achg =  0.; alife = .10000E+16; abn =  0;}
        else if(geid == 49) {name=TString("HE3");          amass = .2809E+01; achg =  2.; alife = .10000E+16; abn =  3;}
        else if(geid == 50) {name=TString("Cerenkov");     amass = .0000E+00; achg =  0.; alife = .10000E+16; abn =  0;}
         
        if(amass < 0) { cout << "Unknown particle " << geid << endl; }
        
        apt             = sqrt(px * px + py * py);
        aptot           = sqrt(apt * apt + pz * pz);
        float etot      = sqrt(amass * amass + aptot * aptot);
        arap            = 0.5 * TMath::Log( (etot + pz) / (etot - pz) );

        if(pz != 0.0 && aptot > 0.001) { 
             TVector3 v_p(px, py, pz);
             aeta = v_p.PseudoRapidity(); 
        } 
        else {
             aeta = 0;
        }
}
