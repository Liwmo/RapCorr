//-----------------------------------------------------------
//
//	just a little program to test my urqmd files at PDSF
//	w.j. llope, 04/21/2016
//	
//-----------------------------------------------------------

void UrQMD_simpleread(){

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
	for (Long64_t jentry = 0; jentry < nentries; jentry++) {
		nb = chain->GetEntry(jentry);
		if (jentry < 10){
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
		for (int itrk = 0; itrk < mpart; itrk++){
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
