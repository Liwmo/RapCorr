void GetInfo(int geid, float px, float py, float pz, 
                TString &name, float &amass, float &achg, float &alife,
                float &aeta, float &arap, float &aphi, 
                float &aptot, float &apt, float &abn);

float amass, achg, alife, aeta, arap, aphi, aptot, apt, abn; 
TString aname;
GetInfo(igid[itrk], gpx[itrk], gpy[itrk], gpz[itrk], 
        aname, amass, achg, alife, aeta, arap, aphi, aptot, apt, abn);

void GetInfo(int geid, float px, float py, float pz, 
        TString &name, float &amass, float &achg, float &alife,
        float &aeta, float &arap, float &aphi, 
        float &aptot, float &apt, float &abn ) {

        amass=achg=alife=-9;
        // GEANT ID            PARTICLE                       MASS         CHARGE        LIFETIME     BARYON #?
        if(geid== 1)      {name=TString("GAMMA");        amass=.0000E+00; achg=  0.; alife=.10000E+16; abn= 0;}
        else if(geid== 2) {name=TString("POSITRON");     amass=.5110E-03; achg=  1.; alife=.10000E+16; abn= 0;}
        else if(geid== 3) {name=TString("ELECTRON");     amass=.5110E-03; achg= -1.; alife=.10000E+16; abn= 0;}
        else if(geid== 4) {name=TString("NEUTRINO");     amass=.0000E+00; achg=  0.; alife=.10000E+16; abn= 0;}
        else if(geid== 5) {name=TString("MUON+");        amass=.1057E+00; achg=  1.; alife=.21970E-05; abn= 0;}
        else if(geid== 6) {name=TString("MUON-");        amass=.1057E+00; achg= -1.; alife=.21970E-05; abn= 0;}
        else if(geid== 7) {name=TString("PION0");        amass=.1350E+00; achg=  0.; alife=.84000E-16; abn= 0;}
        else if(geid== 8) {name=TString("PION+");        amass=.1396E+00; achg=  1.; alife=.26030E-07; abn= 0;}
        else if(geid== 9) {name=TString("PION-");        amass=.1396E+00; achg= -1.; alife=.26030E-07; abn= 0;}
        else if(geid==10) {name=TString("KAON0LONG");    amass=.4977E+00; achg=  0.; alife=.51700E-07; abn= 0;}
        else if(geid==11) {name=TString("KAON+");        amass=.4937E+00; achg=  1.; alife=.12370E-07; abn= 0;}
        else if(geid==12) {name=TString("KAON-");        amass=.4937E+00; achg= -1.; alife=.12370E-07; abn= 0;}
        else if(geid==13) {name=TString("NEUTRON");      amass=.9396E+00; achg=  0.; alife=.88700E+03; abn= 1;}
        else if(geid==14) {name=TString("PROTON");       amass=.9383E+00; achg=  1.; alife=.10000E+16; abn= 1;}
        else if(geid==15) {name=TString("ANTIPROTON");   amass=.9383E+00; achg= -1.; alife=.10000E+16; abn=-1;}
        else if(geid==16) {name=TString("KAON 0 SHORT"); amass=.4977E+00; achg=  0.; alife=.89260E-10; abn= 0;}
        else if(geid==17) {name=TString("ETA");          amass=.5475E+00; achg=  0.; alife=.54850E-18; abn= 0;}
        else if(geid==18) {name=TString("LAMBDA");       amass=.1116E+01; achg=  0.; alife=.26320E-09; abn= 1;}
        else if(geid==19) {name=TString("SIGMA+");       amass=.1189E+01; achg=  1.; alife=.79900E-10; abn= 1;}
        else if(geid==20) {name=TString("SIGMA0");       amass=.1193E+01; achg=  0.; alife=.74000E-19; abn= 1;}
        else if(geid==21) {name=TString("SIGMA-");       amass=.1197E+01; achg= -1.; alife=.14790E-09; abn= 1;}
        else if(geid==22) {name=TString("XI0");          amass=.1315E+01; achg=  0.; alife=.29000E-09; abn= 1;}
        else if(geid==23) {name=TString("XI-");          amass=.1321E+01; achg= -1.; alife=.16390E-09; abn= 1;}
        else if(geid==24) {name=TString("OMEGA-");       amass=.1672E+01; achg= -1.; alife=.82200E-10; abn= 1;}
        else if(geid==25) {name=TString("ANTINEUTRON");  amass=.9396E+00; achg=  0.; alife=.88700E+03; abn=-1;}
        else if(geid==26) {name=TString("ANTILAMBDA");   amass=.1116E+01; achg=  0.; alife=.26320E-09; abn=-1;}
        else if(geid==27) {name=TString("ANTISIGMA-");   amass=.1189E+01; achg= -1.; alife=.79900E-10; abn=-1;}
        else if(geid==28) {name=TString("ANTISIGMA0");   amass=.1193E+01; achg=  0.; alife=.74000E-19; abn=-1;}
        else if(geid==29) {name=TString("ANTISIGMA+");   amass=.1197E+01; achg=  1.; alife=.14790E-09; abn=-1;}
        else if(geid==30) {name=TString("ANTIXI0");      amass=.1315E+01; achg=  0.; alife=.29000E-09; abn=-1;}
        else if(geid==31) {name=TString("ANTIXI+");      amass=.1321E+01; achg=  1.; alife=.16390E-09; abn=-1;}
        else if(geid==32) {name=TString("ANTIOMEGA+");   amass=.1672E+01; achg=  1.; alife=.82200E-10; abn=-1;}
        else if(geid==45) {name=TString("DEUTERON");     amass=.1876E+01; achg=  1.; alife=.10000E+16; abn= 2;}
        else if(geid==46) {name=TString("TRITON");       amass=.2809E+01; achg=  1.; alife=.10000E+16; abn= 3;}
        else if(geid==47) {name=TString("ALPHA");        amass=.3727E+01; achg=  2.; alife=.10000E+16; abn= 4;}
        else if(geid==48) {name=TString("GEANTINO");     amass=.0000E+00; achg=  0.; alife=.10000E+16; abn= 0;}
        else if(geid==49) {name=TString("HE3");          amass=.2809E+01; achg=  2.; alife=.10000E+16; abn= 3;}
        else if(geid==50) {name=TString("Cerenkov");     amass=.0000E+00; achg=  0.; alife=.10000E+16; abn= 0;}
         
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
