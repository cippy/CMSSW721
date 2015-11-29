#define zlljets_Axe_noSkim_light_cxx
#include "EmanTreeAnalysis.h"
//C or C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip> //for input/output manipulators
//ROOT header files
#include <TAxis.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TVirtualFitter.h>
//my headers
#include "functionsForAnalysis.h"
#include "myClasses.h"

using namespace std;
using namespace myAnalyzerTEman;

#ifdef zlljets_Axe_noSkim_light_cxx

zlljets_Axe_noSkim_light::zlljets_Axe_noSkim_light(TTree *tree) : edimarcoTree_v2(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);

}

#endif

void zlljets_Axe_noSkim_light::loop(const char* configFileName)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   fChain->SetBranchStatus("weight",1);   // includes k-factor

   fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
   //fChain->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (pt > 20 + everything else)
   //fChain->SetBranchStatus("nTau18V",1);
   fChain->SetBranchStatus("nTauClean18V",1);

   fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
   fChain->SetBranchStatus("jetclean1",1);      // 1 if jet is cleaned, 0 otherwise
   fChain->SetBranchStatus("jetclean2",1);      // 1 if jet is cleaned, 0 otherwise
   fChain->SetBranchStatus("nJetClean30",1);    // # of jet with pt > 30 & eta < 2.5 and cleaning for against muons misidentified as PFjets   
   fChain->SetBranchStatus("JetClean_pt",1);  
   fChain->SetBranchStatus("JetClean_eta",1);  
   // fChain->SetBranchStatus("nJet",1);         // # of jets with pt > 25 && |eta| < 2.5
   // fChain->SetBranchStatus("nJet30",1);         // # of jets with pt > 30 && |eta| < 2.4
   // fChain->SetBranchStatus("nJet30a",1);       // # of jets with pt > 30 && |eta| < 4.7 
   // fChain->SetBranchStatus("Jet_pt",1);  
   // fChain->SetBranchStatus("Jet_eta",1);  
 
   fChain->SetBranchStatus("nLepGood",1);
   fChain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus
   fChain->SetBranchStatus("LepGood_pt",1);
   fChain->SetBranchStatus("LepGood_eta",1);
   fChain->SetBranchStatus("LepGood_phi",1);   
   fChain->SetBranchStatus("LepGood_mass",1);
   //fChain->SetBranchStatus("LepGood_charge",1);
   fChain->SetBranchStatus("LepGood_tightId",1);
   fChain->SetBranchStatus("LepGood_relIso04",1);
   // fChain->SetBranchStatus("ngenLep",1);         // not needed, using GenPart to study generator level quantities
   // fChain->SetBranchStatus("genLep_pdgId",1);
   // fChain->SetBranchStatus("genLep_pt",1);
   // fChain->SetBranchStatus("genLep_eta",1);
   // fChain->SetBranchStatus("genLep_phi",1);
   fChain->SetBranchStatus("mZ1",1);  // best m(ll) SF/OS
  
   fChain->SetBranchStatus("nGenPart",1);
   fChain->SetBranchStatus("GenPart_pdgId",1);
   fChain->SetBranchStatus("GenPart_motherId",1);
   fChain->SetBranchStatus("GenPart_pt",1);
   fChain->SetBranchStatus("GenPart_eta",1);
   fChain->SetBranchStatus("GenPart_phi",1);
   fChain->SetBranchStatus("GenPart_mass",1);
   fChain->SetBranchStatus("GenPart_motherIndex",1);

   fChain->SetBranchStatus("met_pt",1);
   fChain->SetBranchStatus("met_eta",1);
   fChain->SetBranchStatus("met_phi",1);

   fChain->SetBranchStatus("metNoMu_pt",1);
   // fChain->SetBranchStatus("metNoMu_eta",1);
   // fChain->SetBranchStatus("metNoMu_phi",1);

   //fChain->SetBranchStatus("nVert",1);  // number of good vertices

   // the following two branches are used for spring15_25ns samples
   fChain->SetBranchStatus("genWeight",1);
   fChain->SetBranchStatus("xsec",1);

   char ROOT_FNAME[100];
   char TXT_FNAME[100];
   char TEX_FNAME[100];
   char FLAVOUR[10];                   // e.g. "ele", "mu"
   char LL_FLAVOUR[10];             // e.g. "ee", "mumu"
   char CONTROL_SAMPLE[10];   // e.g. "Z-->ee"

   Double_t LUMI;
   Int_t NJETS;
   Double_t J1PT;
   Double_t J1ETA;
   Double_t J2PT;
   Double_t J2ETA;
   Double_t J1J2DPHI;
   Int_t LEP_PDG_ID;
   Double_t LEP1PT;
   Double_t LEP2PT;
   Double_t LEP1ETA;
   Double_t LEP2ETA;
   Double_t DILEPMASS_LOW;
   Double_t DILEPMASS_UP;
   Double_t LEP_ISO_04;
   Double_t GENLEP1PT;
   Double_t GENLEP2PT;
   Double_t GENLEP1ETA;
   Double_t GENLEP2ETA;
   Double_t GEN_ZMASS_LOW;
   Double_t GEN_ZMASS_UP;
   Int_t TAU_VETO_FLAG;
   // Int_t HLT_FLAG;
   // Double_t HLT_LEP1PT;
   // Double_t HLT_LEP2PT;
   // Double_t HLT_LEP1ETA;
   // Double_t HLT_LEP2ETA;
   // Int_t NVTXS;                           // # of points for study of u_par and u_perp vs # of reconstructed vertices nvtx
   // Int_t FIRST_NVTX;                    // starting number of vertices for met study   
   // Double_t METNOLEP_START;
   string FILENAME_BASE;

   ifstream inputFile(configFileName);

   if (inputFile.is_open()) {

     Double_t value;
     string name;
     string parameterName;
     string parameterType;

     mySpaces(cout,2);
     cout << "Printing content of " << configFileName << " file" << endl;
     mySpaces(cout,1);

     while (inputFile >> parameterType ) {

       if (parameterType == "NUMBER") {

	 inputFile >> parameterName >> value;
	 cout << setw(20) << parameterName << setw(7) << value << endl;

	 if (parameterName == "LUMI") LUMI = value;
	 else if (parameterName == "NJETS") NJETS = value;
	 else if (parameterName == "J1PT") J1PT = value;
	 else if (parameterName == "J1ETA") J1ETA = value;
	 else if (parameterName == "J2PT") J2PT = value;
	 else if (parameterName == "J2ETA") J2ETA = value;
	 else if (parameterName == "J1J2DPHI") J1J2DPHI = value;
	 else if (parameterName == "LEP_PDG_ID") LEP_PDG_ID = value;
	 else if (parameterName == "LEP1PT") LEP1PT = value;
	 else if (parameterName == "LEP2PT") LEP2PT = value;
	 else if (parameterName == "LEP1ETA") LEP1ETA = value;
	 else if (parameterName == "LEP2ETA") LEP2ETA = value;
	 else if (parameterName == "DILEPMASS_LOW") DILEPMASS_LOW = value;
	 else if (parameterName == "DILEPMASS_UP") DILEPMASS_UP = value;
	 else if (parameterName == "LEP_ISO_04") LEP_ISO_04 = value;
	 else if (parameterName == "GENLEP1PT") GENLEP1PT = value;
	 else if (parameterName == "GENLEP2PT") GENLEP2PT = value;
	 else if (parameterName == "GENLEP1ETA") GENLEP1ETA = value;
	 else if (parameterName == "GENLEP2ETA") GENLEP2ETA = value;
	 else if (parameterName == "GEN_ZMASS_LOW") GEN_ZMASS_LOW = value;
	 else if (parameterName == "GEN_ZMASS_UP") GEN_ZMASS_UP = value;
	 else if (parameterName == "TAU_VETO_FLAG") TAU_VETO_FLAG = value;

       } else if (parameterType == "STRING") {
	 
	 inputFile >> parameterName >> name;
	 cout << right << setw(20) << parameterName << "  " << left << name << endl;
	 if (parameterName == "FILENAME_BASE") FILENAME_BASE = name; 

       }

     }
     
     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << endl;
     exit(EXIT_FAILURE);

   }

   //Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   Double_t metBinEdges[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 750., 850., 1000.};
   Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;

   Double_t HTbinEdges[] = {100,200,400,600,50000};
   Int_t nHTbins = (sizeof(HTbinEdges)/sizeof(Double_t)) - 1;
  
   // selections for monojet selection (it also includes veto on muons or electrons depending on the sample
   selection jet1C("jet1C",Form("jet1pt > %4.0lf",(Double_t)J1PT),Form("nJetClean30 >= 1 && JetClean1_pt > %4.0lf && abs(JetClean1_eta) < %1.1lf && jetclean1 > 0.5",(Double_t)J1PT,J1ETA));
   selection jjdphiC("jjdphiC",Form("jjdphi < %1.1lf",J1J2DPHI),Form("only if njets = %i",NJETS));
   selection njetsC("njets","nJetClean30 <= 2");
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");
   selection tauLooseVetoC("tauLooseVetoC","tau veto");
   // additional selections for control sample
   selection oppChargeLeptonsC("oppChargeLeptonsC","OS/SF leptons");
   selection invMassC("invMassC",Form("mass in [%3.0lf,%3.0lf]",DILEPMASS_LOW,DILEPMASS_UP));
   // following selections are set differently in the next "if" statements depending on the lepton flavour 
   selection lepLooseVetoC;
   selection twoLeptonsC;
   selection twoLepLooseC;
   selection lep1tightIdIso04C;
   selection twoLepTightC;
   selection lep1ptC;
   selection lep2ptC;
   selection lep1etaC;  
   selection lep2etaC;
   selection genLepC; 
   // the following are only for electrons
   selection lep2tightIdIso04C;

   //TVector3 met, ele;    // ele is any electron to compute MetNoEle, for muons it's not needed because it's already in the tree
   TVector2 met, ele; 

   // following indices refer to the leading pair of OS/SF in the list of LepGood. They are initialized with 0 and 1 by default, but can be set with function
   // myGetPairIndexInArray (see functionsForAnalysis.cc for reference). 
   // When myGetPairIndexInArray() is called, the index of "correct" particles will be used. If they are not found (e.g. a pair of OS/SF is mismeasured as 2 mu+), 
   // indices are set as 0 and 1 (and following selection asking lep[0] and lep[1] to be OS or whatever will fail).
   Int_t firstIndex = 0;
   Int_t secondIndex = 1;

   Int_t firstIndexGen = 0;
   Int_t secondIndexGen = 1;
   Int_t recoLepFound_flag = 0;
   Int_t genLepFound_flag = 0;

   // following 2 variable are used for acceptance and efficiency selection, define below in the loop: if selection is passed they are set to 1, otherwise they are set to 0
   Int_t acceptanceSelectionDef = 0;
   Int_t efficiencySelectionDef = 0;

   Float_t *ptr_nLepLoose = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches
   Float_t *ptr_nLep10V = NULL;   

   Float_t *ptr_metNoLepPt = NULL;       // only needed for muons, it will point to the branches with the metNoMu_pt, then metNoLepPt = *ptr_metNoLepPt (metNoLepPt defined below)

   Float_t nLepLoose = 0.0;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such
   Float_t nLep10V = 0.0;
   Double_t metNoLepPt = 0.0;        // this variable will be assigned with *ptr_metNoLepPt, where the pointer will point to the branch metNoMu_pt for mu, and with a hand-defined variable for e

   // following vector are needed to compute reco-gen matching to be included in the efficiency computation. 
   //I must add it if I want to use this result to estimate Nzvv/Nzll as Ratio_BR/Axe, because it is included in Nzmm 
   Int_t efficiencyWithRecogenMatch_flag = 1;  // if it is 1, recoGen Match efficiency will be included in the efficiency definition
   Int_t recogenMatch_isPassed = 0;   // if efficiencyWithRecogenMatch_flag = 1, recogenMatch_isPassed will be set to to 0 or 1 and used to asses in the efficiency selection is passed

   TLorentzVector l1gen, l2gen;     // gen level  l1,l2  (Z->(l1 l2)
   TLorentzVector l1reco, l2reco;  // reco level

   Double_t ZgenMass = 0.0;
   Double_t currentWeight = -1.0;
   Int_t htbin = -1;
   Int_t Z_index = -1;

   strcpy(ROOT_FNAME,(FILENAME_BASE + ".root").c_str());
   strcpy(TXT_FNAME,(FILENAME_BASE + ".txt").c_str());
   strcpy(TEX_FNAME,(FILENAME_BASE + ".tex").c_str());

   if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...
     
     strcpy(FLAVOUR,"mu");
     strcpy(LL_FLAVOUR,"mumu");
     strcpy(CONTROL_SAMPLE,"Z-->mumu");
         
     ptr_nLepLoose = &nMu10V;                      // ask 2 muons
     ptr_nLep10V = &nEle10V;                         // veto on electrons
     ptr_metNoLepPt = &metNoMu_pt;               // for muons  get this variable from the tree 

     lepLooseVetoC.set("eLooseVetoC","electrons veto");
     // twoLeptonsC.set("twomuonsC","muons");
     // twoLepLooseC.set("twomuLooseC","2 loose muons");
     // lep1tightIdIso04C.set("mu1tightIdIso04C","leading muon tight","tight ID + relIso04 (as Emanuele)");
     // twoLepTightC.set("twomuTightC","2 tight muons");
     // lep1ptC.set("mu1ptC",Form("mu1pt > %3.0lf",LEP1PT),"leading muon pt");
     // lep2ptC.set("mu2ptC",Form("mu2pt > %3.0lf",LEP2PT),"trailing muon pt");
     // lep1etaC.set("mu1etaC",Form("|mu1eta| < %1.1lf",LEP1ETA),"leading muon eta");  
     // lep2etaC.set("mu2etaC",Form("|mu2eta| < %1.1lf",LEP2ETA),"trailing muon eta");
     genLepC.set("genMuonsC","muons generated");     

   } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

     strcpy(FLAVOUR,"ele");
     strcpy(LL_FLAVOUR,"ee");
     strcpy(CONTROL_SAMPLE,"Z-->ee");

     ptr_nLepLoose = &nEle10V;                      // ask 2 electrons
     ptr_nLep10V = &nMu10V;                         // veto on muons   

     lepLooseVetoC.set("muLooseVetoC","muons veto");
     // twoLeptonsC.set("twoelectronsC","electrons");
     // twoLepLooseC.set("twoeleLooseC","2 loose electrons");
     // lep1tightIdIso04C.set("ele1tightIdIso04C","leading electron tight","tight ID + relIso04 (as Emanuele)");
     // twoLepTightC.set("twoeleTightC","2 tight electrons");
     // lep1ptC.set("ele1ptC",Form("ele1pt > %3.0lf",LEP1PT),"leading electron pt");
     // lep2ptC.set("ele2ptC",Form("ele2pt > %3.0lf",LEP2PT),"trailing electron pt");
     // lep1etaC.set("ele1etaC",Form("|ele1eta| < %1.1lf",LEP1ETA),"leading electron eta");  
     // lep2etaC.set("ele2etaC",Form("|ele2eta| < %1.1lf",LEP2ETA),"trailing electron eta");
     genLepC.set("genElectronsC","electrons generated");     

     // the following are only for electrons
     lep2tightIdIso04C.set("ele2tightIdIso04C","trailing electron tight","tight ID + relIso04 (as Emanuele)");

   }
                       
   selection acceptanceC("acceptanceC","acceptance cuts");   // defined in the loop
   selection efficiencyC("efficiencyC","efficiency cuts");   // defined in the loop

   selection::checkMaskLength();
   selection::printActiveSelections(cout);

   UInt_t maskMonoJetSelection = njetsC.get2ToId() + jet1C.get2ToId() + jjdphiC.get2ToId() +
                                   lepLooseVetoC.get2ToId() + gammaLooseVetoC.get2ToId();

   if ( TAU_VETO_FLAG ) maskMonoJetSelection += tauLooseVetoC.get2ToId();

   mask *lep_acc_eff[nMetBins];
   for ( Int_t i = 0; i < nMetBins; i++) {
     lep_acc_eff[i] = new mask;
     lep_acc_eff[i]->setName(Form("%s_acc_eff:  %3.0lf < met < %3.0lf",FLAVOUR,metBinEdges[i], metBinEdges[i+1]));
     lep_acc_eff[i]->append(genLepC.get2ToId());       // to select only Z->ee or Z->mumu 
     lep_acc_eff[i]->append(maskMonoJetSelection);
     lep_acc_eff[i]->append(acceptanceC.get2ToId());
     lep_acc_eff[i]->append(efficiencyC.get2ToId());     // I can choose whether it includes the reco-gen matching efficiency
   }

   cout << "Opening file " <<ROOT_FNAME<< endl;

   TFile *rootFile = new TFile(ROOT_FNAME,"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout<<"Error: file \""<<ROOT_FNAME<<"\" was not opened."<<endl;
     exit(EXIT_FAILURE);
   }

   // the following are needed only for spring15_25ns
   Double_t SUMWEIGHTS;   // ==============  To be initialized with proper value to compute event weight in MC ========================
   vector<Double_t> sumWeightVector;
   vector<Int_t> eventsInSubsamples;
   Int_t eventCounter = 0;

   Int_t using_spring15_25ns_sample_flag = 0;
   if (FILENAME_BASE.find("spring15_25ns") != std::string::npos) {
     using_spring15_25ns_sample_flag = 1;    
     string suffix = "DYJetsToLL";
     cout << "Using spring15_25ns samples" << endl;
     mySumWeight_filler_spring15_25ns(suffix, sumWeightVector);  // this function fills the vector with the proper values of sumWeight depending on the sample
     myEventsInSubsamples_filler_spring15_25ns(suffix, eventsInSubsamples); 
   }

   // the following two variables coincide if no HLT is applied
   Double_t nTotalWeightedEvents = 0.0;       // total events (including weights) considering HLT (i.e. # of events passing HLT selection, if any)
   //Double_t nTotalWeightedEventsNoHLT = 0.0;   // total events (including weights) without considering the HLT  

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   TH1D* HmonojetW = new TH1D("HmonojetW","",nMetBins,metBinEdges);
   TH1D* HaccW = new TH1D("HaccW","",nMetBins,metBinEdges);
   TH1D* HeffW = new TH1D("HeffW","",nMetBins,metBinEdges);
   TH1D* HacceffW = new TH1D("HacceffW","",nMetBins,metBinEdges);

   TH1D* Hacc = new TH1D("Hacc","",nMetBins,metBinEdges);
   TH1D* Heff = new TH1D("Heff","",nMetBins,metBinEdges);
   TH1D* Hacceff = new TH1D("Hacceff","",nMetBins,metBinEdges);

   TH1D* HaccDivNnoW = new TH1D("HaccDivNnoW","",nMetBins,metBinEdges);
   TH1D* HeffDivNnoW = new TH1D("HeffDivNnoW","",nMetBins,metBinEdges);
   TH1D* HacceffDivNnoW = new TH1D("HacceffDivNnoW","",nMetBins,metBinEdges);

   TEfficiency *TE_EffFillW = new TEfficiency("TE_EffFillW","",nMetBins,metBinEdges);
   TEfficiency *TE_AccFillW = new TEfficiency("TE_AccFillW","",nMetBins,metBinEdges);
   TEfficiency *TE_AccEffFillW = new TEfficiency("TE_AccEffFillW","",nMetBins,metBinEdges);

   TE_EffFillW->SetUseWeightedEvents();
   TE_AccFillW->SetUseWeightedEvents();
   TE_AccEffFillW->SetUseWeightedEvents();

   TH1D *HmonoJetSel_noW_HTbin[nHTbins];;
   TH1D *Hacc_noW_HTbin[nHTbins];
   TH1D *Heff_noW_HTbin[nHTbins];
   TH1D *Hacceff_noW_HTbin[nHTbins];

   for (Int_t i = 0; i < nHTbins; i++) {
     HmonoJetSel_noW_HTbin[i] = new TH1D(Form("HmonoJetSel_noW_HTbin%3.0lfto%3.0lf",HTbinEdges[i],HTbinEdges[i+1]),"",nMetBins,metBinEdges);
     Hacc_noW_HTbin[i] = new TH1D(Form("Hacc_noW_HTbin%3.0lfto%3.0lf",HTbinEdges[i],HTbinEdges[i+1]),"",nMetBins,metBinEdges);
     Heff_noW_HTbin[i] = new TH1D(Form("Heff_noW_HTbin%3.0lfto%3.0lf",HTbinEdges[i],HTbinEdges[i+1]),"",nMetBins,metBinEdges);
     Hacceff_noW_HTbin[i] = new TH1D(Form("Hacceff_noW_HTbin%3.0lfto%3.0lf",HTbinEdges[i],HTbinEdges[i+1]),"",nMetBins,metBinEdges);
   }   

   // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
   TH1D *HmetBinEdges = new TH1D("HmetBinEdges","bin edges for met distributions",nMetBins+1,0.0,nMetBins+1);
   for (Int_t i = 0; i <= nMetBins; i++) {
     HmetBinEdges->SetBinContent(i+1,metBinEdges[i]);
   }

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"zlljets_Axe_noSkim_light::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;

   for (Int_t jentry=0; jentry<nentries; jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     if ((jentry % 500000) == 0) cout << jentry << endl;

     UInt_t eventMask = 0; 
     Double_t newwgt;

     if (using_spring15_25ns_sample_flag) {

	 if (jentry == eventCounter) {
	   eventCounter += eventsInSubsamples[htbin+1];  //note that htbin starts from -1 so that, when I do "htbin++" for the first time, it is 0
	   SUMWEIGHTS = sumWeightVector[htbin+1];
	   htbin++;  // here I use it to access histogram with given HT bin
	   cout << endl;
	   cout << "entry = " << jentry << ":   " ;
	   cout << "htbin = " << htbin+1 << "  --->  ";   // it will print 1, 2, 3 ... but as an index it would be 0, 1, 2 ...
	   cout << "sumWeight = " << SUMWEIGHTS << endl;
	   cout << endl;	   
	 }

	 newwgt = 1000 * LUMI * xsec * genWeight / SUMWEIGHTS; 

     } else {  // this part was the old part of the programme, before introduction of spring15 samples

       newwgt = weight * LUMI;

       if (jentry == 0) {
	 currentWeight = weight;
	 htbin = 0;
       }
       if (currentWeight != weight) {  // when weight changes, it means we are entering new HT bin ( weights are the same within the same HT bin, and are ordered)
	 currentWeight = weight;
	 htbin++;
       }

     }

     nTotalWeightedEvents += newwgt;  // counting events with weights

     nLepLoose = *ptr_nLepLoose;          
     nLep10V = *ptr_nLep10V;

     // Z_PDGID = 23   
     genLepFound_flag = myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, LEP_PDG_ID, 23, firstIndexGen, secondIndexGen, Z_index, GenPart_motherIndex); 

     if (genLepFound_flag) {

       ZgenMass = GenPart_mass[Z_index];

     } else continue;  // if I look for Z->mumu events, then I must skip events where Z doesn't go into muons
     
     recoLepFound_flag = myGetPairIndexInArray(LEP_PDG_ID, nLepGood, LepGood_pdgId, firstIndex, secondIndex);       

     //=============================================

     //enter this part if 2 OS/SF leptons were found among gen and reco particles. Now checking compatibilities between pairs
     // e.g. l1gen = e+, l2gen = e- ; l1reco = e+, l2reco = e- (but the charge order might not coincide)
     // now we require a DeltaR cut between them to assess that lreco comes from lgen
     // since 2 OS/SF were found to get inside here, if !(l1gen->l1reco && l2gen->l2reco) then for sure l1gen->l2reco && l2gen->l1reco
       
     if (genLepFound_flag == 1 && recoLepFound_flag == 1 && efficiencyWithRecogenMatch_flag == 1) {       

       l1gen.SetPtEtaPhiM(GenPart_pt[firstIndexGen],GenPart_eta[firstIndexGen],GenPart_phi[firstIndexGen],GenPart_mass[firstIndexGen]);
       l2gen.SetPtEtaPhiM(GenPart_pt[secondIndexGen],GenPart_eta[secondIndexGen],GenPart_phi[secondIndexGen],GenPart_mass[secondIndexGen]);
       l1reco.SetPtEtaPhiM(LepGood_pt[firstIndex],LepGood_eta[firstIndex],LepGood_phi[firstIndex],LepGood_mass[firstIndex]);
       l2reco.SetPtEtaPhiM(LepGood_pt[secondIndex],LepGood_eta[secondIndex],LepGood_phi[secondIndex],LepGood_mass[secondIndex]);

       Double_t DeltaR_lreco_lgen_pair1 = 0.0;
       Double_t DeltaR_lreco_lgen_pair2 = 0.0;
       
       if(LepGood_pdgId[firstIndex] == GenPart_pdgId[firstIndexGen] && LepGood_pdgId[secondIndex] == GenPart_pdgId[secondIndexGen]) {
	 
	 DeltaR_lreco_lgen_pair1 = l1reco.DeltaR(l1gen);
	 DeltaR_lreco_lgen_pair2 = l2reco.DeltaR(l2gen);

       } else {
	 
	 DeltaR_lreco_lgen_pair1 = l1reco.DeltaR(l2gen);
	 DeltaR_lreco_lgen_pair2 = l2reco.DeltaR(l1gen);
	 
       }
       
       if (DeltaR_lreco_lgen_pair1 < 0.1 && DeltaR_lreco_lgen_pair2 < 0.1) recogenMatch_isPassed = 1;
       else recogenMatch_isPassed = 0;

     }

     //=============================================

     if (fabs(LEP_PDG_ID) == 13) { 

       metNoLepPt = *ptr_metNoLepPt;        // casting might not be necessary: promoting float to double should be ok

       if ( genLepFound_flag && (GenPart_pt[firstIndexGen] > GENLEP1PT) && (GenPart_pt[secondIndexGen] > GENLEP2PT) && ( fabs(GenPart_eta[firstIndexGen]) < GENLEP1ETA) && ( fabs(GenPart_eta[secondIndexGen]) < GENLEP2ETA) && (ZgenMass > GEN_ZMASS_LOW) && (ZgenMass < GEN_ZMASS_UP) )  acceptanceSelectionDef = 1;
       else acceptanceSelectionDef = 0;

       if (recoLepFound_flag && (nLepLoose == 2) && (LepGood_tightId[firstIndex] >0.5) && (LepGood_relIso04[firstIndex] < LEP_ISO_04)) {

	 if (efficiencyWithRecogenMatch_flag == 1) {
	   efficiencySelectionDef = ( (recogenMatch_isPassed == 1) ? 1 : 0);
	 } else efficiencySelectionDef = 1;

       } else efficiencySelectionDef = 0;

     } else if (fabs(LEP_PDG_ID) == 11) { 

       //met.SetPtEtaPhi(met_pt,met_eta,met_phi);
       met.SetMagPhi(met_pt,met_phi);

       for (Int_t i = 0; i < nLepGood; i++) {
	 if (fabs(LepGood_pdgId[i]) == LEP_PDG_ID) {
	   //ele.SetPtEtaPhi(LepGood_pt[i],LepGood_eta[i],LepGood_phi[i]);
	   ele.SetMagPhi(LepGood_pt[i],LepGood_phi[i]);
	   met += ele;
	 }
       }

       // metNoLep vector created summing real met vector and vector sum of all electrons, but then I require exactly 2 loose leptons for the efficiency
       // In principle I should only add the Z: for muons it adds every muon, but then the usual selection for control samples would require exactly 2 muons
       // while for signal we would veto on muons
       //metNoLepPt = met.Pt();  // for electrons we define components by hand, for muons we used the variable in the tree to form the vector
       metNoLepPt = met.Mod();

       if ( genLepFound_flag && (GenPart_pt[firstIndexGen] > GENLEP1PT) && (GenPart_pt[secondIndexGen] > GENLEP2PT) &&
       	    ( fabs(GenPart_eta[firstIndexGen]) < GENLEP1ETA) && ( fabs(GenPart_eta[secondIndexGen]) < GENLEP2ETA) &&
       	    (ZgenMass > GEN_ZMASS_LOW) && (ZgenMass < GEN_ZMASS_UP) ) acceptanceSelectionDef = 1;
       else acceptanceSelectionDef = 0;

       if ( recoLepFound_flag && (nLepLoose == 2) && (LepGood_tightId[firstIndex] >0.5) && (LepGood_tightId[secondIndex] >0.5) &&
       	    (LepGood_relIso04[firstIndex] < LEP_ISO_04 ) && (LepGood_relIso04[secondIndex] < LEP_ISO_04 ) ) {

	 if (efficiencyWithRecogenMatch_flag == 1) {
	   efficiencySelectionDef = ( (recogenMatch_isPassed == 1) ? 1 : 0);
	 } else efficiencySelectionDef = 1;

       } else efficiencySelectionDef = 0;

     }

     // beginning of eventMask building

     eventMask += genLepC.addToMask( genLepFound_flag ); 
     eventMask += acceptanceC.addToMask( acceptanceSelectionDef );
     eventMask += efficiencyC.addToMask( efficiencySelectionDef );

     eventMask += jet1C.addToMask(nJetClean30 >= 1 && JetClean_pt[0] > J1PT && fabs(JetClean_eta[0] < J1ETA && jetclean1 > 0.5));
     eventMask += jjdphiC.addToMask( nJetClean30 == 1 || (nJetClean30 >= NJETS && fabs(dphijj) < J1J2DPHI && jetclean2 > 0.5));
     eventMask += njetsC.addToMask(nJetClean30 <= NJETS);
     eventMask += lepLooseVetoC.addToMask(nLep10V == 0);
     eventMask += tauLooseVetoC.addToMask(nTauClean18V == 0);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     
     // the following make sense only if recoLepFound_flag == 1 (i.e. flag is true), which means that fabs(LepGood_pdgId[firstIndex/secondIndex]) == LEP_PDG_ID) is 
     // true
     // also, 2 OS/SF leptons are present
     // if (recoLepFound_flag) {

     //   eventMask += oppChargeLeptonsC.addToMask( 1);
     //   eventMask += twoLeptonsC.addToMask(1);
     //   eventMask += twoLepLooseC.addToMask(nLepLoose == 2);
     //   eventMask += lep1ptC.addToMask((LepGood_pt[firstIndex] > LEP1PT)); 
     //   eventMask += lep1etaC.addToMask( (fabs(LepGood_eta[firstIndex]) < LEP1ETA) );
     //   eventMask += lep2ptC.addToMask((LepGood_pt[secondIndex] > LEP2PT) );
     //   eventMask += lep2etaC.addToMask((fabs(LepGood_eta[secondIndex]) < LEP2ETA) );
     //   eventMask += invMassC.addToMask((mZ1 > DILEPMASS_LOW) && (mZ1 < DILEPMASS_UP));     
     //   eventMask += lep1tightIdIso04C.addToMask((LepGood_tightId[firstIndex] == 1) && (LepGood_relIso04[firstIndex] < LEP_ISO_04 ) );
     //   if (fabs(LEP_PDG_ID) == 11) { 
     // 	 eventMask += lep2tightIdIso04C.addToMask((LepGood_tightId[secondIndex] == 1) && (LepGood_relIso04[secondIndex] < LEP_ISO_04 ));
     //   }

     // }

     // end of eventMask building
     // now entering analysis in bins of met

     // here I should not require exactly 2 loose leptons (I would to make sure the met computed summing all leptons reduces to met computed by only summing 
     // the Z to real met) because I bias the initial number of events
     // the requirement of having 2 OS leptons is included in genLepFound_flag and recoLepFound_flag

     if ( genLepFound_flag && /*(nLepLoose == 2) && */ ((eventMask & maskMonoJetSelection) == maskMonoJetSelection) ) {

       HmonoJetSel_noW_HTbin[htbin]->Fill(metNoLepPt);
       HmonojetW->Fill(metNoLepPt,newwgt);

       if (acceptanceSelectionDef) {

	 TE_AccFillW->FillWeighted(true,newwgt,metNoLepPt);

	 Hacc_noW_HTbin[htbin]->Fill(metNoLepPt);
	 HaccW->Fill(metNoLepPt,newwgt);

	 if (efficiencySelectionDef) {

	   TE_EffFillW->FillWeighted(true,newwgt,metNoLepPt);
	   Heff_noW_HTbin[htbin]->Fill(metNoLepPt);
	   HeffW->Fill(metNoLepPt,newwgt);

	 } else TE_EffFillW->FillWeighted(false,newwgt,metNoLepPt);

       } else TE_AccFillW->FillWeighted(false,newwgt,metNoLepPt);
	 
       if (acceptanceSelectionDef && efficiencySelectionDef) TE_AccEffFillW->FillWeighted(true,newwgt,metNoLepPt);
       else TE_AccEffFillW->FillWeighted(false,newwgt,metNoLepPt);

     }

     if ((metNoLepPt > metBinEdges[0]) && (metNoLepPt < metBinEdges[nMetBins])) {

       Int_t bin = myGetBin(metNoLepPt,metBinEdges,nMetBins);
       lep_acc_eff[bin]->countEvents(eventMask,newwgt);       
	 
     }                      // end of    if ((metNoLepPt > metBinEdges[0]) && (metNoLepPt < metBinEdges[nMetBins])) 
       
   }                        // end of loop on entries
 
   cout<<"creating file '"<<TXT_FNAME<<"' ..."<<endl;
   ofstream myfile(TXT_FNAME,ios::out);

   if ( !myfile.is_open() ) {

     cout<<"Error: unable to open file "<<TXT_FNAME<<" !"<<endl;
     exit(EXIT_FAILURE);
     
   }

   //opening inputFile named configFileName again to save content in myfile named TXT_FNAME

   inputFile.open(configFileName);

   if (inputFile.is_open()) {
     
     mySpaces(myfile,2);
     cout << "Saving content of " << configFileName << " file in "<< TXT_FNAME << endl;
     myfile << "Content of " << configFileName << endl;
     mySpaces(myfile,1);

     Double_t value;
     string name;
     string parameterName;
     string parameterType;

     while (inputFile >> parameterType ) {

       if (parameterType == "NUMBER") {

	 inputFile >> parameterName >> value;
	 myfile << setw(20) << parameterName << setw(7) << value << endl;

       } else if (parameterType == "STRING") {
	 
	 inputFile >> parameterName >> name;
	 myfile << right << setw(20) << parameterName << "  " << left << name << endl;

       }

     }
     
     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << " to save content in "<< TXT_FNAME << endl;
     exit(EXIT_FAILURE);

   }

   mySpaces(myfile,3);

   TH1D *HmonoJetSelSumHTbin = new TH1D("HmonoJetSelSumHTbin","",nMetBins,metBinEdges);
   TH1D *HaccSumHTbin = new TH1D("HaccSumHTbin","",nMetBins,metBinEdges);
   TH1D *HeffSumHTbin = new TH1D("HeffSumHTbin","",nMetBins,metBinEdges);
   TH1D *HacceffSumHTbin = new TH1D("HacceffSumHTbin","",nMetBins,metBinEdges);

   for ( Int_t j = 0; j < nHTbins; j++ ) {

     HmonoJetSelSumHTbin->Add(HmonoJetSel_noW_HTbin[j]);   // THIS SUMMATION IS NOT CORRECT AND LEADS TO UNCORRECT RESULTS. 
     HaccSumHTbin->Add(Hacc_noW_HTbin[j]);
     HeffSumHTbin->Add(Heff_noW_HTbin[j]);

   }

   Int_t nEvtMonoJetSelPass = -1;
   Int_t nEvtAccSelPass = -1;
   //Int_t nEvtEffSelPass = -1;

   // using [0] element to find step, all elements are equivalent for this purpose
   Int_t stepMonojetSelection_In_lepAccEff = lep_acc_eff[0]->whichStepHas(maskMonoJetSelection);
   Int_t stepAcceptance_In_lepAccEff = lep_acc_eff[0]->whichStepHas(acceptanceC.get2ToId());
   Int_t stepEfficiency_In_lepAccEff = lep_acc_eff[0]->whichStepHas(efficiencyC.get2ToId());
   // cout<<"step: MJ     acc     eff"<<endl;
   // cout<<stepMonojetSelection_In_lepAccEff<<stepAcceptance_In_lepAccEff<<stepEfficiency_In_lepAccEff<<endl;
   Double_t acc, eff, accStatErr, effStatErr, acceff, acceffStatErr;

   mySpaces(cout,2);
   mySpaces(myfile,2);
   cout << "Printing acceptance and efficiency." << endl;
   cout << "MET [GeV]     acc     acc_unc     eff     eff_unc" <<endl;
   myfile << "MET [GeV]     acc     acc_unc     eff     eff_unc" <<endl;

   TH1D *HevtPassMonoJetSel = new TH1D("HevtPassMonoJetSel","events passing monojet selection for A#times#epsilon",nMetBins,metBinEdges);
   TH1D *HevtPassAccSel = new TH1D("HevtPassAccSel","events passing acceptance selection for A#times#epsilon",nMetBins,metBinEdges);
   TH1D *HevtPassEffSel = new TH1D("HevtPassEffSel","events passing efficiency selection for A#times#epsilon",nMetBins,metBinEdges);

   for (Int_t i = 0; i < nMetBins; i++) {
  
     HevtPassMonoJetSel->SetBinContent(i+1,lep_acc_eff[i]->getEvents(stepMonojetSelection_In_lepAccEff));
     HevtPassMonoJetSel->SetBinError(  i+1,lep_acc_eff[i]->getEventsErr(stepMonojetSelection_In_lepAccEff));
     HevtPassAccSel->SetBinContent(    i+1,lep_acc_eff[i]->getEvents(stepAcceptance_In_lepAccEff));
     HevtPassAccSel->SetBinError(      i+1,lep_acc_eff[i]->getEventsErr(stepAcceptance_In_lepAccEff));
     HevtPassEffSel->SetBinContent(    i+1,lep_acc_eff[i]->getEvents(stepEfficiency_In_lepAccEff));
     HevtPassEffSel->SetBinError(      i+1,lep_acc_eff[i]->getEventsErr(stepEfficiency_In_lepAccEff));

     nEvtMonoJetSelPass = HmonoJetSelSumHTbin->GetBinContent(i+1);
     nEvtAccSelPass = HaccSumHTbin->GetBinContent(i+1);
     //nEvtEffSelPass = HeffSumHTbin->GetBinContent(i+1);

     acc = lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff];
     eff = lep_acc_eff[i]->nEvents[stepEfficiency_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff];
     Hacc->SetBinContent(i+1,acc);
     Heff->SetBinContent(i+1,eff);
     HaccDivNnoW->SetBinContent(i+1,acc);
     HeffDivNnoW->SetBinContent(i+1,eff);
     accStatErr = sqrt(acc * (1 - acc) / lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff]);
     Hacc->SetBinError(i+1,accStatErr);
     accStatErr = sqrt(acc * (1 - acc) / nEvtMonoJetSelPass);
     HaccDivNnoW->SetBinError(i+1,accStatErr);
     effStatErr = sqrt(eff * (1 - eff) / lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff]);
     Heff->SetBinError(i+1,effStatErr);
     effStatErr = sqrt(eff * (1 - eff) / nEvtAccSelPass);
     HeffDivNnoW->SetBinError(i+1,effStatErr);

     cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acc<<" "<<accStatErr<<" "<<eff<<" "<<effStatErr<<endl;
     myfile<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acc<<" "<<accStatErr<<" "<<eff<<" "<<effStatErr<<endl;

   }

   HacceffW->Divide(HeffW,HmonojetW,1,1,"B");
   HeffW->Divide(HeffW,HaccW,1,1,"B");
   HaccW->Divide(HaccW,HmonojetW,1,1,"B");

   for (Int_t j = 0; j < nHTbins;j++ ) {

     mySpaces(cout,2);
     mySpaces(myfile,2);

     Hacceff_noW_HTbin[j]->Divide(Heff_noW_HTbin[j],HmonoJetSel_noW_HTbin[j],1,1,"B");  // use binomial statistics
     Heff_noW_HTbin[j]->Divide(Heff_noW_HTbin[j],Hacc_noW_HTbin[j],1,1,"B");
     Hacc_noW_HTbin[j]->Divide(Hacc_noW_HTbin[j],HmonoJetSel_noW_HTbin[j],1,1,"B");

     cout << "Printing acceptance and efficiency: HT in [" << HTbinEdges[j] << "," << HTbinEdges[j+1] << "]" << endl;
     cout << "MET [GeV]     acc     acc_unc     eff     eff_unc" <<endl;
     myfile << "HT in [" << HTbinEdges[j] << "," << HTbinEdges[j+1] << "]" << endl;
     myfile << "MET [GeV]     acc     acc_unc     eff     eff_unc" <<endl;

     for  (Int_t i = 0; i < nMetBins; i++) {

       cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
       cout<< Hacc_noW_HTbin[j]->GetBinContent(i+1) << " "<< Hacc_noW_HTbin[j]->GetBinError(i+1) << " ";
       cout<< Heff_noW_HTbin[j]->GetBinContent(i+1) << " "<< Heff_noW_HTbin[j]->GetBinError(i+1) << endl;
       myfile<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
       myfile<< Hacc_noW_HTbin[j]->GetBinContent(i+1) << " "<< Hacc_noW_HTbin[j]->GetBinError(i+1) << " ";
       myfile<< Heff_noW_HTbin[j]->GetBinContent(i+1) << " "<< Heff_noW_HTbin[j]->GetBinError(i+1) << endl;

     }

   }

   HacceffSumHTbin->Divide(HeffSumHTbin,HmonoJetSelSumHTbin,1,1,"B"); 
   HeffSumHTbin->Divide(HeffSumHTbin,HaccSumHTbin,1,1,"B");
   HaccSumHTbin->Divide(HaccSumHTbin,HmonoJetSelSumHTbin,1,1,"B");


   mySpaces(cout,2);
   mySpaces(myfile,2);
   cout << "Printing acceptance * efficiency" << endl;
   cout << "MET [GeV]     acc*eff     acc*eff_unc" <<endl;
   myfile << "MET [GeV]     acc*eff     acc*eff_unc" <<endl;

   for (Int_t i = 0; i < nMetBins; i++) {
     // do not merge with previous loop: I want to print them after the previous loop because I might copy and paste this output to make acc * eff table
     
     nEvtMonoJetSelPass = HmonoJetSelSumHTbin->GetBinContent(i+1);
     nEvtAccSelPass = HaccSumHTbin->GetBinContent(i+1);
     //nEvtEffSelPass = HeffSumHTbin->GetBinContent(i+1);

     acceff = lep_acc_eff[i]->nEvents[stepEfficiency_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff];
     Hacceff->SetBinContent(i+1,acceff);
     HacceffDivNnoW->SetBinContent(i+1,acceff);
     acceffStatErr = sqrt(acceff * (1 - acceff) / lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff]);
     Hacceff->SetBinError(i+1,acceffStatErr);
     acceffStatErr = sqrt(acceff * (1 - acceff) / nEvtMonoJetSelPass);
     HacceffDivNnoW->SetBinError(i+1,acceffStatErr);

     cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acceff<<" "<<acceffStatErr<<endl;
     myfile<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acceff<<" "<<acceffStatErr<<endl;

   }


   for (Int_t j = 0; j < nHTbins;j++ ) {
     // do not merge with previous loop: I want to print them after the previous loop because I might copy and paste this output to make acc * eff table
  
     mySpaces(cout,2);
     mySpaces(myfile,2);
   
     cout << "Printing acceptance * efficiency: HT in [" << HTbinEdges[j] << "," << HTbinEdges[j+1] << "]" << endl;
     cout << "MET [GeV]     acc*eff     acc*eff_unc" <<endl;
     myfile << "HT in [" << HTbinEdges[j] << "," << HTbinEdges[j+1] << "]" << endl;
     myfile << "MET [GeV]     acc*eff     acc*eff_unc" <<endl;

     for (Int_t i = 0; i < nMetBins; i++) {

       cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
       cout<< Hacceff_noW_HTbin[j]->GetBinContent(i+1) << " "<< Hacceff_noW_HTbin[j]->GetBinError(i+1) << endl;
       myfile<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
       myfile<< Hacceff_noW_HTbin[j]->GetBinContent(i+1) << " "<< Hacceff_noW_HTbin[j]->GetBinError(i+1) << endl;

     } 

   }

   myfile.close();

   TEfficiency *Acceptance = new TEfficiency(*HevtPassAccSel,*HevtPassMonoJetSel);
   TEfficiency *Efficiency = new TEfficiency(*HevtPassEffSel,*HevtPassAccSel);
   TEfficiency *AccEff = new TEfficiency(*HevtPassEffSel,*HevtPassMonoJetSel);

   Acceptance->SetUseWeightedEvents();
   Efficiency->SetUseWeightedEvents();
   AccEff->SetUseWeightedEvents();

   TGraphAsymmErrors *grAE_Acc = Acceptance->CreateGraph();
   TGraphAsymmErrors *grAE_Eff = Efficiency->CreateGraph();
   TGraphAsymmErrors *grAE_AccEff = AccEff->CreateGraph();

   grAE_Acc->Write("grAE_Acc");
   grAE_Eff->Write("grAE_Eff");
   grAE_AccEff->Write("grAE_AccEff");

   TGraphAsymmErrors *grAE_AccFillW = TE_AccFillW->CreateGraph();
   TGraphAsymmErrors *grAE_EffFillW = TE_EffFillW->CreateGraph();
   TGraphAsymmErrors *grAE_AccEffFillW = TE_AccEffFillW->CreateGraph();

   grAE_AccFillW->Write("grAE_AccFillW");
   grAE_EffFillW->Write("grAE_EffFillW");
   grAE_AccEffFillW->Write("grAE_AccEffFillW");

   rootFile->Write();

   rootFile->Close();
   delete rootFile;

   //creating a .tex file to build tables with data
   FILE *fp;
   fp = fopen(TEX_FNAME,"w");

   if ( fp == NULL)  cout<<"Error: '"<<TEX_FNAME<<"' not opened"<<endl;
   else {

     cout<<"creating file '"<<TEX_FNAME<<"' ..."<<endl;
     myAddDefaultPackages(fp,TEX_FNAME);
     fprintf(fp,"\\begin{document}\n");
     fprintf(fp,"\n");
     string commentInTable;           
     commentInTable = "Note that cuts on second jet are applied only if a second jet exists with $p_t$ > 30\\,GeV.";
     fprintf(fp,"\\end{document}\n");      
     fclose(fp);

   }

   // end of tex file

}



