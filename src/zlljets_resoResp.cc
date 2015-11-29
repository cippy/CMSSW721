#define zlljets_resoResp_cxx
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

#ifdef zlljets_resoResp_cxx

zlljets_resoResp::zlljets_resoResp(TTree *tree) : edimarcoTree_v2(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);

}

#endif

void zlljets_resoResp::loop(const char* configFileName, const Int_t ISDATA_FLAG, const Int_t unweighted_event_flag)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   fChain->SetBranchStatus("weight",1);   // includes k-factor
   //fChain->SetBranchStatus("LHEorigWeight",1); // contains negative values: the weight in the event is weight*LHEorigWeight

   fChain->SetBranchStatus("genWeight",1); 

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

   if (!ISDATA_FLAG) {
     fChain->SetBranchStatus("nGenPart",1);
     fChain->SetBranchStatus("GenPart_pdgId",1);
     fChain->SetBranchStatus("GenPart_motherId",1);
     fChain->SetBranchStatus("GenPart_pt",1);
     fChain->SetBranchStatus("GenPart_eta",1);
     fChain->SetBranchStatus("GenPart_phi",1);
     fChain->SetBranchStatus("GenPart_mass",1);
     fChain->SetBranchStatus("GenPart_motherIndex",1);

     fChain->SetBranchStatus("vtxW",1);   // weight to have better agreement between data and MC
     fChain->SetBranchStatus("xsec",1);   // weight to have better agreement between data and MC
     fChain->SetBranchStatus("vtxWeight",1);   // weight to have better agreement between data and MC
   }

   fChain->SetBranchStatus("met_pt",1);
   //fChain->SetBranchStatus("met_eta",1);
   fChain->SetBranchStatus("met_phi",1);

   fChain->SetBranchStatus("metNoMu_pt",1);
   //fChain->SetBranchStatus("metNoMu_eta",1);
   fChain->SetBranchStatus("metNoMu_phi",1);

   fChain->SetBranchStatus("nVert",1);  // number of good vertices 

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
   Double_t RATIO_BR_ZINV_ZLL;
   Double_t UNC_RATIO_BR_ZINV_ZLL;
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
   Int_t HLT_FLAG;
   Double_t HLT_LEP1PT;
   Double_t HLT_LEP2PT;
   Double_t HLT_LEP1ETA;
   Double_t HLT_LEP2ETA;
   Int_t NVTXS;                           // # of points for study of u_par and u_perp vs # of reconstructed vertices nvtx
   Int_t FIRST_NVTX;                    // starting number of vertices for met study   
   Double_t METNOLEP_START;
   Int_t JETS_SELECTION_RESORESP_FLAG;
   Int_t PHOTON_VETO_RESORESP_FLAG;
   Int_t LEPTON_VETO_RESORESP_FLAG;
   Double_t XSEC_OVER_NPROCESSED;
   Double_t SUMWEIGHTS;
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
	 cout << right << setw(20) << parameterName << "  " << left << value << endl;

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
	 else if (parameterName == "TAU_VETO_FLAG") TAU_VETO_FLAG = value;
	 else if (parameterName == "HLT_FLAG") HLT_FLAG = value;
	 else if (parameterName == "HLT_LEP1PT") HLT_LEP1PT = value;
	 else if (parameterName == "HLT_LEP2PT") HLT_LEP2PT = value;
	 else if (parameterName == "HLT_LEP1ETA") HLT_LEP1ETA = value;
	 else if (parameterName == "HLT_LEP2ETA") HLT_LEP2ETA = value;
	 else if (parameterName == "NVTXS") NVTXS = value;
	 else if (parameterName == "FIRST_NVTX") FIRST_NVTX = value;
	 else if (parameterName == "METNOLEP_START") METNOLEP_START = value;
	 else if (parameterName == "JETS_SELECTION_RESORESP_FLAG") JETS_SELECTION_RESORESP_FLAG = value;
	 else if (parameterName == "PHOTON_VETO_RESORESP_FLAG") PHOTON_VETO_RESORESP_FLAG = value;
	 else if (parameterName == "LEPTON_VETO_RESORESP_FLAG") LEPTON_VETO_RESORESP_FLAG = value;

	 if (!ISDATA_FLAG) {

	 if (parameterName == "GENLEP1PT") GENLEP1PT = value;
	 else if (parameterName == "GENLEP2PT") GENLEP2PT = value;
	 else if (parameterName == "GENLEP1ETA") GENLEP1ETA = value;
	 else if (parameterName == "GENLEP2ETA") GENLEP2ETA = value;
	 else if (parameterName == "GEN_ZMASS_LOW") GEN_ZMASS_LOW = value;
	 else if (parameterName == "GEN_ZMASS_UP") GEN_ZMASS_UP = value;
	 else if (parameterName == "XSEC_OVER_NPROCESSED") XSEC_OVER_NPROCESSED = value;
	 else if (parameterName == "SUMWEIGHTS") SUMWEIGHTS = value;

	 }
	 

       } else if (parameterType == "STRING") {
	 
	 inputFile >> parameterName >> name;
	 cout << right << setw(20) << parameterName << "  " << left << name << endl;
	 if (parameterName == "FILENAME_BASE") {

	   FILENAME_BASE = name; 
	   if ( !ISDATA_FLAG && unweighted_event_flag) FILENAME_BASE += "_weq1";  // if using unit weight, add _weq1 to filename (weq1 means weight = 1)

	 }

       }

     }
     
     mySpaces(cout,2);

     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << endl;
     exit(EXIT_FAILURE);

   }

   //Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   Double_t metBinEdges[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 750., 850., 1000.};
   Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;

   // vector<Double_t> metCut;
   // metCut.push_back(250);
   // metCut.push_back(300);
   // metCut.push_back(350);
   // metCut.push_back(400);
   // metCut.push_back(500);

   //selection njetsC("njetsC",Form("njets <= %i",NJETS),"pt > 30; |eta| < 4.7");   // using nJet30a
   //selection njetsEmanC("njetsEmanC","njets","1 or 2 jets, cleaning, pt > 30, |eta| < 2.4");       // using nJet30
   //selection njetsEmanC("njetsEmanC","njets","1 or 2 jets, cleaning, pt > 30, |eta| < 2.5");       // using nJetClen30
   //selection jet1ptC("jet1ptC",Form("jet1pt > %4.0lf",(Double_t)J1PT));
   //selection jet1etaC("jet1etaC",Form("|jet1eta| < %2.1lf",J1ETA));
   //selection jet2etaC("jet2etaC",Form("|jet2eta| < %2.1lf",J2ETA),Form("only if njets = %i",NJETS));
  
   // selections for monojet selection (it also includes veto on muons or electrons depending on the sample
   selection jet1C("jet1C",Form("jet1pt > %4.0lf",(Double_t)J1PT),Form("nJetClean >= 1 && JetClean1_pt > %4.0lf && abs(JetClean1_eta) < %1.1lf && jetclean1 > 0.5",(Double_t)J1PT,J1ETA));
   selection jjdphiC("jjdphiC",Form("jjdphi < %1.1lf",J1J2DPHI),Form("only if njets = %i",NJETS));
   selection njetsC("njets","nJetClean30 <= 2");
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");
   selection tauLooseVetoC;
   if (TAU_VETO_FLAG) tauLooseVetoC.set("tauLooseVetoC","tau veto");
   // additional selections for control sample
   selection oppChargeLeptonsC("oppChargeLeptonsC","OS/SF leptons");
   selection invMassC("invMassC",Form("mass in [%3.0lf,%3.0lf]",DILEPMASS_LOW,DILEPMASS_UP));
   // following selections are set differently in the next "if" statements depending on the lepton flavour 
   //selection metNoLepC[metCut.size()];
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
   selection metNoLepStartC;
   selection HLTlepC;
   // the following are only for electrons
   selection lep2tightIdIso04C;

   // TVector3 metNoLepTV3;   // metNoLep 3D vector, TV3 is to make clear that it's a 3D vector and not a TLorentzVector
   // TVector3 ele; 
   // ele is any electron to compute MetNoEle, for muons it's not needed because it's already in the tree
   //TVector3 met, eleVectorSum;

   // the reason to use TVector3 instead of TVector2 (which would be faster) is that I need to compute dphi between metNoLep and the Z vector (TLorentVector).
   // To do that I use the TVector3::DeltaPhi(const TVector3&) method passing TLorentzVector::Vect() which returns the 3D vector from a Lorentz Vector.
   // there isn't a method giving just the transverse 2D vector (although I could do it myself) from a TLorentVector (but it does exist a TVector3::XYvector()). 
   // Maybe I will change because TVector3 is much slower than TVector2

   TVector2 metNoLepTV, ele;

   TLorentzVector l1gen, l2gen, Zgen;     // gen level Z and l1,l2  (Z->(l1 l2)
   TLorentzVector l1reco, l2reco, Zreco;

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
   //Int_t recoGenMatchDR_flag = 0;
   Int_t genTauFound_flag = 0;
   Int_t Z_index = 0; 

   Double_t nTotalWeightedEvents = 0.0;     
   //Double_t nEventsAfterMatchRecoGen = 0.0;
   Int_t HLT_passed_flag = 1;          // some computations (for e) require a trigger preselection, while other don't. The former will be done if the flag is set to 1
                                                       // it's set to 1 because if the trigger selection is not applied every event must be considered to be a "good" event having passed all preselections
                                                       // Actually in this code the trigger is necessary, but I keep it like this nonetheless.

   // following 2 variable are used for acceptance and efficiency selection, define below in the loop: if selection is passed they are set to 1, otherwise they are set to 0
   // Int_t acceptanceSelectionDef = 0;
   // Int_t efficiencySelectionDef = 0;

   Float_t *ptr_nLepLoose = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches
   Float_t *ptr_nLep10V = NULL;   

   Float_t *ptr_metNoLepPt = NULL;       // only needed for muons, it will point to the branches with the metNoMu_pt, then metNoLepPt = *ptr_metNoLepPt (metNoLepPt defined below)
   //Float_t *ptr_metNoLepEta = NULL; 
   Float_t *ptr_metNoLepPhi = NULL;  

   Float_t nLepLoose = 0.0;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such
   Float_t nLep10V = 0.0;

   Double_t metNoLepPt = 0.0;        // this variable will be assigned with *ptr_metNoLepPt, where the pointer will point to the branch metNoMu_pt for mu, and with a hand-defined variable for e
   //Double_t metNoLepEta = 0.0;
   Double_t metNoLepPhi = 0.0;   // same story as above

   Int_t using_phys14_sample_flag = 0;
   if (FILENAME_BASE.find("phys14") != std::string::npos) {
     using_phys14_sample_flag = 1;    
     cout << "Using phys14 samples" << endl;
   }

   Int_t using_spring15_sample_flag = 0;
   if (FILENAME_BASE.find("spring15") != std::string::npos) {
     using_spring15_sample_flag = 1;    
     cout << "Using spring15 samples" << endl;
   }

   if ( !ISDATA_FLAG && unweighted_event_flag) cout << "Warning: no weight applied to events (w = 1)" << endl;  // if MC with unit weight, make user know

   // if using sample spring15, need to use vtxW to get same Nvtx distribution as seen in data. For older trees it's not used

   // the following flag is needed to enable search for Z->ll at generator level. For MC samples different from DYJetsToLL I must not require 2 gen leptons from Z
   Int_t using_zlljets_MCsample_flag = 0;
   if ( !ISDATA_FLAG && ( (FILENAME_BASE.find("zmumujets") != std::string::npos) || (FILENAME_BASE.find("zeejets") != std::string::npos) ) )  using_zlljets_MCsample_flag = 1; 

   Int_t using_ztautaujets_MCsample_flag = 0;
   if ( !ISDATA_FLAG && (FILENAME_BASE.find("ztautaujets") != std::string::npos) ) using_ztautaujets_MCsample_flag = 1; 

   if (ISDATA_FLAG) {
     strcpy(ROOT_FNAME,(FILENAME_BASE + "_DATA.root").c_str());
     strcpy(TXT_FNAME,(FILENAME_BASE + "_DATA.txt").c_str());
     strcpy(TEX_FNAME,(FILENAME_BASE + "_DATA.tex").c_str());
   } else {
     strcpy(ROOT_FNAME,(FILENAME_BASE + ".root").c_str());
     strcpy(TXT_FNAME,(FILENAME_BASE + ".txt").c_str());
     strcpy(TEX_FNAME,(FILENAME_BASE + ".tex").c_str());
   }

   if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...
     
     strcpy(FLAVOUR,"mu");
     strcpy(LL_FLAVOUR,"mumu");
     strcpy(CONTROL_SAMPLE,"Z-->mumu");
         
     ptr_nLepLoose = &nMu10V;                      // ask 2 muons
     ptr_nLep10V = &nEle10V;                         // veto on electrons
     ptr_metNoLepPt = &metNoMu_pt;               // for muons  get this variable from the tree 
     //ptr_metNoLepEta = &metNoMu_eta;               // for muons  get this variable from the tree 
     ptr_metNoLepPhi = &metNoMu_phi;         // for muons  get this variable from the tree

     // for (Int_t i = 0; i < metCut.size(); i++) {
     //   metNoLepC[i].set(Form("metNoMuC[%i]",i),Form("metNoMu > %3.0lf",metCut.at(i)));
     // }

     lepLooseVetoC.set("eLooseVetoC","electrons veto");
     twoLeptonsC.set("twomuonsC","muons");
     twoLepLooseC.set("twomuLooseC","2 loose muons");
     lep1tightIdIso04C.set("mu1tightIdIso04C","leading muon tight","tight ID + relIso04 (as Emanuele)");
     twoLepTightC.set("twomuTightC","2 tight muons");
     lep1ptC.set("mu1ptC",Form("mu1pt > %3.0lf",LEP1PT),"leading muon pt");
     lep2ptC.set("mu2ptC",Form("mu2pt > %3.0lf",LEP2PT),"trailing muon pt");
     lep1etaC.set("mu1etaC",Form("|mu1eta| < %1.1lf",LEP1ETA),"leading muon eta");  
     lep2etaC.set("mu2etaC",Form("|mu2eta| < %1.1lf",LEP2ETA),"trailing muon eta");
      if (!ISDATA_FLAG && using_zlljets_MCsample_flag) genLepC.set("genMuonsC","muons generated");     
     metNoLepStartC.set("metNoMu200C",Form("metNoMu > %2.0lf",METNOLEP_START));
     HLTlepC.set("HLTmuonC","HLT for muons");
     lep2tightIdIso04C.set("mu2tightIdIso04C","trailing muon tight","tight ID + relIso04 (as Emanuele)");


   } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

     strcpy(FLAVOUR,"ele");
     strcpy(LL_FLAVOUR,"ee");
     strcpy(CONTROL_SAMPLE,"Z-->ee");

     ptr_nLepLoose = &nEle10V;                      // ask 2 electrons
     ptr_nLep10V = &nMu10V;                         // veto on muons   

     // for (Int_t i = 0; i < metCut.size(); i++) {
     //   metNoLepC[i].set(Form("metNoEleC[%i]",i),Form("metNoEle > %3.0lf",metCut.at(i)));
     // }
     
     lepLooseVetoC.set("muLooseVetoC","muons veto");
     twoLeptonsC.set("twoelectronsC","electrons");
     twoLepLooseC.set("twoeleLooseC","2 loose electrons");
     lep1tightIdIso04C.set("ele1tightIdIso04C","leading electron tight","tight ID + relIso04 (as Emanuele)");
     twoLepTightC.set("twoeleTightC","2 tight electrons");
     lep1ptC.set("ele1ptC",Form("ele1pt > %3.0lf",LEP1PT),"leading electron pt");
     lep2ptC.set("ele2ptC",Form("ele2pt > %3.0lf",LEP2PT),"trailing electron pt");
     lep1etaC.set("ele1etaC",Form("|ele1eta| < %1.1lf",LEP1ETA),"leading electron eta");  
     lep2etaC.set("ele2etaC",Form("|ele2eta| < %1.1lf",LEP2ETA),"trailing electron eta");
     if (!ISDATA_FLAG && using_zlljets_MCsample_flag) genLepC.set("genElectronsC","electrons generated");     
     metNoLepStartC.set("metNoEle200C",Form("metNoEle > %2.0lf",METNOLEP_START));
     HLTlepC.set("HLTelectronC","HLT for electrons");
     lep2tightIdIso04C.set("ele2tightIdIso04C","trailing electron tight","tight ID + relIso04 (as Emanuele)");

   }

   selection recoGenLepMatchC;
   if (!ISDATA_FLAG && using_zlljets_MCsample_flag) {

     if (fabs(LEP_PDG_ID) == 13) recoGenLepMatchC.set("recoGenMuMatchC","reco-gen match (DR = 0.1)","only for zlljets: looks for matching of reco and gen particles");      
     else if (fabs(LEP_PDG_ID) == 11) recoGenLepMatchC.set("recoGenEleMatchC","reco-gen match (DR = 0.1)","only for zlljets: looks for matching of reco and gen particles");    
  
   }

   selection genTauC;
   if (!ISDATA_FLAG && using_ztautaujets_MCsample_flag) genTauC.set("genTauC","taus generated");                       
   //selection acceptanceC("acceptanceC","acceptance cuts");
   //selection efficiencyC("efficiencyC","efficiency cuts");

   selection::checkMaskLength();
   selection::printActiveSelections(cout); 

   UInt_t maskJetsSelection = njetsC.get2ToId() + jet1C.get2ToId() + jjdphiC.get2ToId();

   UInt_t maskMonoJetSelection = maskJetsSelection + lepLooseVetoC.get2ToId() + gammaLooseVetoC.get2ToId();

   if ( TAU_VETO_FLAG ) maskMonoJetSelection += tauLooseVetoC.get2ToId();

   UInt_t maskTightTag;   // holds cuts for lepton tight selection, which is different between muons and electrons

   //mask zlljetsControlSample(Form("%s control sample with selection flow as Emanuele's",CONTROL_SAMPLE));

   mask zlljetsControlSample(Form("%s control sample (%s gen if DYJetsToLL MC) with selection flow as Emanuele's",CONTROL_SAMPLE,FLAVOUR));
   if (!ISDATA_FLAG) {
     if (using_zlljets_MCsample_flag) zlljetsControlSample.append(genLepC.get2ToId());
     if (using_ztautaujets_MCsample_flag) zlljetsControlSample.append(genTauC.get2ToId());
   }
   zlljetsControlSample.append(HLTlepC.get2ToId());

   // mask tautaubkgInZll(Form("tau tau background in %s control sample",CONTROL_SAMPLE));
   // tautaubkgInZll.append(genTausC.get2ToId());

   mask resoAndResponse("selection for resolution and response");
   if (!ISDATA_FLAG) {
     if (using_zlljets_MCsample_flag) resoAndResponse.append(genLepC.get2ToId());
     if (using_ztautaujets_MCsample_flag) resoAndResponse.append(genTauC.get2ToId());
   }
   resoAndResponse.append(HLTlepC.get2ToId());
   resoAndResponse.append(oppChargeLeptonsC.get2ToId());
   resoAndResponse.append(twoLepLooseC.get2ToId());
   

   // if (HLT_FLAG) {

   //   //zlljetsControlSample.append(HLTlepC.get2ToId());     
   //   //tautaubkgInZll.append(HLTlepC.get2ToId());

   // }

   if (METNOLEP_START) {
       
     zlljetsControlSample.append(metNoLepStartC.get2ToId());
     //resoAndResponse.append(metNoLepStartC.get2ToId());

   }

   if (fabs(LEP_PDG_ID) == 13) {  

     maskTightTag = lep1tightIdIso04C.get2ToId() + /* lep2tightIdIso04C.get2ToId() +*/lep1ptC.get2ToId() + lep2ptC.get2ToId() + lep1etaC.get2ToId() + lep2etaC.get2ToId();  ;  // for now tight requirements on pt and eta are already included in the loose condition because they coincide (not true for electrons)
   
     // zlljetsControlSample.append(metNoLepStartC.get2ToId());
     // zlljetsControlSample.append(twoLepLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
     // zlljetsControlSample.append(twoLeptonsC.get2ToId());
     // zlljetsControlSample.append(maskTightTag);
     // zlljetsControlSample.append(invMassC.get2ToId());
         
     //zlljetsControlSample.append(metNoLepStartC.get2ToId());
     zlljetsControlSample.append(oppChargeLeptonsC.get2ToId());
     zlljetsControlSample.append(twoLepLooseC.get2ToId());
     zlljetsControlSample.append(twoLeptonsC.get2ToId());
     zlljetsControlSample.append(maskTightTag); 
     zlljetsControlSample.append(invMassC.get2ToId());
   
     // tautaubkgInZll.append(metNoLepStartC.get2ToId());
     // tautaubkgInZll.append(twoLepLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
     // tautaubkgInZll.append(twoLeptonsC.get2ToId());
     // tautaubkgInZll.append(maskTightTag);
     // tautaubkgInZll.append(invMassC.get2ToId());

     resoAndResponse.append(maskTightTag + lep2tightIdIso04C.get2ToId());
     resoAndResponse.append(invMassC.get2ToId());

   } else if (fabs(LEP_PDG_ID) == 11) {  

     maskTightTag = lep1tightIdIso04C.get2ToId() + lep2tightIdIso04C.get2ToId() + lep1ptC.get2ToId() + lep2ptC.get2ToId() + lep1etaC.get2ToId() + lep2etaC.get2ToId();

     // zlljetsControlSample.append(metNoLepStartC.get2ToId());
     // zlljetsControlSample.append(oppChargeLeptonsC.get2ToId()); // skip loose requirement because I wil ask the tight one for both
     // zlljetsControlSample.append(twoLeptonsC.get2ToId());
     // zlljetsControlSample.append(maskTightTag);
     // zlljetsControlSample.append(invMassC.get2ToId());
       
     //zlljetsControlSample.append(metNoLepStartC.get2ToId());
     zlljetsControlSample.append(oppChargeLeptonsC.get2ToId());
     zlljetsControlSample.append(twoLepLooseC.get2ToId());
     zlljetsControlSample.append(twoLeptonsC.get2ToId());
     zlljetsControlSample.append(maskTightTag);
     zlljetsControlSample.append(invMassC.get2ToId());
   
     // tautaubkgInZll.append(metNoLepStartC.get2ToId());
     // tautaubkgInZll.append(oppChargeLeptonsC.get2ToId());
     // tautaubkgInZll.append(twoLeptonsC.get2ToId());
     // tautaubkgInZll.append(maskTightTag);
     // tautaubkgInZll.append(invMassC.get2ToId());
     
     resoAndResponse.append(lep1ptC.get2ToId() + lep2ptC.get2ToId() + lep1etaC.get2ToId() + lep2etaC.get2ToId());
     resoAndResponse.append(lep1tightIdIso04C.get2ToId() + lep2tightIdIso04C.get2ToId());
     resoAndResponse.append(invMassC.get2ToId());

   }
 
   // zlljetsControlSample.append(jet1C.get2ToId());
   // zlljetsControlSample.append(jjdphiC.get2ToId());
   // zlljetsControlSample.append(njetsC.get2ToId());
   // zlljetsControlSample.append(lepLooseVetoC.get2ToId());
   // zlljetsControlSample.append(gammaLooseVetoC.get2ToId());
   
   zlljetsControlSample.append(jet1C.get2ToId());
   zlljetsControlSample.append(jjdphiC.get2ToId());
   zlljetsControlSample.append(njetsC.get2ToId());
   zlljetsControlSample.append(lepLooseVetoC.get2ToId());
   zlljetsControlSample.append(gammaLooseVetoC.get2ToId());

   // tautaubkgInZll.append(jet1C.get2ToId());
   // tautaubkgInZll.append(jjdphiC.get2ToId());
   // tautaubkgInZll.append(njetsC.get2ToId());
   // tautaubkgInZll.append(lepLooseVetoC.get2ToId());
   // tautaubkgInZll.append(gammaLooseVetoC.get2ToId());

   if (JETS_SELECTION_RESORESP_FLAG) resoAndResponse.append(maskJetsSelection);
   else resoAndResponse.append(njetsC.get2ToId());      // I always consider at least < =2 jets (it may also be 0 jets)
   if (PHOTON_VETO_RESORESP_FLAG) resoAndResponse.append(gammaLooseVetoC.get2ToId());
   if (LEPTON_VETO_RESORESP_FLAG) resoAndResponse.append(lepLooseVetoC.get2ToId());

   if (TAU_VETO_FLAG) {
     // zlljetsControlSample.append(tauLooseVetoC.get2ToId());
     zlljetsControlSample.append(tauLooseVetoC.get2ToId());
     resoAndResponse.append(tauLooseVetoC.get2ToId());
   }

   if (!ISDATA_FLAG && using_zlljets_MCsample_flag) {
     zlljetsControlSample.append(recoGenLepMatchC.get2ToId());
     resoAndResponse.append(recoGenLepMatchC.get2ToId());
   }

   cout << "Opening file " <<ROOT_FNAME<< endl;

   TFile *rootFile = new TFile(ROOT_FNAME,"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout<<"Error: file \""<<ROOT_FNAME<<"\" was not opened."<<endl;
     exit(EXIT_FAILURE);
   }
 

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   //TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   //Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};       

   Float_t invMassBinWidth = 1.0;  // invariant mass histogram's bin width in GeV
   Int_t NinvMassBins = (DILEPMASS_UP - DILEPMASS_LOW) / invMassBinWidth;

   //TH1D *HzlljetsYieldsMetBin = new TH1D("HzlljetsYieldsMetBin",Form("yields of %s control sample in bins of met;#slash{E}_{T};# of events",CONTROL_SAMPLE),nMetBins,metBinEdges);
   TH1D *HzlljetsYieldsMetBin = new TH1D("HzlljetsYieldsMetBin",Form("yields of %s control sample (%s gen if DY MC) in bins of met; #slash{E}_{T};# of events",CONTROL_SAMPLE,CONTROL_SAMPLE),nMetBins,metBinEdges);
   //TH1D *HzlljetsYieldsMetBinGenTau = new TH1D("HzlljetsYieldsMetBinGenTau",Form("yields of %s control sample (Z->#tau#tau gen) in bins of met; #slash{E}_{T};# of events",CONTROL_SAMPLE),nMetBins,metBinEdges);
   
   TH1D *HinvMass = new TH1D("HinvMass","",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);    // for MC it's done on Z->mumu or Z->ee at gen level
   TH1D *HvtxDistribution = new TH1D("HvtxDistribution","",40,-0.5,39.5);   
   TH1D *HnjetsDistribution = new TH1D("HnjetsDistribution","njets using nJetClean30",10,-0.5,9.5);   
   TH1D *Hj1j2dphiDistribution = new TH1D("Hj1j2dphiDistribution","",30,0.0,3.0);
   TH1D *Hjet1etaDistribution = new TH1D("Hjet1etaDistribution","",60,-3.0,3.0);
   TH1D *Hjet2etaDistribution = new TH1D("Hjet2etaDistribution","",60,-3.0,3.0);
   TH1D *HmetNoLepDistribution;
   TH1D *HzptDistribution;
   TH1D *Hjet1ptDistribution;  
   TH1D *Hjet2ptDistribution;
   
   TH1D *HinvMass_NoJetCuts = new TH1D("HinvMass_NoJetCuts","",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);  
   TH1D *HvtxDistribution_NoJetCuts = new TH1D("HvtxDistribution_NoJetCuts","",40,-0.5,39.5);
   TH1D *HnjetsDistribution_NoJetCuts = new TH1D("HnjetsDistribution_NoJetCuts","njets using nJetClean30",10,-0.5,9.5);
   TH1D *HmetNoLepDistribution_NoJetCuts;
   TH1D *HzptDistribution_NoJetCuts;
   TH1D *Hjet1ptDistribution_NoJetCuts;

   if (using_phys14_sample_flag) {
     HmetNoLepDistribution = new TH1D("HmetNoLepDistribution","",(Int_t)(1000.0-METNOLEP_START)/10.0,METNOLEP_START,1000.0);
     HzptDistribution = new TH1D("HzptDistribution","",200,0.0,1000.0);
     Hjet1ptDistribution = new TH1D("Hjet1ptDistribution","",97,30,1000); 
     Hjet2ptDistribution = new TH1D("Hjet2ptDistribution","",97,30,1000);
     HmetNoLepDistribution_NoJetCuts = new TH1D("HmetNoLepDistribution_NoJetCuts","",(Int_t)(1000.0-METNOLEP_START)/10.0,METNOLEP_START,1000.0);
     HzptDistribution_NoJetCuts = new TH1D("HzptDistribution_NoJetCuts","",200,0.0,1000.0);
     Hjet1ptDistribution_NoJetCuts = new TH1D("Hjet1ptDistribution_NoJetCuts","",97,30,1000);
   } else if (using_spring15_sample_flag) {
     HmetNoLepDistribution = new TH1D("HmetNoLepDistribution","",60,METNOLEP_START,METNOLEP_START+600);
     HzptDistribution = new TH1D("HzptDistribution","",80,0,400);
     Hjet1ptDistribution = new TH1D("Hjet1ptDistribution","",60,J1PT,J1PT+600); 
     Hjet2ptDistribution = new TH1D("Hjet2ptDistribution","",60,J2PT,J2PT+600);
     HmetNoLepDistribution_NoJetCuts = new TH1D("HmetNoLepDistribution_NoJetCuts","",60,METNOLEP_START,METNOLEP_START+600);
     HzptDistribution_NoJetCuts = new TH1D("HzptDistribution_NoJetCuts","",80,0,400);
     Hjet1ptDistribution_NoJetCuts = new TH1D("Hjet1ptDistribution_NoJetCuts","",60,J1PT,J1PT+600);
   }

   TH1D *HZtoLLRecoPt = new TH1D("HZtoLLRecoPt","",101,0.,1010);   // end at 1010 because I will put the overflow in the last bin
   // the previous histogram is differen from HzptDistribution because the binning is different
   TH1D *HZtoLLGenPt ;
   TH1D *HZtoLLPt_RecoGenRatio;                    // this is the histogram with reco/gen
   TH1D *HZtoLLPt_RecoGenRatio_pdf;             // histogram of reco/gen distribution function
   TH1D *HZtoLLPt_RecoGenRatio_pdf_ZpT600ToInf;   // when we have more data, it will be useful to see comparison between ee and mumu at high ZpT  

   if (!ISDATA_FLAG) {
     HZtoLLGenPt = new TH1D("HZtoLLGenPt","",101,0.,1010);
     HZtoLLPt_RecoGenRatio = new TH1D("HZtoLLPt_RecoGenRatio","",101,0.,1010.);
     HZtoLLPt_RecoGenRatio_pdf = new TH1D("HZtoLLPt_RecoGenRatio_pdf","",100,0.5,1.5);
     HZtoLLPt_RecoGenRatio_pdf_ZpT600ToInf = new TH1D("HZtoLLPt_RecoGenRatio_pdf_ZpT600ToInf","",100,0.5,1.5);  // will be useful when we get to this range in data
   }

   //TH1D *HinvMass[nMetBins];
   TH1D *HzlljetsInvMassMetBinGenLep[nMetBins];
   //TH1D *HzlljetsInvMassMetBinGenTau[nMetBins];
   TH1D *HZtoLLRecoPt_MetBin[nMetBins];
   TH1D *HZtoLLGenPt_MetBin[nMetBins];
   TH1D *HZtoLLPt_RecoGenRatio_MetBin[nMetBins];
   TH1D *HZtoLLPt_RecoGenRatio_pdf_MetBin[nMetBins];

   for (Int_t i = 0; i < nMetBins; i++) {

     //HinvMass[i] = new TH1D(Form("HinvMass_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);
     HzlljetsInvMassMetBinGenLep[i] = new TH1D(Form("HzlljetsInvMassMetBinGenLep_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);
     //HzlljetsInvMassMetBinGenTau[i] = new TH1D(Form("HzlljetsInvMassMetBinGenTau_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);
     HZtoLLRecoPt_MetBin[i] = new TH1D(Form("HZtoLLRecoPt_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",101,0.,1010.);

     if (!ISDATA_FLAG) {
       HZtoLLGenPt_MetBin[i] = new TH1D(Form("HZtoLLGenPt_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",101,0.,1010.);
       HZtoLLPt_RecoGenRatio_MetBin[i] = new TH1D(Form("HZtoLLPt_RecoGenRatio_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",101,0.,1010.);
       HZtoLLPt_RecoGenRatio_pdf_MetBin[i] = new TH1D(Form("HZtoLLPt_RecoGenRatio_pdf_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",100,0.5,1.5);
     }

   } 

   // following 2 histograms hold the inclusive distribution of parallel and ortogonal component of recoil (projected along the ZpT direction)
   TH1D *H_uPerp_Distribution = new TH1D("H_uPerp_Distribution","",50,-200,200);
   TH1D *H_uParMinusZpT_Distribution = new TH1D("H_uParMinusZpT_Distribution","",50,-200,200);

   TH1D *H_uPerp_VS_Nvtx[NVTXS];
   TH1D *H_uParMinusZpT_VS_Nvtx[NVTXS]; 
   TH1D *H_uParMinusZpT_VS_Nvtx_lowZpT[NVTXS];
 
   for (Int_t i = 0; i < NVTXS; i++) {

     H_uPerp_VS_Nvtx[i] = new TH1D(Form("H_uPerp_VS_Nvtx_nvtx%i",FIRST_NVTX+i),"",80,-200,200);  // 5 GeV bins 
     H_uParMinusZpT_VS_Nvtx[i] = new TH1D(Form("H_uParMinusZpT_VS_Nvtx_nvtx%i",FIRST_NVTX+i),"ZpT in [250,500]GeV",80,-200,200);  // 5 GeV bins
     H_uParMinusZpT_VS_Nvtx_lowZpT[i] = new TH1D(Form("H_uParMinusZpT_VS_Nvtx_lowZpT_nvtx%i",FIRST_NVTX+i),"ZpT in [50,250]GeV",80,-200,200);  // 5 GeV bins

   }

   // zpt bin edges for respose studies (which are independent on the CS selection, they are done on the whole sample)
   //Double_t ZptBinEdges[] = {250., 260., 270., 280., 290., 310., 330., 350., 370., 390., 410., 430., 450., 470., 500., 530., 560, 600., 640., 700., 800.};
   //Double_t ZptBinEdges[] = {250., 260., 270., 280., 290., 310., 330., 350., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   //Double_t ZptBinEdges[] = {10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 310., 330., 350., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   //Double_t ZptBinEdges[] = {20., 40., 60., 80., 100., 120., 140., 160., 180., 200., 220., 240., 260., 280., 300., 320., 340., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   // Double_t ZptBinEdgesMC[] = {1., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 220., 240., 260., 280., 300., 320., 340., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   // Double_t ZptBinEdgesMC[] = {1., 10., 20., 30., 40., 60., 80., 100., 120., 140., 180., 220., 260., 300., 340., 380.};
   //Double_t ZptBinEdgesDATA[] = {1., 10., 20., 30., 40., 60., 80., 100., 120., 140., 180., 220.};
   

   Double_t *ZptBinEdges = NULL;
   Int_t nBinsForResponse = 0;   // # of bins for analysis as a function of ZpT
   Double_t ZptBinEdges_Spring15[] = {1., 10., 20., 30., 40., 60., 80., 100., 120., 140., 180., 220.}; 
   Double_t ZptBinEdges_Phys14[] = {250., 260., 270., 280., 290., 310., 330., 350., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   Double_t ZptBinEdges_Phys14_noSkim[] = {1., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 310., 330., 350., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   //Int_t nBinsForResponse = sizeof(ZptBinEdges)/sizeof(ZptBinEdges[0]) - 1;  //number of bins is n-1 where n is the number of ZptBinEdges's elements
   //Int_t nBinsForResponse_0jets = 0;  // for the response curve in events with nJetClean30 = 0

   // It's better to have only one binning, so that it's possible to make ratios (t's not relevant that the mean ZpT is the same for data and MC, but it's importanto to have 
   // corresponding bins). data binedges must be a subset of those for MC

   if (using_spring15_sample_flag) {

     ZptBinEdges = ZptBinEdges_Spring15;
     nBinsForResponse = sizeof(ZptBinEdges_Spring15)/sizeof(ZptBinEdges_Spring15[0]) - 1;  //number of bins is n-1 where n is the number of ZptBinEdges's elements
     //nBinsForResponse_0jets = 6; //use first bins, up to 60 GeV

   } else if (using_phys14_sample_flag) {

     if (METNOLEP_START == 0) {
       ZptBinEdges = ZptBinEdges_Phys14_noSkim;
       nBinsForResponse = sizeof(ZptBinEdges_Phys14_noSkim)/sizeof(ZptBinEdges_Phys14_noSkim[0]) - 1;  //number of bins is n-1 where n is the number of ZptBinEdges's elements
     } else {
       ZptBinEdges = ZptBinEdges_Phys14;
       nBinsForResponse = sizeof(ZptBinEdges_Phys14)/sizeof(ZptBinEdges_Phys14[0]) - 1;  //number of bins is n-1 where n is the number of ZptBinEdges's elements
     }

   }

   TH1D *H_uPerp_VS_ZpT[nBinsForResponse];  
   TH1D *H_uParMinusZpT_VS_ZpT[nBinsForResponse];    // actually it will be (u_par-ZpT)
   TH1D *H_uPar_ZpT_ratio[nBinsForResponse];         // for the response curve: we will compute response in many ways. Here we use <uPar/ZpT>
   TH1D *H_uParMinusZpT_ZpT_ratio[nBinsForResponse]; // here we use <(uPar-ZpT)/ZpT> and will add back 1, so we get the same as above
   TH1D *HZptBinned[nBinsForResponse];
   
   //TH1D *H_uPar_ZpT_ratio_0jets[nBinsForResponse_0jets];  // for the response curve in events with nJetClean30 = 0

   //the following histograms will give the distribution of met|| / wzpt. The mean value will be used to create the response curve, that is (<met|| / wzpt>) vs wzpt
   // for each point, wzpt will be taken as the average wzpt in the range considered
 
   for (Int_t i = 0; i < nBinsForResponse; i++) {   

     //HZptBinned[i] are histograms with 5 bins in the range given by ZptBinEdges[i] and ZptBinEdges[i+1]
     // the mean wzpt in each bin will be computed as the histogram's mean
     HZptBinned[i] = new TH1D(Form("HZptBinned_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",5,ZptBinEdges[i],ZptBinEdges[i+1]); 
     // in the following histogram , range must include negative value: I saw that for low ZoT this distribution tends to be flat, thus if range goes from 0 to 2 (as it was before) the 
     //mean for ZpT tending to 0 will be 1 and not 0 as we would expect.
     H_uPar_ZpT_ratio[i] = new TH1D(Form("H_uPar_ZpT_ratio_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",350,-7.0,7.0); 
     H_uParMinusZpT_ZpT_ratio[i] = new TH1D(Form("H_uParMinusZpT_ZpT_ratio_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",350,-7.0,7.0);
     H_uPerp_VS_ZpT[i] = new TH1D(Form("H_uPerp_VS_ZpT_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",40,-200,200); 
     H_uParMinusZpT_VS_ZpT[i] = new TH1D(Form("H_uParMinusZpT_VS_ZpT_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",40,-200,200); 
     //if ( i < nBinsForResponse_0jets) H_uPar_ZpT_ratio_0jets[i] = new TH1D(Form("H_uPar_ZpT_ratio_0jets_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",350,-7.0,7.0); 

   }

   // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
   TH1D *HmetBinEdges = new TH1D("HmetBinEdges","bin edges for met distributions",nMetBins+1,0.0,nMetBins+1);
   for (Int_t i = 0; i <= nMetBins; i++) {
     HmetBinEdges->SetBinContent(i+1,metBinEdges[i]);
   }

   TH1D *HZptBinEdges = new TH1D("HZptBinEdges","bin edges for ZpT distributions",nBinsForResponse+1,0.0,nBinsForResponse+1);
   for (Int_t i = 0; i <= nBinsForResponse; i++) {
     HZptBinEdges->SetBinContent(i+1,ZptBinEdges[i]);
   }

   TH1D *HnvtxBins = new TH1D("HnvtxBins","# of vertices for studies of variables as a function of nvtx",NVTXS,0.0,NVTXS);
   for (Int_t i = 0; i < NVTXS; i++) {              // watch out: differently from above, i < NVTXS, not <=, because if NVTXS = 3 I need 3 points, not 4
     HnvtxBins->SetBinContent(i+1,FIRST_NVTX+i);
   }

   // deciding  what is the event weight
   Double_t newwgt;

   if (ISDATA_FLAG || unweighted_event_flag) newwgt = 1.0;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"zlljets_resoResp::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;

   for (Int_t jentry=0; jentry<nentries; jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     if (jentry%500000 == 0) cout << jentry << endl;

     UInt_t eventMask = 0; 

     if(!ISDATA_FLAG && !unweighted_event_flag) {

       if (using_spring15_sample_flag) newwgt = 1000 * LUMI * vtxWeight  * xsec * genWeight / SUMWEIGHTS;  
       //if (using_spring15_sample_flag) newwgt = 1000 * LUMI * vtxW  * xsec * genWeight / SUMWEIGHTS;    // 1000 is because LUMI is in fb^-1 and xsec is in pb
       // old wrong one:     newwgt = LUMI * vtxW * weight * LHEorigWeight; 
       else if (using_phys14_sample_flag) newwgt = LUMI * weight;   // for older trees (backward compatibility)
       else newwgt = LUMI * weight;   // for older trees (backward compatibility)

       

     }

     nTotalWeightedEvents += newwgt;  // counting events with weights

     nLepLoose = *ptr_nLepLoose;          
     nLep10V = *ptr_nLep10V;

     //Double_t ZgenMass;        // not used for now
     Double_t ZtoLLGenPt = 0;    // filled below (only if running on MC DYJetsToLL)
     Double_t ZtoLLRecoPt = 0;   // filled below

     // genLepFound_flag is used when analysing DYJetsToLL in MC fo Z->mumu or Z->ee. For other MC samples it's not used.

     if (!ISDATA_FLAG) {

       if (using_zlljets_MCsample_flag) {

	 genLepFound_flag = myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, LEP_PDG_ID, 23, firstIndexGen, secondIndexGen, Z_index, GenPart_motherIndex); 
	 //if (!genLepFound_flag) continue;  // if not found gen ee or mumu for MC DYJetsToLL ( l = mu or e) skip the event. This makes things faster
	 if (genLepFound_flag) {
	   eventMask += genLepC.addToMask( genLepFound_flag );
	   l1gen.SetPtEtaPhiM(GenPart_pt[firstIndexGen],GenPart_eta[firstIndexGen],GenPart_phi[firstIndexGen],GenPart_mass[firstIndexGen]);
	   l2gen.SetPtEtaPhiM(GenPart_pt[secondIndexGen],GenPart_eta[secondIndexGen],GenPart_phi[secondIndexGen],GenPart_mass[secondIndexGen]);
	   Zgen = l1gen + l2gen;
	   ZtoLLGenPt = Zgen.Pt();                             
	 }

       } else if (using_ztautaujets_MCsample_flag) {

	 genTauFound_flag = myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, 15, 23);
	 //if (!genTauFound_flag) continue;  // if not found gen tautau for MC DYJetsToLL ( l = tau) skip the event. This makes things faster
	 eventMask += genTauC.addToMask( genTauFound_flag );

       }

     }

     recoLepFound_flag = myGetPairIndexInArray(LEP_PDG_ID, nLepGood, LepGood_pdgId, firstIndex, secondIndex);  

     if (recoLepFound_flag) {
       l1reco.SetPtEtaPhiM(LepGood_pt[firstIndex],LepGood_eta[firstIndex],LepGood_phi[firstIndex],LepGood_mass[firstIndex]);
       l2reco.SetPtEtaPhiM(LepGood_pt[secondIndex],LepGood_eta[secondIndex],LepGood_phi[secondIndex],LepGood_mass[secondIndex]);
       Zreco = l1reco + l2reco;
       ZtoLLRecoPt = Zreco.Pt();
     }

     if (fabs(LEP_PDG_ID) == 13) { 

       if ( HLT_FLAG ) {

	 // use the dimuon trigger, not the metNoLep trigger
       	 if ( recoLepFound_flag && (fabs(LepGood_eta[firstIndex]) < HLT_LEP1ETA) && (fabs(LepGood_eta[secondIndex]) < HLT_LEP2ETA) && 
       	      (LepGood_pt[firstIndex] > HLT_LEP1PT) && (LepGood_pt[secondIndex] > HLT_LEP2PT) ) HLT_passed_flag = 1; 	 
       	 else HLT_passed_flag = 0; //continue;

       }  // end of   if ( HLT_FLAG )

       metNoLepPt = *ptr_metNoLepPt;       
       //metNoLepEta = *ptr_metNoLepEta; 
       metNoLepPhi = *ptr_metNoLepPhi; 
       //metNoLepTV3.SetPtEtaPhi(metNoLepPt,metNoLepEta,metNoLepPhi);   // will use this 3D vector below
       metNoLepTV.SetMagPhi(metNoLepPt,metNoLepPhi);

     } else if (fabs(LEP_PDG_ID) == 11) { 

       if ( HLT_FLAG ) {

       	 if ( recoLepFound_flag && (LepGood_tightId[firstIndex] > 0.5) && (LepGood_tightId[secondIndex]  > 0.5) && 
       	      (fabs(LepGood_eta[firstIndex]) < HLT_LEP1ETA) && (fabs(LepGood_eta[secondIndex]) < HLT_LEP2ETA) && 
       	      (LepGood_pt[firstIndex] > HLT_LEP1PT) && (LepGood_pt[secondIndex] > HLT_LEP2PT) ) HLT_passed_flag = 1; 	 
	 else HLT_passed_flag = 0;  //continue;

       }  // end of   if ( HLT_FLAG )

       // metNoLepTV3.SetPtEtaPhi(met_pt,met_eta,met_phi);
       // // summing just electrons from Z if found
       // ele.SetPtEtaPhi(LepGood_pt[firstIndex],LepGood_eta[firstIndex],LepGood_phi[firstIndex]);
       // metNoLepTV3 += ele;
       // ele.SetPtEtaPhi(LepGood_pt[secondIndex],LepGood_eta[secondIndex],LepGood_phi[secondIndex]);
       // metNoLepTV3 += ele;

       metNoLepTV.SetMagPhi(met_pt,met_phi);
       // summing just electrons from Z if found
       if (recoLepFound_flag) {
	 ele.SetMagPhi(LepGood_pt[firstIndex],LepGood_phi[firstIndex]);
	 metNoLepTV += ele;
	 ele.SetMagPhi(LepGood_pt[secondIndex],LepGood_phi[secondIndex]);
	 metNoLepTV += ele;
       }

       metNoLepPt = metNoLepTV.Mod();

     }

     // beginning of eventMask building

     
     // genLepC added to mask above if ISDATA_FLAG == false (in order not to repeat here the check) 
     
     eventMask += jet1C.addToMask(nJetClean30 >= 1 && JetClean_pt[0] > J1PT && fabs(JetClean_eta[0] < J1ETA && jetclean1 > 0.5));  //could skip cut on eta 
     eventMask += jjdphiC.addToMask( nJetClean30 == 1 || (nJetClean30 >= NJETS && fabs(dphijj) < J1J2DPHI && jetclean2 > 0.5));
     eventMask += njetsC.addToMask(nJetClean30 <= NJETS);
     eventMask += lepLooseVetoC.addToMask(nLep10V == 0);
     eventMask += tauLooseVetoC.addToMask(nTauClean18V == 0);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     eventMask += metNoLepStartC.addToMask(metNoLepPt > METNOLEP_START);

     // for (Int_t i = 0; i <  metCut.size(); i++) {
     //   eventMask += metNoLepC[i].addToMask(metNoLepPt > metCut[i]);
     // }
     
     // the following make sense only if recoLepFound_flag == 1 (i.e. flag is true), which means that fabs(LepGood_pdgId[firstIndex/secondIndex]) == LEP_PDG_ID) is 
     // true
     // also, 2 OS/SF leptons are present
     if (recoLepFound_flag) {

       eventMask += HLTlepC.addToMask(HLT_passed_flag);     
       eventMask += oppChargeLeptonsC.addToMask( 1);
       eventMask += twoLeptonsC.addToMask(1);
       eventMask += twoLepLooseC.addToMask(nLepLoose == 2);
       eventMask += lep1ptC.addToMask((LepGood_pt[firstIndex] > LEP1PT)); 
       eventMask += lep1etaC.addToMask( (fabs(LepGood_eta[firstIndex]) < LEP1ETA) );
       eventMask += lep2ptC.addToMask((LepGood_pt[secondIndex] > LEP2PT) );
       eventMask += lep2etaC.addToMask((fabs(LepGood_eta[secondIndex]) < LEP2ETA) );
       eventMask += invMassC.addToMask((mZ1 > DILEPMASS_LOW) && (mZ1 < DILEPMASS_UP));     
       eventMask += lep1tightIdIso04C.addToMask((LepGood_tightId[firstIndex] > 0.5 ) && (LepGood_relIso04[firstIndex] < LEP_ISO_04 ) );
       eventMask += lep2tightIdIso04C.addToMask((LepGood_tightId[secondIndex] > 0.5) && (LepGood_relIso04[secondIndex] < LEP_ISO_04 ) );
      
     }

     // end of eventMask building

     // test matching of reco and gen lep for DY MC 
     if (!ISDATA_FLAG && using_zlljets_MCsample_flag) {

       //enter this part if 2 OS/SF leptons were found among gen and reco particles. Now checking compatibilities between pairs
       // e.g. l1gen = e+, l2gen = e- ; l1reco = e+, l2reco = e- (but the charge order might not coincide)
       // now we require a DeltaR cut between them to assess that lreco comes from lgen
       // since 2 OS/SF were found to get inside here, if !(l1gen->l1reco && l2gen->l2reco) then for sure l1gen->l2reco && l2gen->l1reco
       
       if (genLepFound_flag && recoLepFound_flag) {       

	 Double_t DeltaR_lreco_lgen_pair1 = 0.0;
	 Double_t DeltaR_lreco_lgen_pair2 = 0.0;
       
	 if(LepGood_pdgId[firstIndex] == GenPart_pdgId[firstIndexGen] && LepGood_pdgId[secondIndex] == GenPart_pdgId[secondIndexGen]) {
	 
	   DeltaR_lreco_lgen_pair1 = l1reco.DeltaR(l1gen);
	   DeltaR_lreco_lgen_pair2 = l2reco.DeltaR(l2gen);

	 } else {
	 
	   DeltaR_lreco_lgen_pair1 = l1reco.DeltaR(l2gen);
	   DeltaR_lreco_lgen_pair2 = l2reco.DeltaR(l1gen);
	 
	 }
       
	 if (DeltaR_lreco_lgen_pair1 < 0.1 && DeltaR_lreco_lgen_pair2 < 0.1) eventMask += recoGenLepMatchC.addToMask(1);
	 else eventMask += recoGenLepMatchC.addToMask(0);

       }

     }

     zlljetsControlSample.countEvents(eventMask,newwgt);
     resoAndResponse.countEvents(eventMask, newwgt);

     if ( ((eventMask & zlljetsControlSample.globalMask.back()) == zlljetsControlSample.globalMask.back()) ) {
       
       // this histogram holds the final yields in bins of MET
	 HzlljetsYieldsMetBin->Fill(metNoLepPt,newwgt);
	 
	 HinvMass->Fill(mZ1,newwgt);
	 HmetNoLepDistribution->Fill(metNoLepPt,newwgt);
	 HzptDistribution->Fill(ZtoLLRecoPt,newwgt);
	 HvtxDistribution->Fill(nVert,newwgt);
	 HnjetsDistribution->Fill(nJetClean30,newwgt);
	 Hjet1etaDistribution->Fill(JetClean_eta[0],newwgt);
	 Hjet1ptDistribution->Fill(JetClean_pt[0],newwgt);
	 if (nJetClean30 == 2) {
	   Hj1j2dphiDistribution->Fill(dphijj,newwgt);
	   Hjet2etaDistribution->Fill(JetClean_eta[1],newwgt);
	   Hjet2ptDistribution->Fill(JetClean_pt[1],newwgt);
	 }

     }
     
     // if ( ((eventMask & tautaubkgInZll.globalMask.back()) == tautaubkgInZll.globalMask.back()) ) {  
     // 	 // this histogram holds the final yields in bins of MET
     // 	 HzlljetsYieldsMetBinGenTau->Fill(metNoLepPt,newwgt);  
     // }
     
     if ( ((eventMask & resoAndResponse.globalMask.back()) == resoAndResponse.globalMask.back()) ) {  
       
       if (!ISDATA_FLAG && using_zlljets_MCsample_flag) {

	 HZtoLLGenPt->Fill(ZtoLLGenPt,newwgt);
	 if (ZtoLLGenPt != 0) {

	   HZtoLLPt_RecoGenRatio_pdf->Fill(ZtoLLRecoPt/ZtoLLGenPt,newwgt);
	   if (ZtoLLRecoPt > 600) HZtoLLPt_RecoGenRatio_pdf_ZpT600ToInf->Fill(ZtoLLRecoPt/ZtoLLGenPt,newwgt);

	 }

       }

       HZtoLLRecoPt->Fill(ZtoLLRecoPt,newwgt);	 

       HinvMass_NoJetCuts->Fill(mZ1,newwgt);
       HmetNoLepDistribution_NoJetCuts->Fill(metNoLepPt,newwgt);
       HzptDistribution_NoJetCuts->Fill(ZtoLLRecoPt,newwgt);
       HvtxDistribution_NoJetCuts->Fill(nVert,newwgt);
       HnjetsDistribution_NoJetCuts->Fill(nJetClean30,newwgt);      
       if (nJetClean30 > 0) Hjet1ptDistribution_NoJetCuts->Fill(JetClean_pt[0],newwgt);

       // following is done if two OS leptons are found (otherwise there would be no Z)
		
       //metNoLepTV3.SetPtEtaPhi((Double_t)metNoLepPt,(Double_t)metNoLepEta,(Double_t)metNoLepPhi);  // already initialized above
       //Double_t dphiMetNoLepZ = metNoLepTV3.DeltaPhi(Zreco.Vect());
       TVector3 Zreco3D = Zreco.Vect();
       Double_t dphiMetNoLepZ = metNoLepTV.DeltaPhi(Zreco3D.XYvector());

       Double_t u_par = metNoLepPt * TMath::Cos(dphiMetNoLepZ);  // actually u_par is minus this quantity, but then I do u_par-ZpT instead of u_par+ZpT
       Double_t u_perp = metNoLepPt * TMath::Sin(dphiMetNoLepZ);
       Double_t uparMinusZrecoPt = u_par - ZtoLLRecoPt;

       H_uPerp_Distribution->Fill(u_perp,newwgt);
       H_uParMinusZpT_Distribution->Fill(uparMinusZrecoPt,newwgt);

       if (ZtoLLRecoPt > ZptBinEdges[0]) {  

	 Int_t nvtxBin = nVert - FIRST_NVTX;
	 Int_t lastnvtx = NVTXS + FIRST_NVTX;

	 if ((nvtxBin >= 0) && (nVert < lastnvtx)) {

	   H_uPerp_VS_Nvtx[nvtxBin]->Fill(u_perp,newwgt);
	     
	   if (ZtoLLRecoPt < 250 ) {

	     H_uParMinusZpT_VS_Nvtx_lowZpT[nvtxBin]->Fill(uparMinusZrecoPt,newwgt);
	     
	   } else if (ZtoLLRecoPt < 500) {                       // (met||-wzpt) distribution's width depends on Zpt, thus I use this range

	     H_uParMinusZpT_VS_Nvtx[nvtxBin]->Fill(uparMinusZrecoPt,newwgt);
	 
	   }       
    
	 }  // end of   if ((nvtxBin >= 0) && (nVert < lastnvtx))

	 /**************************************************/
	 // computing met responses
	 /**************************************************/

	 // first of all I make sure that wzpt is in the appropriate range
	 if ( ZtoLLRecoPt < ZptBinEdges[nBinsForResponse] ) {

	   Int_t respBin = myGetBin(ZtoLLRecoPt,ZptBinEdges,nBinsForResponse);
	   //cout<<"bin = "<<bin<<endl;
	   HZptBinned[respBin]->Fill(ZtoLLRecoPt,newwgt);        
	   H_uPar_ZpT_ratio[respBin]->Fill(u_par/ZtoLLRecoPt,newwgt);          //the mean value of this histogram is the response
	   H_uParMinusZpT_ZpT_ratio[respBin]->Fill(uparMinusZrecoPt/ZtoLLRecoPt,newwgt);  //the mean value of this histogram +1 is the response
	   H_uPerp_VS_ZpT[respBin]->Fill(u_perp,newwgt);
	   H_uParMinusZpT_VS_ZpT[respBin]->Fill(uparMinusZrecoPt,newwgt);
	   //if (ZtoLLRecoPt < ZptBinEdges[nBinsForResponse_0jets]) H_uPar_ZpT_ratio_0jets[respBin]->Fill(u_par/ZtoLLRecoPt,newwgt);

	 }

       }            // end of if (ZtoLLRecoPt > ZptBinEdges[0])

     }
	

     // now entering analysis in bins of met

     if ((metNoLepPt > metBinEdges[0]) && (metNoLepPt < metBinEdges[nMetBins])) {

       Int_t bin = myGetBin(metNoLepPt,metBinEdges,nMetBins);
       
       // if ((eventMask & zlljetsControlSample.globalMask.back()) == zlljetsControlSample.globalMask.back()) {
       //   // this histogram holds the invariant mass distribution (one for each met bin)
       //   HinvMass[bin]->Fill(mZ1,newwgt);   
       // }

       if ( ((eventMask & zlljetsControlSample.globalMask.back()) == zlljetsControlSample.globalMask.back()) ) { 
 
	 HzlljetsInvMassMetBinGenLep[bin]->Fill(mZ1,newwgt); 

       }

       if (((eventMask & resoAndResponse.globalMask.back()) == resoAndResponse.globalMask.back())) {

	 if (!ISDATA_FLAG && using_zlljets_MCsample_flag) {

	   HZtoLLRecoPt_MetBin[bin]->Fill(ZtoLLRecoPt,newwgt);
	   HZtoLLGenPt_MetBin[bin]->Fill(ZtoLLGenPt,newwgt);
	   if (ZtoLLGenPt != 0) HZtoLLPt_RecoGenRatio_pdf_MetBin[bin]->Fill(ZtoLLRecoPt/ZtoLLGenPt,newwgt);

	 } else HZtoLLRecoPt_MetBin[bin]->Fill(ZtoLLRecoPt,newwgt);  // if running on data just do this

       } 
	 
     }                      // end of    if ((metNoLepPt > metBinEdges[0]) && (metNoLepPt < metBinEdges[nMetBins])) 
       
   }                        // end of loop on entries

   
   /************************************/
   //                    MET|| & MET_|_ VS NVTX & ZpT
   /************************************/

   //resolution vs nvtx

   Double_t xValues[NVTXS];
   Double_t yValues[NVTXS];
   Double_t yValuesErr[NVTXS];

   for (Int_t i = 0; i < NVTXS; i++) {
     xValues[i] = i + FIRST_NVTX;
     yValues[i] = H_uParMinusZpT_VS_Nvtx_lowZpT[i]->GetRMS();
     yValuesErr[i] = H_uParMinusZpT_VS_Nvtx_lowZpT[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepParZvsNvtx_lowZpT = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMetNoLepParZvsNvtx_lowZpT->SetTitle(Form("resolution || from histogram's RMS, ZpT in [%2.0lf,250] GeV",ZptBinEdges[0]));
   GresolutionMetNoLepParZvsNvtx_lowZpT->Draw("AP");
   GresolutionMetNoLepParZvsNvtx_lowZpT->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsNvtx_lowZpT->GetXaxis()->SetTitle("nvtx");
   GresolutionMetNoLepParZvsNvtx_lowZpT->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsNvtx_lowZpT->GetYaxis()->SetTitleOffset(1.4); 
   GresolutionMetNoLepParZvsNvtx_lowZpT->SetName("gr_resolution_uPar_vs_Nvtx_lowZpT");
   GresolutionMetNoLepParZvsNvtx_lowZpT->Write();

   for (Int_t i = 0; i < NVTXS; i++) {
     yValues[i] = H_uParMinusZpT_VS_Nvtx[i]->GetRMS();
     yValuesErr[i] = H_uParMinusZpT_VS_Nvtx[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepParZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMetNoLepParZvsNvtx->SetTitle("resolution || from histogram's RMS, ZpT in [250,500] GeV");
   GresolutionMetNoLepParZvsNvtx->Draw("AP");
   GresolutionMetNoLepParZvsNvtx->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsNvtx->GetXaxis()->SetTitle("nvtx");
   GresolutionMetNoLepParZvsNvtx->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   GresolutionMetNoLepParZvsNvtx->SetName("gr_resolution_uPar_vs_Nvtx");
   GresolutionMetNoLepParZvsNvtx->Write();

   for (Int_t i = 0; i < NVTXS; i++) {
     yValues[i] = H_uPerp_VS_Nvtx[i]->GetRMS();
     yValuesErr[i] = H_uPerp_VS_Nvtx[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepOrtZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMetNoLepOrtZvsNvtx->SetTitle("resolution _|_ from histogram's RMS");
   GresolutionMetNoLepOrtZvsNvtx->Draw("AP");
   GresolutionMetNoLepOrtZvsNvtx->SetMarkerStyle(7);
   GresolutionMetNoLepOrtZvsNvtx->GetXaxis()->SetTitle("nvtx");
   GresolutionMetNoLepOrtZvsNvtx->GetYaxis()->SetTitle("#sigma (u_#perp ) [GeV]");
   GresolutionMetNoLepOrtZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   GresolutionMetNoLepOrtZvsNvtx->SetName("gr_resolution_uPerp_vs_Nvtx");
   GresolutionMetNoLepOrtZvsNvtx->Write();

   // resolution vs ZpT
   
   Double_t meanZpt[nBinsForResponse];
   Double_t meanZptErr[nBinsForResponse];
   Double_t resoMetNoLepParZvsZpt[nBinsForResponse];
   Double_t resoMetNoLepParZvsZptErr[nBinsForResponse];
   Double_t resoMetNoLepOrtZvsZpt[nBinsForResponse];
   Double_t resoMetNoLepOrtZvsZptErr[nBinsForResponse];

   for (Int_t i = 0; i < nBinsForResponse; i++) {

     meanZpt[i] = HZptBinned[i]->GetMean();
     meanZptErr[i] = HZptBinned[i]->GetMeanError();
     resoMetNoLepParZvsZpt[i] = H_uParMinusZpT_VS_ZpT[i]->GetRMS();
     resoMetNoLepParZvsZptErr[i] = H_uParMinusZpT_VS_ZpT[i]->GetRMSError();
     resoMetNoLepOrtZvsZpt[i] = H_uPerp_VS_ZpT[i]->GetRMS();
     resoMetNoLepOrtZvsZptErr[i] = H_uPerp_VS_ZpT[i]->GetRMSError();

   }

   TGraphErrors *GresolutionMetNoLepParZvsZpt = new TGraphErrors(nBinsForResponse,meanZpt,resoMetNoLepParZvsZpt,meanZptErr,resoMetNoLepParZvsZptErr);
   GresolutionMetNoLepParZvsZpt->SetTitle("resolution || from histogram's RMS");
   GresolutionMetNoLepParZvsZpt->Draw("AP");
   GresolutionMetNoLepParZvsZpt->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsZpt->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMetNoLepParZvsZpt->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsZpt->GetYaxis()->SetTitleOffset(1.2); 
   GresolutionMetNoLepParZvsZpt->SetName("gr_resolution_uPar_vs_ZpT");
   GresolutionMetNoLepParZvsZpt->Write();

   TGraphErrors *GresolutionMetNoLepOrtZvsZpt = new TGraphErrors(nBinsForResponse,meanZpt,resoMetNoLepOrtZvsZpt,meanZptErr,resoMetNoLepOrtZvsZptErr);
   GresolutionMetNoLepOrtZvsZpt->SetTitle("resolution _|_ from histogram's RMS");
   GresolutionMetNoLepOrtZvsZpt->Draw("AP");
   GresolutionMetNoLepOrtZvsZpt->SetMarkerStyle(7);
   GresolutionMetNoLepOrtZvsZpt->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMetNoLepOrtZvsZpt->GetYaxis()->SetTitle("#sigma (u_#perp ) [GeV]");
   GresolutionMetNoLepOrtZvsZpt->GetYaxis()->SetTitleOffset(1.2); 
   GresolutionMetNoLepOrtZvsZpt->SetName("gr_resolution_uPerp_vs_ZpT");
   GresolutionMetNoLepOrtZvsZpt->Write();

   // response curve: it can be defined in many ways

   Double_t response[nBinsForResponse];              // <uPar/ZpT>  this is somewhat biased
   Double_t responseErr[nBinsForResponse];

   Double_t response_gausFit[nBinsForResponse];      // <Upar-ZpT>/<ZpT>  + 1: numerator from a gaussian fit (we use gaussian mean since we look for the peak)   
   Double_t responseErr_gausFit[nBinsForResponse];
   TFitResultPtr ptrGausFit;
   Double_t meanUparMinusZpt_gausFit[nBinsForResponse];       // will contain <u_par - ZpT> taken from a gaussian fit
   Double_t meanUparMinusZptErr_gausFit[nBinsForResponse];    // will contain uncertainty on <u_par - ZpT> taken from a gaussian fit

   Double_t response_gausFit_bis[nBinsForResponse];    // <Upar-ZpT/ZpT>  + 1: mean value from a gaussian fit (we use gaussian mean since we look for the peak)   
   Double_t responseErr_gausFit_bis[nBinsForResponse];
   TFitResultPtr ptrGausFit_bis;
   Double_t mean_UparMinusZpt_ZpT_ratio_gausFit_bis[nBinsForResponse];       // will contain <u_par - ZpT> taken from a gaussian fit
   Double_t mean_UparMinusZpt_ZpT_ratioErr_gausFit_bis[nBinsForResponse];    
  
   // using <A/B> is different from using <A>/<B>. In general <f(x,y)> != f(<x>,<y>). The two relation coincide if the variance of x and y can be neglected or
   // the second partial derivatives of f wrt x and y are small enough (note that here is also uPar = uPar(ZpT) so that we actually have 1 independent variable)

   for (Int_t i = 0; i < nBinsForResponse; i++) {    

     response[i] = H_uPar_ZpT_ratio[i]->GetMean();
     responseErr[i] = H_uPar_ZpT_ratio[i]->GetMeanError();
     //cout<<i<<" meanZpt = "<<meanZpt[i]<<" +/- "<<meanZptErr[i]<<"    response = "<<response[i]<<" +/- "<<responseErr[i]<<endl;

     // using "(<u_par -ZpT> / <ZpT>) + 1" to compute response (adding 1 because that ratio is centered around 0). <.> is the mean value
     // do the fit between -3.5 RMS and + 3.5 RMS (found to be good enough) with a gaussian to get mean
     // option WL use loglikelihood fit with weighted histograms (better if there are empty bins)
     // Q is quiet mode (minimum printing on stdout). V prints everything. Default option is between Q and V
     // S is necessary to pass object and access to fit parameter

     Double_t tmpRMS;                                                   // temporary variable with distribution's RMS
     if (H_uParMinusZpT_VS_ZpT[i]->GetEntries() <= 10 ) {                // if empty histogram (including underflows and overflows) no fit is done
       meanUparMinusZpt_gausFit[i] = H_uParMinusZpT_VS_ZpT[i]->GetMean();                               // actually the fit has no sense with few points
       meanUparMinusZptErr_gausFit[i] = H_uParMinusZpT_VS_ZpT[i]->GetMeanError();
     } else {
       tmpRMS = H_uParMinusZpT_VS_ZpT[i]->GetRMS();
       ptrGausFit = H_uParMinusZpT_VS_ZpT[i]->Fit("gaus","Q S","",-3.5*tmpRMS,3.5*tmpRMS);
       meanUparMinusZpt_gausFit[i] = ptrGausFit->Parameter(1);      // 1 is the mean (0 and 2 are normalization and sigma of gaussian)
       meanUparMinusZptErr_gausFit[i] = ptrGausFit->ParError(1);
     }
     response_gausFit[i] = (meanUparMinusZpt_gausFit[i] / meanZpt[i]);    // in a second moment, adding 1 to response, otherwise it would be centered around 0
     responseErr_gausFit[i] = response_gausFit[i] * sqrt( (meanUparMinusZptErr_gausFit[i]*meanUparMinusZptErr_gausFit[i] / (meanUparMinusZpt_gausFit[i]*meanUparMinusZpt_gausFit[i])) + (meanZptErr[i]*meanZptErr[i] / (meanZpt[i]*meanZpt[i])));     // for the uncertainty, using response BEFORE adding 1 
     response_gausFit[i] += 1.;
     // uncertainty computed assuming uncorrelated numerator and denominator

     // now using "(<u_par -ZpT / ZpT>) + 1" to compute response (adding 1 because that ratio is centered around 0)

     if ( H_uParMinusZpT_ZpT_ratio[i]->GetEntries() <= 10 ) {                 // if empty histogram (including underflows and overflows) no fit is done
       mean_UparMinusZpt_ZpT_ratio_gausFit_bis[i] = H_uParMinusZpT_ZpT_ratio[i]->GetMean();                     // actually the fit has no sense with few points
       mean_UparMinusZpt_ZpT_ratioErr_gausFit_bis[i] = H_uParMinusZpT_ZpT_ratio[i]->GetMeanError();
     } else {
       tmpRMS = H_uParMinusZpT_ZpT_ratio[i]->GetRMS();
       ptrGausFit_bis = H_uParMinusZpT_ZpT_ratio[i]->Fit("gaus","Q S","",-3.5*tmpRMS,3.5*tmpRMS);
       mean_UparMinusZpt_ZpT_ratio_gausFit_bis[i] =  ptrGausFit_bis->Parameter(1);  // 1 is the mean (0 and 2 are normalization and sigma of gaussian)
       mean_UparMinusZpt_ZpT_ratioErr_gausFit_bis[i] =  ptrGausFit_bis->ParError(1);
     }     
     response_gausFit_bis[i] = mean_UparMinusZpt_ZpT_ratio_gausFit_bis[i] + 1.;  
     responseErr_gausFit_bis[i] = mean_UparMinusZpt_ZpT_ratioErr_gausFit_bis[i]; // for the uncertainty, using response BEFORE adding 1 

   }

   TGraphErrors *GresponseCurve = new TGraphErrors(nBinsForResponse,meanZpt,response,meanZptErr,responseErr);
   GresponseCurve->SetTitle("response curve");
   GresponseCurve->Draw("AP");
   GresponseCurve->SetMarkerStyle(7);    // 7 is a medium dot
   GresponseCurve->GetXaxis()->SetTitle("ZpT [GeV]");
   GresponseCurve->GetYaxis()->SetTitle(" < u_{||} / ZpT >");
   GresponseCurve->GetYaxis()->SetRangeUser(0.0, 1.1);
   GresponseCurve->GetYaxis()->SetTitleOffset(1.4); 
   GresponseCurve->SetName("gr_responseCurve");
   GresponseCurve->Write();

   TGraphErrors *GresponseCurve_gausFit = new TGraphErrors(nBinsForResponse,meanZpt,response_gausFit,meanZptErr,responseErr_gausFit);
   GresponseCurve_gausFit->SetTitle("response curve as <uPar-ZpT>/<ZpT> +1; numerator from mean of gaussian fit (mode of distribution)");
   GresponseCurve_gausFit->Draw("AP");
   GresponseCurve_gausFit->SetMarkerStyle(7);    // 7 is a medium dot
   GresponseCurve_gausFit->GetXaxis()->SetTitle("ZpT [GeV]");
   GresponseCurve_gausFit->GetYaxis()->SetTitle(" < u_{||} - ZpT > / < ZpT > + 1");
   GresponseCurve_gausFit->GetYaxis()->SetRangeUser(0.0, 1.1);
   GresponseCurve_gausFit->GetYaxis()->SetTitleOffset(1.4); 
   GresponseCurve_gausFit->SetName("gr_responseCurve_gausFit");
   GresponseCurve_gausFit->Write();

   TGraphErrors *GresponseCurve_gausFit_bis = new TGraphErrors(nBinsForResponse,meanZpt,response_gausFit_bis,meanZptErr,responseErr_gausFit_bis);
   GresponseCurve_gausFit_bis->SetTitle("response curve as <uPar-ZpT/ZpT> +1, mean from mean gaussian fit (mode of the distribution)");
   GresponseCurve_gausFit_bis->Draw("AP");
   GresponseCurve_gausFit_bis->SetMarkerStyle(7);    // 7 is a medium dot
   GresponseCurve_gausFit_bis->GetXaxis()->SetTitle("ZpT [GeV]");
   GresponseCurve_gausFit_bis->GetYaxis()->SetTitle(" < u_{||} - ZpT / ZpT > + 1");
   GresponseCurve_gausFit_bis->GetYaxis()->SetRangeUser(0.0, 1.1);
   GresponseCurve_gausFit_bis->GetYaxis()->SetTitleOffset(1.4); 
   GresponseCurve_gausFit_bis->SetName("gr_responseCurve_gausFit_bis");
   GresponseCurve_gausFit_bis->Write();

   // Double_t response_0jets[nBinsForResponse_0jets];
   // Double_t responseErr_0jets[nBinsForResponse_0jets];
   // Double_t meanZpt_0jets[nBinsForResponse_0jets];
   // Double_t meanZptErr_0jets[nBinsForResponse_0jets];

   // for (Int_t i = 0; i < nBinsForResponse_0jets; i++) {
   //   meanZpt_0jets[i] = HZptBinned[i]->GetMean();
   //   meanZptErr_0jets[i] = HZptBinned[i]->GetMeanError();
   //   response_0jets[i] = H_uPar_ZpT_ratio_0jets[i]->GetMean();
   //   responseErr_0jets[i] = H_uPar_ZpT_ratio_0jets[i]->GetMeanError();
   //   //cout<<i<<" meanZpt = "<<meanZpt[i]<<" +/- "<<meanZptErr[i]<<"    response = "<<response[i]<<" +/- "<<responseErr[i]<<endl;
   // }

   // TGraphErrors *GresponseCurve_0jets = new TGraphErrors(nBinsForResponse_0jets,meanZpt_0jets,response_0jets,meanZptErr_0jets,responseErr_0jets);
   // GresponseCurve_0jets->SetTitle("response curve");
   // GresponseCurve_0jets->Draw("AP");
   // GresponseCurve_0jets->SetMarkerStyle(7);    // 7 is a medium dot
   // GresponseCurve_0jets->GetXaxis()->SetTitle("ZpT [GeV]");
   // GresponseCurve_0jets->GetYaxis()->SetTitle(" < u_{||} / ZpT >");
   // GresponseCurve_0jets->GetYaxis()->SetRangeUser(0.0, 1.1);
   // GresponseCurve_0jets->GetYaxis()->SetTitleOffset(1.4); 
   // GresponseCurve_0jets->SetName("gr_responseCurve_0jets");
   // GresponseCurve_0jets->Write();

   // correcting resolution of uPar for the response
  
   Int_t nPoints = GresolutionMetNoLepParZvsZpt->GetN();
   TH1D* Hyresponse = new TH1D("Hyresponse","",nPoints,0,nPoints);
   TH1D* Hyresopar = new TH1D("Hyresopar","",nPoints,0,nPoints);
   Double_t xresopar;  // just to use TGraph::GetPoint()
   Double_t yresopar[nPoints];  
   Double_t yresoparErr[nPoints];
   Double_t yresponse[nPoints];

   for (Int_t i = 0; i < nPoints; i++) {

     GresponseCurve_gausFit_bis->GetPoint(i,xresopar,yresponse[i]);
     Hyresponse->SetBinContent(i+1, *(yresponse + i));
     Hyresponse->SetBinError(i+1, GresponseCurve_gausFit_bis->GetErrorY(i));
     GresolutionMetNoLepParZvsZpt->GetPoint(i,xresopar,yresopar[i]);
     Hyresopar->SetBinContent(i+1, *(yresopar + i));
     Hyresopar->SetBinError(i+1, GresolutionMetNoLepParZvsZpt->GetErrorY(i));

   }

   if ( !Hyresopar->Divide(Hyresponse) ) cout << " Error in Hyresopar->Divide(Hyresponse) " << endl; 
   else {

     for (Int_t i = 0; i < nPoints; i++) {

       *(yresopar + i) = Hyresopar->GetBinContent(i+1);
       *(yresoparErr + i) = Hyresopar->GetBinError(i+1);

     }

   }

   TGraphErrors *GresolutionMetNoLepParZvsZpt_Corrected = new TGraphErrors(nBinsForResponse,meanZpt,yresopar,meanZptErr,yresoparErr);
   GresolutionMetNoLepParZvsZpt_Corrected->SetTitle("resolution || from histogram's RMS, corrected for response");
   GresolutionMetNoLepParZvsZpt_Corrected->Draw("AP");
   GresolutionMetNoLepParZvsZpt_Corrected->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsZpt_Corrected->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMetNoLepParZvsZpt_Corrected->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsZpt_Corrected->GetYaxis()->SetTitleOffset(1.2); 
   GresolutionMetNoLepParZvsZpt_Corrected->SetName("gr_resolution_uPar_vs_ZpT_corrected");
   GresolutionMetNoLepParZvsZpt_Corrected->Write();

   delete Hyresponse;
   delete Hyresopar;

   // end of TGraphs

   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zlljetsControlSample);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &resoAndResponse);
 

   mySpaces(cout,2);
   myPrintYieldsMetBinInStream(cout, HzlljetsYieldsMetBin, metBinEdges, nMetBins);
 
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

   mySpaces(myfile,2);
   if (!ISDATA_FLAG && unweighted_event_flag) myfile << "======   Using unweighted events (w = 1)   ======" << endl;
   mySpaces(myfile,3);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zlljetsControlSample);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &resoAndResponse);
   mySpaces(myfile,2);
   myPrintYieldsMetBinInStream(myfile, HzlljetsYieldsMetBin, metBinEdges, nMetBins);

   // TH1D *HevtPassMonoJetSel = new TH1D("HevtPassMonoJetSel","events passing monojet selection for A#times#epsilon",nMetBins,metBinEdges);
   // TH1D *HevtPassAccSel = new TH1D("HevtPassAccSel","events passing acceptance selection for A#times#epsilon",nMetBins,metBinEdges);
   // TH1D *HevtPassEffSel = new TH1D("HevtPassEffSel","events passing efficiency selection for A#times#epsilon",nMetBins,metBinEdges);

   // // using [0] element to find step, all elements are equivalent for this purpose
   // Int_t stepMonojetSelection_In_lepAccEff = lep_acc_eff[0]->whichStepHas(maskMonoJetSelection);
   // Int_t stepAcceptance_In_lepAccEff = lep_acc_eff[0]->whichStepHas(acceptanceC.get2ToId());
   // Int_t stepEfficiency_In_lepAccEff = lep_acc_eff[0]->whichStepHas(efficiencyC.get2ToId());
   // // cout<<"step: MJ     acc     eff"<<endl;
   // // cout<<stepMonojetSelection_In_lepAccEff<<stepAcceptance_In_lepAccEff<<stepEfficiency_In_lepAccEff<<endl;
   // Double_t acc, eff, accStatErr, effStatErr, acceff, acceffStatErr;

   // mySpaces(cout,2);
   // mySpaces(myfile,2);
   // cout << "Printing acceptance and efficiency." << endl;
   // cout << "MET [GeV]     acc     acc_unc     eff     eff_unc" <<endl;
   // myfile << "MET [GeV]     acc     acc_unc     eff     eff_unc" <<endl;

   // for (Int_t i = 0; i < nMetBins; i++) {
  
   //   HevtPassMonoJetSel->SetBinContent(i+1,lep_acc_eff[i]->getEvents(stepMonojetSelection_In_lepAccEff));
   //   HevtPassMonoJetSel->SetBinError(  i+1,lep_acc_eff[i]->getEventsErr(stepMonojetSelection_In_lepAccEff));
   //   HevtPassAccSel->SetBinContent(    i+1,lep_acc_eff[i]->getEvents(stepAcceptance_In_lepAccEff));
   //   HevtPassAccSel->SetBinError(      i+1,lep_acc_eff[i]->getEventsErr(stepAcceptance_In_lepAccEff));
   //   HevtPassEffSel->SetBinContent(    i+1,lep_acc_eff[i]->getEvents(stepEfficiency_In_lepAccEff));
   //   HevtPassEffSel->SetBinError(      i+1,lep_acc_eff[i]->getEventsErr(stepEfficiency_In_lepAccEff));

   //   acc = lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff];
   //   eff = lep_acc_eff[i]->nEvents[stepEfficiency_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff];
   //   accStatErr = sqrt(acc * (1 - acc) / lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff]);
   //   effStatErr = sqrt(eff * (1 - eff) / lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff]);
   //   Hacc->SetBinContent(i+1,acc);
   //   Hacc->SetBinError(i+1,accStatErr);
   //   Heff->SetBinContent(i+1,eff);
   //   Heff->SetBinError(i+1,effStatErr);

   //   cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acc<<" "<<accStatErr<<" "<<eff<<" "<<effStatErr<<endl;
   //   myfile<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acc<<" "<<accStatErr<<" "<<eff<<" "<<effStatErr<<endl;

   // }

   // mySpaces(cout,2);
   // mySpaces(myfile,2);
   // cout << "Printing acceptance * efficiency" << endl;
   // cout << "MET [GeV]     acc*eff     acc*eff_unc" <<endl;
   // myfile << "MET [GeV]     acc*eff     acc*eff_unc" <<endl;

   // for (Int_t i = 0; i < nMetBins; i++) {
   //   // do not merge with previous loop: I want to print them after the previous loop because I might copy and paste this output to make acc * eff table
     
   //   acceff = lep_acc_eff[i]->nEvents[stepEfficiency_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff];
   //   acceffStatErr = sqrt(acceff * (1 - acceff) / lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff]);
   //   Hacceff->SetBinContent(i+1,acceff);
   //   Hacceff->SetBinError(i+1,acceffStatErr);

   //   cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acceff<<" "<<acceffStatErr<<endl;
   //   myfile<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acceff<<" "<<acceffStatErr<<endl;
  
   // }

   myfile.close();

   // TEfficiency *Acceptance = new TEfficiency(*HevtPassAccSel,*HevtPassMonoJetSel);
   // TEfficiency *Efficiency = new TEfficiency(*HevtPassEffSel,*HevtPassAccSel);
   // TEfficiency *AccTimesEff = new TEfficiency(*HevtPassEffSel,*HevtPassMonoJetSel);

   // Acceptance->SetUseWeightedEvents();
   // Efficiency->SetUseWeightedEvents();
   // AccTimesEff->SetUseWeightedEvents();

   // TGraphAsymmErrors *grAE_Acc = Acceptance->CreateGraph();
   // TGraphAsymmErrors *grAE_Eff = Efficiency->CreateGraph();
   // TGraphAsymmErrors *grAE_AccTimesEff = AccTimesEff->CreateGraph();

   // grAE_Acc->Write("grAE_Acc");
   // grAE_Eff->Write("grAE_Eff");
   // grAE_AccTimesEff->Write("grAE_AccTimesEff");

   // now get Z(inv) estimate as N_Zvv = N_Zll * R / (A*e), R being BR(Zvv)/BR(Zll) where l is either mu or e (R ~ 6)
   // TH1D *H_BR_ratio = new TH1D("H_BR_ratio",Form("BR(Z#nu#nu/BR(%s)",CONTROL_SAMPLE),nMetBins,metBinEdges);

   // for(Int_t i = 0; i <= nMetBins; i++) {
   //   H_BR_ratio->SetBinContent(i,RATIO_BR_ZINV_ZLL);
   //   H_BR_ratio->SetBinError(i,UNC_RATIO_BR_ZINV_ZLL);
   // }

   // HzvvEstimate->Multiply(HzlljetsYieldsMetBin,H_BR_ratio);
   // HzvvEstimate->Divide(Hacceff);
   // delete H_BR_ratio; //no need to save it

   // I add overflow bin's content in the last bin for all histograms where that is needed
   // for those histogram filled with Divide() method, it's not done as long as it was already done on the histograms given as
   // argument to the Divide() method
   myAddOverflowInLastBin(HZtoLLRecoPt);
   if (!ISDATA_FLAG) {
     myAddOverflowInLastBin(HZtoLLGenPt);  
     HZtoLLPt_RecoGenRatio->Divide(HZtoLLRecoPt,HZtoLLGenPt);
   }

   for (Int_t i = 0; i < nMetBins; i++) {
     myAddOverflowInLastBin(HZtoLLRecoPt_MetBin[i]);
   }

   
   if (!ISDATA_FLAG) {
     for (Int_t i = 0; i < nMetBins; i++) {
       myAddOverflowInLastBin(HZtoLLGenPt_MetBin[i]);
       HZtoLLPt_RecoGenRatio_MetBin[i]->Divide(HZtoLLRecoPt_MetBin[i],HZtoLLGenPt_MetBin[i]);  
     }
   }

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
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &zlljetsControlSample,commentInTable);
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &resoAndResponse);
     fprintf(fp,"\\end{document}\n");      
     fclose(fp);

   }

   // end of tex file

}



