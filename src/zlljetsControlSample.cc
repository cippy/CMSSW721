#define zlljetsControlSample_cxx
#include "EmanTreeAnalysis.h"
//#include "AnalysisDarkMatter.h" // already included in EmanTreeAnalysis.h
//C or C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream ; to read array of numbers from a line in a file
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

#ifdef zlljetsControlSample_cxx

zlljetsControlSample::zlljetsControlSample(TTree *tree, const char* inputSuffix) : AnalysisDarkMatter(tree) {
  //cout <<"check in constructor "<<endl;
  suffix = inputSuffix;  // it is the sample name (e.g. QCD, ZJetsToNuNu ecc...)
  //edimarcoTree_v3::Init(tree);
  AnalysisDarkMatter::Init(tree);  // could also be just Init(tree)

}

#endif

void zlljetsControlSample::loop(const char* configFileName, const Int_t ISDATA_FLAG, const Int_t unweighted_event_flag, vector< Double_t > &yRow, vector< Double_t > &eRow, vector< Double_t > &uncRow)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   //fChain->SetBranchStatus("LHEorigWeight",1); // contains negative values: the weight in the event is weight*LHEorigWeight

   //fChain->SetBranchStatus("genWeight",1); 

   fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
   fChain->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (isolation included)
   fChain->SetBranchStatus("nEle20T",1);  // # of electrons passing tight selection (isolation included)
   //fChain->SetBranchStatus("nTau18V",1);
   fChain->SetBranchStatus("nTauClean18V",1);

   fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
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

     //fChain->SetBranchStatus("xsec",1);

     fChain->SetBranchStatus("vtxWeight",1);   // weight to have better agreement between data and MC (added in tree from 06/10/15)
   }

   fChain->SetBranchStatus("met_pt",1);
   //fChain->SetBranchStatus("met_eta",1);
   fChain->SetBranchStatus("met_phi",1);

   fChain->SetBranchStatus("metNoMu_pt",1);
   //fChain->SetBranchStatus("metNoMu_eta",1);
   fChain->SetBranchStatus("metNoMu_phi",1);

   fChain->SetBranchStatus("nVert",1);  // number of good vertices 

   // met filters to be used (the config file has a parameter saying whether they should be used or not)
   fChain->SetBranchStatus("cscfilter",1);
   fChain->SetBranchStatus("ecalfilter",1);
   fChain->SetBranchStatus("hbheFilterNew25ns",1);
   fChain->SetBranchStatus("hbheFilterIso",1);

   //added on November 2015. These are new variables (except for weight, which has just changed in the definition)
   fChain->SetBranchStatus("nBTag15",1);  // for b-jet veto
   fChain->SetBranchStatus("dphijm",1);          // flag for dphi minimum between met and any of the jets in the event (using only the first four jets
   fChain->SetBranchStatus("weight",1);   // modified since 17 November 2015: now it includes the whol weight, e.g. 1000*xsec*genWeight ...
   fChain->SetBranchStatus("events_ntot",1);    // equivalent to SUMWEIGHTS for samples before 17 November 2015
   fChain->SetBranchStatus("JetClean_leadClean",1); // has new cleaning on energy fractions (added on 17 November 2015) 

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
   Int_t HLT_FLAG;
   Double_t HLT_LEP1PT;
   Double_t HLT_LEP2PT;
   Double_t HLT_LEP1ETA;
   Double_t HLT_LEP2ETA;
   Double_t METNOLEP_START;
   Int_t MET_FILTERS_FLAG;
   Double_t JMET_DPHI_MIN;
   // Double_t XSEC_OVER_NPROCESSED;
   // Double_t SUMWEIGHTS;
   string FILENAME_BASE;
   string DIRECTORY_TO_SAVE_FILES;
   string DIRECTORY_NAME;

   vector<Double_t> metBinEdgesVector;  // filled with values in file named configFileName

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
	 else if (parameterName == "METNOLEP_START") METNOLEP_START = value;
	 else if (parameterName == "MET_FILTERS_FLAG") MET_FILTERS_FLAG = value;
	 else if (parameterName == "JMET_DPHI_MIN") JMET_DPHI_MIN = value;

	 if (!ISDATA_FLAG) {

	 if (parameterName == "GENLEP1PT") GENLEP1PT = value;
	 else if (parameterName == "GENLEP2PT") GENLEP2PT = value;
	 else if (parameterName == "GENLEP1ETA") GENLEP1ETA = value;
	 else if (parameterName == "GENLEP2ETA") GENLEP2ETA = value;
	 else if (parameterName == "GEN_ZMASS_LOW") GEN_ZMASS_LOW = value;
	 else if (parameterName == "GEN_ZMASS_UP") GEN_ZMASS_UP = value;
	 // else if (parameterName == "XSEC_OVER_NPROCESSED") XSEC_OVER_NPROCESSED = value;
	 // else if (parameterName == "SUMWEIGHTS") SUMWEIGHTS = value;

	 }
	 

       } else if (parameterType == "STRING") {
	 
	 inputFile >> parameterName >> name;
	 cout << right << setw(20) << parameterName << "  " << left << name << endl;
	 if (parameterName == "FILENAME_BASE") {

	   FILENAME_BASE = name; 
	   if ( !ISDATA_FLAG && unweighted_event_flag) FILENAME_BASE += "_weq1";  // if using unit weight, add _weq1 to filename (weq1 means weight = 1)

	 }

	 if (parameterName == "DIRECTORY_PATH") {  // name of directory where files are saved

	  DIRECTORY_TO_SAVE_FILES = name;
	  //std::cout << "Files will be saved in '" << name << "' ." <<std::endl;

	} 

	if (parameterName == "DIRECTORY_NAME") {  // name of directory where files are saved

	  DIRECTORY_NAME = name;
	  //std::cout << "Files will be saved in directory named '" << name << "' ." <<std::endl;

	} 

       } else if (parameterType == "ARRAY_NUM") {

	 inputFile >> parameterName;

	 if (parameterName == "MET_BIN_EDGES") { 

	   cout << right << setw(20) << parameterName << "  ";
	   string stringvalues;
	   getline(inputFile, stringvalues);    // read whole line starting from current position (i.e. without reading ARRAY_NUM)
	   istringstream iss(stringvalues);
	   Double_t num;

	   while(iss >> num) {
	    
	     metBinEdgesVector.push_back(num);
	     cout << metBinEdgesVector.back() << " ";

	   }

	   cout << endl;

	 }

       }

     }
     
     mySpaces(cout,2);

     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << endl;
     exit(EXIT_FAILURE);

   }

   string outputFolder =  DIRECTORY_TO_SAVE_FILES + DIRECTORY_NAME + "/";

   //Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   // Double_t metBinEdges[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 750., 850., 1000.};
   // Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;
   Int_t nMetBins = (metBinEdgesVector.size()) - 1;

   cout << "MetBinEdges: [ ";
   for(Int_t i = 0; i <= nMetBins; i++) {
     if (i != nMetBins) cout << metBinEdgesVector[i] << ", ";
     else cout << metBinEdgesVector[i] << "]" << endl;
   }
  
   // monojet selection (muon or electron veto set afterwards, depending on the sample)
   selection metFiltersC("metFiltersC","met filters","cscfilter, ecalfilter, hbheFilterNew25ns, hbheFilterIso");
   selection jet1C("jet1C","jet1",Form("nJetClean >= 1 && JetClean1_pt > %4.0lf && abs(JetClean1_eta) < %1.1lf",(Double_t)J1PT,J1ETA));
   selection jetMetDphiMinC("jetMetDphiMinC",Form("min[dphi(jets,MET)] > %1.1lf",JMET_DPHI_MIN),"minimum dphi between jets and MET (using only the first 4 jets)");
   selection jetNoiseCleaningC("jetNoiseCleaningC","noise cleaning","energy fractions (only for jet1): CH > 0.1; NH < 0.8");
   selection bjetVetoC("bjetVetoC","b-jets veto");
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");
   selection tauLooseVetoC;
   if (TAU_VETO_FLAG) tauLooseVetoC.set("tauLooseVetoC","tau veto");
   // additional selections for control sample
   selection oppChargeLeptonsC("oppChargeLeptonsC","OS/SF leptons");
   selection invMassC("invMassC",Form("mass in [%3.0lf,%3.0lf]",DILEPMASS_LOW,DILEPMASS_UP));
   // following selections are set differently in the next "if" statements depending on the lepton flavour 
   //selection metNoLepC[metCut.size()];
   selection lepLooseVetoC;
   // selection twoLeptonsC;
   selection twoLepLooseC;
   selection tightLepC;
   // selection lep1tightIdIso04C;
   // selection twoLepTightC;
   // selection lep1ptC;
   // selection lep2ptC;
   // selection lep1etaC;  
   //selection lep2etaC;
   selection genLepC;  
   selection metNoLepStartC;
   selection HLTlepC;
   selection recoGenLepMatchC;
   //selection lep2tightIdIso04C;

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

   // Double_t SUMWEIGHTS;
   // vector<Double_t> sumWeightVector;  // used for spring15_25ns
   // vector<Int_t> eventsInSubsamples;

   Double_t nTotalWeightedEvents = 0.0;     
   //Double_t nEventsAfterMatchRecoGen = 0.0;
   Int_t HLT_passed_flag = 1;          // some computations (for e) require a trigger preselection, while others don't. The former will be done if the flag is set to 1
                                                       // it's set to 1 because if the trigger selection is not applied every event must be considered to be a "good" event having passed all preselections
                                                       // Actually in this code the trigger is not always necessary, but I keep it like this nonetheless.


   Float_t *ptr_nLepLoose = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches
   Float_t *ptr_nLep10V = NULL;   

   Float_t *ptr_nLepTight = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches 

   Float_t *ptr_metNoLepPt = NULL;       // only needed for muons, it will point to the branches with the metNoMu_pt, then metNoLepPt = *ptr_metNoLepPt (metNoLepPt defined below)
   //Float_t *ptr_metNoLepEta = NULL; 
   Float_t *ptr_metNoLepPhi = NULL;  

   Float_t nLepLoose = 0.0;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such
   Float_t nLep10V = 0.0;

   Float_t nLepTight = 0.0;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such

   Double_t metNoLepPt = 0.0;        // this variable will be assigned with *ptr_metNoLepPt, where the pointer will point to the branch metNoMu_pt for mu, and with a hand-defined variable for e
   //Double_t metNoLepEta = 0.0;
   Double_t metNoLepPhi = 0.0;   // same story as above

   if ( !ISDATA_FLAG && unweighted_event_flag) cout << "Warning: no weight applied to events (w = 1)" << endl;  // if MC with unit weight, make user know

   // if using sample spring15, need to use vtxW to get same Nvtx distribution as seen in data. For older trees it's not used

   // the following flag is needed to enable search for Z->ll at generator level. For MC samples different from DYJetsToLL I must not require 2 gen leptons from Z
   Int_t using_zlljets_MCsample_flag = 0;
   if ( !ISDATA_FLAG && ( suffix == "DYJetsToLL" )  )  using_zlljets_MCsample_flag = 1; 

   Int_t using_ztautaujets_MCsample_flag = 0;
   if ( !ISDATA_FLAG && ( suffix == "ZJetsToTauTau" )) using_ztautaujets_MCsample_flag = 1; 

   if (ISDATA_FLAG) {
     strcpy(ROOT_FNAME,(FILENAME_BASE + "_DATA.root").c_str());
     strcpy(TXT_FNAME,(FILENAME_BASE + "_DATA.txt").c_str());
     strcpy(TEX_FNAME,(FILENAME_BASE + "_DATA.tex").c_str());
   } else {
     strcpy(ROOT_FNAME,(FILENAME_BASE + "_" + suffix + ".root").c_str());
     strcpy(TXT_FNAME,(FILENAME_BASE + "_" + suffix + ".txt").c_str());
     strcpy(TEX_FNAME,(FILENAME_BASE + "_" + suffix + ".tex").c_str());
   }

   if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...
     
     strcpy(FLAVOUR,"mu");
     strcpy(LL_FLAVOUR,"mumu");
     strcpy(CONTROL_SAMPLE,"Z-->mumu");
         
     ptr_nLepLoose = &nMu10V;                      // ask 2 muons
     ptr_nLep10V = &nEle10V;                         // veto on electrons
     ptr_nLepTight = &nMu20T;                                            
     ptr_metNoLepPt = &metNoMu_pt;               // for muons  get this variable from the tree 
     //ptr_metNoLepEta = &metNoMu_eta;               // for muons  get this variable from the tree 
     ptr_metNoLepPhi = &metNoMu_phi;         // for muons  get this variable from the tree

     lepLooseVetoC.set("eLooseVetoC","electrons veto");
     //twoLeptonsC.set("twomuonsC","muons");
     twoLepLooseC.set("twomuLooseC","2 loose muons");
     tightLepC.set("tightMuC",">0 tight muons");
     // lep1tightIdIso04C.set("mu1tightIdIso04C","leading muon tight","tight ID + relIso04 (as Emanuele)");
     // twoLepTightC.set("twomuTightC","2 tight muons");
     // lep1ptC.set("mu1ptC",Form("mu1pt > %3.0lf",LEP1PT),"leading muon pt");
     // lep2ptC.set("mu2ptC",Form("mu2pt > %3.0lf",LEP2PT),"trailing muon pt");
     // lep1etaC.set("mu1etaC",Form("|mu1eta| < %1.1lf",LEP1ETA),"leading muon eta");  
     // lep2etaC.set("mu2etaC",Form("|mu2eta| < %1.1lf",LEP2ETA),"trailing muon eta");
     if (!ISDATA_FLAG && using_zlljets_MCsample_flag) {
       genLepC.set("genMuonsC","muons generated");     
       recoGenLepMatchC.set("recoGenMuMatchC","reco-gen match (DR = 0.1)","only for zlljets: looks for matching of reco and gen particles");      
     }

     metNoLepStartC.set("metNoMu200C",Form("metNoMu > %2.0lf",METNOLEP_START));
     HLTlepC.set("HLTmuonC","HLT for muons");
     // lep2tightIdIso04C.set("mu2tightIdIso04C","trailing muon tight","tight ID + relIso04 (as Emanuele)");


   } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

     strcpy(FLAVOUR,"ele");
     strcpy(LL_FLAVOUR,"ee");
     strcpy(CONTROL_SAMPLE,"Z-->ee");

     ptr_nLepLoose = &nEle10V;                      // ask 2 electrons
     ptr_nLep10V = &nMu10V;                         // veto on muons   
     ptr_nLepTight = &nEle20T;

     lepLooseVetoC.set("muLooseVetoC","muons veto");
     //twoLeptonsC.set("twoelectronsC","electrons");
     twoLepLooseC.set("twoeleLooseC","2 loose electrons");
     tightLepC.set("tightEleC",">0 tight electrons");
     // lep1tightIdIso04C.set("ele1tightIdIso04C","leading electron tight","tight ID + relIso04 (as Emanuele)");
     // twoLepTightC.set("twoeleTightC","2 tight electrons");
     // lep1ptC.set("ele1ptC",Form("ele1pt > %3.0lf",LEP1PT),"leading electron pt");
     // lep2ptC.set("ele2ptC",Form("ele2pt > %3.0lf",LEP2PT),"trailing electron pt");
     // lep1etaC.set("ele1etaC",Form("|ele1eta| < %1.1lf",LEP1ETA),"leading electron eta");  
     // lep2etaC.set("ele2etaC",Form("|ele2eta| < %1.1lf",LEP2ETA),"trailing electron eta");
     if (!ISDATA_FLAG && using_zlljets_MCsample_flag) {
       genLepC.set("genElectronsC","electrons generated");     
       recoGenLepMatchC.set("recoGenEleMatchC","reco-gen match (DR = 0.1)","only for zlljets: looks for matching of reco and gen particles");    
     }

     metNoLepStartC.set("metNoEle200C",Form("metNoEle > %2.0lf",METNOLEP_START));
     HLTlepC.set("HLTelectronC","HLT for electrons");
     // lep2tightIdIso04C.set("ele2tightIdIso04C","trailing electron tight","tight ID + relIso04 (as Emanuele)");

   }

   selection genTauC;
   if (!ISDATA_FLAG && using_ztautaujets_MCsample_flag) genTauC.set("genTauC","taus generated");                       

   selection::checkMaskLength();
   selection::printActiveSelections(cout); 

   UInt_t maskTightTag;   // holds cuts for lepton tight selection, which is different between muons and electrons

   mask zlljetsControlSample(Form("%s control sample (%s gen if DYJetsToLL MC) with selection flow as Emanuele's",CONTROL_SAMPLE,FLAVOUR));
   if (!ISDATA_FLAG) {
     if (using_zlljets_MCsample_flag) zlljetsControlSample.append(genLepC.get2ToId());
     if (using_ztautaujets_MCsample_flag) zlljetsControlSample.append(genTauC.get2ToId());
   }

   if (MET_FILTERS_FLAG != 0) zlljetsControlSample.append(metFiltersC.get2ToId());

   if ( HLT_FLAG != 0 ) zlljetsControlSample.append(HLTlepC.get2ToId());

   if (fabs(LEP_PDG_ID) == 13) {  

     // maskTightTag = lep1tightIdIso04C.get2ToId() + /* lep2tightIdIso04C.get2ToId() +*/lep1ptC.get2ToId() + lep2ptC.get2ToId() + lep1etaC.get2ToId() + lep2etaC.get2ToId();  ;  // for now tight requirements on pt and eta are already included in the loose condition because they coincide (not true for electrons)

     zlljetsControlSample.append(twoLepLooseC.get2ToId());
     zlljetsControlSample.append(tightLepC.get2ToId());
     zlljetsControlSample.append(oppChargeLeptonsC.get2ToId());     
     zlljetsControlSample.append(invMassC.get2ToId());

   } else if (fabs(LEP_PDG_ID) == 11) {  

     //maskTightTag = lep1tightIdIso04C.get2ToId() /*+ lep2tightIdIso04C.get2ToId() */+ lep1ptC.get2ToId() + lep2ptC.get2ToId() + lep1etaC.get2ToId() + lep2etaC.get2ToId();

     zlljetsControlSample.append(twoLepLooseC.get2ToId());
     zlljetsControlSample.append(tightLepC.get2ToId());
     zlljetsControlSample.append(oppChargeLeptonsC.get2ToId());
     zlljetsControlSample.append(invMassC.get2ToId());

   }
   
   zlljetsControlSample.append(lepLooseVetoC.get2ToId());
   if (TAU_VETO_FLAG) zlljetsControlSample.append(tauLooseVetoC.get2ToId());
   zlljetsControlSample.append(gammaLooseVetoC.get2ToId());
   zlljetsControlSample.append(bjetVetoC.get2ToId());
   if (METNOLEP_START != 0) zlljetsControlSample.append(metNoLepStartC.get2ToId());
   zlljetsControlSample.append(jet1C.get2ToId());
   zlljetsControlSample.append(jetNoiseCleaningC.get2ToId());
   zlljetsControlSample.append(jetMetDphiMinC.get2ToId());
   
   
  
   if (!ISDATA_FLAG && using_zlljets_MCsample_flag) zlljetsControlSample.append(recoGenLepMatchC.get2ToId());


   cout << "Opening file " <<ROOT_FNAME<< " in folder " << outputFolder << endl;

   TFile *rootFile = new TFile((outputFolder + ROOT_FNAME).c_str(),"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout<<"Error: file \""<<outputFolder + ROOT_FNAME<<"\" was not opened."<<endl;
     exit(EXIT_FAILURE);
   }
 

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   //TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   //Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};       

   Float_t invMassBinWidth = 1.0;  // invariant mass histogram's bin width in GeV
   Int_t NinvMassBins = (DILEPMASS_UP - DILEPMASS_LOW) / invMassBinWidth;

   TH1D *HzlljetsYieldsMetBin = new TH1D("HYieldsMetBin",Form("yields of %s control sample (%s gen if DY MC) in bins of met; #slash{E}_{T};# of events",CONTROL_SAMPLE,CONTROL_SAMPLE),nMetBins,metBinEdgesVector.data());
   
   TH1D *HinvMass = new TH1D("HinvMass","",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);    // for MC it's done on Z->mumu or Z->ee at gen level
   TH1D *HvtxDistribution = new TH1D("HvtxDistribution","",40,-0.5,39.5);   
   TH1D *HnjetsDistribution = new TH1D("HnjetsDistribution","njets using nJetClean30",10,-0.5,9.5);   
   TH1D *Hj1j2dphiDistribution = new TH1D("Hj1j2dphiDistribution","",30,0.0,3.0);
   TH1D *HjetMetDphiMinDistribution = new TH1D("HjetMetDphiMinDistribution","",30,0.0,3.2);
   TH1D *Hjet1etaDistribution = new TH1D("Hjet1etaDistribution","",60,-3.0,3.0);
   TH1D *Hjet2etaDistribution = new TH1D("Hjet2etaDistribution","",60,-3.0,3.0);
   TH1D *HmetNoLepDistribution = new TH1D("HmetNoLepDistribution","",100,0.0,1000.0);
   TH1D *HzptDistribution = new TH1D("HzptDistribution","",200,0.0,1000.0);
   TH1D *Hjet1ptDistribution = new TH1D("Hjet1ptDistribution","",100,0.0,1000.0); 
   TH1D *Hjet2ptDistribution = new TH1D("Hjet2ptDistribution","",100,0.0,1000.0);
     
   TH1D *HzlljetsInvMassMetBinGenLep[nMetBins];

   for (Int_t i = 0; i < nMetBins; i++) {

     HzlljetsInvMassMetBinGenLep[i] = new TH1D(Form("HzlljetsInvMassMetBinGenLep_met%2.0lfTo%2.0lf",metBinEdgesVector[i],metBinEdgesVector[i+1]),"",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);

   } 

   // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
   TH1D *HmetBinEdges = new TH1D("HmetBinEdges","bin edges for met distributions",nMetBins+1,0.0,nMetBins+1);
   for (Int_t i = 0; i <= nMetBins; i++) {
     HmetBinEdges->SetBinContent(i+1,metBinEdgesVector[i]);
   }

   // deciding  what is the event weight
   Double_t newwgt; 

   if (ISDATA_FLAG || unweighted_event_flag) newwgt = 1.0;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"zlljetsControlSample::loop()"<<endl;
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

       newwgt = LUMI * weight /*/ events_ntot*/;  // starting from 17 November, "events_ntot" substitutes SUMWEIGHT and is already present in the trees. Same for weight, which is now defined as "1000 * xsec * genWeight" (1000*xsec is the cross section in fb, since xsec is in pb.)
       // I found out that division by events_ntot was already included in weight definition

     }

     nTotalWeightedEvents += newwgt;  // counting events with weights

     nLepLoose = *ptr_nLepLoose;          
     nLep10V = *ptr_nLep10V;
     nLepTight = *ptr_nLepTight;

     //Double_t ZgenMass;        // not used for now
     Double_t ZtoLLGenPt = 0;    // filled below (only if running on MC DYJetsToLL)
     Double_t ZtoLLRecoPt = 0;   // filled below

     // genLepFound_flag is used when analysing DYJetsToLL in MC fo Z->mumu or Z->ee. For other MC samples it's not used.

     if (!ISDATA_FLAG) {

       if (using_zlljets_MCsample_flag) {

	 //genLepFound_flag = 1;
	 genLepFound_flag = myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, LEP_PDG_ID, 23, firstIndexGen, secondIndexGen, Z_index, GenPart_motherIndex); 
	 //if (!genLepFound_flag) continue;  // if not found gen ee or mumu for MC DYJetsToLL ( l = mu or e) skip the event. This makes things faster
	 if (genLepFound_flag != 0) {
	   eventMask += genLepC.addToMask(1);
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

     //cout << "LepGood_pdgId[0]; LepGood_pdgId[1] --> " << LepGood_pdgId[0] << "   " << LepGood_pdgId[1] << endl;

     // recoLepFound_flag = myGetPairIndexInArray(LEP_PDG_ID, nLepGood, LepGood_pdgId, firstIndex, secondIndex);  

     // look if the first two good leptons are OS/SF (flavour depending on config file)
     // the check for 2 OS/SF leptons in LepGood list is just useful to speed up things, but would not be necessary if other variables or computations that require this condition are only used at the end of selection (which already includes 2 OS/SF condition).
     // thus, we evaluate this right now and use this pice of information for later use

     if ( (fabs(LepGood_pdgId[firstIndex]) ==  LEP_PDG_ID) && (LepGood_pdgId[firstIndex] == (-1) * LepGood_pdgId[secondIndex]) ) recoLepFound_flag = 1;
     else recoLepFound_flag = 0;

     if (recoLepFound_flag) {
       l1reco.SetPtEtaPhiM(LepGood_pt[firstIndex],LepGood_eta[firstIndex],LepGood_phi[firstIndex],LepGood_mass[firstIndex]);
       l2reco.SetPtEtaPhiM(LepGood_pt[secondIndex],LepGood_eta[secondIndex],LepGood_phi[secondIndex],LepGood_mass[secondIndex]);
       Zreco = l1reco + l2reco;
       ZtoLLRecoPt = Zreco.Pt();
     }

     if (fabs(LEP_PDG_ID) == 13) { 

       if ( HLT_FLAG != 0) {

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

       if ( HLT_FLAG != 0 ) {

       	 if ( recoLepFound_flag && (LepGood_tightId[firstIndex] > 0.5) && (LepGood_tightId[secondIndex]  > 0.5) && 
       	      (fabs(LepGood_eta[firstIndex]) < HLT_LEP1ETA) && (fabs(LepGood_eta[secondIndex]) < HLT_LEP2ETA) && 
       	      (LepGood_pt[firstIndex] > HLT_LEP1PT) && (LepGood_pt[secondIndex] > HLT_LEP2PT) ) HLT_passed_flag = 1; 	 
	 else HLT_passed_flag = 0;  //continue;

       }  // end of   if ( HLT_FLAG )

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
     
     eventMask += jet1C.addToMask(nJetClean30 >= 1 && JetClean_pt[0] > J1PT /*&& fabs(JetClean_eta[0]) < J1ETA*/);
     eventMask += jetMetDphiMinC.addToMask(fabs(dphijm > JMET_DPHI_MIN));
     eventMask += jetNoiseCleaningC.addToMask(JetClean_leadClean[0] > 0.5);
     eventMask += bjetVetoC.addToMask(nBTag15 == 0);
     eventMask += lepLooseVetoC.addToMask(nLep10V == 0);
     eventMask += tauLooseVetoC.addToMask(nTauClean18V == 0);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     eventMask += metNoLepStartC.addToMask(metNoLepPt > METNOLEP_START);
     eventMask += metFiltersC.addToMask(cscfilter == 1 && ecalfilter == 1 && hbheFilterNew25ns == 1 && hbheFilterIso == 1);  

     // the following make sense only if recoLepFound_flag == 1 (i.e. flag is true), which means that fabs(LepGood_pdgId[firstIndex/secondIndex]) == LEP_PDG_ID) is 
     // true
     // also, 2 OS/SF leptons are present

     if (recoLepFound_flag == 1) {
       eventMask += HLTlepC.addToMask(HLT_passed_flag);     
       eventMask += oppChargeLeptonsC.addToMask(1); // included in recoLepFound_flag togehter with correct flavour 
       eventMask += twoLepLooseC.addToMask(((Int_t) nLepLoose) == 2);
       eventMask += tightLepC.addToMask(nLepTight > 0);
       // eventMask += lep1ptC.addToMask((LepGood_pt[firstIndex] > LEP1PT)); 
       // eventMask += lep1etaC.addToMask( (fabs(LepGood_eta[firstIndex]) < LEP1ETA) );
       // eventMask += lep2ptC.addToMask((LepGood_pt[secondIndex] > LEP2PT) );
       // eventMask += lep2etaC.addToMask((fabs(LepGood_eta[secondIndex]) < LEP2ETA) );
       eventMask += invMassC.addToMask((mZ1 > DILEPMASS_LOW) && (mZ1 < DILEPMASS_UP));     
       // eventMask += lep1tightIdIso04C.addToMask((LepGood_tightId[firstIndex] > 0.5 ) && (LepGood_relIso04[firstIndex] < LEP_ISO_04 ) );
       // eventMask += lep2tightIdIso04C.addToMask((LepGood_tightId[secondIndex] > 0.5) && (LepGood_relIso04[secondIndex] < LEP_ISO_04 ) );
       
     }

     // end of eventMask building

     // test matching of reco and gen lep for DY MC 
     if (!ISDATA_FLAG && using_zlljets_MCsample_flag) {
    
       // now we require a DeltaR cut between gen-level and reco-level leptons (mu or e) from Z in Drell-Yan events to assess that lreco comes from lgen
       
       // since genLepFound_flag && recoLepFound_flag require 2 OS/SF leptons (with correct flavour) to be found in their lists, we have two possible cases:
       // 1) the first gen lepton and the first reco lepton have same charge
       // 2) the first gen lepton and the first reco lepton have opposite charge
       // if the first case is not satisfied, then the second will automatically be

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
	 HjetMetDphiMinDistribution->Fill(dphijm,newwgt);
	 if (nJetClean30 >= 2) {
	   Hj1j2dphiDistribution->Fill(dphijj,newwgt);
	   Hjet2etaDistribution->Fill(JetClean_eta[1],newwgt);
	   Hjet2ptDistribution->Fill(JetClean_pt[1],newwgt);
	 }

     }
	

     // now entering analysis in bins of met

     if ((metNoLepPt > metBinEdgesVector[0]) && (metNoLepPt < metBinEdgesVector[nMetBins])) {

       Int_t bin = myGetBin(metNoLepPt,metBinEdgesVector.data(),nMetBins);
       
       // if ((eventMask & zlljetsControlSample.globalMask.back()) == zlljetsControlSample.globalMask.back()) {
       //   // this histogram holds the invariant mass distribution (one for each met bin)
       //   HinvMass[bin]->Fill(mZ1,newwgt);   
       // }

       if ( ((eventMask & zlljetsControlSample.globalMask.back()) == zlljetsControlSample.globalMask.back()) ) { 
 
	 HzlljetsInvMassMetBinGenLep[bin]->Fill(mZ1,newwgt); 

       }
	 
     }                      // end of    if ((metNoLepPt > metBinEdges[0]) && (metNoLepPt < metBinEdges[nMetBins])) 
       
   }                        // end of loop on entries

   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zlljetsControlSample);

   mySpaces(cout,2);
   myPrintYieldsMetBinInStream(cout, HzlljetsYieldsMetBin, metBinEdgesVector.data(), nMetBins);
 
   cout<<"creating file '"<<TXT_FNAME<<"' in folder "<< outputFolder <<" ..."<<endl;

   ofstream myfile((outputFolder + TXT_FNAME).c_str(),ios::out);

   if ( !myfile.is_open() ) {

     cout<<"Error: unable to open file "<<TXT_FNAME<<" !"<<endl;
     exit(EXIT_FAILURE);
     
   }

   //opening inputFile named configFileName again to save content in myfile named TXT_FNAME
   /*
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
   */

   if (!ISDATA_FLAG && unweighted_event_flag) myfile << "======   Using unweighted events (w = 1)   ======" << endl;
   mySpaces(myfile,3);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zlljetsControlSample);
   mySpaces(myfile,2);
   myPrintYieldsMetBinInStream(myfile, HzlljetsYieldsMetBin, metBinEdgesVector.data(), nMetBins);

   myfile.close();

   // filling with yields and efficiency: I will use efficiency with respect to total and not to previous step, but I could make this choice in the config file

   // entry point
   if (using_zlljets_MCsample_flag == 1 || using_ztautaujets_MCsample_flag == 1) {
     yRow.push_back(zlljetsControlSample.nEvents[0]); // [0] refers to genLep, which is the first selection in these cases
     eRow.push_back(1.0000);
     uncRow.push_back(sqrt(zlljetsControlSample.nEvents[0]));
   } else {
     yRow.push_back(nTotalWeightedEvents); // [0] 
     eRow.push_back(1.0000);
     uncRow.push_back(sqrt(nTotalWeightedEvents));
   }

   vector<Int_t> selStep;   //array to store index of step to form selection flow (might want to consider two or more steps together and not separated)
   //first step is the preselection before OS condition: doing like this because it might change (there can be or not MetNoLep cut, but I want the last step)

   // in case a step wold be present for some sample but not for others (e.g. the RecoGen match done only in Zll MC), the step is referred to as -1 and the corresponding values are set to -1, so that, when printing the table, yields will be filled with " / / " which means " uneffected" (because that step was not done)

   if (MET_FILTERS_FLAG != 0) selStep.push_back(zlljetsControlSample.whichStepHas(metFiltersC.get2ToId()));
   selStep.push_back(zlljetsControlSample.whichStepHas(twoLepLooseC.get2ToId()));
   selStep.push_back(zlljetsControlSample.whichStepHas(tightLepC.get2ToId()));
   selStep.push_back(zlljetsControlSample.whichStepHas(oppChargeLeptonsC.get2ToId()));  
   selStep.push_back(zlljetsControlSample.whichStepHas(invMassC.get2ToId()));
   selStep.push_back(zlljetsControlSample.whichStepHas(lepLooseVetoC.get2ToId()));
   if (TAU_VETO_FLAG) selStep.push_back(zlljetsControlSample.whichStepHas(tauLooseVetoC.get2ToId()));
   selStep.push_back(zlljetsControlSample.whichStepHas(gammaLooseVetoC.get2ToId()));
   selStep.push_back(zlljetsControlSample.whichStepHas(bjetVetoC.get2ToId()));
   if (METNOLEP_START != 0) selStep.push_back(zlljetsControlSample.whichStepHas(metNoLepStartC.get2ToId()));
   selStep.push_back(zlljetsControlSample.whichStepHas(jet1C.get2ToId()));
   selStep.push_back(zlljetsControlSample.whichStepHas(jetNoiseCleaningC.get2ToId()));
   selStep.push_back(zlljetsControlSample.whichStepHas(jetMetDphiMinC.get2ToId()));
   if (!ISDATA_FLAG && using_zlljets_MCsample_flag)  selStep.push_back(zlljetsControlSample.whichStepHas(recoGenLepMatchC.get2ToId()));
   else selStep.push_back(selStep.back());  // in this case copy the previous entry

   for(Int_t i = 0; i < selStep.size(); i++) {
     
     yRow.push_back(zlljetsControlSample.nEvents[selStep[i]]);
     uncRow.push_back(sqrt(yRow.back()));
     if (i == 0) eRow.push_back(zlljetsControlSample.nEvents[selStep[i]]/nTotalWeightedEvents);
     else if( (i != 0) && (zlljetsControlSample.nEvents[selStep[i]-1] == 0) ) eRow.push_back(1.0000);
     else eRow.push_back(zlljetsControlSample.nEvents[selStep[i]]/zlljetsControlSample.nEvents[selStep[i]-1]);
   
   }

   // fill last bin with overflow 
   myAddOverflowInLastBin(HmetNoLepDistribution);
   myAddOverflowInLastBin(HzptDistribution);
   myAddOverflowInLastBin(Hjet1ptDistribution);
   myAddOverflowInLastBin(Hjet2ptDistribution);

   rootFile->Write();

   rootFile->Close();
   delete rootFile;

   //creating a .tex file to build tables with data
   FILE *fp;
   fp = fopen((outputFolder + TEX_FNAME).c_str(),"w");

   if ( fp == NULL)  cout<<"Error: '"<<TEX_FNAME<<"' not opened"<<endl;
   else {

     cout<<"creating file '"<<TEX_FNAME<<" in folder " << outputFolder << "' ..."<<endl;
     myAddDefaultPackages(fp,TEX_FNAME);
     fprintf(fp,"\\begin{document}\n");
     fprintf(fp,"\n");
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &zlljetsControlSample);
     fprintf(fp,"\\end{document}\n");      
     fclose(fp);

   }

   // end of tex file

}




