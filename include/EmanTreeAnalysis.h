#ifndef EmanTree_Analysis_h
#define EmanTree_Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <cstdlib>
#include <vector>

#include "AnalysisDarkMatter.h"

namespace myAnalyzerTEman {

  class zlljets_Axe_noSkim_light : public AnalysisDarkMatter /*,public edimarcoTreeFriend*/ {
  public:

    // only computes acceptance and efficiency in order to be faster
    zlljets_Axe_noSkim_light(TTree *tree);
    virtual ~zlljets_Axe_noSkim_light() { std::cout<<"~zlljets_Axe_noSkim_light() called"<<std::endl; }
  
    void loop(const char* configFileName);

  };

  class zlljets_resoResp : public AnalysisDarkMatter /*,public edimarcoTreeFriend*/ {
  public:

    zlljets_resoResp(TTree *tree);
    virtual ~zlljets_resoResp() { std::cout<<"~zlljets_resoResp() called"<<std::endl; }
  
    void loop(const char* configFileName, const Int_t ISDATA_FLAG, const Int_t unweighted_event_flag);

  };

  class zlljetsControlSample : public AnalysisDarkMatter {
  public:

    zlljetsControlSample(TTree *tree);  
    virtual ~zlljetsControlSample() { std::cout<<"~zlljetsControlSample() called"<<std::endl; }
  
    void setSelections(); 
    void setMask();
    void setHistograms();
    void setNumberParameterValue(const std::string, const Double_t);  // set some numeric variables if they are in config file
    void setControlSampleSpecificParameter();
    void setVarFromConfigFile();
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);

    selection oppChargeLeptonsC;
    selection invMassC;
   // selection twoLeptonsC;
    selection lepLooseVetoC;  //if Z->mumu, veto on electrons and viceversa  // could use muon or electron veto selection in AnalysisDarkMatter, but then I should define a selection* pointing to the correct one to be used (for now I define a new selection)
   selection twoLepLooseC;   //if Z->mumu, select loose muons (if Z->ee select electrons)
   selection tightLepC;
   // selection lep1tightIdIso04C;
   // selection twoLepTightC;
   // selection lep1ptC;
   // selection lep2ptC;
   // selection lep1etaC;  
   //selection lep2etaC;
   selection genLepC;  
   selection recoGenLepMatchC;
   selection genTauC;

   Int_t LEP_PDG_ID;           // choose electrons or muons
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
   Double_t HLT_LEP1PT;
   Double_t HLT_LEP2PT;
   Double_t HLT_LEP1ETA;
   Double_t HLT_LEP2ETA;

   // control samples specific parameters
   char FLAVOUR[10];                   // e.g. "ele", "mu"
   char LL_FLAVOUR[10];             // e.g. "ee", "mumu"
   char CONTROL_SAMPLE[10];   // e.g. "Z-->ee"
   Double_t invMassBinWidth;  // invariant mass histogram's bin width in GeV
   Int_t NinvMassBins;            // number of bins (computed given the range and invMassBinWidth
   // the following flag is needed to enable search for Z->ll at generator level. For MC samples different from DYJetsToLL I must not require 2 gen leptons from Z
   // unless it is Z->tautau, in which case I start from generated taus and apply selection (tau can produce muon or electron)
   Int_t using_zlljets_MCsample_flag; 
   Int_t using_ztautaujets_MCsample_flag;
    
   TH1D *HinvMass;
   TH1D *HzptDistribution;    

  };


  class wlnujetsControlSample : public AnalysisDarkMatter {
  public:

    wlnujetsControlSample(TTree *tree);  
    virtual ~wlnujetsControlSample() { std::cout<<"~wlnujetsControlSample() called"<<std::endl; }
  
    void setSelections(); 
    void setMask();
    void setHistograms();
    void setNumberParameterValue(const std::string, const Double_t);  // set some numeric variables if they are in config file
    void setControlSampleSpecificParameter();
    void setVarFromConfigFile();
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);

    selection lepLooseVetoC;  //if W->munu, veto on electrons and viceversa  // could use muon or electron veto selection in AnalysisDarkMatter, but then I should define a selection* pointing to the correct one to be used (for now I define a new selection)
   selection oneLepLooseC;   //if W->munu, select loose muon (if W->enu select electron)
   selection tightLepC;
   // selection lep1tightIdIso04C;
   // selection twoLepTightC;
   // selection lep1ptC;
   // selection lep2ptC;
   // selection lep1etaC;  
   //selection lep2etaC;
   selection genLepC;  
   selection recoGenLepMatchC;
   selection genTauC;

   Int_t LEP_PDG_ID;           // choose electrons or muons
   Double_t LEP1PT;
   //Double_t LEP2PT;
   Double_t LEP1ETA;
   //Double_t LEP2ETA;
   Double_t LEP_ISO_04;
   Double_t GENLEP1PT;
   //Double_t GENLEP2PT;
   Double_t GENLEP1ETA;
   //Double_t GENLEP2ETA;
   Double_t HLT_LEP1PT;
   //Double_t HLT_LEP2PT;
   Double_t HLT_LEP1ETA;
   //Double_t HLT_LEP2ETA;

   // control samples specific parameters
   char FLAVOUR[10];                   // e.g. "ele", "mu"
   char LL_FLAVOUR[10];             // e.g. "ee", "mumu"
   char CONTROL_SAMPLE[10];   // e.g. "Z-->ee"
      
   // the following flag is needed to enable search for W->lnu at generator level. For MC samples different from WJetsToLNu I must not require 2 gen leptons from Z
   // unless it is Z->tautau, in which case I start from generated taus and apply selection (tau can produce muon or electron)
   Int_t using_wlnujets_MCsample_flag; 
   Int_t using_wtaunujets_MCsample_flag;
        
   TH1D *HtransverseMass;

  };

  class monojet_SignalRegion : public AnalysisDarkMatter {
  public:

    monojet_SignalRegion(TTree *tree); 
    virtual ~monojet_SignalRegion() { std::cout<<"~monojet_SignalRegion() called"<<std::endl; }
  
    void setSelections(); 
    void setMask();
    void setHistograms();
    void setNumberParameterValue(const std::string, const Double_t);
    void setVarFromConfigFile();
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);

  };

  class zlljets_metResoResp : public AnalysisDarkMatter {
  public:

    zlljets_metResoResp(TTree *tree); 
    virtual ~zlljets_metResoResp() { std::cout<<"~zlljets_metResoResp() called"<<std::endl; }
  
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
  
  };



}
#endif




