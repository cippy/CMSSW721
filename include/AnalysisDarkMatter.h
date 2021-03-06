#ifndef ANALYSIS_DARK_MATTER_h
#define ANALYSIS_DARK_MATTER_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>

#include "edimarcoTree_v3.h"
//#include "functionsForAalysis.h"
#include "myClasses.h"

class AnalysisDarkMatter : public edimarcoTree_v3 {
 public:

  AnalysisDarkMatter(TTree *tree);

  virtual ~AnalysisDarkMatter() { std::cout<<"~AnalysisDarkMatter() called"<<std::endl; }
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree);

  //common selections (between signal and control regions) are declared here
  selection HLTC;
  selection metFiltersC;
  selection metNoLepC;
  selection jet1C;
  selection jetMetDphiMinC;
  selection jetNoiseCleaningC;
  selection bjetVetoC;
  selection muonLooseVetoC; 
  selection electronLooseVetoC; 
  selection tauLooseVetoC;
  selection gammaLooseVetoC;

  mask analysisMask;
  selectionManager analysisSelectionManager;

  virtual void setBasicConf(const char* inputSuffix, const std::string inputUncertainty, const char* inputConfigFileName, const Int_t inputIsDataFlag, const Int_t inputUnweightedEeventFlag, const Int_t inputHasSFfriendFlag);
  virtual void setCalibEleFlag();
  virtual void setNumberParameterValue(const std::string, const Double_t);
  virtual void setVarFromConfigFile();
  virtual void setSelections();

  char ROOT_FNAME[100];
  char TXT_FNAME[100];
  char TEX_FNAME[100];

  Double_t LUMI;
  //Int_t NJETS;
  Double_t J1PT;
  Double_t J1ETA;
  //Double_t J2PT;
  //Double_t J2ETA;
  //Double_t J1J2DPHI;
  Int_t TAU_VETO_FLAG;
  Int_t HLT_FLAG;                  // usage depends on specific analysis
  Double_t METNOLEP_START;
  Int_t MET_FILTERS_FLAG;
  Double_t JMET_DPHI_MIN;
  std::string FILENAME_BASE;
  std::string DIRECTORY_TO_SAVE_FILES;
  std::string DIRECTORY_NAME;

  std::string outputFolder;

  std::vector<Double_t> metBinEdgesVector;  // filled with values in file named configFileName

  //obsolete, substituted by selectionManager class
  std::vector<Int_t> selStep;  //set when setting the mask analysisMask 
  //array to store index of step to form selection flow (might want to consider two or more steps together and not separated)
   // in case a step wold be present for some sample but not for others (e.g. the RecoGen match done only in Zll MC), the step is referred to as -1 and the corresponding values are set to -1, so that, when printing the table, yields will be filled with " / / " which means " uneffected" (because that step was not done)

  Int_t nMetBins;

  // the following can be set by setBasicConf(), but are initialized in construcor
  std::string suffix;   // it is the sample name (e.g. QCD, ZJetsToNuNu ecc...)
  std::string uncertainty; //uncertainty on yields (from MC, Poisson, X%)
  char* configFileName;
  Int_t ISDATA_FLAG;
  Int_t unweighted_event_flag;
  Int_t hasSFfriend_flag;  //tells if sfFriend are present

  Int_t calibEle_flag; // tells if using calibrated electron properties or traditional ones
    
  //root histograms: these are common among all analysis (signal ad control region)
  TH1D *HYieldsMetBin = NULL;
  TH1D *HvtxDistribution = NULL;   
  TH1D *HnjetsDistribution = NULL;   
  TH1D *Hj1j2dphiDistribution = NULL;
  TH1D *HjetMetDphiMinDistribution = NULL;
  TH1D *Hjet1etaDistribution = NULL;
  TH1D *Hjet2etaDistribution = NULL;
  TH1D *HmetNoLepDistribution = NULL;
  TH1D *Hjet1ptDistribution = NULL;
  TH1D *Hjet2ptDistribution = NULL;
  TH1D *HmetBinEdges = NULL;

  virtual void setHistograms();

};

#endif
