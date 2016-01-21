#ifndef ANALYSIS_DARK_MATTER_h
#define ANALYSIS_DARK_MATTER_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>

#include "edimarcoTree_v3.h"
//#include "functionsForAnalysis.h"
#include "myClasses.h"

class AnalysisDarkMatter : public edimarcoTree_v3 {
 public:

  AnalysisDarkMatter(TTree *tree);

  virtual ~AnalysisDarkMatter() { std::cout<<"~AnalysisDarkMatter() called"<<std::endl; }
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree);

  //common selections (between signal and control regions) are declared here
  selection metFiltersC;
  selection metNoMuC;
  selection jet1C;
  selection jetMetDphiMinC;
  selection jetNoiseCleaningC;
  selection bjetVetoC;
  selection muonLooseVetoC; 
  selection electronLooseVetoC; 
  selection tauLooseVetoC;
  selection gammaLooseVetoC;

  virtual void SetVarFromConfigFile();
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
  // Int_t HLT_FLAG;                  // not needed: monojet default trigger paths should be 100% efficient wrt offline selection
  Double_t METNOLEP_START;
  Int_t MET_FILTERS_FLAG;
  Double_t JMET_DPHI_MIN;
  std::string FILENAME_BASE;
  std::string DIRECTORY_TO_SAVE_FILES;
  std::string DIRECTORY_NAME;

  std::string outputFolder;

  std::vector<Double_t> metBinEdgesVector;  // filled with values in file named configFileName
  Int_t nMetBins;

  std::string suffix;
  std::string uncertainty;
  char* configFileName;
  Int_t ISDATA_FLAG;
  Int_t unweighted_event_flag;
    

};

#endif
