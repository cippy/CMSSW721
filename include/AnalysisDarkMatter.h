#ifndef ANALYSIS_DARK_MATTER_h
#define ANALYSIS_DARK_MATTER_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <cstdlib>
#include <vector>

#include "edimarcoTree_v3.h"

class AnalysisDarkMatter : public edimarcoTree_v3 {
 public:

  AnalysisDarkMatter(TTree *tree);

  virtual ~AnalysisDarkMatter() { std::cout<<"~AnalysisDarkMatter() called"<<std::endl; }
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree);

};

#endif
