#ifndef EmanTree_Analysis_h
#define EmanTree_Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <cstdlib>
#include <vector>

#include "edimarcoTree_v2.h"

namespace myAnalyzerTEman {

  class AnalysisDarkMatter : public edimarcoTree_v2 {
  public:

  AnalysisDarkMatter(TTree *tree) : edimarcoTree_v2(tree) {      
      Init(tree);
    }

    virtual ~AnalysisDarkMatter() { std::cout<<"~AnalysisDarkMatter() called"<<std::endl; }
    
  };

  class zlljets_Axe_noSkim_light : public edimarcoTree_v2 /*,public edimarcoTreeFriend*/ {
  public:

    // only computes acceptance and efficiency in order to be faster
    zlljets_Axe_noSkim_light(TTree *tree);
    virtual ~zlljets_Axe_noSkim_light() { std::cout<<"~zlljets_Axe_noSkim_light() called"<<std::endl; }
  
    void loop(const char* configFileName);

  };

  class zlljets_resoResp : public edimarcoTree_v2 /*,public edimarcoTreeFriend*/ {
  public:

    zlljets_resoResp(TTree *tree);
    virtual ~zlljets_resoResp() { std::cout<<"~zlljets_resoResp() called"<<std::endl; }
  
    void loop(const char* configFileName, const Int_t ISDATA_FLAG, const Int_t unweighted_event_flag);

  };

  class zlljetsControlSample : public edimarcoTree_v2 {
  public:

    zlljetsControlSample(TTree *tree, const char*);  // the char* is the sample's name (e.g. QCD, Znunu, ecc...)
    virtual ~zlljetsControlSample() { std::cout<<"~zlljetsControlSample() called"<<std::endl; }
  
    void loop(const char* configFileName, const Int_t ISDATA_FLAG, const Int_t unweighted_event_flag, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    std::string suffix;
  };

  class monojet_SignalRegion : public edimarcoTree_v2 {
  public:

    monojet_SignalRegion(TTree *tree, const char*);  // the char* is the sample's name (e.g. QCD, Znunu, ecc...)
    virtual ~monojet_SignalRegion() { std::cout<<"~monojet_SignalRegion() called"<<std::endl; }
  
    void loop(const char* configFileName, const Int_t ISDATA_FLAG, const Int_t unweighted_event_flag, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    std::string suffix;
  };

  class zlljets_metResoResp : public edimarcoTree_v2 {
  public:

    zlljets_metResoResp(TTree *tree, const char*);  // the char* is the sample's name (e.g. QCD, Znunu, ecc...)
    virtual ~zlljets_metResoResp() { std::cout<<"~zlljets_metResoResp() called"<<std::endl; }
  
    void loop(const char* configFileName, const Int_t ISDATA_FLAG, const Int_t unweighted_event_flag, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    std::string suffix;
  };



}
#endif




