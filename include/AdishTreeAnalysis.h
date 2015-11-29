#ifndef AdishTree_Analysis_h
#define AdishTree_Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <cstdlib>

#include "adishTree.h"

namespace myAnalyzerTAdish {

  class znunujetsAnaTAdish : public adishTree {
  public:

    znunujetsAnaTAdish(TTree *tree);
    virtual ~znunujetsAnaTAdish() { std::cout<<"~znunujetsAnaTAdish() called"<<std::endl; }
  
    void loop();

  };

  class zmumujetsAnaTAdish : public adishTree {
  public:

    zmumujetsAnaTAdish(TTree *tree);
    virtual ~zmumujetsAnaTAdish() { std::cout<<"~zmumujetsAnaTAdish() called"<<std::endl; }
  
    void loop();

  };

  class zlljets_resoResp_forAdish : public adishTree {
  public:

    zlljets_resoResp_forAdish(TTree *tree);
    virtual ~zlljets_resoResp_forAdish() { std::cout<<"~zlljets_resoResp_forAdish() called"<<std::endl; }
  
    void loop(const char* configFileName, const Int_t ISDATA_FLAG, const Int_t unweighted_event_flag);

  };

}
#endif




