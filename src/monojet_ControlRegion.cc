#define monojet_ControlRegion_cxx
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

#ifdef monojet_ControlRegion_cxx

monojet_ControlRegion::monojet_ControlRegion(TTree *tree) : AnalysisDarkMatter(tree) {
  //cout <<"check in constructor "<<endl;
  //edimarcoTree_v3::Init(tree);
  // suffix = "";
  // uncertainty = "";
  // configFileName = NULL;
  // ISDATA_FLAG = 0;
  // unweighted_event_flag = 0;
  // hasSFfriend_flag = 0;
  //AnalysisDarkMatter::Init(tree);  // could also be just Init(tree)
  // maybe this is useless because the constructor already calls Init

}

#endif

//===============================================

// void monojet_ControlRegion::Init(TTree *tree) {
//   AnalysisDarkMatter::Init(tree);
// } 

//===============================================

void monojet_ControlRegion::setSelections() {

  AnalysisDarkMatter::setSelections();

  if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...

    lepLooseVetoC.set("ele veto","electrons veto");
    if (METNOLEP_START != 0) metNoLepC.set(Form("recoil > %2.0lf",METNOLEP_START),Form("metNoMu > %2.0lf",METNOLEP_START));

  } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

    lepLooseVetoC.set("muon veto","muons veto");
    if (METNOLEP_START != 0) metNoLepC.set(Form("recoil > %2.0lf",METNOLEP_START),Form("metNoEle > %2.0lf",METNOLEP_START));
  }

  selection::checkMaskLength();
  // following is directly called in the derived class
  // selection::printActiveSelections(cout); 

}

//===============================================

void monojet_ControlRegion::setHistograms() {

  AnalysisDarkMatter::setHistograms();
  
  Hlep1ptDistribution = new TH1D("Hlep1ptDistribution","",200,0.0,1000.0);
  Hlep1etaDistribution = new TH1D("Hlep1etaDistribution","",100,-5.0,5.0);

}

//===============================================

void monojet_ControlRegion::setNumberParameterValue(const std::string parameterName, const Double_t value) {

  AnalysisDarkMatter::setNumberParameterValue(parameterName, value);

  if (parameterName == "LEP_PDG_ID") LEP_PDG_ID = (value < 0) ? (-0.5 + value) : (0.5 + value);
  else if (parameterName == "LEP1PT") LEP1PT = value;
  else if (parameterName == "LEP_ISO_04") LEP_ISO_04 = value;
  else if (parameterName == "LEP1ETA") LEP1ETA = value;
  else if (parameterName == "HLT_LEP1PT") HLT_LEP1PT = value;
  else if (parameterName == "HLT_LEP1ETA") HLT_LEP1ETA = value;

  if (!ISDATA_FLAG) {

    if (parameterName == "GENLEP_TAG") GENLEP_TAG = (value < 0) ? (-0.5 + value) : (0.5 + value);
    else if (parameterName == "GENLEP1PT") GENLEP1PT = value;
    else if (parameterName == "GENLEP1ETA") GENLEP1ETA = value;

  }

}

//===============================================

void monojet_ControlRegion::setControlSampleSpecificParameter() {

  // no common specific parameter among control samples for now
  
}

//===============================================

void monojet_ControlRegion::setVarFromConfigFile() {

  AnalysisDarkMatter::setVarFromConfigFile();
  //monojet_ControlRegion::setControlSampleSpecificParameter();

}

//===============================================





