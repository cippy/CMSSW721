#define AnalysisDarkMatter_cxx
#include "edimarcoTree_v3.h"
#include "AnalysisDarkMatter.h"
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

using namespace std;

#ifdef AnalysisDarkMatter_cxx

AnalysisDarkMatter::AnalysisDarkMatter(TTree *tree) : edimarcoTree_v3(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);
}

Int_t AnalysisDarkMatter::GetEntry(Long64_t entry) {
  edimarcoTree_v2::GetEntry(entry);
}

Long64_t AnalysisDarkMatter::LoadTree(Long64_t entry) {
  edimarcoTree_v2::LoadTree(entry);
}

void AnalysisDarkMatter::Init(TTree *tree) {
  edimarcoTree_v2::Init(tree);
} 

#endif
