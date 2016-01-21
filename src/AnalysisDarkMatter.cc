#define AnalysisDarkMatter_cxx
#include "edimarcoTree_v3.h"
#include "AnalysisDarkMatter.h"
//#include "functionsForAnalysis.h"
//#include "myClasses.h"
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

//===============================================

AnalysisDarkMatter::AnalysisDarkMatter(TTree *tree) : edimarcoTree_v3(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);
}

//===============================================

Int_t AnalysisDarkMatter::GetEntry(Long64_t entry) {
  edimarcoTree_v3::GetEntry(entry);
}

//===============================================

Long64_t AnalysisDarkMatter::LoadTree(Long64_t entry) {
  edimarcoTree_v3::LoadTree(entry);
}

//===============================================

void AnalysisDarkMatter::Init(TTree *tree) {
  edimarcoTree_v3::Init(tree);
} 

//===============================================

void AnalysisDarkMatter::SetVarFromConfigFile() {
  
  ifstream inputFile(configFileName);

  if (inputFile.is_open()) {

    Double_t value;
    string name;
    string parameterName;
    string parameterType;

    //mySpaces(cout,2);
    cout <<endl; cout <<endl;
    cout << "Printing content of " << configFileName << " file" << endl;
    //mySpaces(cout,1);
    cout <<endl;

    while (inputFile >> parameterType ) {

      if (parameterType == "NUMBER") {

	inputFile >> parameterName >> value;
	cout << right << setw(20) << parameterName << "  " << left << value << endl;

	if (parameterName == "LUMI") LUMI = value;
	//else if (parameterName == "NJETS") NJETS = value;
	else if (parameterName == "J1PT") J1PT = value;
	else if (parameterName == "J1ETA") J1ETA = value;
	//else if (parameterName == "J2PT") J2PT = value;
	//else if (parameterName == "J2ETA") J2ETA = value;
	//else if (parameterName == "J1J2DPHI") J1J2DPHI = value;
	else if (parameterName == "TAU_VETO_FLAG") TAU_VETO_FLAG = value;
	else if (parameterName == "METNOLEP_START") METNOLEP_START = value;
	else if (parameterName == "MET_FILTERS_FLAG") MET_FILTERS_FLAG = value;
	else if (parameterName == "JMET_DPHI_MIN") JMET_DPHI_MIN = value;

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
     
    //mySpaces(cout,2);
    cout <<endl; cout <<endl;

    inputFile.close();
                                                                                                                         
  } else {

    cout << "Error: could not open file " << configFileName << endl;
    exit(EXIT_FAILURE);
    
  }
  
  outputFolder =  DIRECTORY_TO_SAVE_FILES + DIRECTORY_NAME + "/";

   //Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   // Double_t metBinEdges[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 750., 850., 1000.};
   // Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;
   nMetBins = (metBinEdgesVector.size()) - 1;

   cout << "MetBinEdges: [ ";
   for(Int_t i = 0; i <= nMetBins; i++) {
     if (i != nMetBins) cout << metBinEdgesVector[i] << ", ";
     else cout << metBinEdgesVector[i] << "]" << endl;
   }

}

//===============================================

void AnalysisDarkMatter::setSelections() {

    metFiltersC.set("metFiltersC","met filters","cscfilter, ecalfilter, hbheFilterNew25ns, hbheFilterIso");
    metNoMuC.set("metNoMuC",Form("metNoMu > %4.0lf",METNOLEP_START),"first cut on met");
    jet1C.set("jet1C","jet1",Form("nJetClean >= 1 && JetClean1_pt > %4.0lf",(Double_t)J1PT));
    jetMetDphiMinC.set("jetMetDphiMinC",Form("min[dphi(jets,MET)] > %1.1lf",JMET_DPHI_MIN),"minimum dphi between jets and MET (using only the first 4 jets)");
    jetNoiseCleaningC.set("jetNoiseCleaningC","noise cleaning","energy fractions (only for jet1): CH > 0.1; NH < 0.8");
    bjetVetoC.set("bjetVetoC","b-jets veto");
    muonLooseVetoC.set("muonLooseVetoC","muons veto");
    electronLooseVetoC.set("electronLooseVetoC","electrons veto");
    if (TAU_VETO_FLAG) tauLooseVetoC.set("tauLooseVetoC","tau veto");
    gammaLooseVetoC.set("gammaLooseVetoC","photons veto");

  }

#endif
