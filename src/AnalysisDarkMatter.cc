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

void AnalysisDarkMatter::setBasicConf(const char* inputSuffix, const string inputUncertainty, const char* inputConfigFileName, const Int_t inputIsDataFlag, const Int_t inputUnweightedEeventFlag, const Int_t inputHasSFfriendFlag) {

  suffix = inputSuffix;  // it is the sample name (e.g. QCD, ZJetsToNuNu ecc...)
  uncertainty = inputUncertainty; //sample uncertainty (poisson, MC, X%),
  configFileName = (char*) inputConfigFileName;
  ISDATA_FLAG = inputIsDataFlag;
  unweighted_event_flag = inputUnweightedEeventFlag;
  hasSFfriend_flag = inputHasSFfriendFlag;

}

//===============================================

void AnalysisDarkMatter::setCalibEleFlag() {

  calibEle_flag = 1;

}


//===============================================

void AnalysisDarkMatter::setNumberParameterValue(const string parameterName, const Double_t value) {

  if (parameterName == "LUMI") LUMI = value;
  //else if (parameterName == "NJETS") NJETS = value;
  else if (parameterName == "J1PT") J1PT = value;
  else if (parameterName == "J1ETA") J1ETA = value;
  //else if (parameterName == "J2PT") J2PT = value;
  //else if (parameterName == "J2ETA") J2ETA = value;
  //else if (parameterName == "J1J2DPHI") J1J2DPHI = value;
  else if (parameterName == "TAU_VETO_FLAG") TAU_VETO_FLAG = value;
  else if (parameterName == "HLT_FLAG") HLT_FLAG = value;
  else if (parameterName == "METNOLEP_START") METNOLEP_START = value;
  else if (parameterName == "MET_FILTERS_FLAG") MET_FILTERS_FLAG = value;
  else if (parameterName == "JMET_DPHI_MIN") JMET_DPHI_MIN = value;

}

//===============================================

void AnalysisDarkMatter::setVarFromConfigFile() {
  
  if (configFileName == NULL) {
    cout << "Error: configFileName points to 'NULL'. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

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

	setNumberParameterValue(parameterName, value);

      } else if (parameterType == "STRING") {
	 
	inputFile >> parameterName >> name;
	cout << right << setw(20) << parameterName << "  " << left << name << endl;

	if (parameterName == "FILENAME_BASE") {

	  FILENAME_BASE = name; 

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
  
  if (calibEle_flag == 1) outputFolder =  DIRECTORY_TO_SAVE_FILES + DIRECTORY_NAME + "_CalibEle/";
  else if (unweighted_event_flag == 1) outputFolder =  DIRECTORY_TO_SAVE_FILES + DIRECTORY_NAME + "_weq1/";
  else outputFolder =  DIRECTORY_TO_SAVE_FILES + DIRECTORY_NAME + "/";

   //Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   // Double_t metBinEdges[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 750., 850., 1000.};
   // Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;
   nMetBins = (metBinEdgesVector.size()) - 1;

   cout << "MetBinEdges: [ ";
   for(Int_t i = 0; i <= nMetBins; i++) {
     if (i != nMetBins) cout << metBinEdgesVector[i] << ", ";
     else cout << metBinEdgesVector[i] << "]" << endl;
   }

   if ( !ISDATA_FLAG && unweighted_event_flag) cout << "Warning: no weight applied to events (w = 1)" << endl;  // if MC with unit weight, make user know

   if (ISDATA_FLAG) {
     strcpy(ROOT_FNAME,(FILENAME_BASE + "_DATA.root").c_str());
     strcpy(TXT_FNAME,(FILENAME_BASE + "_DATA.txt").c_str());
     strcpy(TEX_FNAME,(FILENAME_BASE + "_DATA.tex").c_str());
   } else {
     strcpy(ROOT_FNAME,(FILENAME_BASE + "_" + suffix + ".root").c_str());
     strcpy(TXT_FNAME,(FILENAME_BASE + "_" + suffix + ".txt").c_str());
     strcpy(TEX_FNAME,(FILENAME_BASE + "_" + suffix + ".tex").c_str());
   }

}

//===============================================

void AnalysisDarkMatter::setSelections() {

  metFiltersC.set("met filters","met filters","cscfilter, ecalfilter, hbheFilterNew25ns, hbheFilterIso, Flag_eeBadScFilter");
  //metNoMuC.set("metNoMuC",Form("metNoMu > %4.0lf",METNOLEP_START),"first cut on met");
  jet1C.set("jet1",Form("jet1pt > %3.0lf",J1PT),Form("nJetClean >= 1 && JetClean1_pt > %4.0lf",(Double_t)J1PT));
  jetMetDphiMinC.set("dphiMin(j,MET)",Form("min[dphi(j,MET)] > %1.1lf",JMET_DPHI_MIN),"minimum dphi between jets and MET (using only the first 4 jets)");
  jetNoiseCleaningC.set("jet1 cleaning","noise cleaning","energy fractions (only for jet1): CH > 0.1; NH < 0.8");
  bjetVetoC.set("bjet veto","b-jets veto");
  // muonLooseVetoC.set("muonLooseVetoC","muons veto");    
  // electronLooseVetoC.set("electronLooseVetoC","electrons veto");
  if (HLT_FLAG != 0) HLTC.set("trigger","trigger");
  if (TAU_VETO_FLAG) tauLooseVetoC.set("tau veto","tau veto");
  gammaLooseVetoC.set("photon veto","photons veto");

}

//===============================================

void AnalysisDarkMatter::setHistograms() {

  HYieldsMetBin = new TH1D("HYieldsMetBin","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
   
  HvtxDistribution = new TH1D("HvtxDistribution","",40,-0.5,39.5);   
  HnjetsDistribution = new TH1D("HnjetsDistribution","njets using nJetClean30",10,-0.5,9.5);   
  Hj1j2dphiDistribution = new TH1D("Hj1j2dphiDistribution","",30,0.0,3.0);
  HjetMetDphiMinDistribution = new TH1D("HjetMetDphiMinDistribution","",30,0.0,3.2);
  Hjet1etaDistribution = new TH1D("Hjet1etaDistribution","",60,-3.0,3.0);
  Hjet2etaDistribution = new TH1D("Hjet2etaDistribution","",60,-3.0,3.0);
  HmetNoLepDistribution = new TH1D("HmetNoLepDistribution","",100,0.0,1000.0);
  Hjet1ptDistribution = new TH1D("Hjet1ptDistribution","",100,0.0,1000.0); 
  Hjet2ptDistribution = new TH1D("Hjet2ptDistribution","",100,0.0,1000.0);

  // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
  HmetBinEdges = new TH1D("HmetBinEdges","bin edges for met distributions",nMetBins+1,0.0,nMetBins+1);
  for (Int_t i = 0; i <= nMetBins; i++) {
    HmetBinEdges->SetBinContent(i+1,metBinEdgesVector[i]);
  }

}

#endif
