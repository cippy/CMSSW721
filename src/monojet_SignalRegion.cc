#define monojet_SignalRegion_cxx
#include "EmanTreeAnalysis.h"
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

#ifdef monojet_SignalRegion_cxx

monojet_SignalRegion::monojet_SignalRegion(TTree *tree) : AnalysisDarkMatter(tree) {
  //cout <<"check in constructor "<<endl;
  // suffix = "";
  // uncertainty = "";
  // configFileName = NULL;
  // ISDATA_FLAG = 0;
  // unweighted_event_flag = 0;
  // hasSFfriend_flag = 0;
  // Init(tree);

}

#endif

//===============================================

void monojet_SignalRegion::setSelections() {

  AnalysisDarkMatter::setSelections();

  metNoLepC.set(Form("recoil > %2.0lf",METNOLEP_START),Form("metNoMu > %4.0lf",METNOLEP_START),"first cut on met");
  muonLooseVetoC.set("muon veto","muons veto");
  electronLooseVetoC.set("ele veto","electrons veto");

  selection::checkMaskLength();
  selection::printActiveSelections(cout); 

}

//===============================================

void monojet_SignalRegion::setMask() {

  analysisMask.setName("monojet signal selection");

  if (HLT_FLAG != 0) analysisMask.append(HLTC.get2ToId());
  if (MET_FILTERS_FLAG != 0) analysisMask.append(metFiltersC.get2ToId());
  analysisMask.append(muonLooseVetoC.get2ToId());
  analysisMask.append(electronLooseVetoC.get2ToId());
  if (TAU_VETO_FLAG) analysisMask.append(tauLooseVetoC.get2ToId());
  analysisMask.append(gammaLooseVetoC.get2ToId());
  analysisMask.append(bjetVetoC.get2ToId());
  if (METNOLEP_START != 0) analysisMask.append(metNoLepC.get2ToId());  
  analysisMask.append(jet1C.get2ToId());
  analysisMask.append(jetNoiseCleaningC.get2ToId());
  analysisMask.append(jetMetDphiMinC.get2ToId());

  analysisSelectionManager.SetMaskPointer(&analysisMask);

  if (HLT_FLAG != 0) analysisSelectionManager.append(&HLTC);
  if (MET_FILTERS_FLAG != 0) analysisSelectionManager.append(&metFiltersC);
  analysisSelectionManager.append(&muonLooseVetoC);
  analysisSelectionManager.append(&electronLooseVetoC);
  if (TAU_VETO_FLAG) analysisSelectionManager.append(&tauLooseVetoC);
  analysisSelectionManager.append(&gammaLooseVetoC);
  analysisSelectionManager.append(&bjetVetoC);
  analysisSelectionManager.append(&metNoLepC); 
  analysisSelectionManager.append(&jet1C);
  analysisSelectionManager.append(&jetNoiseCleaningC);
  analysisSelectionManager.append(&jetMetDphiMinC);

}

//===============================================

void monojet_SignalRegion::setHistograms() {

  AnalysisDarkMatter::setHistograms();

  // following histograms are really needed only when using Zvv or Wlv samples, so they are not instantiated in AnalysisDarkMatter::setHistograms();

  if (suffix == "ZJetsToNuNu" || suffix == "WJetsToLNu") {

    HYieldsMetBin_qcdRenScaleUp = new TH1D("HYieldsMetBin_qcdRenScaleUp","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
    HYieldsMetBin_qcdRenScaleDown = new TH1D("HYieldsMetBin_qcdRenScaleDown","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
    HYieldsMetBin_qcdFacScaleUp = new TH1D("HYieldsMetBin_qcdFacScaleUp","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
    HYieldsMetBin_qcdFacScaleDown = new TH1D("HYieldsMetBin_qcdFacScaleDown","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
    HYieldsMetBin_qcdPdfUp = new TH1D("HYieldsMetBin_qcdPdfUp","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
    HYieldsMetBin_qcdPdfDown = new TH1D("HYieldsMetBin_qcdPdfDown","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
    HYieldsMetBin_ewkUp = new TH1D("HYieldsMetBin_ewkUp","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
    HYieldsMetBin_ewkDown = new TH1D("HYieldsMetBin_ewkDown","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
    
    HSyst_qcdRenScale = new TH1D("HSyst_qcdRenScale","systematic uncertainty for QCD renormalization scale",nMetBins,metBinEdgesVector.data());
    HSyst_qcdFacScale = new TH1D("HSyst_qcdFacScale","systematic uncertainty for QCD factorization scale",nMetBins,metBinEdgesVector.data());
    HSyst_qcdPdf = new TH1D("HSyst_qcdPdf","systematic uncertainty for QCD due to PDF uncertainty",nMetBins,metBinEdgesVector.data());
    HSyst_ewk = new TH1D("HSyst_ewk","systematic uncertainty for EWK",nMetBins,metBinEdgesVector.data());
    HSyst_total = new TH1D("HSyst_total","total systematic uncertainty (sum in quadrature of all single systematics)",nMetBins,metBinEdgesVector.data());

  }

  //histograms specific to monojet selection (but not to control regions) must be set here as in AnalysisDarkMatter::setHistograms()

}

//===============================================

void monojet_SignalRegion::setNumberParameterValue(const std::string parameterName, const Double_t value) {
	 	 
  AnalysisDarkMatter::setNumberParameterValue(parameterName, value);

  //parameters specific to monojet selection (but not to control regions) must be set here as in AnalysisDarkMatter::setNumberParameterValue()

}

//===============================================

void monojet_SignalRegion::setVarFromConfigFile() {

  AnalysisDarkMatter::setVarFromConfigFile();

}

//===============================================

Double_t monojet_SignalRegion::computeEventWeight() {

  if (ISDATA_FLAG || unweighted_event_flag) return 1.0;
  else {
    // sf_nlo_weight = (*ptr_sf_nlo_QCD) * (*ptr_sf_nlo_EWK);
    // return LUMI * weight * vtxWeight * SF_BTag * sf_nlo_weight; //SF_BTag is in evVarFriend, not sfFriend
    return LUMI * weight * vtxWeight * SF_BTag * SF_NLO_QCD * SF_NLO_EWK; //SF_BTag is in evVarFriend, not sfFriend
  }

}

//===============================================

void monojet_SignalRegion::loop(vector< Double_t > &yRow, vector< Double_t > &eRow, vector< Double_t > &uncRow)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   //fChain->SetBranchStatus("LHEorigWeight",1); // contains negative values: the weight in the event is weight*LHEorigWeight

   //fChain->SetBranchStatus("genWeight",1); 

   fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
   //fChain->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (pt > 20 + everything else)
   //fChain->SetBranchStatus("nTau18V",1);
   fChain->SetBranchStatus("nTauClean18V",1);

   fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
   fChain->SetBranchStatus("nJetClean30",1);    // # of jet with pt > 30 & eta < 2.5 and cleaning for against muons misidentified as PFjets   
   fChain->SetBranchStatus("JetClean_pt",1);  
   fChain->SetBranchStatus("JetClean_eta",1);  
   
   fChain->SetBranchStatus("nLepGood",1);
   fChain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus
   fChain->SetBranchStatus("LepGood_pt",1);
   fChain->SetBranchStatus("LepGood_eta",1);
   fChain->SetBranchStatus("LepGood_phi",1);   
   fChain->SetBranchStatus("LepGood_mass",1);
   //fChain->SetBranchStatus("LepGood_charge",1);
   fChain->SetBranchStatus("LepGood_tightId",1);
   fChain->SetBranchStatus("LepGood_relIso04",1);
   // fChain->SetBranchStatus("ngenLep",1);         // not needed, using GenPart to study generator level quantities
   // fChain->SetBranchStatus("genLep_pdgId",1);
   // fChain->SetBranchStatus("genLep_pt",1);
   // fChain->SetBranchStatus("genLep_eta",1);
   // fChain->SetBranchStatus("genLep_phi",1);
   fChain->SetBranchStatus("mZ1",1);  // best m(ll) SF/OS

   //fChain->SetBranchStatus("met_pt",1);
   //fChain->SetBranchStatus("met_eta",1);
   //fChain->SetBranchStatus("met_phi",1);

   fChain->SetBranchStatus("metNoMu_pt",1);
   //fChain->SetBranchStatus("metNoMu_eta",1);
   fChain->SetBranchStatus("metNoMu_phi",1);

   fChain->SetBranchStatus("nVert",1);  // number of good vertices 
   fChain->SetBranchStatus("HLT_MonoJetMetNoMuMHT90",1);
   fChain->SetBranchStatus("HLT_MonoJetMetNoMuMHT120",1);

   // met filters to be used (the config file has a parameter saying whether they should be used or not)
   fChain->SetBranchStatus("cscfilter",1);
   fChain->SetBranchStatus("ecalfilter",1);
   fChain->SetBranchStatus("hbheFilterNew25ns",1);
   fChain->SetBranchStatus("hbheFilterIso",1);
   fChain->SetBranchStatus("Flag_eeBadScFilter",1);

   //added on November 2015. These are new variables (except for weight, which has just changed in the definition)
   fChain->SetBranchStatus("nBTag15",1);  // for b-jet veto
   fChain->SetBranchStatus("dphijm",1);          // flag for dphi minimum between met and any of the jets in the event (using only the first four jets
   fChain->SetBranchStatus("weight",1);   // modified since 17 November 2015: now it includes the whol weight, e.g. 1000*xsec*genWeight ...
   fChain->SetBranchStatus("events_ntot",1);    // equivalent to SUMWEIGHTS for samples before 17 November 2015
   fChain->SetBranchStatus("JetClean_leadClean",1); // has new cleaning on energy fractions (added on 17 November 2015) 

   //added on 23/01/2016
   fChain->SetBranchStatus("nEle40T",1);
   fChain->SetBranchStatus("nCalibEle",1);
   fChain->SetBranchStatus("CalibEle_pt",1);
   fChain->SetBranchStatus("CalibEle_energy",1);
   fChain->SetBranchStatus("CalibEle_eta",1);
   fChain->SetBranchStatus("CalibEle_phi",1);
   fChain->SetBranchStatus("CalibEle_mass",1);

   if (!ISDATA_FLAG) {
     // fChain->SetBranchStatus("nGenPart",1);
     // fChain->SetBranchStatus("GenPart_pdgId",1);
     // fChain->SetBranchStatus("GenPart_motherId",1);
     // fChain->SetBranchStatus("GenPart_pt",1);
     // fChain->SetBranchStatus("GenPart_eta",1);
     // fChain->SetBranchStatus("GenPart_phi",1);
     // fChain->SetBranchStatus("GenPart_mass",1);
     // fChain->SetBranchStatus("GenPart_motherIndex",1);

     //fChain->SetBranchStatus("vtxW",1);   // weight to have better agreement between data and MC (will not be always used)  ** NOW OBSOLETE **
     fChain->SetBranchStatus("vtxWeight",1);   // weight to have better agreement between data and MC (added the 10th of October, substituting vtxW in the new set of trees, keeping both for backward compatibility)
     //fChain->SetBranchStatus("xsec",1);  
 
     //added on 23/01/2016
     fChain->SetBranchStatus("SF_BTag",1);
     //variables in sfFriend tree
     fChain->SetBranchStatus("SF_trig1lep",1);
     fChain->SetBranchStatus("SF_trigmetnomu",1);
     fChain->SetBranchStatus("SF_LepTightLoose",1);
     fChain->SetBranchStatus("SF_LepTight",1);
     fChain->SetBranchStatus("SF_NLO",1);
     fChain->SetBranchStatus("SF_NLO_QCD",1);
     fChain->SetBranchStatus("SF_NLO_QCD_renScaleUp",1);
     fChain->SetBranchStatus("SF_NLO_QCD_renScaleDown",1);
     fChain->SetBranchStatus("SF_NLO_QCD_facScaleUp",1);
     fChain->SetBranchStatus("SF_NLO_QCD_facScaleDown",1);
     fChain->SetBranchStatus("SF_NLO_QCD_pdfUp",1);
     fChain->SetBranchStatus("SF_NLO_QCD_pdfDown",1);
     fChain->SetBranchStatus("SF_NLO_EWK",1);
     fChain->SetBranchStatus("SF_NLO_EWK_up",1);
     fChain->SetBranchStatus("SF_NLO_EWK_down",1);

   }

   //cout << "CHECK1 " << endl;
   setVarFromConfigFile();
   //cout << "CHECK2 " << endl;
   setSelections();
   //cout << "CHECK3 " << endl;
   setMask();
   //cout << "CHECK4 " << endl;
   //set_SF_NLO_pointers(sf_nlo, ptr_sf_nlo_QCD, ptr_sf_nlo_EWK);

   cout << "Opening file " <<ROOT_FNAME<< " in folder " << outputFolder << endl;

   TFile *rootFile = new TFile((outputFolder + ROOT_FNAME).c_str(),"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout<<"Error: file \""<<outputFolder + ROOT_FNAME<<"\" was not opened."<<endl;
     exit(EXIT_FAILURE);
   }
 

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   //TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   setHistograms();

   // Double_t nTotalWeightedEvents = 0.0;     
   // // deciding  what is the event weight
   // Double_t newwgt;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"monojet_SignalRegion::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;

   for (Int_t jentry=0; jentry<nentries; jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     UInt_t eventMask = 0; 

     if (jentry%500000 == 0) {
       cout << "entry: " << jentry << endl;
     }

     // if (ISDATA_FLAG || unweighted_event_flag) newwgt = 1.0;
     // else {

     //   newwgt = LUMI * vtxWeight * weight * SF_trigmetnomu * SF_BTag;/* / events_ntot*/;  // starting from 17 November, "events_ntot" substitutes SUMWEIGHT and is already present in the trees. Same for weight, which is now defined as "1000 * xsec * genWeight" (1000*xsec is the cross section in fb, since xsec is in pb.)

     // }

     newwgt = computeEventWeight();
     nTotalWeightedEvents += newwgt;  // counting events with weights

     // beginning of eventMask building
     
     //eventMask += jet1C.addToMask(nJetClean30 >= 1 && JetClean_pt[0] > J1PT && fabs(JetClean_eta[0] < J1ETA && jetclean1 > 0.5));  //could skip cut on eta 
     eventMask += jet1C.addToMask(nJetClean30 >= 1 && JetClean_pt[0] > J1PT /* && fabs(JetClean_eta[0]) < J1ETA*/);
     eventMask += jetMetDphiMinC.addToMask(fabs(dphijm > JMET_DPHI_MIN));
     eventMask += jetNoiseCleaningC.addToMask(JetClean_leadClean[0] > 0.5);
     eventMask += bjetVetoC.addToMask(nBTag15 < 0.5);
     eventMask += muonLooseVetoC.addToMask(nMu10V < 0.5);
     eventMask += electronLooseVetoC.addToMask(nEle10V < 0.5);
     eventMask += tauLooseVetoC.addToMask(nTauClean18V < 0.5);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V < 0.5);
     eventMask += metNoLepC.addToMask(metNoMu_pt > METNOLEP_START);
     eventMask += metFiltersC.addToMask(cscfilter == 1 && ecalfilter == 1 && hbheFilterNew25ns == 1 && hbheFilterIso == 1 && Flag_eeBadScFilter > 0.5);
     eventMask += HLTC.addToMask(HLT_MonoJetMetNoMuMHT90 > 0.5 || HLT_MonoJetMetNoMuMHT120 > 0.5); //HLT_* variables are stored as float, so using "== 1" might yield unexpected results
     
     // end of eventMask building

     analysisMask.countEvents(eventMask,newwgt);

     if ( ((eventMask & analysisMask.globalMask.back()) == analysisMask.globalMask.back()) ) {
       
       // this histogram holds the final yields in bins of MET
	 HYieldsMetBin->Fill(metNoMu_pt,newwgt);
	 
	 HmetNoLepDistribution->Fill(metNoMu_pt,newwgt);
	 HvtxDistribution->Fill(nVert,newwgt);
	 HnjetsDistribution->Fill(nJetClean30,newwgt);
	 Hjet1etaDistribution->Fill(JetClean_eta[0],newwgt);
	 Hjet1ptDistribution->Fill(JetClean_pt[0],newwgt);
	 HjetMetDphiMinDistribution->Fill(dphijm,newwgt);
	 if (nJetClean30 >= 2) {
	   Hj1j2dphiDistribution->Fill(dphijj,newwgt);
	   Hjet2etaDistribution->Fill(JetClean_eta[1],newwgt);
	   Hjet2ptDistribution->Fill(JetClean_pt[1],newwgt);
	 }

	 if (suffix == "ZJetsToNuNu" || suffix == "WJetsToLNu") {

	   HYieldsMetBin_qcdRenScaleUp->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_renScaleUp/ SF_NLO_QCD));
	   HYieldsMetBin_qcdRenScaleDown->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_renScaleDown/ SF_NLO_QCD));
	   HYieldsMetBin_qcdFacScaleUp->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_facScaleUp/ SF_NLO_QCD));
	   HYieldsMetBin_qcdFacScaleDown->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_facScaleDown/ SF_NLO_QCD));
	   HYieldsMetBin_qcdPdfUp->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_pdfUp/ SF_NLO_QCD));
	   HYieldsMetBin_qcdPdfDown->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_pdfDown/ SF_NLO_QCD));
	   HYieldsMetBin_ewkUp->Fill(metNoMu_pt,(newwgt * SF_NLO_EWK_up/ SF_NLO_EWK));
	   HYieldsMetBin_ewkDown->Fill(metNoMu_pt,(newwgt * SF_NLO_EWK_down/ SF_NLO_EWK));

	 }


     }
       
   }                        // end of loop on entries

   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &analysisMask);

   mySpaces(cout,2);
   myPrintYieldsMetBinInStream(cout, HYieldsMetBin, metBinEdgesVector.data(), nMetBins);
 
   cout<<"creating file '"<<TXT_FNAME<<"' in folder "<< outputFolder <<" ..."<<endl;
   ofstream myfile((outputFolder + TXT_FNAME).c_str(),ios::out);

   if ( !myfile.is_open() ) {

     cout<<"Error: unable to open file "<<TXT_FNAME<<" !"<<endl;
     exit(EXIT_FAILURE);
     
   }

  
   if (!ISDATA_FLAG && unweighted_event_flag) myfile << "======   Using unweighted events (w = 1)   ======" << endl;
   mySpaces(myfile,3);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &analysisMask);
   mySpaces(myfile,2);
   myPrintYieldsMetBinInStream(myfile, HYieldsMetBin, metBinEdgesVector.data(), nMetBins);

   myfile.close();

   // filling with yields and efficiency: I will use efficiency with respect to total and not to previous step, but I could make this choice in the config file

   // entry point
   // yRow.push_back(nTotalWeightedEvents);
   // eRow.push_back(1.0000);
   // uncRow.push_back(sqrt(nTotalWeightedEvents)); //should use a kind of myGetUncertainty function, but I don't save sum of newwgt^2 so I can't use MC uncertainty

   // for(Int_t i = 0; i < selStep.size(); i++) {
  
   //   yRow.push_back(analysisMask.nEvents[selStep[i]]);
   //   uncRow.push_back(myGetUncertainty(&analysisMask, selStep[i], uncertainty));
   //   if (i == 0) eRow.push_back(analysisMask.nEvents[selStep[i]]/nTotalWeightedEvents);
   //   else if( (i != 0) && (analysisMask.nEvents[selStep[i]-1] == 0) ) eRow.push_back(1.0000);
   //   else eRow.push_back(analysisMask.nEvents[selStep[i]]/analysisMask.nEvents[selStep[i]-1]);

   // }

   
   // entry point  
   yRow.push_back(nTotalWeightedEvents); // [0] 
   eRow.push_back(1.0000);
   uncRow.push_back(sqrt(nTotalWeightedEvents)); //should use a kind of myGetUncertainty function, but I don't save sum of newwgt^2 so I can't use MC uncertainty
  

   for(Int_t i = 0; i < analysisSelectionManager.getVectorSize(); i++) {
     
     yRow.push_back(analysisMask.nEvents[analysisSelectionManager.getStepIndex(i)]);
     //uncRow.push_back(sqrt(yRow.back()));
     uncRow.push_back(myGetUncertainty(&analysisMask, analysisSelectionManager.getStepIndex(i), uncertainty));
     if (i == 0) eRow.push_back(analysisMask.nEvents[analysisSelectionManager.getStepIndex(i)]/nTotalWeightedEvents);
     else if( (i != 0) && (analysisMask.nEvents[analysisSelectionManager.getStepIndex(i)-1] == 0) ) eRow.push_back(1.0000);
     else eRow.push_back(analysisMask.nEvents[analysisSelectionManager.getStepIndex(i)]/analysisMask.nEvents[analysisSelectionManager.getStepIndex(i)-1]);
   
   }


   // filling last bin with overflow
   myAddOverflowInLastBin(HYieldsMetBin);
   myAddOverflowInLastBin(HmetNoLepDistribution);
   myAddOverflowInLastBin(Hjet1ptDistribution);
   myAddOverflowInLastBin(Hjet2ptDistribution);

   if (suffix == "ZJetsToNuNu" || suffix == "WJetsToLNu") {

     createSystematicsHistogram();
     // following is in function above
     /*
     myAddOverflowInLastBin(HYieldsMetBin_qcdRenScaleUp);
     myAddOverflowInLastBin(HYieldsMetBin_qcdRenScaleDown);
     myAddOverflowInLastBin(HYieldsMetBin_qcdFacScaleUp);
     myAddOverflowInLastBin(HYieldsMetBin_qcdFacScaleDown);
     myAddOverflowInLastBin(HYieldsMetBin_qcdPdfUp);
     myAddOverflowInLastBin(HYieldsMetBin_qcdPdfDown);
     myAddOverflowInLastBin(HYieldsMetBin_ewkUp);
     myAddOverflowInLastBin(HYieldsMetBin_ewkDown);

     // computing systematic uncertainties and saving them as histograms.
     myBuildSystematicsHistogram(HSyst_qcdRenScale, HYieldsMetBin, HYieldsMetBin_qcdRenScaleUp, HYieldsMetBin_qcdRenScaleDown);
     myBuildSystematicsHistogram(HSyst_qcdFacScale, HYieldsMetBin, HYieldsMetBin_qcdFacScaleUp, HYieldsMetBin_qcdFacScaleDown);
     myBuildSystematicsHistogram(HSyst_qcdPdf, HYieldsMetBin, HYieldsMetBin_qcdPdfUp, HYieldsMetBin_qcdPdfDown);
     myBuildSystematicsHistogram(HSyst_ewk, HYieldsMetBin, HYieldsMetBin_ewkUp, HYieldsMetBin_ewkDown);

     // define an empty histogram to sum uncertainties in a clean way
     TH1D *Htmp = new TH1D("Htmp","",nMetBins,metBinEdgesVector.data());
     vector<TH1D*> hptr;
     hptr.push_back(HSyst_qcdRenScale);
     hptr.push_back(HSyst_qcdFacScale);
     hptr.push_back(HSyst_qcdPdf);
     hptr.push_back(HSyst_ewk);
     
     for (Int_t i = 0; i < hptr.size(); i++) {
       Htmp->Multiply(hptr[i],hptr[i]); // square of bin content for each single systematic histogram
       HSyst_total->Add(Htmp);             // adding the squares
     }

     for (Int_t i = 0; i <= (HSyst_total->GetNbinsX() + 1); i++) {  // computing square root of each bin's content (from underflow to overflow bin, but they should be empty)
       HSyst_total->SetBinContent(i, sqrt(HSyst_total->GetBinContent(i)));
     }

     delete Htmp;
     */   

   }

   rootFile->Write();

   rootFile->Close();
   delete rootFile;

   //creating a .tex file to build tables with data
   myCreateTexTable(TEX_FNAME, outputFolder, LUMI,nTotalWeightedEvents, &analysisMask);

   // end of tex file

   mySpaces(cout,2);

}




