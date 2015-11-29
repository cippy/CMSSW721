#define zlljets_resoResp_forAdish_cxx
#include "AdishTreeAnalysis.h"
//C or C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
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
using namespace myAnalyzerTAdish;

#ifdef zlljets_resoResp_forAdish_cxx

zlljets_resoResp_forAdish::zlljets_resoResp_forAdish(TTree *tree) : adishTree(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);

}

#endif

void zlljets_resoResp_forAdish::loop(const char* configFileName, const Int_t ISDATA_FLAG, const Int_t unweighted_event_flag)
{

  if (fChain == 0) return;

  fChain->SetBranchStatus("*",0);  

  fChain->SetBranchStatus("wgt",1);
  fChain->SetBranchStatus("xsec",1);
  if (!ISDATA_FLAG) fChain->SetBranchStatus("wgtsum",1);
  fChain->SetBranchStatus("kfact",1);

  fChain->SetBranchStatus("nvtx",1);
  fChain->SetBranchStatus("zmass",1);
  fChain->SetBranchStatus("zpt",1);
  fChain->SetBranchStatus("zeta",1);
  fChain->SetBranchStatus("zphi",1);
  fChain->SetBranchStatus("pfmet",1);
  fChain->SetBranchStatus("pfmetphi",1);
  fChain->SetBranchStatus("metnohf",1);
  fChain->SetBranchStatus("metnohfphi",1);
  fChain->SetBranchStatus("hmet",1);
  fChain->SetBranchStatus("hmetphi",1);
  fChain->SetBranchStatus("amet",1);
  fChain->SetBranchStatus("ametphi",1);
  fChain->SetBranchStatus("bmet",1);
  fChain->SetBranchStatus("bmetphi",1);
  fChain->SetBranchStatus("cmet",1);
  fChain->SetBranchStatus("cmetphi",1);
  fChain->SetBranchStatus("emet",1);
  fChain->SetBranchStatus("emetphi",1);
  fChain->SetBranchStatus("mmet",1);
  fChain->SetBranchStatus("mmetphi",1);
  fChain->SetBranchStatus("pmet",1);
  fChain->SetBranchStatus("pmetphi",1); 
  fChain->SetBranchStatus("njets",1);
  fChain->SetBranchStatus("signaljetpt",1);
  fChain->SetBranchStatus("secondjetpt",1);
  
  char ROOT_FNAME[100];
  char TXT_FNAME[100];
  char FLAVOUR[10];                   // e.g. "ele", "mu"
  char LL_FLAVOUR[10];             // e.g. "ee", "mumu"
  char CONTROL_SAMPLE[10];   // e.g. "Z-->ee"

  Double_t LUMI;
  Double_t J1PT;
  Double_t J1ETA;
  Double_t J2PT;
  Double_t J2ETA;
  Double_t J1J2DPHI;
  Int_t LEP_PDG_ID;
  Double_t DILEPMASS_LOW;
  Double_t DILEPMASS_UP;
  Int_t NVTXS;                           // # of points for study of u_par and u_perp vs # of reconstructed vertices nvtx
  Int_t FIRST_NVTX;
  Double_t METNOLEP_START;
  //string MET_TYPE_NAME;
  string FILENAME_BASE;

  ifstream inputFile(configFileName);

  if (inputFile.is_open()) {

    Double_t value;
    string name;
    string parameterName;
    string parameterType;
    vector<Double_t> parameterValue;

    mySpaces(cout,2);
    cout << "Printing content of " << configFileName << " file" << endl;
    mySpaces(cout,1);

    while (inputFile >> parameterType ) {

      if (parameterType == "NUMBER") {

	inputFile >> parameterName >> value;
	parameterValue.push_back(value);
	cout << setw(20) << parameterName << setw(7) << value << endl;

      } else if (parameterType == "STRING") {
	 
	inputFile >> parameterName >> name;
	cout << right << setw(20) << parameterName << "  " << left << name << endl;
	if (parameterName == "FILENAME_BASE") {

	  FILENAME_BASE = name + "_Adish"; // to distinguish from files done with Emanuele's tree 
	  if ( !ISDATA_FLAG && unweighted_event_flag) FILENAME_BASE += "_weq1";  // if using unit weight, add _weq1 to filename (weq1 means weight = 1)

	}

	//if (parameterName == "MET_TYPE_NAME") MET_TYPE_NAME = name;

      }

    }

    LUMI = parameterValue[0];
    J1PT = parameterValue[1];
    J1ETA = parameterValue[2];
    J2PT = parameterValue[3];
    J2ETA = parameterValue[4];
    J1J2DPHI = parameterValue[5];
    LEP_PDG_ID = (Int_t) parameterValue[6];
    DILEPMASS_LOW = parameterValue[7];
    DILEPMASS_UP = parameterValue[8];
    NVTXS = (Int_t) parameterValue[9];
    FIRST_NVTX = (Int_t) parameterValue[10];
    METNOLEP_START = parameterValue[11];
    mySpaces(cout,2);

    inputFile.close();
                                                                                                                         
  } else {

    cout << "Error: could not open file " << configFileName << endl;
    exit(EXIT_FAILURE);

  }

  TVector2 metNoLepTV;
  //TLorentzVector Zreco;

  Double_t *metpt_ptr = NULL;
  Double_t *metphi_ptr = NULL;

  //if ( !(std::strcmp("pfmet",MET_TYPE_NAME.c_str())) ) {

    metpt_ptr = &pmet;
    metphi_ptr = &pmetphi;
    FILENAME_BASE += "_pmet";

    //}

    //FILENAME_BASE += "_allmet";

  Double_t metpt;
  Double_t metphi;
  Double_t metNoLepPt;

  //Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
  Double_t metBinEdges[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 750., 850., 1000.};
  Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;

  Int_t using_spring15_sample_flag = 0;
  if (FILENAME_BASE.find("spring15") != std::string::npos) {
    using_spring15_sample_flag = 1;    
    cout << "Using spring15 samples" << endl;
  }

  if ( !ISDATA_FLAG && unweighted_event_flag) cout << "Warning: no weight applied to events (w = 1)" << endl;  // if MC with unit weight, make user know

  Int_t using_zlljets_MCsample_flag = 0;
  if ( !ISDATA_FLAG && ( (FILENAME_BASE.find("zmumujets") != std::string::npos) || (FILENAME_BASE.find("zeejets") != std::string::npos) ) )  using_zlljets_MCsample_flag = 1; 

  if (ISDATA_FLAG) {
    strcpy(ROOT_FNAME,(FILENAME_BASE + "_DATA.root").c_str());
    strcpy(TXT_FNAME,(FILENAME_BASE + "_DATA.txt").c_str());
  } else {
    strcpy(ROOT_FNAME,(FILENAME_BASE + ".root").c_str());
    strcpy(TXT_FNAME,(FILENAME_BASE + ".txt").c_str());
  }

  if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...
     
    strcpy(FLAVOUR,"mu");
    strcpy(LL_FLAVOUR,"mumu");
    strcpy(CONTROL_SAMPLE,"Z-->mumu");

  } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

    strcpy(FLAVOUR,"ele");
    strcpy(LL_FLAVOUR,"ee");
    strcpy(CONTROL_SAMPLE,"Z-->ee");

  }

  TFile *rootFilePUrwt = NULL;     // the very next lines are for PU reweighting when using MC
  TH1F *HPUrwt = NULL;
  if (!ISDATA_FLAG) {

    cout << "Opening file purwtForAdish.root (with weights for PU reweighting)"<< endl;

    rootFilePUrwt = new TFile("purwtForAdish.root","READ");
    if (!rootFilePUrwt || !rootFilePUrwt->IsOpen()) {
      cout<<"Error: file \"purwtForAdish.root\" was not opened. End of programme."<<endl;
      exit(EXIT_FAILURE);
    }
    HPUrwt = (TH1F*)rootFilePUrwt->Get("puhist");
    if (!HPUrwt) {
      cout << "Error: could not get 'puhist' from 'purwtForAdish.root'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
  } 

  cout << "Opening file " <<ROOT_FNAME<< endl;

  TFile *rootFile = new TFile(ROOT_FNAME,"RECREATE");
  if (!rootFile || !rootFile->IsOpen()) {
    cout<<"Error: file \""<<ROOT_FNAME<<"\" was not opened."<<endl;
    exit(EXIT_FAILURE);
  }
  

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
  //TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
  TVirtualFitter::SetDefaultFitter("Minuit");

  //Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};       

  Float_t invMassBinWidth = 1.0;  // invariant mass histogram's bin width in GeV
  Int_t NinvMassBins = (DILEPMASS_UP - DILEPMASS_LOW) / invMassBinWidth;

  TH1D *HzlljetsYieldsMetBinGenLep = new TH1D("HzlljetsYieldsMetBinGenLep",Form("yields of %s control sample (%s gen if MC) in bins of met; #slash{E}_{T};# of events",CONTROL_SAMPLE,CONTROL_SAMPLE),nMetBins,metBinEdges);
  //TH1D *HzlljetsYieldsMetBinGenTau = new TH1D("HzlljetsYieldsMetBinGenTau",Form("yields of %s control sample (Z->#tau#tau gen) in bins of met; #slash{E}_{T};# of events",CONTROL_SAMPLE),nMetBins,metBinEdges);
   
  TH1D *HmetNoLepDistribution = new TH1D("HmetNoLepDistribution","",60,METNOLEP_START,METNOLEP_START+600);
  TH1D *HzptDistribution = new TH1D("HzptDistribution","",80,0,400);
  TH1D *Hjet1ptDistribution = new TH1D("Hjet1ptDistribution","",60,J1PT,J1PT+600);
  TH1D *Hjet2ptDistribution = new TH1D("Hjet2ptDistribution","",60,J2PT,J2PT+600);
  TH1D *HinvMass = new TH1D("HinvMass","",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);    // for MC it's done on Z->mumu or Z->ee at gen level
  TH1D *HvtxDistribution = new TH1D("HvtxDistribution","",40,-0.5,39.5);
  TH1D *HnjetsDistributions = new TH1D("HnjetsDistribution","njets using nJetClean30",10,-0.5,9.5);

  TH1D *HZtoLLRecoPt = new TH1D("HZtoLLRecoPt","",101,0.,1010);   // end at 1010 because I will put the overflow in the last bin
  // the previous histogram is differen from HzptDistribution because the binning is different

  TH1D *HzlljetsInvMassMetBinGenLep[nMetBins];
  TH1D *HZtoLLRecoPt_MetBin[nMetBins];

  for (Int_t i = 0; i < nMetBins; i++) {

    HzlljetsInvMassMetBinGenLep[i] = new TH1D(Form("HzlljetsInvMassMetBinGenLep_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);
    HZtoLLRecoPt_MetBin[i] = new TH1D(Form("HZtoLLRecoPt_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",101,0.,1010.);
    
  }

  // following 2 histograms hold the inclusive distribution of parallel and ortogonal component of recoil (projected along the ZpT direction)
   TH1D *H_uPerp_Distribution = new TH1D("H_uPerp_Distribution","",50,-200,200);
   TH1D *H_uParMinusZpT_Distribution = new TH1D("H_uParMinusZpT_Distribution","",50,-200,200);

   TH1D *H_uPerp_VS_Nvtx[NVTXS];
   TH1D *H_uParMinusZpT_VS_Nvtx[NVTXS]; 
   TH1D *H_uParMinusZpT_VS_Nvtx_lowZpT[NVTXS];
 
   for (Int_t i = 0; i < NVTXS; i++) {

     H_uPerp_VS_Nvtx[i] = new TH1D(Form("H_uPerp_VS_Nvtx_nvtx%i",FIRST_NVTX+i),"",80,-200,200);  // 5 GeV bins 
     H_uParMinusZpT_VS_Nvtx[i] = new TH1D(Form("H_uParMinusZpT_VS_Nvtx_nvtx%i",FIRST_NVTX+i),"ZpT in [250,500]GeV",80,-200,200);  // 5 GeV bins
     H_uParMinusZpT_VS_Nvtx_lowZpT[i] = new TH1D(Form("H_uParMinusZpT_VS_Nvtx_lowZpT_nvtx%i",FIRST_NVTX+i),"ZpT in [50,250]GeV",80,-200,200);  // 5 GeV bins

   }

   // zpt bin edges for respose studies (which are independent on the CS selection, they are done on the whole sample)
   //Double_t ZptBinEdges[] = {250., 260., 270., 280., 290., 310., 330., 350., 370., 390., 410., 430., 450., 470., 500., 530., 560, 600., 640., 700., 800.};
   //Double_t ZptBinEdges[] = {250., 260., 270., 280., 290., 310., 330., 350., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   //Double_t ZptBinEdges[] = {10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 310., 330., 350., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   //Double_t ZptBinEdges[] = {20., 40., 60., 80., 100., 120., 140., 160., 180., 200., 220., 240., 260., 280., 300., 320., 340., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   Double_t *ZptBinEdges = NULL;
   Double_t ZptBinEdgesMC[] = {1., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 220., 240., 260., 280., 300., 320., 340., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   Double_t ZptBinEdgesDATA[] = {1., 10., 20., 40., 60., 80., 100., 120., 140., 170., 200., 250.};
   Int_t nBinsForResponse = 0;   // # of bins for analysis as a function of ZpT

   if (ISDATA_FLAG) {

     ZptBinEdges = ZptBinEdgesDATA;
     nBinsForResponse = sizeof(ZptBinEdgesDATA)/sizeof(Double_t) - 1;  //number of bins is n-1 where n is the number of ZptBinEdges's elements

   } else {

     ZptBinEdges = ZptBinEdgesMC;
     nBinsForResponse = sizeof(ZptBinEdgesMC)/sizeof(Double_t) - 1;  //number of bins is n-1 where n is the number of ZptBinEdges's elements

   }

   TH1D *H_uPerp_VS_ZpT[nBinsForResponse];  
   TH1D *H_uParMinusZpT_VS_ZpT[nBinsForResponse];    // actually it will be (u_par-ZpT)
   TH1D *H_uPar_ZpT_ratio[nBinsForResponse];         // for the response curve: we will compute response in many ways. Here we use <uPar/ZpT>
   TH1D *H_uParMinusZpT_ZpT_ratio[nBinsForResponse]; // here we use <(uPar-ZpT)/ZpT> and will add back 1, so we get the same as above
   TH1D *HZptBinned[nBinsForResponse];

   //the following histograms will give the distribution of met|| / wzpt. The mean value will be used to create the response curve, that is (<met|| / wzpt>) vs wzpt
   // for each point, wzpt will be taken as the average wzpt in the range considered
 
   for (Int_t i = 0; i < nBinsForResponse; i++) {   

     //HZptBinned[i] are histograms with 5 bins in the range given by ZptBinEdges[i] and ZptBinEdges[i+1]
     // the mean wzpt in each bin will be computed as the histogram's mean
     HZptBinned[i] = new TH1D(Form("HZptBinned_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",5,ZptBinEdges[i],ZptBinEdges[i+1]); 
     // in the following histogram , range must include negative value: I saw that for low ZoT this distribution tends to be flat, thus if range goes from 0 to 2 (as it was before) the 
     //mean for ZpT tending to 0 will be 1 and not 0 as we would expect.
     H_uPar_ZpT_ratio[i] = new TH1D(Form("H_uPar_ZpT_ratio_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",350,-7.0,7.0); 
     H_uParMinusZpT_ZpT_ratio[i] = new TH1D(Form("H_uParMinusZpT_ZpT_ratio_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",350,-7.0,7.0);
     H_uPerp_VS_ZpT[i] = new TH1D(Form("H_uPerp_VS_ZpT_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",40,-200,200); 
     H_uParMinusZpT_VS_ZpT[i] = new TH1D(Form("H_uParMinusZpT_VS_ZpT_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",40,-200,200); 

   }
   
   // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
   TH1D *HmetBinEdges = new TH1D("HmetBinEdges","bin edges for met distributions",nMetBins+1,0.0,nMetBins+1);
   for (Int_t i = 0; i <= nMetBins; i++) {
     HmetBinEdges->SetBinContent(i+1,metBinEdges[i]);
   }

   TH1D *HZptBinEdges = new TH1D("HZptBinEdges","bin edges for ZpT distributions",nBinsForResponse+1,0.0,nBinsForResponse+1);
   for (Int_t i = 0; i <= nBinsForResponse; i++) {
     HZptBinEdges->SetBinContent(i+1,ZptBinEdges[i]);
   }

   TH1D *HnvtxBins = new TH1D("HnvtxBins","# of vertices for studies of variables as a function of nvtx",NVTXS,0.0,NVTXS);
   for (Int_t i = 0; i < NVTXS; i++) {              // watch out: differently from above, i < NVTXS, not <=, because if NVTXS = 3 I need 3 points, not 4
     HnvtxBins->SetBinContent(i+1,FIRST_NVTX+i);
   }
   
   // deciding  what is the event weight
   Double_t newwgt;

   if (ISDATA_FLAG || unweighted_event_flag) newwgt = 1.0;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"zlljets_resoResp::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;

   for (Int_t jentry=0; jentry<nentries; jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     if (jentry%500000 == 0) cout<<jentry<<endl;

     if(!ISDATA_FLAG && !unweighted_event_flag) {

       Double_t purwt = 0.0;
       if (nvtx <= 49) purwt = HPUrwt->GetBinContent(nvtx);     // purwt is the reweighting for PU
       newwgt = LUMI * wgt * xsec * purwt / wgtsum;
       // if (jentry%50000 == 0) {      // debugging
       // 	 cout<<newwgt<<endl;
       // 	 cout<<wgt<<endl;
       // 	 cout<<xsec<<endl;
       // 	 cout<<purwt<<endl;
       // 	 cout<<wgtsum<<endl;
       // }

     }

     TVector2 Zreco2D;
     TVector2 metPart;
     Zreco2D.SetMagPhi(zpt,zphi);

     //Zreco.SetPtEtaPhiM(zpt,zeta,zphi,zmass);
     //metNoLepTV.SetMagPhi(metpt,metphi);

     // metNoLepTV.SetMagPhi(0,0);
     // if (!ISDATA_FLAG) metPart.SetMagPhi(hmet/1.1,hmetphi);
     // else metPart.SetMagPhi(hmet,hmetphi);
     // metNoLepTV += metPart;
     // if (!ISDATA_FLAG) metPart.SetMagPhi(amet/1.6,ametphi);
     // else metPart.SetMagPhi(amet,ametphi);
     // metNoLepTV += metPart;
     // metPart.SetMagPhi(bmet,bmetphi); metNoLepTV += metPart;
     // metPart.SetMagPhi(cmet,cmetphi); metNoLepTV += metPart;
     // metPart.SetMagPhi(emet,emetphi); metNoLepTV += metPart;
     // metPart.SetMagPhi(mmet,mmetphi); metNoLepTV += metPart;
     // metPart.SetMagPhi(pmet,pmetphi); metNoLepTV += metPart;
     
     metpt = *metpt_ptr;
     metphi = *metphi_ptr;
     metNoLepTV.SetMagPhi(metpt,metphi);

     metNoLepTV += Zreco2D; // adding Z to met
     metNoLepPt = metNoLepTV.Mod();

     HzlljetsYieldsMetBinGenLep->Fill(metNoLepPt,newwgt);
     HZtoLLRecoPt->Fill(zpt,newwgt);	 

     HinvMass->Fill(zmass,newwgt);
     HmetNoLepDistribution->Fill(metNoLepPt,newwgt);
     HzptDistribution->Fill(zpt,newwgt);
     Hjet1ptDistribution->Fill(signaljetpt,newwgt);
     Hjet2ptDistribution->Fill(secondjetpt,newwgt);
     HvtxDistribution->Fill(nvtx,newwgt);
     HnjetsDistributions->Fill(njets,newwgt);    

     //TVector3 Zreco3D = Zreco.Vect();
     Double_t dphiMetNoLepZ = metNoLepTV.DeltaPhi(Zreco2D);

     Double_t u_par = metNoLepPt * TMath::Cos(dphiMetNoLepZ);  // actually u_par is minus this quantity, but then I do u_par-ZpT instead of u_par+ZpT
     Double_t u_perp = metNoLepPt * TMath::Sin(dphiMetNoLepZ);
     Double_t uparMinusZrecoPt = u_par - zpt;

     H_uPerp_Distribution->Fill(u_perp,newwgt);
     H_uParMinusZpT_Distribution->Fill(uparMinusZrecoPt,newwgt);

     if (zpt > ZptBinEdges[0]) {  

       Int_t nvtxBin = nvtx - FIRST_NVTX;
       Int_t lastnvtx = NVTXS + FIRST_NVTX;

       if ((nvtxBin >= 0) && (nvtx < lastnvtx)) {

	 H_uPerp_VS_Nvtx[nvtxBin]->Fill(u_perp,newwgt);
	     
	 if (zpt < 250 ) {

	   H_uParMinusZpT_VS_Nvtx_lowZpT[nvtxBin]->Fill(uparMinusZrecoPt,newwgt);
	     
	 } else if (zpt < 500) {                       // (met||-wzpt) distribution's width depends on Zpt, thus I use this range

	   H_uParMinusZpT_VS_Nvtx[nvtxBin]->Fill(uparMinusZrecoPt,newwgt);
	 
	 }       
    
       }  // end of   if ((nvtxBin >= 0) && (nvtx < lastnvtx))

       /**************************************************/
       // computing met responses
       /**************************************************/

       // first of all I make sure that wzpt is in the appropriate range
       if ( zpt < ZptBinEdges[nBinsForResponse] ) {

	 Int_t respBin = myGetBin(zpt,ZptBinEdges,nBinsForResponse);
	 //cout<<"bin = "<<bin<<endl;
	 HZptBinned[respBin]->Fill(zpt,newwgt);        
	 H_uPar_ZpT_ratio[respBin]->Fill(u_par/zpt,newwgt);          //the mean value of this histogram is the response
	 H_uParMinusZpT_ZpT_ratio[respBin]->Fill(uparMinusZrecoPt/zpt,newwgt);  //the mean value of this histogram +1 is the response
	 H_uPerp_VS_ZpT[respBin]->Fill(u_perp,newwgt);
	 H_uParMinusZpT_VS_ZpT[respBin]->Fill(uparMinusZrecoPt,newwgt);

       }

     }            // end of if (zpt > ZptBinEdges[0])

     // now entering analysis in bins of met

     if ((metNoLepPt > metBinEdges[0]) && (metNoLepPt < metBinEdges[nMetBins])) {

       Int_t bin = myGetBin(metNoLepPt,metBinEdges,nMetBins);
       
       HzlljetsInvMassMetBinGenLep[bin]->Fill(zmass,newwgt); 
       HZtoLLRecoPt_MetBin[bin]->Fill(zpt,newwgt);

     }



   }  // end of loop on entries 
	
   // if (!ISDATA_FLAG) {
   //   delete rootFilePUrwt;  // no need of them anymore and I don't want troubles with writing in the output file
   //   delete HPUrwt;
   // }

   // cout <<"check" << endl;

   /************************************/
   //                    MET|| & MET_|_ VS NVTX & ZpT
   /************************************/

   //resolution vs nvtx

   Double_t xValues[NVTXS];
   Double_t yValues[NVTXS];
   Double_t yValuesErr[NVTXS];

   for (Int_t i = 0; i < NVTXS; i++) {
     xValues[i] = i + FIRST_NVTX;
     yValues[i] = H_uParMinusZpT_VS_Nvtx_lowZpT[i]->GetRMS();
     yValuesErr[i] = H_uParMinusZpT_VS_Nvtx_lowZpT[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepParZvsNvtx_lowZpT = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMetNoLepParZvsNvtx_lowZpT->SetTitle(Form("resolution || from histogram's RMS, ZpT in [%2.0lf,250] GeV",ZptBinEdges[0]));
   GresolutionMetNoLepParZvsNvtx_lowZpT->Draw("AP");
   GresolutionMetNoLepParZvsNvtx_lowZpT->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsNvtx_lowZpT->GetXaxis()->SetTitle("nvtx");
   GresolutionMetNoLepParZvsNvtx_lowZpT->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsNvtx_lowZpT->GetYaxis()->SetTitleOffset(1.4); 
   GresolutionMetNoLepParZvsNvtx_lowZpT->SetName("gr_resolution_uPar_vs_Nvtx_lowZpT");
   GresolutionMetNoLepParZvsNvtx_lowZpT->Write();

   for (Int_t i = 0; i < NVTXS; i++) {
     yValues[i] = H_uParMinusZpT_VS_Nvtx[i]->GetRMS();
     yValuesErr[i] = H_uParMinusZpT_VS_Nvtx[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepParZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMetNoLepParZvsNvtx->SetTitle("resolution || from histogram's RMS, ZpT in [250,500] GeV");
   GresolutionMetNoLepParZvsNvtx->Draw("AP");
   GresolutionMetNoLepParZvsNvtx->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsNvtx->GetXaxis()->SetTitle("nvtx");
   GresolutionMetNoLepParZvsNvtx->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   GresolutionMetNoLepParZvsNvtx->SetName("gr_resolution_uPar_vs_Nvtx");
   GresolutionMetNoLepParZvsNvtx->Write();

   for (Int_t i = 0; i < NVTXS; i++) {
     yValues[i] = H_uPerp_VS_Nvtx[i]->GetRMS();
     yValuesErr[i] = H_uPerp_VS_Nvtx[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepOrtZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMetNoLepOrtZvsNvtx->SetTitle("resolution _|_ from histogram's RMS");
   GresolutionMetNoLepOrtZvsNvtx->Draw("AP");
   GresolutionMetNoLepOrtZvsNvtx->SetMarkerStyle(7);
   GresolutionMetNoLepOrtZvsNvtx->GetXaxis()->SetTitle("nvtx");
   GresolutionMetNoLepOrtZvsNvtx->GetYaxis()->SetTitle("#sigma (u_#perp ) [GeV]");
   GresolutionMetNoLepOrtZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   GresolutionMetNoLepOrtZvsNvtx->SetName("gr_resolution_uPerp_vs_Nvtx");
   GresolutionMetNoLepOrtZvsNvtx->Write();

   // resolution vs ZpT
   
   Double_t meanZpt[nBinsForResponse];
   Double_t meanZptErr[nBinsForResponse];
   Double_t resoMetNoLepParZvsZpt[nBinsForResponse];
   Double_t resoMetNoLepParZvsZptErr[nBinsForResponse];
   Double_t resoMetNoLepOrtZvsZpt[nBinsForResponse];
   Double_t resoMetNoLepOrtZvsZptErr[nBinsForResponse];

   for (Int_t i = 0; i < nBinsForResponse; i++) {

     meanZpt[i] = HZptBinned[i]->GetMean();
     meanZptErr[i] = HZptBinned[i]->GetMeanError();
     resoMetNoLepParZvsZpt[i] = H_uParMinusZpT_VS_ZpT[i]->GetRMS();
     resoMetNoLepParZvsZptErr[i] = H_uParMinusZpT_VS_ZpT[i]->GetRMSError();
     resoMetNoLepOrtZvsZpt[i] = H_uPerp_VS_ZpT[i]->GetRMS();
     resoMetNoLepOrtZvsZptErr[i] = H_uPerp_VS_ZpT[i]->GetRMSError();

   }

   TGraphErrors *GresolutionMetNoLepParZvsZpt = new TGraphErrors(nBinsForResponse,meanZpt,resoMetNoLepParZvsZpt,meanZptErr,resoMetNoLepParZvsZptErr);
   GresolutionMetNoLepParZvsZpt->SetTitle("resolution || from histogram's RMS");
   GresolutionMetNoLepParZvsZpt->Draw("AP");
   GresolutionMetNoLepParZvsZpt->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsZpt->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMetNoLepParZvsZpt->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsZpt->GetYaxis()->SetTitleOffset(1.2); 
   GresolutionMetNoLepParZvsZpt->SetName("gr_resolution_uPar_vs_ZpT");
   GresolutionMetNoLepParZvsZpt->Write();

   TGraphErrors *GresolutionMetNoLepOrtZvsZpt = new TGraphErrors(nBinsForResponse,meanZpt,resoMetNoLepOrtZvsZpt,meanZptErr,resoMetNoLepOrtZvsZptErr);
   GresolutionMetNoLepOrtZvsZpt->SetTitle("resolution _|_ from histogram's RMS");
   GresolutionMetNoLepOrtZvsZpt->Draw("AP");
   GresolutionMetNoLepOrtZvsZpt->SetMarkerStyle(7);
   GresolutionMetNoLepOrtZvsZpt->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMetNoLepOrtZvsZpt->GetYaxis()->SetTitle("#sigma (u_#perp ) [GeV]");
   GresolutionMetNoLepOrtZvsZpt->GetYaxis()->SetTitleOffset(1.2); 
   GresolutionMetNoLepOrtZvsZpt->SetName("gr_resolution_uPerp_vs_ZpT");
   GresolutionMetNoLepOrtZvsZpt->Write();

   // response curve: it can be defined in many ways

   Double_t response[nBinsForResponse];              // <uPar/ZpT>  this is somewhat biased
   Double_t responseErr[nBinsForResponse];

   Double_t response_gausFit[nBinsForResponse];      // <Upar-ZpT>/<ZpT>  + 1: numerator from a gaussian fit (we use gaussian mean since we look for the peak)   
   Double_t responseErr_gausFit[nBinsForResponse];
   TFitResultPtr ptrGausFit;
   Double_t meanUparMinusZpt_gausFit[nBinsForResponse];       // will contain <u_par - ZpT> taken from a gaussian fit
   Double_t meanUparMinusZptErr_gausFit[nBinsForResponse];    // will contain uncertainty on <u_par - ZpT> taken from a gaussian fit

   Double_t response_gausFit_bis[nBinsForResponse];    // <Upar-ZpT/ZpT>  + 1: mean value from a gaussian fit (we use gaussian mean since we look for the peak)   
   Double_t responseErr_gausFit_bis[nBinsForResponse];
   TFitResultPtr ptrGausFit_bis;
   Double_t mean_UparMinusZpt_ZpT_ratio_gausFit_bis[nBinsForResponse];       // will contain <u_par - ZpT> taken from a gaussian fit
   Double_t mean_UparMinusZpt_ZpT_ratioErr_gausFit_bis[nBinsForResponse];    
  
   // using <A/B> is different from using <A>/<B>. In general <f(x,y)> != f(<x>,<y>). The two relation coincide if the variance of x and y can be neglected or
   // the second partial derivatives of f wrt x and y are small enough (note that here is also uPar = uPar(ZpT) so that we actually have 1 independent variable)

   for (Int_t i = 0; i < nBinsForResponse; i++) {    

     response[i] = H_uPar_ZpT_ratio[i]->GetMean();
     responseErr[i] = H_uPar_ZpT_ratio[i]->GetMeanError();
     //cout<<i<<" meanZpt = "<<meanZpt[i]<<" +/- "<<meanZptErr[i]<<"    response = "<<response[i]<<" +/- "<<responseErr[i]<<endl;

     // using "(<u_par -ZpT> / <ZpT>) + 1" to compute response (adding 1 because that ratio is centered around 0). <.> is the mean value
     // do the fit between -3.5 RMS and + 3.5 RMS (found to be good enough) with a gaussian to get mean
     // option WL use loglikelihood fit with weighted histograms (better if there are empty bins)
     // Q is quiet mode (minimum printing on stdout). V prints everything. Default option is between Q and V
     // S is necessary to pass object and access to fit parameter

     Double_t tmpRMS;                       // temporary variable with distribution's RMS
     Int_t fitStatus = 0;  // when 0 the fit was successful                                            
     if (H_uParMinusZpT_VS_ZpT[i]->GetEntries() <= 10 ) {                // if empty histogram (including underflows and overflows) no fit is done
       meanUparMinusZpt_gausFit[i] = H_uParMinusZpT_VS_ZpT[i]->GetMean();                               // actually the fit has no sense with few points
       meanUparMinusZptErr_gausFit[i] = H_uParMinusZpT_VS_ZpT[i]->GetMeanError();
     } else {
       tmpRMS = H_uParMinusZpT_VS_ZpT[i]->GetRMS();
       ptrGausFit = H_uParMinusZpT_VS_ZpT[i]->Fit("gaus","Q S","",-3.5*tmpRMS,3.5*tmpRMS);
       fitStatus = ptrGausFit;
       if (fitStatus == 0) {
	 meanUparMinusZpt_gausFit[i] = ptrGausFit->Parameter(1);      // 1 is the mean (0 and 2 are normalization and sigma of gaussian)
	 meanUparMinusZptErr_gausFit[i] = ptrGausFit->ParError(1);
       } else {
	 meanUparMinusZpt_gausFit[i] = H_uParMinusZpT_VS_ZpT[i]->GetMean();            // in case of some error in the fit, use TH1::GetMean()
	 meanUparMinusZptErr_gausFit[i] = H_uParMinusZpT_VS_ZpT[i]->GetMeanError();
       }
     }
     response_gausFit[i] = (meanUparMinusZpt_gausFit[i] / meanZpt[i]);    // in a second moment, adding 1 to response, otherwise it would be centered around 0
     responseErr_gausFit[i] = response_gausFit[i] * sqrt( (meanUparMinusZptErr_gausFit[i]*meanUparMinusZptErr_gausFit[i] / (meanUparMinusZpt_gausFit[i]*meanUparMinusZpt_gausFit[i])) + (meanZptErr[i]*meanZptErr[i] / (meanZpt[i]*meanZpt[i])));     // for the uncertainty, using response BEFORE adding 1 
     response_gausFit[i] += 1.;
     // uncertainty computed assuming uncorrelated numerator and denominator

     // now using "(<u_par -ZpT / ZpT>) + 1" to compute response (adding 1 because that ratio is centered around 0)

     if ( H_uParMinusZpT_ZpT_ratio[i]->GetEntries() <= 10 ) {                 // if empty histogram (including underflows and overflows) no fit is done
       mean_UparMinusZpt_ZpT_ratio_gausFit_bis[i] = H_uParMinusZpT_ZpT_ratio[i]->GetMean();                     // actually the fit has no sense with few points
       mean_UparMinusZpt_ZpT_ratioErr_gausFit_bis[i] = H_uParMinusZpT_ZpT_ratio[i]->GetMeanError();
     } else {
       tmpRMS = H_uParMinusZpT_ZpT_ratio[i]->GetRMS();
       ptrGausFit_bis = H_uParMinusZpT_ZpT_ratio[i]->Fit("gaus","Q S","",-3.5*tmpRMS,3.5*tmpRMS);
       fitStatus = ptrGausFit_bis;
       if (fitStatus == 0) {
	 mean_UparMinusZpt_ZpT_ratio_gausFit_bis[i] =  ptrGausFit_bis->Parameter(1);  // 1 is the mean (0 and 2 are normalization and sigma of gaussian)
	 mean_UparMinusZpt_ZpT_ratioErr_gausFit_bis[i] =  ptrGausFit_bis->ParError(1);
       } else {
	 mean_UparMinusZpt_ZpT_ratio_gausFit_bis[i] = H_uParMinusZpT_ZpT_ratio[i]->GetMean();        // in case of some error in the fit, use TH1::GetMean()
	 mean_UparMinusZpt_ZpT_ratioErr_gausFit_bis[i] = H_uParMinusZpT_ZpT_ratio[i]->GetMeanError();
       }
     }     
     response_gausFit_bis[i] = mean_UparMinusZpt_ZpT_ratio_gausFit_bis[i] + 1.;  
     responseErr_gausFit_bis[i] = mean_UparMinusZpt_ZpT_ratioErr_gausFit_bis[i]; // for the uncertainty, using response BEFORE adding 1 

   }

   TGraphErrors *GresponseCurve = new TGraphErrors(nBinsForResponse,meanZpt,response,meanZptErr,responseErr);
   GresponseCurve->SetTitle("response curve");
   GresponseCurve->Draw("AP");
   GresponseCurve->SetMarkerStyle(7);    // 7 is a medium dot
   GresponseCurve->GetXaxis()->SetTitle("ZpT [GeV]");
   GresponseCurve->GetYaxis()->SetTitle(" < u_{||} / ZpT >");
   GresponseCurve->GetYaxis()->SetRangeUser(0.0, 1.1);
   GresponseCurve->GetYaxis()->SetTitleOffset(1.4); 
   GresponseCurve->SetName("gr_responseCurve");
   GresponseCurve->Write();

   TGraphErrors *GresponseCurve_gausFit = new TGraphErrors(nBinsForResponse,meanZpt,response_gausFit,meanZptErr,responseErr_gausFit);
   GresponseCurve_gausFit->SetTitle("response curve as <uPar-ZpT>/<ZpT> +1; numerator from mean of gaussian fit (mode of distribution)");
   GresponseCurve_gausFit->Draw("AP");
   GresponseCurve_gausFit->SetMarkerStyle(7);    // 7 is a medium dot
   GresponseCurve_gausFit->GetXaxis()->SetTitle("ZpT [GeV]");
   GresponseCurve_gausFit->GetYaxis()->SetTitle(" < u_{||} - ZpT > / < ZpT > + 1");
   GresponseCurve_gausFit->GetYaxis()->SetRangeUser(0.0, 1.1);
   GresponseCurve_gausFit->GetYaxis()->SetTitleOffset(1.4); 
   GresponseCurve_gausFit->SetName("gr_responseCurve_gausFit");
   GresponseCurve_gausFit->Write();

   TGraphErrors *GresponseCurve_gausFit_bis = new TGraphErrors(nBinsForResponse,meanZpt,response_gausFit_bis,meanZptErr,responseErr_gausFit_bis);
   GresponseCurve_gausFit_bis->SetTitle("response curve as <uPar-ZpT/ZpT> +1, mean from mean gaussian fit (mode of the distribution)");
   GresponseCurve_gausFit_bis->Draw("AP");
   GresponseCurve_gausFit_bis->SetMarkerStyle(7);    // 7 is a medium dot
   GresponseCurve_gausFit_bis->GetXaxis()->SetTitle("ZpT [GeV]");
   GresponseCurve_gausFit_bis->GetYaxis()->SetTitle(" < u_{||} - ZpT / ZpT > + 1");
   GresponseCurve_gausFit_bis->GetYaxis()->SetRangeUser(0.0, 1.1);
   GresponseCurve_gausFit_bis->GetYaxis()->SetTitleOffset(1.4); 
   GresponseCurve_gausFit_bis->SetName("gr_responseCurve_gausFit_bis");
   GresponseCurve_gausFit_bis->Write();

   // correcting resolution of uPar for the response
  
   Int_t nPoints = GresolutionMetNoLepParZvsZpt->GetN();
   TH1D* Hyresponse = new TH1D("Hyresponse","",nPoints,0,nPoints);
   TH1D* Hyresopar = new TH1D("Hyresopar","",nPoints,0,nPoints);
   Double_t xresopar;  // just to use TGraph::GetPoint()
   Double_t yresopar[nPoints];  
   Double_t yresoparErr[nPoints];
   Double_t yresponse[nPoints];

   for (Int_t i = 0; i < nPoints; i++) {

     GresponseCurve_gausFit_bis->GetPoint(i,xresopar,yresponse[i]);
     Hyresponse->SetBinContent(i+1, *(yresponse + i));
     Hyresponse->SetBinError(i+1, GresponseCurve_gausFit_bis->GetErrorY(i));
     GresolutionMetNoLepParZvsZpt->GetPoint(i,xresopar,yresopar[i]);
     Hyresopar->SetBinContent(i+1, *(yresopar + i));
     Hyresopar->SetBinError(i+1, GresolutionMetNoLepParZvsZpt->GetErrorY(i));

   }

   if ( !Hyresopar->Divide(Hyresponse) ) cout << " Error in Hyresopar->Divide(Hyresponse) " << endl; 
   else {

     for (Int_t i = 0; i < nPoints; i++) {

       *(yresopar + i) = Hyresopar->GetBinContent(i+1);
       *(yresoparErr + i) = Hyresopar->GetBinError(i+1);

     }

   }

   TGraphErrors *GresolutionMetNoLepParZvsZpt_Corrected = new TGraphErrors(nBinsForResponse,meanZpt,yresopar,meanZptErr,yresoparErr);
   GresolutionMetNoLepParZvsZpt_Corrected->SetTitle("resolution || from histogram's RMS, corrected for response");
   GresolutionMetNoLepParZvsZpt_Corrected->Draw("AP");
   GresolutionMetNoLepParZvsZpt_Corrected->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsZpt_Corrected->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMetNoLepParZvsZpt_Corrected->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsZpt_Corrected->GetYaxis()->SetTitleOffset(1.2); 
   GresolutionMetNoLepParZvsZpt_Corrected->SetName("gr_resolution_uPar_vs_ZpT_corrected");
   GresolutionMetNoLepParZvsZpt_Corrected->Write();

   delete Hyresponse;
   delete Hyresopar;

   // end of TGraphs

   mySpaces(cout,2);
   myPrintYieldsMetBinInStream(cout, HzlljetsYieldsMetBinGenLep, metBinEdges, nMetBins);

   cout<<"creating file '"<<TXT_FNAME<<"' ..."<<endl;
   ofstream myfile(TXT_FNAME,ios::out);

   if ( !myfile.is_open() ) {

     cout<<"Error: unable to open file "<<TXT_FNAME<<" !"<<endl;
     exit(EXIT_FAILURE);
     
   }

   //opening inputFile named configFileName again to save content in myfile named TXT_FNAME

   inputFile.open(configFileName);

   if (inputFile.is_open()) {
     
     mySpaces(myfile,2);
     cout << "Saving content of " << configFileName << " file in "<< TXT_FNAME << endl;
     myfile << "Content of " << configFileName << endl;
     mySpaces(myfile,1);

     Double_t value;
     string name;
     string parameterName;
     string parameterType;

     while (inputFile >> parameterType ) {

       if (parameterType == "NUMBER") {

	 inputFile >> parameterName >> value;
	 myfile << setw(20) << parameterName << setw(7) << value << endl;

       } else if (parameterType == "STRING") {
	 
	 inputFile >> parameterName >> name;
	 myfile << right << setw(20) << parameterName << "  " << left << name << endl;

       }

     }
     
     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << " to save content in "<< TXT_FNAME << endl;
     exit(EXIT_FAILURE);

   }

   mySpaces(myfile,2);
   if (!ISDATA_FLAG && unweighted_event_flag) myfile << "======   Using unweighted events (w = 1)   ======" << endl;
   mySpaces(myfile,3);
   myPrintYieldsMetBinInStream(myfile, HzlljetsYieldsMetBinGenLep, metBinEdges, nMetBins);

   myfile.close();

   // I add overflow bin's content in the last bin for all histograms where that is needed
   myAddOverflowInLastBin(HZtoLLRecoPt);

   for (Int_t i = 0; i < nMetBins; i++) {
     myAddOverflowInLastBin(HZtoLLRecoPt_MetBin[i]);
   }

   rootFile->Write();

   rootFile->Close();
   delete rootFile;


}
