#include "TMath.h"
#include <TAxis.h>
#include <TH1D.h>
#include <TF1.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip> //for input/output manipulators
#include <fstream>
#include <vector>
#include <string>

#include "myClasses.h" 
#include "functionsForAnalysis.h"

using namespace std;

void myDashes(ostream &myOutStream = cout) {

    myOutStream<<"-----------------------------------------------------------------------------------"<<endl; 

}


void mySpaces(ostream &myOutStream, Int_t nSpace) {

  for (Int_t i = 0; i < nSpace; i++) {
    myOutStream<<endl; 
  }

}

Double_t myGetPhi(const Double_t x_component, const Double_t y_component) {

  Double_t phi = 0;
  //some useful variables
   Double_t myPi = TMath::Pi(); //just not to call TMath::Pi() every time
   Double_t twoPi = 2*TMath::Pi(); //just not to calculate 2*TMath::Pi() every time

     //while evaluating phi =atan(y/x), we require x!=0. If x=0, then 
     //phi= +/- Pi()/2 with the same sign as the y 
     //for x!=0, since atan(y/x) yields a value in [-Pi()/2, +Pi()/2]
     //we must sum or subtract Pi() depending on the dial

     if (x_component!=0) {

       if ( (x_component<0) && (y_component<0) ) {
         phi = atan(y_component/x_component)-myPi;
       } else if ( (x_component<0) && (y_component>0) ) {
         phi = atan(y_component/x_component)+myPi;
       } else {
         phi = atan(y_component/x_component);
       }

     } else {

      if ( (y_component<0) ) {
	 phi = -myPi/2.;
       } else { 
	 phi = myPi/2.;
      }
     }
     return phi;

}

void myFillOverflowBin(TH1D *histo, const Double_t value) {
  
  //this function put the overflows in the last bin of an histogram

  TAxis *axis = 0;
  axis = histo->GetXaxis();
  if(value>(axis->GetBinUpEdge(axis->GetLast()))) {
    histo->Fill(axis->GetBinCenter(axis->GetLast()));
  } else {
    histo->Fill(value);
  }

}

void myFillOverflowBinW(TH1D *histo, const Double_t value, const Double_t weight) {

  // filling last bin with overflow for histograms with weights

  TAxis *axis = 0;
  axis = histo->GetXaxis();
  if(value > (axis->GetBinUpEdge(axis->GetLast()))) {
    histo->Fill(axis->GetBinCenter(axis->GetLast()),weight);
  } else {
    histo->Fill(value, weight);
  }

}

void myDrawOverflow(const TH1D *oldHisto, const char* myOptions) {

  // this function create a new histogram with one more bin than oldHisto. This additional bin will contain the overflow of oldHisto. 
  // The width of this additional bin is set equal to the width of the last bin of oldHisto. This is not relevant for histograms where the
  // width is always the same for each bin, but care should be taken for histograms whose bin's width depends on the bit
  // By definition, the overflow bin of newHisto is empty
  // The strategy is creating a temporary array htemp which, before being deleted (since it's created with "new" it' better to delete it at the end
  // of the function), will be copied into newHisto

  TH1::SetDefaultSumw2();            // all the following histograms will automatically call TH1::Sumw2() 
  TH1::StatOverflows();                 // enable use of underflows and overflows for statistics computation 

  TAxis *axis = 0;
  axis = oldHisto->GetXaxis();
  Double_t overflowBinWidth = axis->GetBinWidth(axis->GetLast());
  Double_t newUpperEdge = axis->GetBinUpEdge(axis->GetLast()) + overflowBinWidth;
  Double_t newLowerEdge = axis->GetBinLowEdge(axis->GetFirst());
  Int_t newNBins = oldHisto->GetNbinsX() + 1;   //  +1 because I want to add an additional bin for the overflow

  // some checks to keep things under control 
  // cout<<"inside myOverflowInLastBin()"<<endl;
  // cout<<"newNBins = "<<newNBins<<endl;
  // cout<<"newUpperEdge = "<<newUpperEdge<<endl;
  // cout<<"oldHisto->GetBinContent(newNBins) = "<<oldHisto->GetBinContent(newNBins)<<endl;

  TH1D *htemp = new TH1D("htemp","",newNBins,newLowerEdge,newUpperEdge);

  // setting bin content and error
  for (Int_t i = 1; i <= newNBins; i++) {
    htemp->SetBinContent(i,oldHisto->GetBinContent(i));
    htemp->SetBinError(i,oldHisto->GetBinError(i));
  }
  // filling underflow bin 
  htemp->SetBinContent(0,oldHisto->GetBinContent(0));
  htemp->SetBinError(0,oldHisto->GetBinError(0));
  // setting line color and axis titles
  htemp->GetXaxis()->SetTitle(oldHisto->GetXaxis()->GetTitle());
  htemp->GetYaxis()->SetTitle(oldHisto->GetYaxis()->GetTitle());
  htemp->SetLineColor(oldHisto->GetLineColor());
  htemp->SetFillColor(oldHisto->GetFillColor());
  // setting number of entries to effective entries (same as entries for unweighted histograms) + content 
  // of overflow bin (which is not included by Get(Effective)Entries())
  htemp->SetEntries(oldHisto->GetBinContent(newNBins) + oldHisto->GetEffectiveEntries());

  htemp->Draw(myOptions);

}

void myDrawOverflow(const TH1D *oldHisto, const char* myOptions, const Int_t mySetStats = 1) {

  // this function create a new histogram with one more bin than oldHisto. This additional bin will contain the overflow of oldHisto. 
  // The width of this additional bin is set equal to the width of the last bin of oldHisto. This is not relevant for histograms where the
  // width is always the same for each bin, but care should be taken for histograms whose bin's width depends on the bit
  // By definition, the overflow bin of newHisto is empty
  // The strategy is creating a temporary array htemp which, before being deleted (since it's created with "new" it' better to delete it at the end
  // of the function), will be copied into newHisto

  TH1::SetDefaultSumw2();            // all the following histograms will automatically call TH1::Sumw2() 
  TH1::StatOverflows();                 // enable use of underflows and overflows for statistics computation 

  TAxis *axis = 0;
  axis = oldHisto->GetXaxis();
  Double_t overflowBinWidth = axis->GetBinWidth(axis->GetLast());
  Double_t newUpperEdge = axis->GetBinUpEdge(axis->GetLast()) + overflowBinWidth;
  Double_t newLowerEdge = axis->GetBinLowEdge(axis->GetFirst());
  Int_t newNBins = oldHisto->GetNbinsX() + 1;   //  +1 because I want to add an additional bin for the overflow

  // some checks to keep things under control 
  // cout<<"inside myOverflowInLastBin()"<<endl;
  // cout<<"newNBins = "<<newNBins<<endl;
  // cout<<"newUpperEdge = "<<newUpperEdge<<endl;
  // cout<<"oldHisto->GetBinContent(newNBins) = "<<oldHisto->GetBinContent(newNBins)<<endl;

  TH1D *htemp = new TH1D("htemp","",newNBins,newLowerEdge,newUpperEdge);
  // TH1D histo("histo","",newNBins,newLowerEdge,newUpperEdge);
  // TH1D *htemp = &histo;

  // setting bin content and error
  for (Int_t i = 1; i <= newNBins; i++) {
    htemp->SetBinContent(i,oldHisto->GetBinContent(i));
    htemp->SetBinError(i,oldHisto->GetBinError(i));
  }
  // filling underflow bin 
  htemp->SetBinContent(0,oldHisto->GetBinContent(0));
  htemp->SetBinError(0,oldHisto->GetBinError(0));
  // setting line color and axis titles
  htemp->GetXaxis()->SetTitle(oldHisto->GetXaxis()->GetTitle());
  htemp->GetYaxis()->SetTitle(oldHisto->GetYaxis()->GetTitle());
  htemp->SetLineColor(oldHisto->GetLineColor());
  htemp->SetFillColor(oldHisto->GetFillColor());
  // setting number of entries to effective entries (same as entries for unweighted histograms) + content 
  // of overflow bin (which is not included by Get(Effective)Entries())
  htemp->SetEntries(oldHisto->GetBinContent(newNBins) + oldHisto->GetEffectiveEntries());

  if( !mySetStats ) htemp->SetStats(kFALSE);
  htemp->Draw(myOptions);

}

TH1D *myOverflowInLastBin(const TH1D *oldHisto) {

  // this function create a new histogram with one more bin than oldHisto. This additional bin will contain the overflow of oldHisto. 
  // The width of this additional bin is set equal to the width of the last bin of oldHisto. This is not relevant for histograms where the
  // width is always the same for each bin, but care should be taken for histograms whose bin's width depends on the bit
  // By definition, the overflow bin of newHisto is empty
  // The strategy is creating a temporary array htemp which, before being deleted (since it's created with "new" it' better to delete it at the end
  // of the function), will be copied into newHisto

  TH1::SetDefaultSumw2();            // all the following histograms will automatically call TH1::Sumw2() 
  TH1::StatOverflows();                 // enable use of underflows and overflows for statistics computation 

  TAxis *axis = 0;
  axis = oldHisto->GetXaxis();
  Double_t overflowBinWidth = axis->GetBinWidth(axis->GetLast());
  Double_t newUpperEdge = axis->GetBinUpEdge(axis->GetLast()) + overflowBinWidth;
  Double_t newLowerEdge = axis->GetBinLowEdge(axis->GetFirst());
  Int_t newNBins = oldHisto->GetNbinsX() + 1;   //  +1 because I want to add an additional bin for the overflow

  // some checks to keep things under control 
  // cout<<"inside myOverflowInLastBin()"<<endl;
  // cout<<"newNBins = "<<newNBins<<endl;
  // cout<<"newUpperEdge = "<<newUpperEdge<<endl;
  // cout<<"oldHisto->GetBinContent(newNBins) = "<<oldHisto->GetBinContent(newNBins)<<endl;

  static TH1D *htemp1 = new TH1D("htemp1","",newNBins,newLowerEdge,newUpperEdge);

  // setting bin content and error
  for (Int_t i = 1; i <= newNBins; i++) {
    htemp1->SetBinContent(i,oldHisto->GetBinContent(i));
    htemp1->SetBinError(i,oldHisto->GetBinError(i));
  }
  // filling underflow bin 
  htemp1->SetBinContent(0,oldHisto->GetBinContent(0));
  htemp1->SetBinError(0,oldHisto->GetBinError(0));
  // setting line color and axis titles
  htemp1->GetXaxis()->SetTitle(oldHisto->GetXaxis()->GetTitle());
  htemp1->GetYaxis()->SetTitle(oldHisto->GetYaxis()->GetTitle());
  htemp1->SetLineColor(oldHisto->GetLineColor());
  htemp1->SetFillColor(oldHisto->GetFillColor());
  // setting number of entries to effective entries (same as entries for unweighted histograms) + content 
  // of overflow bin (which is not included by Get(Effective)Entries())
  htemp1->SetEntries(oldHisto->GetBinContent(newNBins) + oldHisto->GetEffectiveEntries());

  return htemp1;

}

void myOverflowInLastBin2(TH1D *newHisto, const TH1D *oldHisto) {

  // this function, when called with a histo defined as "TH1D *histo = 0" as newHisto, produces segmentation fault!!!

  // this function create a new histogram with one more bin than oldHisto. This additional bin will contain the overflow of oldHisto. 
  // The width of this additional bin is set equal to the width of the last bin of oldHisto. This is not relevant for histograms where the
  // width is always the same for each bin, but care should be taken for histograms whose bin's width depends on the bit
  // By definition, the overflow bin of newHisto is empty
  // The strategy is creating a temporary array htemp which, before being deleted (since it's created with "new" it' better to delete it at the end
  // of the function), will be copied into newHisto

  TH1::SetDefaultSumw2();            // all the following histograms will automatically call TH1::Sumw2() 
  TH1::StatOverflows();                 // enable use of underflows and overflows for statistics computation 

  TAxis *axis = 0;
  axis = oldHisto->GetXaxis();
  Double_t overflowBinWidth = axis->GetBinWidth(axis->GetLast());
  Double_t newUpperEdge = axis->GetBinUpEdge(axis->GetLast()) + overflowBinWidth;
  Double_t newLowerEdge = axis->GetBinLowEdge(axis->GetFirst());
  Int_t newNBins = oldHisto->GetNbinsX() + 1;   //  +1 because I want to add an additional bin for the overflow

  // some checks to keep things under control 
  // cout<<"inside myOverflowInLastBin()"<<endl;
  // cout<<"newNBins = "<<newNBins<<endl;
  // cout<<"newUpperEdge = "<<newUpperEdge<<endl;
  // cout<<"oldHisto->GetBinContent(newNBins) = "<<oldHisto->GetBinContent(newNBins)<<endl;

  TH1D htemp2("htemp2","",newNBins,newLowerEdge,newUpperEdge);

  // setting bin content and error
  for (Int_t i = 1; i <= newNBins; i++) {
    htemp2.SetBinContent(i,oldHisto->GetBinContent(i));
    htemp2.SetBinError(i,oldHisto->GetBinError(i));
  }
  // filling underflow bin 
  htemp2.SetBinContent(0,oldHisto->GetBinContent(0));
  htemp2.SetBinError(0,oldHisto->GetBinError(0));
  // setting line color and axis titles
  htemp2.GetXaxis()->SetTitle(oldHisto->GetXaxis()->GetTitle());
  htemp2.GetYaxis()->SetTitle(oldHisto->GetYaxis()->GetTitle());
  htemp2.SetLineColor(oldHisto->GetLineColor());
  htemp2.SetFillColor(oldHisto->GetFillColor());
  // setting number of entries to effective entries (same as entries for unweighted histograms) + content 
  // of overflow bin (which is not included by Get(Effective)Entries())
  htemp2.SetEntries(oldHisto->GetBinContent(newNBins) + oldHisto->GetEffectiveEntries());

  *newHisto = htemp2;

}

void myAddOverflowInLastBin(TH1D *h) {

  // to avoid problems regarding memory leak for not deleting htemp in previous function, I sum directly the content of overflow bin in last bin

  Int_t lastBinNumber = h->GetNbinsX();
  Int_t overflowBinNumber = 1 + lastBinNumber;
  Double_t lastBinContent = h->GetBinContent(lastBinNumber);
  Double_t overflowBinContent = h->GetBinContent(overflowBinNumber);
  Double_t lastBinError = h->GetBinError(lastBinNumber);
  Double_t overflowBinError = h->GetBinError(overflowBinNumber);

  // add content of overflow bin in last bin and set error as square root of sum of error squares (with the assumption that they are uncorrelated)
  h->SetBinContent(lastBinNumber, lastBinContent + overflowBinContent);
  h->SetBinError(lastBinNumber, sqrt(lastBinError * lastBinError + overflowBinError * overflowBinError));

}

void myAddUnderflowInFirstBin(TH1D *h) {

  // to avoid problems regarding memory leak for not deleting htemp in previous function, I sum directly the content of overflow bin in last bin

  // Int_t underflowBinNumber = 0;
  Double_t firstBinContent = h->GetBinContent(1);
  Double_t underflowBinContent = h->GetBinContent(0);
  Double_t firstBinError = h->GetBinError(1);
  Double_t underflowBinError = h->GetBinError(0);

  // add content of underflow bin in first bin and set error as square root of sum of error squares (with the assumption that they are uncorrelated)
  h->SetBinContent(1, firstBinContent + underflowBinContent);
  h->SetBinError(1, sqrt(firstBinError * firstBinError + underflowBinError * underflowBinError));

}

//calculation of rejection factors
Double_t myRejectionFactor(Int_t selected, Int_t total) {
 
  return (Double_t) total/selected;

}

Double_t myCrystalBall(double* x, double* par) {

  Double_t xcur = x[0];
  Double_t alpha = par[0];
  Double_t n = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t invAbsAlpha = 1./absAlpha;

  if ( t >= -absAlpha)  {
    return N*exp(-0.5*t*t);
  } else {
    Double_t A = TMath::Power(n*invAbsAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    Double_t B = n*invAbsAlpha - absAlpha;
    return N*A/TMath::Power(B-t,n);
  }

}

Double_t my2sideCrystalBall(double* x, double* par) {

  //a priori we allow for different shape of right and left tail, thus two values of alpha and n 

  Double_t xcur = x[0];
  Double_t alphaL = par[0];
  Double_t nL = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t alphaR = par[5];
  Double_t nR = par[6];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlphaL = fabs((Double_t)alphaL);
  Double_t invAbsAlphaL = 1./absAlphaL;
  Double_t absAlphaR = fabs((Double_t)alphaR);
  Double_t invAbsAlphaR = 1./absAlphaR;

  
  if ( t<-absAlphaL ) {
    //cout<<"checkpoint dscb left"<<endl;
    Double_t AL = TMath::Power(nL*invAbsAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    Double_t BL = nL*invAbsAlphaL - absAlphaL;
    return N*AL*TMath::Power(BL-t,-nL);
  } else if ( t <= absAlphaR )  {
    //cout<<"checkpoint dscb gaussian"<<endl;
    return N*exp(-0.5*t*t);
  } else {
    //cout<<"checkpoint dscb right"<<endl;
    Double_t AR = TMath::Power(nR*invAbsAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    Double_t BR = nR*invAbsAlphaR - absAlphaR;
    return N*AR*TMath::Power(BR+t,-nR);
  }

}

Double_t myResolutionFunction(double* x, double* par) {

  // resolution function dE/E is given by the sum in quadrature of 3 terms: a / sqrt(E) ; b / E ; c
  // thus, sigma(E) is obtained by multiplying dE / E by E  -->  sigma(E) = a*sqrt(E) + b + c*E where + is the sum in quadrature
  // a, b and c are the stochastic, noise and constant term respectively

  Double_t E = x[0];      // energy, pT or whatever
  Double_t a = par[0];
  Double_t b = par[1];
  Double_t c = par[2];

  return sqrt( a*a * E + b*b + c*c * E*E );

}

Double_t myResolutionFunctionNoB(double* x, double* par) {

  // resolution function dE/E is given by the sum in quadrature of 3 terms: a / sqrt(E) ; b / E ; c
  // thus, sigma(E) is obtained by multiplying dE / E by E  -->  sigma(E) = a*sqrt(E) + b + c*E where + is the sum in quadrature
  // a, b and c are the stochastic, noise and constant term respectively

  // the term with b could be neglected at big values of E

  Double_t E = x[0];      // energy, pT or whatever
  Double_t a = par[0];
  Double_t c = par[1];

  return sqrt( a*a * E + c*c * E*E );

}

void myAddDefaultPackages(FILE *fp, const char* filename) {

  if (fp == NULL) {
    cout<<"Error loading default packages in file \""<<filename<<"\".\nThese will not be added to the file."<<endl;
  } else {
    fprintf(fp,"\\documentclass[11pt,a4paper]{article}\n");
    fprintf(fp,"\\usepackage[T1]{fontenc}\n");
    fprintf(fp,"\\usepackage[utf8]{inputenc}\n");
    fprintf(fp,"\\usepackage[italian,english]{babel}\n");
    fprintf(fp,"\\usepackage{graphicx}\n");
    fprintf(fp,"\\usepackage{wrapfig}\n");
    fprintf(fp,"\\usepackage{sidecap}\n");
    fprintf(fp,"\\usepackage{booktabs}\n");
    fprintf(fp,"\\usepackage{amsmath}\n");
    fprintf(fp,"\\usepackage{amsfonts}\n");
    fprintf(fp,"\\usepackage{amssymb}\n");
    fprintf(fp,"\\usepackage{multirow}\n");
    fprintf(fp,"\\usepackage{subfig}\n");
    //fprintf(fp,"\\usepackage{siunitx}\n");
    fprintf(fp,"\n");
  }
  
}

void makeTableTex(FILE *fp, const Double_t lumi, const Double_t nTotalWeightedEvents, const mask *m, const string addString) {

  // begin table
  fprintf(fp,"\\begin{table}\n");
  fprintf(fp,"\\caption{\\emph{Entries after each selection step, normalised to %4.2lf $fb^{-1}$; n is the number of entries after the i-th "
	  "selection step; %s}}\n",lumi,addString.c_str());
  fprintf(fp,"\\[\n");
  fprintf(fp,"\\begin{array}{rcrrr}\n");
  fprintf(fp,"\\toprule\n");
  fprintf(fp,"\\ \\textbf{step}& \\textbf{definition}  & n  &  n_i/n_0 & n_i/n_{i-1}  \\\\ \n");
  fprintf(fp,"\\midrule\n");
  fprintf(fp," %i  & \\text{all entries} & %6.0lf & %1.4lf & %1.4lf \\\\\n",0, nTotalWeightedEvents, 1.0, 1.0);
  fprintf(fp,"\\midrule\n");
  for (Int_t i = 0; i < m->getMaskSize(); i++) {
    if(m->nEvents[i] < 100 ) {
      if (i == 0) {
	fprintf(fp," %i  & & %6.2lf & %1.4lf & %1.4lf \\\\\n",i+1, m->nEvents[i], m->nEvents[i]/nTotalWeightedEvents,m->nEvents[i]/nTotalWeightedEvents);
      } else {
	if(m->nEvents[i-1] == 0) {
	  fprintf(fp," %i  & & %6.2lf & %1.4lf & %1.4lf \\\\\n",i+1,m->nEvents[i],m->nEvents[i]/nTotalWeightedEvents,1.0);
	} else {
	  fprintf(fp," %i  & & %6.2lf & %1.4lf & %1.4lf \\\\\n",i+1,m->nEvents[i],m->nEvents[i]/nTotalWeightedEvents,m->nEvents[i]/m->nEvents[i-1]);
	}
      }
    } else {
      if (i == 0) {
	fprintf(fp," %i  & & %6.0lf & %1.4lf & %1.4lf \\\\\n",i+1, m->nEvents[i], m->nEvents[i]/nTotalWeightedEvents,m->nEvents[i]/nTotalWeightedEvents);
      } else {
	if(m->nEvents[i-1] == 0) {
	  fprintf(fp," %i  & & %6.0lf & %1.4lf & %1.4lf \\\\\n",i+1,m->nEvents[i],m->nEvents[i]/nTotalWeightedEvents,1.0);
	} else {
	  fprintf(fp," %i  & & %6.0lf & %1.4lf & %1.4lf \\\\\n",i+1,m->nEvents[i],m->nEvents[i]/nTotalWeightedEvents,m->nEvents[i]/m->nEvents[i-1]);
	}
      }
    }
    for (Int_t j = 0; j < selection::getNSelections(); j++) {
      if ((m->singleMask[i] >> j) & 1) {
	string str = selection::listOfSelections[j]->getDefinition();
	fprintf(fp,"& \\text{%s} & & & \\\\\n",str.c_str());
      }
    }
    fprintf(fp,"\\midrule\n");
  }
  fprintf(fp,"\\bottomrule\n");
  fprintf(fp,"\\end{array}\n");
  fprintf(fp,"\\]\n");
  fprintf(fp,"\\end{table}\n");
  //end of table
  fprintf(fp,"\n");

}

void makeTableTex(FILE *fp, const Double_t lumi, const Double_t nTotalWeightedEvents, const mask *m) {

  // begin table
  fprintf(fp,"\\begin{table}\n");
  fprintf(fp,"\\caption{\\emph{Entries after each selection step, normalised to %4.2lf $fb^{-1}$; n is the number of entries after the i-th selection step.}}\n",lumi);
  fprintf(fp,"\\[\n");
  fprintf(fp,"\\begin{array}{rcrrr}\n");
  fprintf(fp,"\\toprule\n");
  fprintf(fp,"\\ \\textbf{step}& \\textbf{definition}  & n  & $n\\_i/n\\_0$ & $n\\_i/n\\_{i-1}$ \\\\ \n");
  fprintf(fp,"\\midrule\n");
  fprintf(fp," %i  & \\text{all entries} & %6.0lf & %1.4lf & %1.4lf \\\\\n",0, nTotalWeightedEvents, 1.0, 1.0);
  fprintf(fp,"\\midrule\n");
  for (Int_t i = 0; i < m->getMaskSize(); i++) {
    if(m->nEvents[i] < 100 ) {
      if (i == 0) {
	fprintf(fp," %i  & & %6.2lf & %1.4lf & %1.4lf \\\\\n",i+1, m->nEvents[i], m->nEvents[i]/nTotalWeightedEvents,m->nEvents[i]/nTotalWeightedEvents);
      } else {
	if(m->nEvents[i-1] == 0) {
	  fprintf(fp," %i  & & %6.2lf & %1.4lf & %1.4lf \\\\\n",i+1,m->nEvents[i],m->nEvents[i]/nTotalWeightedEvents,1.0);
	} else {
	  fprintf(fp," %i  & & %6.2lf & %1.4lf & %1.4lf \\\\\n",i+1,m->nEvents[i],m->nEvents[i]/nTotalWeightedEvents,m->nEvents[i]/m->nEvents[i-1]);
	}
      }
    } else {
      if (i == 0) {
	fprintf(fp," %i  & & %6.0lf & %1.4lf & %1.4lf \\\\\n",i+1, m->nEvents[i], m->nEvents[i]/nTotalWeightedEvents,m->nEvents[i]/nTotalWeightedEvents);
      } else {
	if(m->nEvents[i-1] == 0) {
	  fprintf(fp," %i  & & %6.0lf & %1.4lf & %1.4lf \\\\\n",i+1,m->nEvents[i],m->nEvents[i]/nTotalWeightedEvents,1.0);
	} else {
	  fprintf(fp," %i  & & %6.0lf & %1.4lf & %1.4lf \\\\\n",i+1,m->nEvents[i],m->nEvents[i]/nTotalWeightedEvents,m->nEvents[i]/m->nEvents[i-1]);
	}
      }
    }
    for (Int_t j = 0; j < selection::getNSelections(); j++) {
      if ((m->singleMask[i] >> j) & 1) {
	string str = selection::listOfSelections[j]->getDefinition();
	fprintf(fp,"& \\text{%s} & & & \\\\\n",str.c_str());
      }
    }
    fprintf(fp,"\\midrule\n");
  }
  fprintf(fp,"\\bottomrule\n");
  fprintf(fp,"\\end{array}\n");
  fprintf(fp,"\\]\n");
  fprintf(fp,"\\end{table}\n");
  //end of table
  fprintf(fp,"\n");

}


void makeTableTexNoEff(FILE *fp, const Double_t lumi, const Double_t nTotalWeightedEvents, const mask *m, const string addString) {

  // begin table
  fprintf(fp,"\\begin{table}\n");
  fprintf(fp,"\\caption{\\emph{Entries after each selection step, normalised to %4.2lf $fb^{-1}$; n is the number of entries after the i-th "
	  "selection step; %s}}\n",lumi,addString.c_str());
  fprintf(fp,"\\[\n");
  fprintf(fp,"\\begin{array}{rcr}\n");
  fprintf(fp,"\\toprule\n");
  fprintf(fp,"\\ \\textbf{step}& \\textbf{definition}  & n   \\\\ \n");
  fprintf(fp,"\\midrule\n");
  fprintf(fp," %i  & \\text{all entries} & %6.0lf \\\\\n",0, nTotalWeightedEvents);
  fprintf(fp,"\\midrule\n");
  for (Int_t i = 0; i < m->getMaskSize(); i++) {
    if(m->nEvents[i] < 100 ) {
      fprintf(fp," %i  & & %6.2lf \\\\\n",i+1, m->nEvents[i]);          
    } else {  
      fprintf(fp," %i  & & %6.0lf  \\\\\n",i+1, m->nEvents[i]);
    }
    for (Int_t j = 0; j < selection::getNSelections(); j++) {
      if ((m->singleMask[i] >> j) & 1) {
	string str = selection::listOfSelections[j]->getDefinition();
	fprintf(fp,"& \\text{%s} & \\\\\n",str.c_str());
      }
    }
    fprintf(fp,"\\midrule\n");
  }
  fprintf(fp,"\\bottomrule\n");
  fprintf(fp,"\\end{array}\n");
  fprintf(fp,"\\]\n");
  fprintf(fp,"\\end{table}\n");
  //end of table
  fprintf(fp,"\n");

}

void makeTableTexNoEff(FILE *fp, const Double_t lumi, const Double_t nTotalWeightedEvents, const mask *m) {

  // begin table
  fprintf(fp,"\\begin{table}\n");
  fprintf(fp,"\\caption{\\emph{Entries after each selection step, normalised to %4.2lf $fb^{-1}$; n is the number of entries after the i-th selection step.}}\n",lumi);
  fprintf(fp,"\\[\n");
  fprintf(fp,"\\begin{array}{rcr}\n");
  fprintf(fp,"\\toprule\n");
  fprintf(fp,"\\ \\textbf{step}& \\textbf{definition}  & n   \\\\ \n");
  fprintf(fp,"\\midrule\n");
  fprintf(fp," %i  & \\text{all entries} & %6.0lf  \\\\\n",0, nTotalWeightedEvents);
  fprintf(fp,"\\midrule\n");
  for (Int_t i = 0; i < m->getMaskSize(); i++) {
    if(m->nEvents[i] < 100 ) {     
      fprintf(fp," %i  & & %6.2lf \\\\\n",i+1, m->nEvents[i]);      
    } else {
      	fprintf(fp," %i  & & %6.0lf  \\\\\n",i+1, m->nEvents[i]);	
      }
    for (Int_t j = 0; j < selection::getNSelections(); j++) {
      if ((m->singleMask[i] >> j) & 1) {
	string str = selection::listOfSelections[j]->getDefinition();
	fprintf(fp,"& \\text{%s} &  \\\\\n",str.c_str());
      }
    }
    fprintf(fp,"\\midrule\n");
  }
  fprintf(fp,"\\bottomrule\n");
  fprintf(fp,"\\end{array}\n");
  fprintf(fp,"\\]\n");
  fprintf(fp,"\\end{table}\n");
  //end of table
  fprintf(fp,"\n");

}


char myAskToSaveFile(const char* filename) {

  char answer = '\0';

  do {
     cout<<"Save event yields on file \""<<filename<<"\"?(y/n) : ";
     cin>>answer;
       if (cin.fail()) {
	 cout<<"Error: cin failed. Data will not be written on \""<<filename<<"\""<<endl;
	 answer = 'n';
       }
       cout<<"answer is \""<<answer<<"\""<<endl;
   } while ( (answer != 'y') && (answer != 'n') );

  return answer;

}

Int_t myGetBin(Double_t x, Double_t *inputArray, Int_t nBins) {
  
// given an array containing a set of edges defining some bins, this function returns the bin which the value x belongs to. 
// e.g.:  array = {1.2, 3.4, 4.5} has 3 elements and defines 2 bins, so that nBins = 2. If x = 3.7, the function returns 1

  if ( (x < *(inputArray)) || (x > *(inputArray + nBins)) ) {
    cout<<"Warning in function myGetBin(): bin not found! End of programme!"<<endl;
    exit(EXIT_FAILURE);
  }

  Int_t correctBin = 0;
  for (Int_t bin = 0; bin < nBins; bin++) {
    if (x > *(inputArray+bin) )  correctBin = bin;  
  } 

  return correctBin;
 
}


void myPrintEventYields(ostream & myOutStream, const Double_t lumi, const Int_t cutSteps, const Double_t * eventsInStep) {

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*        EVENT YIELDS        *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  //nwentries is the number of events taking all weights into account (cuts due to triggers were not taken into account)
  myOutStream<<"nwentries:  weighted total number of entries = "<<eventsInStep[0]<<endl;
  myOutStream<<"n:               number of events after i-th cut"<<endl;
  myOutStream<<"aR:            absolute ratio = n(i)/n(1)"<<endl;
  myOutStream<<"rR:             relative ratio = n(i)/n(i-1)"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"data normalised to "<<lumi<<" fb^-1"<<endl;
  myOutStream<<setw(3)<<"cut"<<setw(12)<<"n"<<setw(12)<<"aR"<<setw(8)<<"rR"<<endl;
  for (Int_t i = 1; i <= cutSteps; i++) {
    if (i == 1) {	  
      myOutStream<<setw(3)<<i<<setw(12)<<eventsInStep[i]<<fixed<<setprecision(4)<<setw(12)<<1.0<<setw(8)<<1.0<<endl;     
    } else {
      myOutStream<<setw(3)<<i<<setw(12)<<eventsInStep[i]<<fixed<<setprecision(4)<<setw(12)<<eventsInStep[i]/eventsInStep[1]<<
	setw(8)<<eventsInStep[i]/eventsInStep[i-1]<<endl;  
    }
  }
  myOutStream<<"--------------------------------------------------------"<<endl;
  myOutStream<<endl;

}


void myPrintEventYields(ostream & myOutStream, const Double_t lumi, const vector<Double_t> eventsInStep) {

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*        EVENT YIELDS        *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  //nwentries is the number of events taking all weights into account (cuts due to triggers were not taken into account)
  myOutStream<<"nwentries:  weighted total number of entries = "<<eventsInStep[0]<<endl;
  myOutStream<<"n:               number of events after i-th cut"<<endl;
  myOutStream<<"aR:            absolute ratio = n(i)/nwentries"<<endl;
  myOutStream<<"rR:             relative ratio = n(i)/n(i-1)"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"data normalised to "<<lumi<<" fb^-1"<<endl;
  myOutStream<<setw(5)<<"cut"<<setw(12)<<"n"<<setw(12)<<"aR"<<setw(8)<<"rR"<<endl;    
  for (Int_t i = 1; i < eventsInStep.size(); i++) {
    myOutStream<<setw(5)<<i<<setw(12)<<eventsInStep[i]<<fixed<<setprecision(4)<<setw(12)<<eventsInStep[i]/eventsInStep[0]<<
      setw(8)<<eventsInStep[i]/eventsInStep[i-1]<<endl;  
  }
  myOutStream<<"--------------------------------------------------------"<<endl;
  myOutStream<<endl;

}

Int_t myPartGenAlgo(const Int_t nGenParticles, const Int_t* particleId, const Int_t* particleMotherId, const Int_t part_pdgId, const Int_t mother_pdgId ) {

  // This function looks whether there are two leptons with given PDG_ID generated from a decaying particle. It works when using 
  //Emanuele's tree.
  // The algorithm looks for two opposite sign leptons originating from mother in the list of generated particles.
  // The function returns 1 if the search is successful, and 0 otherwise

  // first it looks for first leptons; if it's found, flag gets 1, thus only the "else if" will be evaluated until the second is found

  Int_t found = 0;  // will become 1 if search is successful, that is to say, two OS particles coming from Z are found
  Int_t flag = 0;   // will become 1 when first tau from Z is found
  Int_t flavour = 0;  // when first tau is found, it's set to its pdgID
  Int_t i = 0;      // index of the array

  while (!found && i < nGenParticles) {
    if ( (flag == 0) && (fabs(particleId[i]) == part_pdgId) && (fabs(particleMotherId[i]) == mother_pdgId) ) {
      flag = 1;
      flavour = particleId[i];
    } else if ( (flag == 1) && ((particleId[i] + flavour) == 0) && (fabs(particleMotherId[i]) == mother_pdgId) ) {
      found = 1;
    }
    i++;
  } 
 
  return found; 

}

Int_t myPartGenAlgo(const Int_t nGenParticles, const Int_t* particleId, const Int_t* particleMotherId, const Int_t part_pdgId, const Int_t mother_pdgId, Int_t &firstIndex, Int_t &secondIndex) {

  // This function looks whether there are two leptons with given PDG_ID generated from a decaying particle and takes their indices in 
  //the array. It works when using Emanuele's tree.
  // The algorithm looks for two opposite sign leptons originating from mother in the list of generated particles.
  // The function returns 1 if the search is successful, and 0 otherwise

  // first it looks for first leptons; if it's found, flag gets 1, thus only the else if will be evaluated until the secon is found

  Int_t found = 0;  // will become 1 if search is successful, that is to say, two OS particles coming from Z are found
  Int_t flag = 0;   // will become 1 when first tau from Z is found
  Int_t flavour = 0;  // when first tau is found, it's set to its pdgID
  Int_t i = 0;      // index of the array

  //initialize indices with non valid elements in the list
  firstIndex = nGenParticles;
  secondIndex = nGenParticles;

  while (!found && (i < nGenParticles)) {
    if ( !flag && (fabs(particleId[i]) == part_pdgId) && (fabs(particleMotherId[i]) == mother_pdgId) ) {
      flag = 1;
      flavour = particleId[i];
      firstIndex = i;
    } else if ( (flag == 1) && ((particleId[i] + flavour) == 0) && (fabs(particleMotherId[i]) == mother_pdgId) ) {
      found = 1;
      secondIndex = i;
    }
    i++;
  } 
 
  if (!found) {
    // in case one of the two is not found, assign first 2 elements in the list (but function will return 0)
    firstIndex = 0;
    secondIndex = 1;
  }

  return found; 

}

Int_t myPartGenAlgo(const Int_t nGenParticles, const Int_t* particleId, const Int_t* particleMotherId, const Int_t part_pdgId, const Int_t mother_pdgId, Int_t &firstIndex, Int_t &secondIndex, Int_t &motherIndex, const Int_t* particleMotherIndex) {

  // This function looks whether there are two leptons with given PDG_ID generated from a decaying particle and takes their indices in 
  //the array. It works when using Emanuele's tree.
  // It also looks for the index of the mother in the list of particles (for instance a Z)
  // The algorithm looks for two opposite sign leptons originating from mother in the list of generated particles.
  // The function returns 1 if the search is successful, and 0 otherwise

  // first it looks for first leptons; if it's found, flag gets 1, thus only the else if will be evaluated until the secon is found

  Int_t found = 0;  // will become 1 if search is successful, that is to say, two OS particles coming from Z are found
  Int_t flag = 0;   // will become 1 when first tau from Z is found
  Int_t flavour = 0;  // when first tau is found, it's set to its pdgID
  Int_t i = 0;      // index of the array

  //initialize indices with non valid elements in the list
  firstIndex = nGenParticles;
  secondIndex = nGenParticles;

  while (!found && (i < nGenParticles)) {
    if ( !flag && (fabs(particleId[i]) == part_pdgId) && (fabs(particleMotherId[i]) == mother_pdgId) ) {
      flag = 1;
      flavour = particleId[i];
      firstIndex = i;
    } else if ( (flag == 1) && ((particleId[i] + flavour) == 0) && (fabs(particleMotherId[i]) == mother_pdgId) ) {
      found = 1;
      secondIndex = i;
      motherIndex = particleMotherIndex[firstIndex]; 
    }
    i++;
  } 
 
  if (!found) {
    // in case one of the two is not found, assign first 3 elements in the list (but function will return 0)
    firstIndex = 0;
    secondIndex = 1;
    motherIndex = 2;
  }

  return found; 

}

Int_t myGetPairIndexInArray(const Int_t id, const Int_t nLep, const Int_t *lep_pdgId, Int_t &firstIndex, Int_t &secondIndex) {

  // this function aims at finding the index of two SF/OS leptons with |pdg_id| = id in a list of particles. In Emanuele's tree the array is sorted by decreasing pt
  // and usual algorythm takes the two leading ones and apply selections such as OS/SF. Suppose we are looking for electrons: with this function we allow for 
  // events with e+, mu+, e- or even mu+, e+, e-  or  mu+, mu+, e+, e+, e-  to be analyzed (and discarded later in the selection flow, if that is the case).

  // In case there are more than two potentially good leptons (say, 3 electrons) we take the two leading ones as the pair. For example, if we have e+, mu-, e+, e-,
  // then the first lepton, e+, and the last one, e-, are selected as the pair. Of course in principle the pair from a decaying Z might be made of the last two, so we 
  // are making the assumption that the interesting pair is made by the hardest possible leptons: it might be possible to measure the efficiency of such a choice for 
  // a MC sample where we can get information about the generated particles

  // In case one of the two indices is not found (i.e. at least one is equal to nLep), they are set to 0 and 1 to ensure that they are valid indices: consequently, they will 
  // be discarded in the selection asking for a pair of OS/SF leptons of a given flavour

  // However, the function returns 1 if particles are found, and 0 otherwise

  // as a start we assign non valid values to the indices 
  firstIndex = nLep;
  secondIndex = nLep;

  Int_t found = 0;
  Int_t first_pdgId = 0;
  Int_t i = 0;
  Int_t absId = fabs(id);

  if (nLep < 2) return found;

  do {

    if ( fabs(lep_pdgId[i]) == absId ) {
      firstIndex = i;
      first_pdgId = lep_pdgId[i];
    }
    i++;

  } while ( (firstIndex == nLep) && (i < nLep) );
  // when first index is found, we look for the second, starting from firstIndex in the list. It might happen that there's no valid secondIndex. 
  //Note than now i = firstIndex + 1 and it might happen that  i = nLep which would not be valid: thus we now use a "while() ... do" loop instead of a "do ... while".
  while ((i < nLep) && (secondIndex == nLep) ) {

    if ( (lep_pdgId[i] + first_pdgId) == 0 ) {
      secondIndex = i;
      found = 1;
    }

    i++;

  }

  if ( !found ) {
    //cout<<"Warning in myGetPairIndexInArray(): couldn't find two valid indices. 0 and 1 will be used."<<endl;
    firstIndex = 0;
    secondIndex = 1;
  }  
  
  return found;

}


Int_t myGetPartIndex(const Int_t id, const Int_t nPart, const Int_t* partId) {

  //this function return the index of a given particle in the array of particle_pdgId. The first occurrence of the index is returned as the position in the array

  Int_t i = 0;
  Int_t index = nPart;

  do {
    if ( fabs(partId[i]) == id ) index = i;
    i++;
  } while ( (index == nPart) && (i < nPart) );
  if (index < nPart) return index;
  else {
    cout<<" Warning in myGetPartIndex: index not found. End of Programme"<<endl;
    exit(EXIT_FAILURE);
  }

}

void myPrintYieldsMetBin(const TH1D* histo, const Double_t *metBinEdges, const Int_t nMetBins, const char* YieldsFileName) {

  ofstream myfile(YieldsFileName,ios::out);

  if ( !myfile.is_open() ) {
    cout<<"Error: unable to open file "<<YieldsFileName<<" !"<<endl;  
  } else {

   cout<<"Saving yields in bins of met in file \""<<YieldsFileName<<"\" ..."<<endl;
   myfile<<"# "<<histo->GetName()<<" : yields in met bin"<<endl;
   myfile<<"#bin   yield   uncertainty"<<endl;
   myfile<<"<"<<(Int_t)metBinEdges[0]<<" ";
   myfile<<histo->GetBinContent(0)<<" "<<histo->GetBinError(0)<<endl;
   for (Int_t i = 0; i < nMetBins; i++ ) {  //from underflow to overflow bin 
     myfile<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
     myfile<<histo->GetBinContent(i+1)<<" "<<histo->GetBinError(i+1)<<endl;
   }
   myfile<<">"<<(Int_t)metBinEdges[nMetBins]<<" ";
   myfile<<histo->GetBinContent(nMetBins + 1)<<" "<<histo->GetBinError(nMetBins + 1)<<endl;
   myfile.close();
  }

}

void myPrintYieldsMetBin(const TH1D* histo, const Double_t *metBinEdges, const Int_t nMetBins) {

  
   cout<<histo->GetName()<<" : yields in met bin"<<endl;
   cout<<"#bin   yield   uncertainty"<<endl;
   cout<<"<"<<(Int_t)metBinEdges[0]<<" ";
   cout<<histo->GetBinContent(0)<<" "<<histo->GetBinError(0)<<endl;
   for (Int_t i = 0; i < nMetBins; i++ ) {  //from underflow to overflow bin 
     cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
     cout<<histo->GetBinContent(i+1)<<" "<<histo->GetBinError(i+1)<<endl;
   }
   cout<<">"<<(Int_t)metBinEdges[nMetBins]<<" ";
   cout<<histo->GetBinContent(nMetBins + 1)<<" "<<histo->GetBinError(nMetBins + 1)<<endl;


}


void myPrintYieldsMetBinInStream(ostream & myOutStream, const TH1D* histo, const Double_t *metBinEdges, const Int_t nMetBins) {
 
   myOutStream<<histo->GetName()<<" : yields in met bin"<<endl;
   myOutStream<<"#bin   yield   MC uncertainty     sqrt(yield)"<<endl;
   myOutStream<<"<"<<(Int_t)metBinEdges[0]<<" ";
   myOutStream<<histo->GetBinContent(0)<<" "<<histo->GetBinError(0)<<" "<<sqrt(histo->GetBinContent(0))<<endl;
   for (Int_t i = 0; i < nMetBins; i++ ) {  //from underflow to overflow bin 
     myOutStream<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
     myOutStream<<histo->GetBinContent(i+1)<<" "<<histo->GetBinError(i+1)<<" "<<sqrt(histo->GetBinContent(i+1))<<endl;
   }
   myOutStream<<">"<<(Int_t)metBinEdges[nMetBins]<<" ";
   myOutStream<<histo->GetBinContent(nMetBins + 1)<<" "<<histo->GetBinError(nMetBins + 1)<<" "<<sqrt(histo->GetBinContent(nMetBins + 1))<<endl;

}


void myPrintYieldsMetBinInStream(ostream & myOutStream, const TH1D* histo, const std::vector<Double_t> &metBinEdges) {
 
  Int_t nMetBins = metBinEdges.size() - 1; // -1 because metBinEdges has the lower edges of the histogram bins, so that metBinEdges = {0,10,20,30} is to be interpreted as 3 bins from 0 to 30 (values > 30 fill the overflow bin)

  myOutStream<<histo->GetName()<<" : yields in met bin"<<endl;
  myOutStream<<"#bin   yield   MC uncertainty     sqrt(yield)"<<endl;
  myOutStream<<"<"<<(Int_t)metBinEdges[0]<<" ";
  myOutStream<<histo->GetBinContent(0)<<" "<<histo->GetBinError(0)<<" "<<sqrt(histo->GetBinContent(0))<<endl;
  for (Int_t i = 0; i < nMetBins; i++ ) {  //from underflow to overflow bin 
    myOutStream<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
    myOutStream<<histo->GetBinContent(i+1)<<" "<<histo->GetBinError(i+1)<<" "<<sqrt(histo->GetBinContent(i+1))<<endl;
  }
  myOutStream<<">"<<(Int_t)metBinEdges[nMetBins]<<" ";
  myOutStream<<histo->GetBinContent(nMetBins + 1)<<" "<<histo->GetBinError(nMetBins + 1)<<" "<<sqrt(histo->GetBinContent(nMetBins + 1))<<endl;

}


void mySumWeight_filler_spring15_25ns(const std::string suffix,  std::vector<Double_t> & sumWeightVector) {

  // fill vector with sumWeight for samples used in spring15 at 25 ns. The problem was that, for a given sample,  each HT bin had a different value of sumWeight. Now, since I had to
  // merge these bins to run the analysis, I had to keep track of all these values and use the proper value when the i-th bin is being analyzed in the code 
  // trees with skim on metNoMu > 200

  if (suffix.find("DYJetsToLL") != std::string::npos) {
     sumWeightVector.push_back(2625679.0);
     sumWeightVector.push_back(955972.0);
     sumWeightVector.push_back(1048047.0);
     sumWeightVector.push_back(987977.0);
   }

   if (suffix.find("ZJetsToNuNu") != std::string::npos) {
     sumWeightVector.push_back(3937143.0);
     sumWeightVector.push_back(5032927.0);
     sumWeightVector.push_back(1014139.0);
     sumWeightVector.push_back(1015904.0);
   }

   if (suffix.find("WJetsToLNu") != std::string::npos) {  // binning 100-200-400-600-Inf
     sumWeightVector.push_back(10142187.0);
     sumWeightVector.push_back(5231856.0);
     sumWeightVector.push_back(1901705.0);
     sumWeightVector.push_back(1033899.0);
   }

   if (suffix.find("GJets") != std::string::npos) {  //HT200to400 bin missing
     sumWeightVector.push_back(4875882.0);
     sumWeightVector.push_back(2419265.0);
     sumWeightVector.push_back(2544631.0);
   }

   if (suffix.find("Top") != std::string::npos) { // all samples except TTJets_LO
     sumWeightVector.push_back(1000000.0);
     sumWeightVector.push_back(86689822712.6);
     sumWeightVector.push_back(18928638.0);
     sumWeightVector.push_back(3318695.54646);
     sumWeightVector.push_back(1336800136.23);
     sumWeightVector.push_back(995600.0);
   }

   if (suffix.find("QCD") != std::string::npos) {  // HT500to1000 bin missing
     sumWeightVector.push_back(18717349.0);
     sumWeightVector.push_back(20086103.0);
     sumWeightVector.push_back(4963895.0);
   }

   if (suffix.find("Diboson") != std::string::npos) {  // WW, WZ, ZZ
     sumWeightVector.push_back(1930000.0);
     sumWeightVector.push_back(991232.0);
     sumWeightVector.push_back(996168.0);
   }

}
 

void myEventsInSubsamples_filler_spring15_25ns(const std::string suffix,  std::vector<Int_t> & eventsInSubsamples) {

  // fill vector with number of events for samples used in spring15 at 25 ns. The problem was that, for a given sample,  each HT bin had a different value of sumWeight. Now, since I had to
  // merge these bins to run the analysis, I had to keep track of all these values and use the proper value when the i-th bin is being analyzed in the code 
  // trees with skim on metNoMu > 200

  if (suffix.find("DYJetsToLL") != std::string::npos) {
     eventsInSubsamples.push_back(5534);
     eventsInSubsamples.push_back(44507);
     eventsInSubsamples.push_back(98826);
     eventsInSubsamples.push_back(127087);   
   }

   if (suffix.find("ZJetsToNuNu") != std::string::npos) {
     eventsInSubsamples.push_back(24169);
     eventsInSubsamples.push_back(601862);
     eventsInSubsamples.push_back(225671);
     eventsInSubsamples.push_back(296792);    
   }

   if (suffix.find("WJetsToLNu") != std::string::npos) {  // binning 100-200-400-600-Inf
     eventsInSubsamples.push_back(28150);
     eventsInSubsamples.push_back(317221);
     eventsInSubsamples.push_back(249117);
     eventsInSubsamples.push_back(189073);   
   }

   if (suffix.find("GJets") != std::string::npos) {  //HT200to400 bin missing
     eventsInSubsamples.push_back(46);
     eventsInSubsamples.push_back(1756);
     eventsInSubsamples.push_back(7723);
   }

   if (suffix.find("Top") != std::string::npos) { // all samples except TTJets_LO
     eventsInSubsamples.push_back(21274);
     eventsInSubsamples.push_back(948245);
     eventsInSubsamples.push_back(414595);
     eventsInSubsamples.push_back(13074);   
     eventsInSubsamples.push_back(227709);
     eventsInSubsamples.push_back(21192);
   }

   if (suffix.find("QCD") != std::string::npos) {  // HT500to1000 bin missing
     eventsInSubsamples.push_back(45);
     eventsInSubsamples.push_back(141);
     eventsInSubsamples.push_back(22397);
   }

   if (suffix.find("Diboson") != std::string::npos) {  // WW, WZ, ZZ
     eventsInSubsamples.push_back(36070);
     eventsInSubsamples.push_back(9497);
     eventsInSubsamples.push_back(9690);
   }

}


void mySumWeight_filler_spring15_25ns_2lepSkim(const std::string suffix,  std::vector<Double_t> & sumWeightVector) {

  // fill vector with sumWeight for samples used in spring15 at 25 ns. The problem was that, for a given sample,  each HT bin had a different value of sumWeight. Now, since I had to
  // merge these bins to run the analysis, I had to keep track of all these values and use the proper value when the i-th bin is being analyzed in the code 
  // trees with skim on 2 leptons

  if (suffix.find("DYJetsToLL") != std::string::npos) {
     sumWeightVector.push_back(4.52712069254e+11);
   }

   if (suffix.find("WJetsToLNu") != std::string::npos) {  // binning 100-200-400-600-800-1200 (I don't have friend for other bins)
     sumWeightVector.push_back(9571170.0);
     sumWeightVector.push_back(5231856.0);
     sumWeightVector.push_back(1901705.0);
     sumWeightVector.push_back(3984529.0);
     sumWeightVector.push_back(1574633.0);
   }

   if (suffix.find("Top") != std::string::npos) { // all samples available (TToLeptons_sch missing)
     sumWeightVector.push_back(1000000.0);
     sumWeightVector.push_back(89002380545.3);
     sumWeightVector.push_back(4972773.0);
     sumWeightVector.push_back(1336800136.23);
     sumWeightVector.push_back(995600.0);
   }

   if (suffix.find("Diboson") != std::string::npos) {  // WW, WZ, ZZ
     sumWeightVector.push_back(1930000.0);
     sumWeightVector.push_back(991232.0);
     sumWeightVector.push_back(996168.0);
   }

}
 

void myEventsInSubsamples_filler_spring15_25ns_2lepSkim(const std::string suffix,  std::vector<Int_t> & eventsInSubsamples) {

  // fill vector with number of events for samples used in spring15 at 25 ns. The problem was that, for a given sample,  each HT bin had a different value of sumWeight. Now, since I had to
  // merge these bins to run the analysis, I had to keep track of all these values and use the proper value when the i-th bin is being analyzed in the code 
  // trees with skim on 2 leptons

  if (suffix.find("DYJetsToLL") != std::string::npos) {
     eventsInSubsamples.push_back(7770587);  
   }

   if (suffix.find("WJetsToLNu") != std::string::npos) {  // binning 100-200-400-600-800-1200 (I don't have friend for other bins)
     eventsInSubsamples.push_back(43651);
     eventsInSubsamples.push_back(34688);
     eventsInSubsamples.push_back(15648);
     eventsInSubsamples.push_back(35130); 
     eventsInSubsamples.push_back(13870); 
   }

   if (suffix.find("Top") != std::string::npos) { // all samples available (TToLeptons_sch missing)
     eventsInSubsamples.push_back(50771);
     eventsInSubsamples.push_back(2295837);
     eventsInSubsamples.push_back(1864475);
     eventsInSubsamples.push_back(359116); 
     eventsInSubsamples.push_back(50671); 
   }

   if (suffix.find("Diboson") != std::string::npos) {  // WW, WZ, ZZ
     eventsInSubsamples.push_back(515306);
     eventsInSubsamples.push_back(43528);
     eventsInSubsamples.push_back(74472);   
   }

}


Double_t myGetUncertainty(const mask* m, const Int_t mstep, const std::string uncDef) {

  if (uncDef == "" || uncDef == "poisson") return  m->getSqrtEvents(mstep);
  else if (uncDef == "MC") return m->getEventsErr(mstep);
  else if (uncDef.find("%") != std::string::npos) {
    std::size_t pos = uncDef.find("%");  // e.g. uncDef = "20.5%" 
    std::string substr = uncDef.substr(0,pos);  //copy till %, e.g. 20.5
    char *p; 
    Double_t unc = strtod(substr.c_str(),&p) / 100.; //e.g. 20% --> 0.2
    return  unc * m->getEvents(mstep);      // return fraction of yield according to percentage
  }

}
