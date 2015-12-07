#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <math.h>
#include <sys/types.h>  // for mkdir and stat
#include <sys/stat.h>    // for mkdir and stat
#include <unistd.h>      // for stat

#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
#include <TFile.h>
#include "EmanTreeAnalysis.h"
#include "AdishTreeAnalysis.h"

using namespace myAnalyzerTEman;
using namespace myAnalyzerTAdish;

using namespace std;


int main(int argc, char* argv[]) {


//================ Creating chain 

  std::cout << std::endl;

  if (argc < 2) {
    std::cout << "Not enough arguments: launch as -> "; 
    std::cout << argv[0] << " inputfile " <<std::endl;
    std::cout << "inputfile is a configuration file containing thresholds and other parameters" << std::endl;
    exit(EXIT_FAILURE);
  }

  char configFileName[200];
  std::strcpy(configFileName,argv[1]);

  Int_t lepton_PDGID;
  Int_t isdata_flag;
  Int_t tau_veto_flag;
  Double_t metnolep_start;
  Int_t met_filters_flag;
  string treePath;
  string friendTreePath;
  string option = "";

  Int_t unweighted_event_flag = 0;  // this flag tells the user if the MC uses unit weight (using w = 1 is basically for debugging purposes)
  Int_t adishTree_flag = 0;  // tells user I'm using Adish's tree

  // following are for CS analysis
  Int_t signalRegion_flag = 0;
  Int_t controlSample_flag = 0;
  Int_t metResolutionAndResponse_flag = 0;
  string fileWithSamplesPath = "";
  string filename_base = "";
  string directory_to_save_files = "";
  string directory_name = "";
  string outputFolder = "./"; // current directory by default, but it could be set as 'directory_to_save_files + directory_name'

  if (argc > 2 ) {

    for (Int_t i = 2; i < argc; i++) {   // look at all possible options passed

      string thisArgument(argv[i]);

      if (thisArgument  == "-nw" ) {
	
	unweighted_event_flag = 1;    //-nw option stands for "no weight"
	cout << "Option " << thisArgument << " passed: using no event weight (w =1) if allowed" << std::endl;

      }

      if (thisArgument  == "-at" ) {   // "at" means Adish's tree
	     
	adishTree_flag = 1;
	cout << "Option " << thisArgument << " passed: using Adish's trees" << std::endl;

      }

    }

  }

  ifstream inputFile(configFileName);

  if (inputFile.is_open()) {

    Double_t value;
    string name;
    string parameterName;
    string parameterType;

    while (inputFile >> parameterType ) {  // read only first object  here

      if (parameterType == "NUMBER") {

	inputFile >> parameterName >> value;  

	if (parameterName == "LEP_PDG_ID") {

	  lepton_PDGID = (Int_t) value;
	  std::cout << "lepton_pdgID = " << value <<std::endl;

	  if (fabs(lepton_PDGID) == 13) {
	    std::cout << "Analysis of Z->mumu" << std::endl;
	  } else if (fabs(lepton_PDGID) == 11) {
	    std::cout << "Analysis of Z->ee" << std::endl;
	  } else if (fabs(lepton_PDGID) == 12 || fabs(lepton_PDGID) == 14 || fabs(lepton_PDGID) == 16) {
	    std::cout << "Analysis of Z->nunu" << std::endl;
	  }

	} 

	if (parameterName == "ISDATA_FLAG") {

	  isdata_flag = (Int_t) value;

	  if (isdata_flag == 0) std::cout << "Running on MonteCarlo" << std::endl;
	  else std::cout << "Running on data" << std::endl;

	}

	if (parameterName == "TAU_VETO_FLAG") {

	  tau_veto_flag = (Int_t) value;

	  if (tau_veto_flag == 0) std::cout << "Not applying tau veto" << std::endl;
	  else std::cout << "Applying tau veto" << std::endl;

	}

	if (parameterName == "METNOLEP_START") {

	  metnolep_start = value;

	  if (metnolep_start != 0) std::cout << "Cutting on recoil > " << metnolep_start << " GeV" << std::endl;
	  else std::cout << "No cut on recoil" << std::endl;

	}

	if (parameterName == "MET_FILTERS_FLAG") {

	  met_filters_flag = (Int_t) value;

	  if (met_filters_flag != 0) std::cout << "Applying met filters" << std::endl;
	  else std::cout << "Not applying met filters" << std::endl;

	}

      } else if (parameterType == "STRING") {

	inputFile >> parameterName >> name;

	if (parameterName == "TREE_PATH") {

	  treePath = name;
	  std::cout << setw(20) << "tree : " << treePath <<std::endl;

	}  

	if (parameterName == "FRIEND_TREE_PATH") {

	  friendTreePath = name;
	  std::cout << setw(20) << "friend tree : " << friendTreePath <<std::endl;

	}  

	if (parameterName == "OPTION") {

	  option = name;
	  std::cout << setw(20) << "option : " << option <<std::endl;
  

	} 

	if (parameterName == "FILENAME_BASE") {  // base for file's name

	  filename_base = name;
	  std::cout << setw(20) << "filename_base: " << filename_base <<std::endl;

	} 

	if (parameterName == "PATH_TO_SAMPLES_4SR") {  // path to file with signal region's samples and some options

	  signalRegion_flag = 1;
	  fileWithSamplesPath = name;
	  std::cout << "Performing analysis on signal region." <<std::endl;
	  std::cout << setw(20) << "File pointing to samples: " << fileWithSamplesPath<<std::endl;

	} 

	if (parameterName == "PATH_TO_SAMPLES_4CS") {  // path to file with control region's samples and some options

	  controlSample_flag = 1;
	  fileWithSamplesPath = name;
	  std::cout << "Performing analysis on control samples." <<std::endl;
	  std::cout << setw(20) << "File pointing to samples: " << fileWithSamplesPath<<std::endl;

	} 

	if (parameterName == "PATH_TO_SAMPLES_4MET_RESO_RESP") {  // path to file with CS samples and some options

	  metResolutionAndResponse_flag = 1;
	  fileWithSamplesPath = name;
	  std::cout << "Performing MET resolution and response analysis." <<std::endl;
	  std::cout << setw(20) << "File pointing to samples: " << fileWithSamplesPath<<std::endl;

	} 	

	if (parameterName == "DIRECTORY_PATH") {  // name of directory where files are saved

	  directory_to_save_files = name;
	  std::cout << "Files will be saved in '" << name << "' ." <<std::endl;

	} 

	if (parameterName == "DIRECTORY_NAME") {  // name of directory where files are saved

	  directory_name = name;
	  std::cout << "Files will be saved in directory named '" << name << "' ." <<std::endl;

	} 

      }

    }

    inputFile.close();
                                                                                                                         
  } else {

    std::cout << "Error: could not open file " << configFileName << std::endl;
    exit(EXIT_FAILURE);

  }

  // =========  Creating directory where files are saved ============

  if ( (directory_to_save_files != "") && (directory_name != "") ) {

    struct stat st = {0};

    outputFolder = directory_to_save_files + directory_name + "/";
    std::cout << "Creating new directory " << outputFolder << " ... " << std::endl;

    if (stat(outputFolder.c_str(), &st) == -1) {

      if (mkdir(outputFolder.c_str(),0755) == 0) {   // 755 refers to access rights

	std::cout << "Directory was created successfully!" << std::endl; 
    
      } else std::cout << "Error occurred when creating directory!" << std::endl; 
      // error will never occur with stat(): stat is -1 if directory doesn't exist, so it is created by mkdir(), which fails if directory already exists (but in this case the stat() prevents programme from entering and doing mkdir()

    } else std::cout << "Warning: maybe directory already exists" << std::endl;

    // ================== saving content of config file in "report.txt" =====

    string reportFileName = "report.txt";
    ofstream reportFile((outputFolder + reportFileName).c_str(),ios::out);

    if ( !reportFile.is_open() ) {

      cout<<"Error: unable to open file " << reportFileName <<" !"<<endl;
      exit(EXIT_FAILURE);
     
    } else {

      string str;

      inputFile.open(configFileName); //opening inputFile named configFileName again to save content in file named "report.txt"

      if (inputFile.is_open()) {
     
	cout << "Saving content of " << configFileName << " file in "<< reportFileName << " located in " << outputFolder << endl;
	reportFile << "Content of " << configFileName << endl;
	reportFile << endl;

	while (getline(inputFile,str)) {

	  reportFile << str << endl;

	}
     
	inputFile.close();
        
	reportFile << endl;
	reportFile << endl;
                                                                                                             
      } else {

	cout << "Error: could not open file " << configFileName << " to save content in "<< reportFileName << endl;
	exit(EXIT_FAILURE);

      }

      //now also opening inputFile named fileWithSamplesPath to save content in file named "report.txt"

      inputFile.open(fileWithSamplesPath.c_str());

      if (inputFile.is_open()) {

	cout << "Saving content of " << fileWithSamplesPath << " file in "<< reportFileName << " located in " << outputFolder << endl;
	reportFile << "Content of " << fileWithSamplesPath << endl;
	reportFile << endl;

	while (getline(inputFile,str)) {

	  reportFile << str << endl;

	}

	inputFile.close();
  
      } else {

	cout << "Error: could not open file " << fileWithSamplesPath << " to save content in "<< reportFileName << endl;
	exit(EXIT_FAILURE);

      }

      reportFile.close();

    }

  }

  // ==================================================

  if (signalRegion_flag == 1 || controlSample_flag == 1 || metResolutionAndResponse_flag == 1 ) {

    Double_t value;
    string name;
    string parameterName;
    string parameterType;

    vector< Double_t > yieldsRow;
    vector< Double_t > efficiencyRow;
    vector< Double_t > uncertaintyRow;
    Int_t nSample = 0;
    std::vector<std::string> sampleName;
    
    std::vector<std::string> selectionDefinition;

    if (signalRegion_flag == 1) {
    
      selectionDefinition.push_back("entry point");  
      if (met_filters_flag == 1) selectionDefinition.push_back("met filters");
      selectionDefinition.push_back("muon veto");
      selectionDefinition.push_back("electron veto");
      if (tau_veto_flag) selectionDefinition.push_back("tau veto");
      selectionDefinition.push_back("photon veto");
      selectionDefinition.push_back("bjet veto");
      selectionDefinition.push_back("metNoMu > 200");
      selectionDefinition.push_back("jet1pt");
      selectionDefinition.push_back("jet1 cleaning");
      selectionDefinition.push_back("dphiMin(j,Met)");    

    } else if (controlSample_flag == 1) {

      selectionDefinition.push_back("entry point"); // include genLep for DY
      if (met_filters_flag == 1) selectionDefinition.push_back("met filters");
      selectionDefinition.push_back("2lep loose");
      selectionDefinition.push_back(">0 tight lep");
      selectionDefinition.push_back("2lep SF/OS");
      selectionDefinition.push_back("mll");
      if (fabs(lepton_PDGID) == 13) selectionDefinition.push_back("electron veto");
      else if (fabs(lepton_PDGID) == 11) selectionDefinition.push_back("muon veto");
      if (tau_veto_flag) selectionDefinition.push_back("tau veto");
      selectionDefinition.push_back("photon veto");
      if (metnolep_start != 0) selectionDefinition.push_back("recoil > 200");
      selectionDefinition.push_back("bjet veto");
      selectionDefinition.push_back("jet1pt");
      selectionDefinition.push_back("jet1 cleaning");
      selectionDefinition.push_back("dphiMin(j,Met)");        
      selectionDefinition.push_back("lep match");

    } else if (metResolutionAndResponse_flag == 1) {

      selectionDefinition.push_back("entry point");        
      selectionDefinition.push_back("preselection");   // include genLep, HLT, MetNoLep
      selectionDefinition.push_back("2lep SF/OS");
      selectionDefinition.push_back("2lep loose");
      if (fabs(lepton_PDGID) == 13) selectionDefinition.push_back("muons");
      else if (fabs(lepton_PDGID) == 11) selectionDefinition.push_back("electrons");
      selectionDefinition.push_back("tight Tag");
      selectionDefinition.push_back("mll");
      selectionDefinition.push_back("bjet veto");
      selectionDefinition.push_back("jet1pt");
      selectionDefinition.push_back("dphiMin(j,Met)");
      selectionDefinition.push_back("jet1 cleaning");
      if (fabs(lepton_PDGID) == 13) selectionDefinition.push_back("electron veto");
      else if (fabs(lepton_PDGID) == 11) selectionDefinition.push_back("muon veto");
      selectionDefinition.push_back("photon veto");
      if (tau_veto_flag) selectionDefinition.push_back("tau veto");
      selectionDefinition.push_back("lep match");

    }

    ifstream sampleFile(fileWithSamplesPath.c_str());
    Int_t fileEndReached_flag = 0;

    if (sampleFile.is_open()) {

      while (fileEndReached_flag == 0) {

	//============================================/
  
	while ( (sampleFile >> parameterType) && (!(parameterType == "#")) ) {  // read only first object  here: if it is '#' it signals that another part is starting

	  if (parameterType == "#STOP") fileEndReached_flag = 1;   //tells me that the file is ended and I don't need to go on reading it.

	  if (parameterType == "NUMBER") {

	    sampleFile >> parameterName >> value;  

	    if (parameterName == "ISDATA_FLAG") {

	      isdata_flag = (Int_t) value;

	      if (isdata_flag == 0) std::cout << "Running on MonteCarlo" << std::endl;
	      else std::cout << "Running on data" << std::endl;

	    }

	  } else if (parameterType == "STRING") {

	    sampleFile >> parameterName >> name;

	    if (parameterName == "SAMPLE_NAME") {

	      sampleName.push_back(name.c_str());
	      std::cout << setw(20) << "sample name : " << name <<std::endl;

	    }  

	    if (parameterName == "TREE_PATH") {

	      treePath = name;
	      std::cout << setw(20) << "tree : " << treePath <<std::endl;

	    }   

	    if (parameterName == "FRIEND_TREE_PATH") {

	      friendTreePath = name;
	      std::cout << setw(20) << "friend tree : " << friendTreePath <<std::endl;

	    }  

	  }

	}

	if ( !fileEndReached_flag) {

	  std::cout << "Creating chain ..." << std::endl;
	  TChain* chain = new TChain("tree");
	  chain->Add(TString(treePath.c_str()));

	  std::cout << "Adding friend to chain ..." << std::endl;
	  TChain* chFriend = new TChain("mjvars/t");
	  chain->AddFriend("mjvars/t",TString(friendTreePath.c_str()));
	
	  if(!chain) {
	    std::cout << "Error: chain not created. End of programme" << std::endl;
	    exit(EXIT_FAILURE);
	  }
	  std::cout << "analysing sample n " << (nSample+1) << " : " << sampleName[nSample]  << std::endl; 
	  std::cout<<chain->GetEntries()<<std::endl;      
	  //================ Run Analysis
	  //zmumujetsAna tree( chain );
 
	  if (signalRegion_flag == 1) {

	    monojet_SignalRegion tree( chain , sampleName[nSample].c_str());
	    tree.loop(configFileName, isdata_flag, unweighted_event_flag, yieldsRow, efficiencyRow, uncertaintyRow); 

	  } else if (controlSample_flag == 1) {

	    zlljetsControlSample tree( chain , sampleName[nSample].c_str());
	    tree.loop(configFileName, isdata_flag, unweighted_event_flag, yieldsRow, efficiencyRow, uncertaintyRow); 

	  } else if (metResolutionAndResponse_flag == 1) {

	    zlljets_metResoResp tree( chain , sampleName[nSample].c_str());
	    tree.loop(configFileName, isdata_flag, unweighted_event_flag, yieldsRow, efficiencyRow, uncertaintyRow); 

	  }

	  cout << endl;   cout << endl;
	  nSample++;
	  delete chain;
	  delete chFriend;

	}

	//============================================/
 
      }  //end of while(fileEndReanched_flag == 0)

      sampleFile.close();

    } else {   // end of if (sampleFile.is_open())

      cout << "Error: could not open file " << fileWithSamplesPath<< endl;
      exit(EXIT_FAILURE);

    }

    if ( yieldsRow.size() != efficiencyRow.size() ) {
      std::cout << "Warning:  different number of steps for yields and efficiencies." << std::endl; 
      std::cout << "Will use the bigger one." << std::endl; 
    }

    Int_t selectionSize = (yieldsRow.size() >= efficiencyRow.size()) ? (yieldsRow.size()/nSample) : (efficiencyRow.size()/nSample);

    // ==================================

    // before printing the table, I add another row with the sum of all "non data" column. 
    // For the SR, it would be the sum of all backgrounds, regardless they are data-driven or MC estimate
    // For the CR, it should be the sum of MC background for the Z(ll) sample (now for semplicity it is the sum of all MC).

    // suppose I have tre samples A, B, C and a selection with three steps 1, 2, 3. Then, yieldsRow is filled as A123 B123 C123 (9 entries, the letter refers to the specific sample)
    // to build the table, which is created row by row, we read the vector as A1, B1, C1, A2, B2 ... so that the table looks like
    //
    //   sample A     sample B     sample C
    //        A1                 B1                C1
    //        A2                 B2                C2
    //        A3                 B3                C3
    //
    // next to yields there are the efficiency or the uncertainties (could add both but table would get too crowded)
    //
    // now we want to add an additional column with the sum of entries for the different samples
    // so now we add S123 to yieldsRow (S stands for sum), where Si = Ai + Bi +Ci

    for (Int_t i = 0; i < selectionSize; i++) {

      Double_t lastValue = yieldsRow.back();   // keep track of last element to make the efficiency ratio 
      // when i = 0, it is the last value before adding the "non-data" column (C3 in the previous example), but it is not used because the efficiency is automatically set to 1.0
      // for the other values of i, it is the last value (call it a), now we compute the following value (call it b) and compute the efficiency as b/a

      yieldsRow.push_back(0.0);   // adding new element for each selection step
      efficiencyRow.push_back(0.0);
      uncertaintyRow.push_back(0.0);

      for(Int_t j = 0; j < nSample; j++) {  

	// must skip sample with data, if present

	if ( (std::strcmp("data",sampleName[j].c_str())) ) {  //std::strcmp returns 0 when the strings are equal. otherwise it returns a non zero value

	  // adding all values in the same row (i.e. for the same selection step). If an entry is negative (because that step was not considered for that sample), the previous step is summed)
	  Int_t vectorElement = i + j * selectionSize;

	  if (yieldsRow.at(vectorElement) < 0) vectorElement = (i - 1) + j * selectionSize;  
	  // can be negative when that step was not filled for a given sample (e.g. the recoGen match is only for DY sample in Z+jets CS)
	  
	  yieldsRow.back() += yieldsRow.at(vectorElement);  
	  uncertaintyRow.back() += uncertaintyRow.at(vectorElement) * uncertaintyRow.at(vectorElement);   // sum in quadrature of samples' uncertainties  

	}

      }     // end of for(Int_t j = 0; j < nSample; j++)

      uncertaintyRow.back() = sqrt(uncertaintyRow.back());
      if (i == 0) efficiencyRow.back() = 1.0;
      else if ( (i != 0) && ( lastValue == 0 )  ) efficiencyRow.back() = 1.0000;  
      // in the line above, if previous yield is 0, the next is also 0 and the efficiency is set to 1.0  (otherwise it would be of the form 0/0)
      else efficiencyRow.back() = yieldsRow.back()/lastValue;
     
    }

    // ==================================

    selectionSize = (yieldsRow.size() >= efficiencyRow.size()) ? (yieldsRow.size()/(nSample+1)) : (efficiencyRow.size()/(nSample+1));  
    //(nSample+1) because now there is also the sum on MC entries

    FILE* fp;
    string finalFileName = filename_base;
    finalFileName += "_yieldsTable";
    if (unweighted_event_flag) finalFileName += "_weq1";   //means with weights equal to 1 (for debugging purposes)
    finalFileName += ".dat";

    if ( (fp=fopen((outputFolder + finalFileName).c_str(),"w")) == NULL) {

      cout<<"Error: '"<<finalFileName<<"' not opened"<<endl;

    } else {

      cout<<"creating file '"<<finalFileName<<"' to save table with yields in folder " << outputFolder << " ..."<<endl;
      fprintf(fp,"#    step         ");

      for(Int_t i = 0; i <= nSample; i++) {

	if (i == nSample) fprintf(fp,"%-16s ","all non-data");  // last column in file will hold the sum af all MC or backgrounds
	else fprintf(fp,"%-16s ",sampleName[i].c_str());

      }

      fprintf(fp,"\n");

      for (Int_t i = 0; i < selectionSize; i++) {

	fprintf(fp,"%-16s",selectionDefinition[i].c_str());

	for(Int_t j = 0; j <= nSample; j++) {
	  
	  if (j == nSample) {

	    if (yieldsRow.at( i + j * selectionSize) < 0) {

	      string space = "//";
	      fprintf(fp,"%7s ",space.c_str());
	      fprintf(fp,"%5s    ",space.c_str());

	    } else { 

	      if (yieldsRow.at( i + j * selectionSize) < 10) fprintf(fp,"%7.1lf ",yieldsRow.at( i + j * selectionSize));  //	j * selectionSize refers to number for a sample, i refers to the selection step   
	      else fprintf(fp,"%7.0lf ",yieldsRow.at( i + j * selectionSize));

	      if (uncertaintyRow.at( i + j * selectionSize) < 10) fprintf(fp,"%7.1lf   ",uncertaintyRow.at( i + j * selectionSize));
	      else fprintf(fp,"%7.0lf   ",uncertaintyRow.at( i + j * selectionSize));

	    }

	  } else {

	    if (yieldsRow.at( i + j * selectionSize) < 0) {

	      string space = "//";
	      fprintf(fp,"%7s ",space.c_str());
	      fprintf(fp,"%5s    ",space.c_str());

	    } else { 

	      if (yieldsRow.at( i + j * selectionSize) < 10) fprintf(fp,"%7.1lf ",yieldsRow.at( i + j * selectionSize));  //	j * selectionSize refers to number for a sample, i refers to the selection step   
	      else fprintf(fp,"%7.0lf ",yieldsRow.at( i + j * selectionSize));
	      fprintf(fp,"%5.1lf%%   ",(100 * efficiencyRow.at( i + j * selectionSize)));

	    }

	  }

	} 

	fprintf(fp,"\n");

      }

      fclose(fp);

    }

    return 0;

  }  //end of  "if (controlSample_flag == 1 || signalRegion_flag == 1 || ...)"


  //if (!isdata_flag && unweighted_event_flag) std::cout << "Using unweighted events (w = 1) " << std::endl;
  //std::cout<<"Using Emanuele's trees with no skim"<<std::endl;
  std::cout<< std::endl;
  if (adishTree_flag) std::cout<<"Using Adish's trees "<<std::endl;
  else std::cout<<"Using Emanuele's trees "<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain;
  if (adishTree_flag) chain = new TChain("tree/tree");
  else chain = new TChain("tree");
  //chain->Add("/cmshome/ciprianim/edimarcoTree/tree_noSkim/DYJetsToLL_M50/tree.root");
  chain->Add(treePath.c_str());
  TChain* chFriend = new TChain("mjvars/t");
  //chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_noSkim/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");
  if (!adishTree_flag) chain->AddFriend("mjvars/t",friendTreePath.c_str());

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  if (adishTree_flag) {

    zlljets_resoResp_forAdish tree( chain );
    tree.loop(configFileName, isdata_flag, unweighted_event_flag);

  } else {

    if ( !(std::strcmp("axel",option.c_str()))) {

      zlljets_Axe_noSkim_light tree( chain );
      tree.loop(configFileName);

    } else if ( !(std::strcmp("rerel",option.c_str()))) {

      zlljets_resoResp tree( chain );
      tree.loop(configFileName, isdata_flag, unweighted_event_flag);

    } else {

      std::cout << "Option '" << option << "' not found. End of programme" << std::endl;
      return 0;

    }

  }
  
  delete chain;
  delete chFriend;

  return 0;
}

