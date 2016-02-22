#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>      // std::istringstream ; to read array of numbers from a line in a file
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
  Int_t HLT_flag;
  string uncertainty; // say which uncertainty is to be used for yields in a sample. Can be "poisson", "MC", "X%"
  string treePath;
  string friendTreePath;
  string sf_friendTreePath = ""; // friend with scale factor

  Int_t unweighted_event_flag = 0;  // this flag tells the user if the MC uses unit weight (using w = 1 is basically for debugging purposes)

  Int_t using_EmanTreesWithSubSamples_flag = 0; //used if reading trees from Emanuele's area
  Int_t sf_friend_flag = 0; // set to 1 if using sf_friend when reading trees from Emanuele's area
  Int_t calibEle_flag = 0;

  Int_t adishTree_flag = 0;

  Int_t signalRegion_flag = 0;
  Int_t controlSample_flag = 0;
  Int_t metResolutionAndResponse_flag = 0;
  string controlSample_boson = ""; // Z, W, GAMMA
  string fileWithSamplesPath = "";
  string filename_base = "";
  string directory_to_save_files = "";
  string directory_name = "";
  string outputFolder = "./"; // current directory by default, but it could be set as 'directory_to_save_files + directory_name'
  string dirName_suffix = ""; // user can add a suffix without modifying the name of directory in config file
  //string sf_nlo_option = "";

  if (argc > 2 ) {

    cout << "========================================================="<< endl;

    for (Int_t i = 2; i < argc; i++) {   // look at all possible options passed

      string thisArgument(argv[i]);

      if (thisArgument  == "-nw" ) {
	
	unweighted_event_flag = 1;    //-nw option stands for "no weight"
	cout << "Option " << thisArgument << " passed: using no event weight (w =1) if allowed" << std::endl;

      } else if (thisArgument  == "-at" ) {   // "at" means Adish's tree
	     
	adishTree_flag = 1;
	cout << "Option " << thisArgument << " passed: using Adish's trees" << std::endl;

      } else if (thisArgument  == "-calibEle" ) {   // use electron calibrated properties (pT, energy ecc...) instead of traditional ones
	     
	calibEle_flag = 1;
	cout << "Option " << thisArgument << " passed: using calibrated variables for electrons" << std::endl;

      } else if (thisArgument  == "-dns" ) {   // dns stands for directory name suffix: user can add a suffix to directory name without changing the config file
	// this option requires a string just after it
	 
	dirName_suffix = argv[i+1];
	cout << "Option " << thisArgument << " passed: attaching suffix \""<<argv[i+1] <<"\" at the end of directory name" << std::endl;

      } // else if (thisArgument  == "-sfnlo" ) {   // scale factor for NLO xsec for Z and W
      // 	// this option requires a string just after it, which is the name of scale factors in friend trees names sfFriend_{<sample_name>}.root. If those friend trees are not used, this option is irrelevant
	 
      // 	sf_nlo_option = argv[i+1];
      // 	cout << "Option " << thisArgument << " passed: using  "<< argv[i+1]<< " scale factor" << std::endl;

      // }


    }

    cout << "========================================================="<< endl;

  }

  ifstream inputFile(configFileName);

  if (inputFile.is_open()) {

    Double_t value;
    string name;
    string parameterName;
    string parameterType;

    cout << "========================================================="<< endl;

    while (inputFile >> parameterType ) {  // read only first object  here

      if (parameterType == "NUMBER") {

	inputFile >> parameterName >> value;  

	if (parameterName == "LEP_PDG_ID") {

	  lepton_PDGID = (Int_t) value;
	  std::cout << "lepton_pdgID = " << value <<std::endl;

	  if (fabs(lepton_PDGID) == 13) {
	    std::cout << "Analysis of Z/W->muons" << std::endl;
	  } else if (fabs(lepton_PDGID) == 11) {
	    std::cout << "Analysis of Z/W->electrons" << std::endl;
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

	if (parameterName == "HLT_FLAG") {

	  HLT_flag = (Int_t) value;

	  if (HLT_flag != 0) std::cout << "Applying HLT cut" << std::endl;
	  else std::cout << "Not applying HLT cut" << std::endl;

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

	if (parameterName == "SF_FRIEND_TREE_PATH") {

	  sf_friendTreePath = name;
	  std::cout << setw(20) << "sf_friend tree : " << sf_friendTreePath <<std::endl;

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

	if (parameterName == "PATH_TO_SAMPLES_4CR") {  // path to file with control region's samples and some options

	  controlSample_flag = 1;
	  fileWithSamplesPath = name;
	  std::cout << "Performing analysis on control samples." <<std::endl;
	  std::cout << setw(20) << "File pointing to samples: " << fileWithSamplesPath<<std::endl;

	} 

	if (parameterName == "CONTROL_SAMPLE_BOSON") {  // path to file with control region's samples and some options

	  controlSample_boson = name;
	  std::cout << "Performing analysis on " << name << " control samples." <<std::endl;

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

	} 

      }

    }

    if (directory_name != "") {
      if (calibEle_flag == 1 && fabs(lepton_PDGID) == 11) directory_name += "_CalibEle";
      if (unweighted_event_flag == 1) directory_name += "_weq1";
      //if (sf_nlo_option != "") directory_name += ("_" + sf_nlo_option); // keep commented it's a nightmare to have different directories when I just need one more histogram
      if (dirName_suffix != "") directory_name += ("_" + dirName_suffix); // this addition should be left as the last one when forming directory name
      std::cout << "Files will be saved in directory named '" << directory_name  << "' ." <<std::endl;
    }

    cout << "========================================================="<< endl;

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

      cout << "========================================================="<< endl;

      reportFile.close();

    }

  }

  // ==================================================


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
  selectionDefinition.push_back("entry point"); //other steps will be set in the following using the selectionManager object in the analyzer

  ifstream sampleFile(fileWithSamplesPath.c_str());
  Int_t fileEndReached_flag = 0;

  if (sampleFile.is_open()) {

    while (fileEndReached_flag == 0) {

      //std::cout << "CHECK IN READING SAMPLE FILE" << std::endl;
      std::vector<std::string> subSampleNameVector; // if using Emanuele's tree, for each sample (e.g. ZJetsToNuNu) there is more than one subsample (e.g. different HT bins or different processes). This vector holds this subsamples's name to be chained (previously I used to merge and copy them from Emanuele's area to mine, so that I only had one tree or friend if any). For backward compatibility, I keep the possibility to use the merged trees when reading the file

      //============================================/

      cout << "========================================================="<< endl;

      while ( (sampleFile >> parameterType) && (!(parameterType == "#")) ) {  // read only first object  here: if it is '#' it signals that another part is starting, and this while loop ends, and analysis starts.
	//std::cout<<"CHECK 2 " <<std::endl;
	//std::cout << "ParameterType : " << parameterType << std::endl;

	if (parameterType == "#STOP") fileEndReached_flag = 1;   //tells me that the file is ended and I don't need to go on reading it.

	if (parameterType == "NUMBER") {

	  sampleFile >> parameterName >> value;  

	  if (parameterName == "ISDATA_FLAG") {

	    isdata_flag = (Int_t) value;

	    if (isdata_flag == 0) std::cout << "Running on MonteCarlo" << std::endl;
	    else std::cout << "Running on data" << std::endl;

	  }

	  if (parameterName == "SF_FRIEND_FLAG") {

	    sf_friend_flag = (Int_t) value;
	    if (sf_friend_flag == 0) std::cout << "Not using scale factor friend trees (sf_friend)" << std::endl;
	    else std::cout << "Using scale factor friend trees (sf_friend)" << std::endl;

	  }

	} else if (parameterType == "STRING") {

	  sampleFile >> parameterName >> name;

	  //std::cout << "parameterName" << parameterName << std::endl;

	  if (parameterName == "SAMPLE_NAME") {

	    sampleName.push_back(name.c_str());
	    std::cout << setw(20) << "sample name : " << name <<std::endl;

	  }  

	  if (parameterName == "UNCERTAINTY") {

	    uncertainty = name;
	    std::cout << setw(20) << "sample uncetainty : " << name <<std::endl;

	  }  

	  if (parameterName == "TREE_PATH") {

	    treePath = name;
	    std::cout << setw(20) << "tree : " << name <<std::endl;

	  }

	  if (parameterName == "FRIEND_TREE_PATH") {

	    friendTreePath = name;
	    std::cout << setw(20) << "friend tree : " << name <<std::endl;

	  } 

	  if (parameterName == "SF_FRIEND_TREE_PATH") {

	    sf_friendTreePath = name;
	    std::cout << setw(20) << "sf_friend tree : " << name <<std::endl;

	  } 

	} else if(parameterType == "ARRAY_STR") {

	  //std::cout<<"CHECK 0 " <<std::endl;
	  using_EmanTreesWithSubSamples_flag = 1;
	  sampleFile >> parameterName;

	  if (parameterName == "SUB_SAMPLE_NAMES") { 

	    cout << right << setw(20) << parameterName << "  ";
	    string stringvalues;
	    getline(sampleFile, stringvalues);    // read whole line starting from current position (i.e. without reading ARRAY_STR and SAMPLE_NAME)
	    istringstream iss(stringvalues);
	    std::string subSample;

	    while(iss >> subSample) {
	    
	      subSampleNameVector.push_back(subSample);
	      cout << subSampleNameVector.back() << " ";

	    }

	    cout << endl;

	  }

	}

      }

      if ( !fileEndReached_flag) {

	TChain* chain = new TChain("tree");
	TChain* chFriend = new TChain("mjvars/t");
	TChain* chSfFriend = new TChain("sf/t");

	//std::cout<<"CHECK 3 " <<std::endl;

	std::cout << "using_EmanTreesWithSubSamples_flag : " << using_EmanTreesWithSubSamples_flag << std::endl;

	if(using_EmanTreesWithSubSamples_flag == 1) {

	  //std::cout<<"CHECK 1 " <<std::endl;

	  std::cout << "Creating chain ..." << std::endl;

	  for(Int_t i = 0; i < subSampleNameVector.size(); i++) {
	   
	    std::string treeRootFile = treePath + subSampleNameVector[i] + "/treeProducerDarkMatterMonoJet/tree.root"; 
	    std::string friend_treeRootFile = treePath + "friends/evVarFriend_" + subSampleNameVector[i]+ ".root"; 
	    std::string sf_friend_treeRootFile = treePath + "friends/sfFriend_" + subSampleNameVector[i]+ ".root"; 

	    chain->Add(TString(treeRootFile.c_str()));
	    chFriend->Add(TString(friend_treeRootFile.c_str()));
	    if (sf_friend_flag != 0) chSfFriend->Add(TString(sf_friend_treeRootFile.c_str()));

	  }

	  std::cout << "Adding friend to chain ..." << std::endl;	    
	  chain->AddFriend(chFriend);  //adding whole friend chain as friend

	  if (sf_friend_flag != 0) {
			  
	    std::cout << "Adding friend with scale factors to chain ..." << std::endl;
	    chain->AddFriend(chSfFriend);  //adding whole friend chain as friend

	  }

	} else {  // here started from if(using_EmanTreesWithSubSamples_flag == 1)

	  //std::cout<<"CHECK 4 " <<std::endl;

	  std::cout << "Creating chain ..." << std::endl;
	  chain->Add(TString(treePath.c_str()));

	  std::cout << "Adding friend to chain ..." << std::endl;
	  chain->AddFriend("mjvars/t",TString(friendTreePath.c_str()));

	  if (sf_friendTreePath != "") {
	    std::cout << "Adding friend with scale factors to chain ..." << std::endl;	  
	    chain->AddFriend("sf/t",TString(sf_friendTreePath.c_str()));
	  }

	}	 

	if(!chain) {
	  std::cout << "Error: chain not created. End of programme" << std::endl;
	  exit(EXIT_FAILURE);
	}
	std::cout << "analysing sample n " << (nSample+1) << " : " << sampleName[nSample]  << std::endl; 
	std::cout<<chain->GetEntries()<<std::endl;      
	//================ Run Analysis
	//zmumujetsAna tree( chain );
 
	if (signalRegion_flag == 1) {

	  monojet_SignalRegion tree( chain);
	  //cout << " CHECK IN MAIN " << endl;
	  tree.setBasicConf(sampleName[nSample].c_str(), uncertainty, configFileName, isdata_flag, unweighted_event_flag, sf_friend_flag);
	  tree.setDirNameSuffix(dirName_suffix);
	  // tree.set_SF_NLO_name(sf_nlo_option);
	  tree.loop(yieldsRow, efficiencyRow, uncertaintyRow); 
	  tree.analysisSelectionManager.exportDefinition(&selectionDefinition);

	} else if (controlSample_flag == 1) {

	  if (controlSample_boson == "Z") {
	    zlljetsControlSample tree( chain);
	    tree.setBasicConf(sampleName[nSample].c_str(), uncertainty, configFileName, isdata_flag, unweighted_event_flag, sf_friend_flag);
	    tree.setDirNameSuffix(dirName_suffix);
	    // tree.set_SF_NLO_name(sf_nlo_option);
	    if (calibEle_flag == 1 && fabs(lepton_PDGID) == 11) tree.setCalibEleFlag();
	    tree.loop(yieldsRow, efficiencyRow, uncertaintyRow);
	    tree.analysisSelectionManager.exportDefinition(&selectionDefinition);
	  } else if (controlSample_boson == "W") {
	    wlnujetsControlSample tree( chain);
	    tree.setBasicConf(sampleName[nSample].c_str(), uncertainty, configFileName, isdata_flag, unweighted_event_flag, sf_friend_flag);
	    tree.setDirNameSuffix(dirName_suffix);
	    // tree.set_SF_NLO_name(sf_nlo_option);
	    if (calibEle_flag == 1 && fabs(lepton_PDGID) == 11) tree.setCalibEleFlag();
	    tree.loop(yieldsRow, efficiencyRow, uncertaintyRow);
	    tree.analysisSelectionManager.exportDefinition(&selectionDefinition);
	  } 

	} else if (metResolutionAndResponse_flag == 1) {

	  zlljets_metResoResp tree( chain);
	  tree.setBasicConf(sampleName[nSample].c_str(), uncertainty, configFileName, isdata_flag, unweighted_event_flag, sf_friend_flag);
	  tree.loop(yieldsRow, efficiencyRow, uncertaintyRow);
	  tree.analysisSelectionManager.exportDefinition(&selectionDefinition);

	}

	// now that chains are set, reset some variables' value
	using_EmanTreesWithSubSamples_flag = 0;
	sf_friendTreePath = ""; //set to void string, then in the previous loop it will be initialized if present in the samplePathFile (sometimes there are no sf_friend and in that case, the sf_friend is not added as friend: there is an if condition asking whether sf_friendTreePath == "" or not)
	sf_friend_flag = 0; // set to zero so that for the next sample it will be left 0 or set to 1 depending on the samplePathFile
	//sf_friend_flag is used only when reading trees from emanuele's area: it signals that there are sf_friend to add as friend tree (friend with scale factors). We assume that normal friends are always present. In a future, we might always have sf_friend as well and this flag will not be needed any more

	cout << endl;   cout << endl;
	nSample++;
	delete chain;
	delete chFriend;
	delete chSfFriend;

      } // end of  if ( !fileEndReached_flag)

      cout << "========================================================="<< endl;

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
  //if (unweighted_event_flag) finalFileName += "_weq1";   //means with weights equal to 1 (for debugging purposes)  // not needed anymore since files are saved in a directory whose name would already have the _weq1 string in the name
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

}

