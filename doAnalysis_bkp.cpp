// g++  -o doAnalysis `root-config --cflags --glibs` doAnalysis.cpp
// g++  -o doAnalysis doAnalysis.cpp `root-config --cflags --glibs`
// ./doAnalysis WaveForms_BTF/scan.dat

#include "TApplication.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"

#include <cstdio>
#include <cstdlib>
#include <stdint.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream> 
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>

#include "analysis_tools.h"
#include "InitTree_BTF.h"

int main (int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>");

  std::string inputName = std::string(argv[1]);

  // get run number
  char split_char = '/';
  std::vector<std::string> tokens;
  std::istringstream split(inputName);
  for(std::string each; getline(split, each, split_char); tokens.push_back(each));

  split_char = '_';
  std::vector<std::string> tokens_name;
  std::istringstream split_name(tokens.at(1));
  for(std::string each; getline(split_name, each, split_char); tokens_name.push_back(each));

  const int Ch_ref1 = atoi((tokens_name.at(1)).c_str());
  const int Ch_ref2 = atoi((tokens_name.at(3)).c_str());
  const int Ch_1 = atoi((tokens_name.at(5)).c_str());
  const int Ch_2 = atoi((tokens_name.at(7)).c_str());
  const int Ch_3 = atoi((tokens_name.at(9)).c_str());

  std::vector<std::string> nameMCP;
  nameMCP.push_back("MiB1");
  nameMCP.push_back("MiB2");
  nameMCP.push_back("ScB");
  nameMCP.push_back("Planacon");
  nameMCP.push_back("MiB3");
  nameMCP.push_back("RomaX");

  std::vector<std::string> pcMCP;                                                                                                                         
  for(unsigned int ii=0; ii<nameMCP.size(); ++ii) pcMCP.push_back("");
  pcMCP.at(Ch_ref1) = tokens_name.at(2);
  pcMCP.at(Ch_ref2) = tokens_name.at(4);
  pcMCP.at(Ch_1) = tokens_name.at(6);
  pcMCP.at(Ch_2) = tokens_name.at(8);
  pcMCP.at(Ch_3) = tokens_name.at(10);

  int nScan = 0;
  int nFiles = 1;
  ifstream log (argv[1], ios::in);

  TFile* out = TFile::Open("outHistos.root","recreate");  
  out->cd();

  while(log >> nFiles)
    {
      ++nScan;
      std::cout << " nFiles = " << nFiles << std::endl;
      vector<float> digiCh[9];
      float timeCF[9];
      float baseline[9];
      
      int count[5] = {0,0,0,0,0};
      int spare[5] = {0,0,0,0,0};
      int spare2[5] = {0,0,0,0,0};
      int tot_tr1 = 0, tot_tr0 = 0, trig = 0;
      int HV1 = 0, HV2 = 0, HV3 = 0;

      TH1F* hBaseLine[9];
      char h2[10];
      for (int iw=0;iw<9;iw++){
	sprintf (h2,"Bs_Ch%d_Scan%d",iw, nScan);
	hBaseLine[iw] = new TH1F( h2, "", 4096, 0., 4096.);
	hBaseLine[iw] -> SetXTitle ("BaseLine (ADC)");
      }
      TH1F* hTimeCF[9];
      char h3[10];
      for (int iw=0;iw<9;iw++){
	sprintf (h3,"TCF_Chs%d,_Scan%d",iw, nScan);
	hTimeCF[iw] = new TH1F( h3, "", 4096, 0., 4096.);
	hTimeCF[iw] -> SetXTitle ("TimeCF (ADC)");
      }

      TChain* chain = new TChain("eventRawData");
      InitTree(chain);
      for(int iFiles=0; iFiles<nFiles; ++iFiles){
	int id;
	log >> id;
	char id_str[30];
	sprintf(id_str, "WaveForms_BTF/run_IMCP_%d_*.root", id);
	chain->Add(id_str);
	cout << "Reading:  WaveForms_BTF/run_IMCP_" << id << endl;
      }
      // mib2 mib3 plana
      log >> HV1 >> HV2 >> HV3; //---config2

      for(int iEntry=0; iEntry<chain->GetEntries(); ++iEntry){
	for(int iCh=0; iCh<9; ++iCh){
                digiCh[iCh].clear();
	}
            
	chain->GetEntry(iEntry);
	if(evtNumber % 10 == 0){   //---Run<145
	//if(evtNumber % 1 == 0){      //---Run>=145

	  for(int iCh=0; iCh<nAdcChannels; ++iCh){
	    if(adcBoard[iCh] == 1 && adcChannel[iCh] == 0){
	      if(adcData[iCh] < 500) trig = 0;
	      if(adcData[iCh] > 500 && adcData[iCh] < 1500 && adcBoard[iCh]) trig = 1;
	      if(adcData[iCh] > 1500 && adcData[iCh] < 2500 && adcBoard[iCh]) trig = 2;
	    } 
	  }
	  if(trig > 1) continue; 
	  
	  //---Read digitizer samples
	  for(int iSample=0; iSample<nDigiSamples; ++iSample){
	      digiCh[digiChannel[iSample]].push_back(digiSampleValue[iSample]);
	  }
	  for(int iCh=0; iCh<6; iCh++)
	  {
	    baseline[iCh] = SubtractBaseline(5, 80, &digiCh[iCh]);
            if(iCh == 3) 
               for(unsigned int ii = 0; ii < digiCh[iCh].size(); ii++)
                   digiCh[iCh].at(ii) = -1*digiCh[iCh].at(ii);
     
	    hBaseLine[iCh]->Fill(-1.*baseline[iCh]);
	    timeCF[iCh] = TimeConstFrac(200, 275, &digiCh[iCh], 0.5);
	    hTimeCF[iCh]->Fill(timeCF[iCh]);
	  }    
	  if(ComputeIntegral(200, 275, &digiCh[Ch_ref1]) < -130 && ComputeIntegral(200, 275, &digiCh[Ch_ref2]) < -130 && trig==1){
	    ++tot_tr1;
	    if(ComputeIntegral(200, 275, &digiCh[Ch_1]) < -130) count[1] += 1;
	    if(baseline[Ch_1] < -20) spare[1] += 1;
	    if(ComputeIntegral(200, 275, &digiCh[Ch_2]) < -130) count[2] += 1;
	    if(baseline[Ch_2] < -20) spare[2] += 1;
	    if(ComputeIntegral(200, 275, &digiCh[Ch_3]) < -130) count[3] += 1;
	    if(baseline[Ch_3] < -10) spare[3] = 1;
	  }
	  else if(ComputeIntegral(200, 275, &digiCh[Ch_ref1]) >= -130 && ComputeIntegral(200, 275, &digiCh[Ch_ref2]) >= -130 && trig==0){
	    ++tot_tr0;
	    if(ComputeIntegral(200, 275, &digiCh[Ch_1]) < -130) spare2[1] += 1; 
	    if(ComputeIntegral(200, 275, &digiCh[Ch_2]) < -130) spare2[2] += 1; 
	    if(ComputeIntegral(200, 275, &digiCh[Ch_3]) < -130) spare2[3] += 1;
	  }
	} // good events
      }//loop over entries
    
        std::cout << "--------------------------" << std::endl;
        std::cout << "number of events:  " << chain->GetEntries()/10 << std::endl;
        std::cout << "Double:  " << tot_tr1 << std::endl;
        std::cout << "No e- :  " << tot_tr0 << std::endl;
        std::cout << "--------------------------" << std::endl;
        std::cout << "Ch_ref1: " << Ch_ref1 << " Ch_ref2: " << Ch_ref2 << std::endl;
        std::cout << "Measuring Eff for => Ch_1: " << Ch_1 << " Ch_2: " << Ch_2 << " Ch_3: " << Ch_3 << std::endl;
        std::cout << "Ch_1:  " << count[1] << "  " << spare[1] << "  " << spare2[1] << std::endl;
        std::cout << "Ch_2:  " << count[2] << "  " << spare[2] << "  " << spare2[2] << std::endl;
        std::cout << "Ch_3:  " << count[3] << "  " << spare[3] << "  " << spare2[3] << std::endl;
        std::cout << "--------------------------" << std::endl;
	std::cout << "HV1 = " << HV1 << " HV2 = " << HV2 << " HV3 = " << HV3 << std::endl;
        std::cout << "--------------------------" << std::endl;
        double eff1 = ((double)count[1]-(double)spare[1])/(double)tot_tr1;
        double eff2 = ((double)count[2]-(double)spare[2])/(double)tot_tr1;
        double eff3 = ((double)count[3]-(double)spare[3])/(double)tot_tr1;
        std::cout << "Ch_1 eff:        " << eff1 << std::endl;
        std::cout << "Ch_1 e_err:      " << TMath::Sqrt((eff1*(1-eff1))/tot_tr1) << std::endl;
        std::cout << "Ch_2 eff:        " << eff2 << std::endl;
        std::cout << "Ch_2 e_err:      " << TMath::Sqrt((eff2*(1-eff2))/tot_tr1) << std::endl;
        std::cout << "Ch_3 eff:    " << eff3 << std::endl;
        std::cout << "Ch_3 e_err:  " << TMath::Sqrt((eff3*(1-eff3))/tot_tr1) << std::endl;
        std::cout << "--------------------------" << std::endl;
    
	std::ofstream data1(("Data_plateau/plateau_"+nameMCP.at(Ch_1)+"_pc_"+pcMCP.at(Ch_1)+".dat").c_str(), ios::app);
	data1 << HV1 << " " <<  eff1 << " " << 0 << " " << TMath::Sqrt((eff1*(1-eff1))/tot_tr1) << std::endl;
	  
	std::ofstream data2(("Data_plateau/plateau_"+nameMCP.at(Ch_2)+"_pc_"+pcMCP.at(Ch_2)+".dat").c_str(), ios::app);
	data2 << HV2 << " " <<  eff2 << " " << 0 << " " << TMath::Sqrt((eff2*(1-eff2))/tot_tr1) << std::endl;
	  
	std::ofstream data3(("Data_plateau/plateau_"+nameMCP.at(Ch_3)+"_pc_"+pcMCP.at(Ch_3)+".dat").c_str(), ios::app);
	data3 << HV3 << " " << eff3 << " " << 0 << " " << TMath::Sqrt((eff3*(1-eff3))/tot_tr1) << std::endl;


	hBaseLine[Ch_1]->Write();
	hBaseLine[Ch_2]->Write();
	hBaseLine[Ch_3]->Write();

        chain->Delete();
        //delete chain;
    }
  out->Close();
  return 0;
}

        
