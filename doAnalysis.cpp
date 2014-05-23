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
    // get run number
    std::string inputName = std::string(argv[1]);
    char split_char = '/';
    std::vector<std::string> tokens;
    std::istringstream split(inputName);
    for(std::string each; getline(split, each, split_char); 
        tokens.push_back(each));
    
    split_char = '_';
    std::vector<std::string> tokens_name;
    std::istringstream split_name(tokens.at(1));
    for(std::string each; getline(split_name, each, split_char); 
        tokens_name.push_back(each));

    const int Ch_ref1 = atoi((tokens_name.at(1)).c_str());
    const int Ch_ref2 = atoi((tokens_name.at(3)).c_str());
    const int Ch_1 = atoi((tokens_name.at(5)).c_str());
    const int Ch_2 = atoi((tokens_name.at(7)).c_str());
    const int Ch_3 = atoi((tokens_name.at(9)).c_str());

    std::vector<std::string> nameMCP;
    nameMCP.push_back("MiB1");
    nameMCP.push_back("MiB2");
    if(argc > 2) nameMCP.at(2) = "Roma1";
    nameMCP.push_back("ScB");
    nameMCP.push_back("Planacon");
    nameMCP.push_back("MiB3");
    nameMCP.push_back("Roma2");

    std::vector<std::string> pcMCP;                                     
    for(unsigned int ii=0; ii<nameMCP.size(); ++ii) pcMCP.push_back("");
    pcMCP.at(Ch_ref1) = tokens_name.at(2);
    pcMCP.at(Ch_ref2) = tokens_name.at(4);
    pcMCP.at(Ch_1) = tokens_name.at(6);
    pcMCP.at(Ch_2) = tokens_name.at(8);
    pcMCP.at(Ch_3) = tokens_name.at(10);
    pcMCP.at(Ch_3).erase(pcMCP.at(Ch_3).size()-4, pcMCP.at(Ch_3).size()-1);
    
    int nFiles=1;
    //---open output files    
    std::ofstream data1(("Data_plateau/"+tokens_name.at(0)+"_"+nameMCP.at(Ch_1)+"_pc_"+pcMCP.at(Ch_1)+".dat").c_str());
    std::ofstream data2(("Data_plateau/"+tokens_name.at(0)+"_"+nameMCP.at(Ch_2)+"_pc_"+pcMCP.at(Ch_2)+".dat").c_str());
	std::ofstream data3(("Data_plateau/"+tokens_name.at(0)+"_"+nameMCP.at(Ch_3)+"_pc_"+pcMCP.at(Ch_3)+".dat").c_str());
    //---do runs loop
    ifstream log (argv[1], ios::in);
    while(log >> nFiles)
    {
        vector<float> digiCh[9];
        float timeCF[9];
        float baseline[9];
        float intBase[9], intSignal[9];
        int count[5]={0,0,0,0,0}, spare[5]={0,0,0,0,0}, spare2[5]={0,0,0,0,0};
        int tot_tr1=0, tot_tr0=0, trig=0;
        int HV1=0, HV2=0, HV3=0, id=0, goodEvt=1;
        TChain* chain = new TChain("eventRawData");
        InitTree(chain);
        for(int iFiles=0; iFiles<nFiles; iFiles++)
        {
            log >> id;
            char id_str[40];
            sprintf(id_str, "WaveForms_BTF/run_IMCP_%d_*.root", id);
            chain->Add(id_str);
            cout << "Reading:  WaveForms_BTF/run_IMCP_" << id << endl;
        }
        log >> HV1 >> HV2 >> HV3; 
        for(int iEntry=0; iEntry<chain->GetEntries(); iEntry++)
        {
            for(int iCh=0; iCh<9; iCh++)
            {
                digiCh[iCh].clear();
            }
            
            chain->GetEntry(iEntry);
            if(id < 145) goodEvt = 10;
            else goodEvt = 1;
            if(evtNumber % goodEvt == 0)   //---Run<145
            {
                trig = 1;
                for(int iCh=0; iCh<nAdcChannels; iCh++)
                {
                    if(adcData[iCh] > 1500 && adcBoard[iCh] == 1 && adcChannel[iCh] == 0) trig=2;
                    if(adcData[iCh] < 500 && adcBoard[iCh] == 1 && adcChannel[iCh] == 0) trig=0;
                }
                if(trig > 1) continue; 
                //---Read digitizer samples
                for(int iSample=0; iSample<nDigiSamples; iSample++)
                {
                    if(digiChannel[iSample] == 3) 
                        digiCh[digiChannel[iSample]].push_back(-digiSampleValue[iSample]);
                    else
                        digiCh[digiChannel[iSample]].push_back(digiSampleValue[iSample]);
                }
                for(int iCh=0; iCh<6; iCh++)
                {
                    baseline[iCh]=SubtractBaseline(5, 25, &digiCh[iCh]);
                    timeCF[iCh]=TimeConstFrac(30, 500, &digiCh[iCh], 0.5);
                    int t1 = (int)timeCF[iCh]/0.2 - 3;
                    int t2 = (int)timeCF[iCh]/0.2 + 17;
                    if(t1 > 30 && t1 < 1024 && t2 > 30 && t2 < 1024)
                        intSignal[iCh] = ComputeIntegral(t1, t2, &digiCh[iCh]);
                        intBase[iCh] = ComputeIntegral(26, 46, &digiCh[iCh]);
                }    
                if(intSignal[Ch_ref1] < -150 && intSignal[Ch_ref2] < -150 && trig==1) 
                {
                    tot_tr1++;
                    if(intSignal[Ch_1] < -150) count[1]=count[1]+1;
                    if(intBase[Ch_1] < -150) spare[1]=spare[1]+1;
                    if(intSignal[Ch_2] < -70) count[2]=count[2]+1;
                    if(intBase[Ch_2] < -70) spare[2]=spare[2]+1;
                    if(intSignal[Ch_3] < -150) count[3]=count[3]+1;
                    if(intBase[Ch_3] < -150) spare[3]=spare[3]+1;
                }
                else if(intSignal[Ch_ref1] >= -150 && intSignal[Ch_ref2] >= -150 && trig==0) 
                {
                    tot_tr0++;
                    if(intSignal[Ch_1] < -150) spare2[1]=spare2[1]+1; 
                    if(intSignal[Ch_1] < -70) spare2[2]=spare2[2]+1; 
                    if(intSignal[Ch_1] < -150) spare2[3]=spare2[3]+1;
                }
            }
        }
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
        std::cout << "Ch_1 eff:    " << eff1 << std::endl;
        std::cout << "Ch_1 e_err:    " << TMath::Sqrt((eff1*(1-eff1))/tot_tr1) << std::endl;
        std::cout << "Ch_2 eff:    " << eff2 << std::endl;
        std::cout << "Ch_2 e_err:    " << TMath::Sqrt((eff2*(1-eff2))/tot_tr1) << std::endl;
        std::cout << "Ch_3 eff:    " << eff3 << std::endl;
        std::cout << "Ch_3 e_err:  " << TMath::Sqrt((eff3*(1-eff3))/tot_tr1) << std::endl;
        std::cout << "--------------------------" << std::endl;
    
	    data1 << HV1 << " " <<  eff1 << " " << 0 << " " << TMath::Sqrt((eff1*(1-eff1))/tot_tr1) << std::endl;
	    data2 << HV2 << " " <<  eff2 << " " << 0 << " " << TMath::Sqrt((eff2*(1-eff2))/tot_tr1) << std::endl;
	    data3 << HV3 << " " << eff3 << " " << 0 << " " << TMath::Sqrt((eff3*(1-eff3))/tot_tr1) << std::endl;

        chain->Delete();
    }
    data1.close();
    data2.close();
    data3.close();
}

        
