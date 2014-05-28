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
    
    int nFiles=1, iRun=0, goodEvt=0;
    //---usefull histos
    TH1F* chHistoBase[9];
    TH1F* chHistoSignal[9];
    TH1F* timeDiffHisto[9];
            //---histos initialization
        for(int iCh=0; iCh<9; iCh++)
        {
            char h1[30], h2[30], h3[30];
            sprintf(h1, "histoBase%d_%d", iCh, iRun);
            sprintf(h2, "histoSignal%d_%d", iCh, iRun);
            sprintf(h3, "histoTime%d_%d", iCh, iRun);
            chHistoBase[iCh] = new TH1F(h1,h1,2000,-1000,1000);
            chHistoSignal[iCh] = new TH1F(h2, h2,30000,-30000,1000);
            timeDiffHisto[iCh] = new TH1F(h3, h3,4000,-10,10);
        } 
    TFile* out = TFile::Open("outHistos.root","recreate");  
    out->cd();
    
    //---do runs loop
    ifstream log (argv[1], ios::in);
    while(log >> nFiles)
    {
        vector<float> digiCh[9];
        float timeCF[9];
        float baseline[9];
        int count[5]={0,0,0,0,0}, spare[5]={0,0,0,0,0}, spare2[5]={0,0,0,0,0};
        int tot_tr1=0, tot_tr0=0, trig=0;
        int HV1=0, HV2=0, HV3=0;
        //---data chain
        TChain* chain = new TChain("eventRawData");
        InitTree(chain);
        for(int iFiles=0; iFiles<nFiles; iFiles++)
        {
            log >> iRun;
            char iRun_str[30];
            sprintf(iRun_str, "WaveForms_BTF/run_IMCP_%d_*.root", iRun);
            chain->Add(iRun_str);
            cout << "Reading:  WaveForms_BTF/run_IMCP_" << iRun << endl;
        }
        log >> HV1 >> HV2 >> HV3;
        //if(iRun != 68) continue; 
        for(int iEntry=0; iEntry<chain->GetEntries(); iEntry++)
        {
            //---always clear the std::vector !!!
            for(int iCh=0; iCh<9; iCh++)
            {
                digiCh[iCh].clear();
            }
            //---Read the entry
            chain->GetEntry(iEntry);
            //---DAQ bug workaround
            if(iRun < 145) goodEvt = 10;
            else goodEvt = 1;
            if(evtNumber % goodEvt == 0)   
            {
                //---Read SciFront ADC value and set the e- multiplicity 
                //---(default = 1)
                trig = 1;
                for(int iCh=0; iCh<nAdcChannels; iCh++)
                {
                    if(adcData[iCh] > 1500 && adcBoard[iCh] == 1 && adcChannel[iCh] == 0) trig=2;
                    if(adcData[iCh] < 500 && adcBoard[iCh] == 1 && adcChannel[iCh] == 0) trig=0;
                }
                if(trig > 1) continue; 
                //---Read digitizer samples
                for(int iSample=0; iSample<nDigiSamples; iSample++)
                        digiCh[digiChannel[iSample]].push_back(digiSampleValue[iSample]);
                for(int iCh=0; iCh<6; iCh++)
                {
                    baseline[iCh]=SubtractBaseline(5, 25, &digiCh[iCh]);
                    if(iCh == 3)
                    {
                        for(int iSample=0; iSample<digiCh[iCh].size(); iSample++)
                            digiCh[iCh].at(iSample) = -digiCh[iCh].at(iSample);
                    }
                    timeCF[iCh]=TimeConstFrac(30, 500, &digiCh[iCh], 0.5);
                    int t1 = (int)(timeCF[iCh]/0.2) - 3;
                    int t2 = (int)(timeCF[iCh]/0.2) + 17;
                    //---Fill the signal integral histo only if the e- multiplicity is 1
                    if(t1 > 30 && t1 < 1024 && t2 > 30 && t2 < 1024 && trig==1)
                    {
                            chHistoSignal[iCh]->Fill(ComputeIntegral(t1, t2, &digiCh[iCh]));
                    }
                    chHistoBase[iCh]->Fill(ComputeIntegral(26, 46, &digiCh[iCh]));
                    timeDiffHisto[iCh]->Fill(timeCF[iCh]*0.2-timeCF[0]*0.2); 
                }
            }
        }  	
        chain->Delete();
      }
        chHistoBase[Ch_1]->Write();
        chHistoBase[Ch_2]->Write();
        chHistoBase[Ch_3]->Write();
   	    chHistoBase[Ch_ref1]->Write();
        chHistoBase[Ch_ref2]->Write();
   	    chHistoSignal[Ch_1]->Write();
        chHistoSignal[Ch_2]->Write();
        chHistoSignal[Ch_3]->Write();
   	    chHistoSignal[Ch_ref1]->Write();
        chHistoSignal[Ch_ref2]->Write();
        timeDiffHisto[Ch_1]->Write();
        timeDiffHisto[Ch_2]->Write();
        timeDiffHisto[Ch_3]->Write();
   	    timeDiffHisto[Ch_ref1]->Write();
        timeDiffHisto[Ch_ref2]->Write();
    //}

}

        
