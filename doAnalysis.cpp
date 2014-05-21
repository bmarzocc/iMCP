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

#include "analysis_tools.h"
#include "InitTree_BTF.h"

int main (int argc, char** argv)
{
    gROOT->ProcessLine("#include <vector>");
    TChain* chain = new TChain("eventRawData");
    for(int i=1; i<argc; i++)
    {
        chain->Add(argv[1]);
    }
    TH1F* base = new TH1F("base","base",2000,0,2000);
    TH1F* pulse = new TH1F("pulse","pulse",2000,0,2000);
    pulse->SetLineColor(kRed);
    //TApplication* app = new TApplication("app",0,0);
    
    InitTree(chain);
    
    vector<float> digiCh[9];
    float timeCF[9];
    float baseline[9];
    int count[5]={0,0,0,0,0}, spare[5]={0,0,0,0,0}, spare2[5]={0,0,0,0,0};
    int tot_tr1=0, tot_tr0, trig=0;
    
    for(int iEntry=0; iEntry<chain->GetEntries(); iEntry++)
    {
        for(int iCh=0; iCh<9; iCh++)
        {
            digiCh[iCh].clear();
        }
        
        if(iEntry % 1 != 0) 
        {
            chain->GetEntry(iEntry);
        }
        else 
        {
            trig = 1;
            chain->GetEntry(iEntry);
            for(int iCh=0; iCh<nAdcChannels; iCh++)
            {
                if(adcData[iCh] > 1500 && adcBoard[iCh] == 1 && adcChannel[iCh] == 0) trig=2;
                if(adcData[iCh] < 500 && adcBoard[iCh] == 1 && adcChannel[iCh] == 0) trig=0;
            }
            if(trig > 1) continue; 
            for(int iCh=0; iCh<7; iCh++)
            {
                digiCh[iCh].clear();
            }
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
                baseline[iCh]=SubtractBaseline(5, 80, &digiCh[iCh]);
                timeCF[iCh]=TimeConstFrac(200, 275, &digiCh[iCh], 0.5);
            }    
            if(AmpMax(200, 275, &digiCh[0]) < -20 && AmpMax(200, 275, &digiCh[5]) < -20 && trig==1) 
            {
                tot_tr1++;
                if(AmpMax(200, 275, &digiCh[1]) < -20) count[1]=count[1]+1;
                if(baseline[1] < -20) spare[1]=spare[1]+1;
                if(AmpMax(200, 275, &digiCh[4]) < -20) count[2]=count[2]+1;
                if(baseline[4] < -20) spare[2]=spare[2]+1;
                if(AmpMax(200, 275, &digiCh[3]) < -20) count[3]=count[3]+1;
                if(baseline[3] < -20) spare[3]=spare[3]+1;
            }
            else if(AmpMax(200, 275, &digiCh[0]) >= -20 && AmpMax(200, 275, &digiCh[5]) >= -20 && trig==0)
            {
                tot_tr0++;
                if(AmpMax(200, 275, &digiCh[1]) < -20) spare2[1]=spare2[1]+1; 
                if(AmpMax(200, 275, &digiCh[4]) < -20) spare2[2]=spare2[2]+1; 
                if(AmpMax(200, 275, &digiCh[3]) < -20) spare2[3]=spare2[3]+1;
            }
        }
    }

    cout << "--------------------------" << endl;
    cout << "number of events:  " << chain->GetEntries() << endl;
    cout << "Double:  " << tot_tr1 << endl;
    cout << "No e- :  " << tot_tr0 << endl;
    cout << "--------------------------" << endl;
    cout << "MiB2:  " << count[1] << "  " << spare[1] << "  " << spare2[1] << endl;
    cout << "MiB3:  " << count[2] << "  " << spare[2] << "  " << spare2[2] << endl;
    cout << "Planacon:  " << count[3] << "  " << spare[3] << "  " << spare2[3] << endl;
    cout << "--------------------------" << endl;
    float eff1 = ((float)count[1]-(float)spare[1])/(float)tot_tr1;
    float eff2 = ((float)count[2]-(float)spare[2])/(float)tot_tr1;
    float eff3 = ((float)count[3]-(float)spare[3])/(float)tot_tr1;
    cout << "MiB2 eff:        " << eff1 << endl;
    cout << "MiB2 e_err:      " << TMath::Sqrt((eff1*(1-eff1))/tot_tr1) << endl;
    cout << "MiB3 eff:        " << eff2 << endl;
    cout << "MiB3 e_err:      " << TMath::Sqrt((eff2*(1-eff2))/tot_tr1) << endl;
    cout << "Planacon eff:    " << eff3 << endl;
    cout << "Planacon e_err:  " << TMath::Sqrt((eff3*(1-eff3))/tot_tr1) << endl;
    cout << "--------------------------" << endl;
    
//altro
        /*
        TH1F* MiB1 = new TH1F("MiB1","MiB1",1024,0,1023*0.2);
        MiB1->SetLineColor(kRed);
        TH1F* MiB2 = new TH1F("MiB2","MiB2",1024,0,1023*0.2);
        MiB2->SetLineColor(kRed+2);
        TH1F* MiB3 = new TH1F("MiB3","MiB3",1024,0,1023*0.2);
        MiB3->SetLineColor(kRed+4);
        TH1F* Planacon = new TH1F("Planacon","Planacon",1024,0,1023*0.2);
        Planacon->SetLineColor(kBlue);
        TH1F* Roma = new TH1F("Roma","Roma",1024,0,1023*0.2);
        Roma->SetLineColor(kGreen);
        for(int bin=0; bin<1024; bin++)
        {
            MiB1->SetBinContent(bin+1,digiCh[0].at(bin));
            MiB2->SetBinContent(bin+1,digiCh[1].at(bin));
            Planacon->SetBinContent(bin+1,-digiCh[3].at(bin));
            MiB3->SetBinContent(bin+1,digiCh[4].at(bin));
            Roma->SetBinContent(bin+1,digiCh[5].at(bin));
        }        
        MiB1->Draw();
        MiB2->Draw("same");
        MiB3->Draw("same");
        Planacon->Draw("same");
        Roma->Draw("same");   
        app->Run();
        int check;
        cin >> check;*/


}

        
