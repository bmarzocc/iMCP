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
    int nFiles=1;
    
    ifstream log (argv[1], ios::in);
    while(log >> nFiles)
    {
        vector<float> digiCh[9];
        float timeCF[9];
        float baseline[9];
        int count[5]={0,0,0,0,0}, spare[5]={0,0,0,0,0}, spare2[5]={0,0,0,0,0};
        int tot_tr1=0, tot_tr0=0, trig=0;
        int HV1=0, HV2=0, HV3=0;
        TChain* chain = new TChain("eventRawData");
        InitTree(chain);
        for(int iFiles=0; iFiles<nFiles; iFiles++)
        {
            int id;
            log >> id;
            char id_str[30];
            sprintf(id_str, "WaveForms_BTF/run_IMCP_%d_*.root", id);
            chain->Add(id_str);
            cout << "Reading:  WaveForms_BTF/run_IMCP_" << id << endl;
        }
        //log >> HV1 >> HV2 >> HV3; //---config1
        log >> HV3 >> HV1 >> HV2; //---config2
        for(int iEntry=0; iEntry<chain->GetEntries(); iEntry++)
        {
            for(int iCh=0; iCh<9; iCh++)
            {
                digiCh[iCh].clear();
            }
            
            chain->GetEntry(iEntry);
            if(evtNumber % 10 == 0)   //---Run<145
            //if(evtNumber % 1 == 0)      //---Run>=145
            {
                trig = 1;
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
                if(AmpMax(200, 275, &digiCh[0]) < -20 && AmpMax(200, 275, &digiCh[4]) < -20 && trig==1) 
                {
                    tot_tr1++;
                    if(AmpMax(200, 275, &digiCh[1]) < -20) count[1]=count[1]+1;
                    if(baseline[1] < -20) spare[1]=spare[1]+1;
                    if(AmpMax(200, 275, &digiCh[5]) < -20) count[2]=count[2]+1;
                    if(baseline[5] < -20) spare[2]=spare[2]+1;
                    if(AmpMax(200, 275, &digiCh[3]) < -10) count[3]=count[3]+1;
                    if(baseline[3] < -10) spare[3]=spare[3]+1;
                }
                else if(AmpMax(200, 275, &digiCh[0]) >= -20 && AmpMax(200, 275, &digiCh[4]) >= -20 && trig==0)
                {
                    tot_tr0++;
                    if(AmpMax(200, 275, &digiCh[1]) < -20) spare2[1]=spare2[1]+1; 
                    if(AmpMax(200, 275, &digiCh[5]) < -20) spare2[2]=spare2[2]+1; 
                    if(AmpMax(200, 275, &digiCh[3]) < -10) spare2[3]=spare2[3]+1;
                }
            }
        }
    
        cout << "--------------------------" << endl;
        cout << "number of events:  " << chain->GetEntries()/10 << endl;
        cout << "Double:  " << tot_tr1 << endl;
        cout << "No e- :  " << tot_tr0 << endl;
        cout << "--------------------------" << endl;
        cout << "MiB2:  " << count[1] << "  " << spare[1] << "  " << spare2[1] << endl;
        cout << "MiB3:  " << count[2] << "  " << spare[2] << "  " << spare2[2] << endl;
        cout << "Planacon:  " << count[3] << "  " << spare[3] << "  " << spare2[3] << endl;
        cout << "--------------------------" << endl;
        double eff1 = ((double)count[1]-(double)spare[1])/(double)tot_tr1;
        double eff2 = ((double)count[2]-(double)spare[2])/(double)tot_tr1;
        double eff3 = ((double)count[3]-(double)spare[3])/(double)tot_tr1;
        cout << "MiB2 eff:        " << eff1 << endl;
        cout << "MiB2 e_err:      " << TMath::Sqrt((eff1*(1-eff1))/tot_tr1) << endl;
        cout << "MiB3 eff:        " << eff2 << endl;
        cout << "MiB3 e_err:      " << TMath::Sqrt((eff2*(1-eff2))/tot_tr1) << endl;
        cout << "Planacon eff:    " << eff3 << endl;
        cout << "Planacon e_err:  " << TMath::Sqrt((eff3*(1-eff3))/tot_tr1) << endl;
        cout << "--------------------------" << endl;
    
        //---config 1
/*        ofstream MiB2_data ("Data_plateau/plateau_MiB2_pc_off.dat", ios::app);
        MiB2_data << HV1 << " " <<  eff1 << " " << 0 << " " << TMath::Sqrt((eff1*(1-eff1))/tot_tr1) << endl;
        ofstream MiB3_data ("Data_plateau/plateau_MiB3_pc_on.dat", ios::app);
        MiB3_data << HV2 << " " <<  eff2 << " " << 0 << " " << TMath::Sqrt((eff2*(1-eff2))/tot_tr1) << endl;
        ofstream Planacon_data ("Data_plateau/plateau_Planacon_pc_on.dat", ios::app);
        Planacon_data << HV3 << " " << eff3 << " " << 0 << " " << TMath::Sqrt((eff3*(1-eff3))/tot_tr1) << endl;
*/
        //---config 2
        ofstream MiB2_data ("Data_plateau/plateau_Rm1_pc_off.dat", ios::app);
        MiB2_data << HV1 << " " <<  eff1 << " " << 0 << " " << TMath::Sqrt((eff1*(1-eff1))/tot_tr1) << endl;
        ofstream MiB3_data ("Data_plateau/plateau_Rm2_pc_on.dat", ios::app);
        MiB3_data << HV2 << " " <<  eff2 << " " << 0 << " " << TMath::Sqrt((eff2*(1-eff2))/tot_tr1) << endl;
        ofstream Planacon_data ("Data_plateau/plateau_Planacon_pc_off_2.dat", ios::app);
        Planacon_data << HV3 << " " << eff3 << " " << 0 << " " << TMath::Sqrt((eff3*(1-eff3))/tot_tr1) << endl;    

        MiB2_data.close();
        MiB3_data.close();
        Planacon_data.close();
        
        chain->Delete();
        //delete chain;
    }
}

        
