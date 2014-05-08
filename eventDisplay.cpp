/*******************************************************************************
this program is call by wavedump at every trigger and plot the waves from
the active channels
compile with ---> c++ -o eventDisplay eventDisplay.cpp `root-config --cflags --glibs`
*******************************************************************************/

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h> 
#include <stdio.h> 
#include <unistd.h>
#include <cmath> 

using namespace std;

//******************************************************************************
// main

int main (int argc, char** argv)
{
//-----------------Definition---------------------------------------------------
    gROOT->Reset();
    gStyle->SetOptStat("e");
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.3);
    TApplication *myapp = new TApplication("myapp",0,0);
    TLegend *legend = new TLegend(0.85,0.8,0.99,0.95);
    legend->SetFillColor(kWhite);
    TCanvas* display = new TCanvas();
    TH1F* double_coinc = new TH1F("double coincidence","double coincidence",500,-50,50);
    TH1F* double_coinc_2 = new TH1F("selected double coincidence","selected double coincidence",500,-50,50);
    double_coinc_2->SetFillColor(kGreen);
    TH1F* triple_coinc = new TH1F("triple coincidence","triple coincidence ",500,-50,50);
    TH1F* channels[9];
    int iEvent=0, nCh=0, nSample=0, run=1;
    float buffer=0, diff_2=0, diff_3=0;
    vector<float> channels_max;
    
//-----------------Initialization-----------------------------------------------
    cin >> nCh >> nSample;
    for(int iCh=0; iCh<=nCh; iCh++)
    {
        int color = 2+iCh;
        if(iCh == 3) color++; //no yellow line please!
        char tmp_h[30];
        TString name_tmp;
        if(iCh == 0) sprintf(tmp_h, "Trigger");
        else sprintf(tmp_h, "Channel_%d", iCh);
        name_tmp = tmp_h;
        channels[iCh] = new TH1F(name_tmp.Data(), "MCP display", nSample, 0, 0.2*nSample);
        channels[iCh]->SetFillColor(0);
        channels[iCh]->SetLineColor(color);
        channels[iCh]->GetXaxis()->SetTitle("t [ns]");
        channels[iCh]->GetYaxis()->SetTitle("Amp ADC");
        channels[iCh]->GetYaxis()->SetRangeUser(0,4096);
        channels[iCh]->SetStats(0);
        char tmp_l[3];
        if(iCh == 0) sprintf(tmp_l, "Tr0");
        else sprintf(tmp_l, "Ch%d", iCh-1);
        legend->AddEntry(channels[iCh],tmp_l,"L");
            
    }

    display->Divide(2,1);
    display->cd(2)->Divide(1,2);

    cin >> run;
//-----------------Trigger loop-------------------------------------------------
    while(run == 1)
    {
        for(int iCh=0; iCh<=nCh; iCh++)
        {
            for(int iSample=1; iSample<=nSample; iSample++)
            {
                cin >> buffer;
                channels[iCh]->SetBinContent(iSample, buffer);
            }
        }
        //---Display---
        display->cd(1);
        channels[0]->Draw();
        for(int iCh=1; iCh<=nCh; iCh++)
        {
            channels[iCh]->Draw("same");
            if(iCh == 3 && channels[iCh]->GetBinContent(channels[iCh]->GetMaximumBin()) < 1250) continue;
            channels_max.push_back(channels[iCh]->GetMinimumBin()*0.2);
        }
        //---Draw Legend---
        legend->Draw("same");
        //---double coinc
        display->cd(2)->cd(1);
        diff_2 = channels_max.at(1) - channels_max.at(0);
        double_coinc->Fill(diff_2);
        double_coinc->GetXaxis()->SetRange(double_coinc->FindFirstBinAbove(0)-5, 
                                           double_coinc->FindLastBinAbove(0)+5);
        double_coinc->Draw();
        //---triple coinc
        if( channels_max.size() > 2) 
        {   
            diff_3 = (channels_max.at(1) + channels_max.at(0))/2 - channels_max.at(2);
            display->cd(2)->cd(2);
            triple_coinc->Fill(diff_3);
            triple_coinc->GetXaxis()->SetRange(triple_coinc->FindFirstBinAbove(0)-5,
                                               triple_coinc->FindLastBinAbove(0)+5);
            triple_coinc->Draw();
            display->cd(2)->cd(1);
            double_coinc_2->Fill(diff_2);
            double_coinc_2->Draw("same");
        }
        //---interactive canvas
        display->Modified();
        display->Update();
        gSystem->ProcessEvents();
        //---wait for the next trigger
        channels_max.clear();
        cin >> run;
    }
//-----------------Quit---------------------------------------------------------
    return 0;
}
