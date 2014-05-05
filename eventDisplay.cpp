/*******************************************************************************
	this program is call by wavedump at every trigger and plot the waves from
	the active channels
	
	compile with --->  c++ -o eventDisplay eventDisplay.cpp `root-config --cflags --glibs`
*******************************************************************************/

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TString.h"
#include "TApplication.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

//******************************************************************************
// main

int main (int argc, char** argv)
{
//-----------------Definition---------------------------------------------------
	gStyle->SetOptStat("e");
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.3);
	TApplication *myapp = new TApplication("myapp",0,0);
	TCanvas* display = new TCanvas();
	TH1F* double_coinc = new TH1F("double coincidence","double coincidence",2500,0,50);
	TH1F* triple_coinc = new TH1F("triple coincidence","triple coincidence ",2500,0,50);
	TH1F* channels[9];
	int iEvent=0, nCh=0, nSample=0, run=1;
	float buffer=0, diff_2=0, diff_3=0;
	vector<float> channels_max;
	
//-----------------Initialization-----------------------------------------------
	nSample = atoi(argv[1]); 
	nCh = atoi(argv[2]);
	for(int iCh=0; iCh<=nCh; iCh++)
	{
		char tmp[30];
		TString name_tmp;
		if(iCh == 0) sprintf(tmp, "Trigger");
		else sprintf(tmp, "Channel_%d", iCh);
		name_tmp = tmp;
		channels[iCh] = new TH1F(name_tmp.Data(), "MCP display", nSample, 0, 0.2*nSample);
		channels[iCh]->SetFillColor(0);
		channels[iCh]->SetLineColor(iCh+2);
		channels[iCh]->GetXaxis()->SetTitle("t [ns]");
		channels[iCh]->GetYaxis()->SetTitle("Amp ADC");
		channels[iCh]->GetYaxis()->SetRangeUser(0,4096);
		channels[iCh]->SetStats(0);
	} 
	
	display->Divide(2,1);
	display->cd(2)->Divide(1,2);
	
	cin >> run;
//-----------------Trigger loop-------------------------------------------------
	while(run == 1)
	{
		ifstream input ("event_tmp.txt", ios::in);
		input >> iEvent >> nSample >> nCh;
 
		for(int iSample=1; iSample<=nSample; iSample++)
		{
			for(int iCh=0; iCh<=nCh; iCh++)
			{
				input >> buffer;
				channels[iCh]->SetBinContent(iSample, buffer);
			}
		}
		//---Display---	
		display->cd(1);
		channels[0]->Draw();
		for(int iCh=1; iCh<=nCh; iCh++)
		{
			channels_max.push_back(channels[iCh]->GetMaximumBin()*0.2);
			channels[iCh]->Draw("same");
		}
		diff_2 = 12;//channels_max.at(1) - channels_max.at(0);
		diff_3 = 12;//(channels_max.at(1) + channels_max.at(0))/2 - channels_max.at(2);
		//---double coinc
		display->cd(2)->cd(1);
		double_coinc->Fill(diff_2);
		double_coinc->GetXaxis()->SetRange(double_coinc->FindFirstBinAbove(0)-250, 
							   			   double_coinc->FindLastBinAbove(0)+250);
		double_coinc->Draw();
		//---triple coinc
		display->cd(2)->cd(2);
		triple_coinc->Fill(diff_3);
		triple_coinc->GetXaxis()->SetRange(triple_coinc->FindFirstBinAbove(0)-250, 
							   			   triple_coinc->FindLastBinAbove(0)+250);
		triple_coinc->Draw();
		//---interactive canvas
		gPad->Modified();
	    gPad->Update();
		gSystem->ProcessEvents();
		//---wait for the next trigger
		input.close();
		cin >> run;
	}
//-----------------Quit---------------------------------------------------------
	return 0;
}
					
		
	
