//c++ plottersimple_new plottersimple_new.cpp `root-config --cflags --glibs` -o  

/*
	Ricordati che la media sara' circa 10ns, cosi' recuperi la scala
sembra buono:	emailData147_140327184601

*/
#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TApplication.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

TGraphErrors *eff_hv = new TGraphErrors(); 
TGraphErrors *sigma_hv = new TGraphErrors();

TH2F *dtvsamp = new TH2F("dtvsamp","",1000,-0.1,0.,1200,-20.,20.);
TGraph *chi = new TGraph();
TGraph *sigma = new TGraph();
std::vector<float> time0;
std::vector<float> time1;
std::vector<float> time2;

TProfile *p_dtvsamp = new TProfile("dtvsamp","",100,-0.1,0.);//,1200,-20.,20.);

int frac = 40;

//costant fraction function
float trigger (float x1, float y1, float x2, float y2, float x3, float y3, float frazione, float amp, int step, int graf)
{
float ritorno;
float denominatore = (3*(x1*x1 + x2*x2 + x3*x3) - (x1+x2+x3)*(x1+x2+x3));
float m = (3*(x1*y1 + x2*y2 + x3*y3) - (x1+x2+x3)*(y1+y2+y3))/denominatore;
float q = ((y1+y2+y3)*(x1*x1 + x2*x2 + x3*x3) - (x1+x2+x3)*(x1*y1 + x2*y2 + x3*y3))/denominatore;
if (denominatore == 0.) 
{
ritorno = x2;
cout<<"interpolation failed, return tStamp with non changes: step "<<step<<", histo: "<<graf<<", tStamp = "<<x2<<endl;
}
else
ritorno = (amp*frazione/100. - q)/m; 
return ritorno;
}

//costant fraction slope function
float deriv (float x1, float y1, float x2, float y2, float x3, float y3, float frazione, float amp, int step, int graf)
{
float denominatore = (3*(x1*x1 + x2*x2 + x3*x3) - (x1+x2+x3)*(x1+x2+x3));
float m = (3*(x1*y1 + x2*y2 + x3*y3) - (x1+x2+x3)*(y1+y2+y3))/denominatore;
if(denominatore == 0)
cout<<"NAN!!!!!!!!!!!!!!!!!  ==> "<<endl;
return m;
}

//main function
void plottersimple(TString runFolder, float HV, int point)
{
std::cout<<frac<<std::endl;
// General reset
gROOT->Reset();

// Define histograms
TH1F* histo[3];
TH1F *h_base[3];
TH1F *h_coeff[3];
char h1[10];
char hb[10];
char hc[10];
	
	for (int iw=0;iw<3;iw++)
	{
    	sprintf (h1,"Ch%d",iw);
	sprintf (hb,"h_base%d",iw);
	sprintf (hc,"h_coeff%d",iw);
	h_base[iw] = new TH1F(hb,hb,1000,-0.1,0.01);
 	h_coeff[iw] = new TH1F(hc,hc,1000,-0.1,0.01);
   	histo[iw] = new TH1F(h1,h1,2000,0.,100.);	//1000,0.,50. -> all but Run012, which uses 2000,0.,100.
    	histo[iw] -> SetXTitle ("Time (ns)");
    	histo[iw] -> SetYTitle ("Amplitude (V)");
    	histo[iw]->SetLineWidth(2);
    	histo[iw]->SetMaximum(0.15);
    	histo[iw]->SetMinimum(-0.15);
  	}

histo[0]->SetLineColor(kRed+1);   
histo[1]->SetLineColor(kBlue+1);  
histo[2]->SetLineColor(kGreen+1); 

// Tree branches and tree strcuture
float piedistallo[3] = {0.,0.,0.};
float piedrms[3] = {0.,0.,0.};
float coeffmedio[3] = {0.,0.,0.};
float ampMax[3];
float tStamp[3];
float baseline[3]; 
float coeff[3];

TTree *nt = new TTree("nt","nt");

nt->Branch("ampMax",&ampMax,"ampMax[3]/F"); 
nt->Branch("tStamp",&tStamp,"tStamp[3]/F"); 
nt->Branch("baseline",&baseline,"baseline[3]/F"); 
nt->Branch("coeff",&coeff,"coeff[3]/F"); 
nt->Branch("coeffmedio",&coeffmedio,"coeffmedio[3]/F"); 
nt->Branch("piedistallo",&piedistallo,"piedistallo[3]/F"); 
nt->Branch("piedrms",&piedrms,"piedrms[3]/F"); 

// define graphics
TCanvas *c = new TCanvas("c","c");
c ->cd(); 
c->SetGrid();
histo[0] -> Draw("");

TLegend* leg = new TLegend(0.88,0.65,0.98,0.85);
leg->SetFillColor(0);
leg->SetTextFont(41);

	for (int iw=0;iw<3;iw++) 
	leg->AddEntry(histo[iw],histo[iw]->GetTitle(),"l");
  
leg->Draw();

// dump the waveform file list into a temporary file
//  TString command = "ls ../WaveForms/" + runFolder + "/*Ch1* > input.tmp"; 
TString command = "ls " + runFolder + "/*Ch1* > input.tmp"; 
gSystem -> Exec(command); 

// loop over the file and read data
ifstream inFile("input.tmp");
std::string line;
int ifile=0;
int isave=0;
int iwin =0;

// define file for wf
TString fileWF = runFolder + ".root";// "_wf_save.root"; 
TFile *fw = new TFile(fileWF,"RECREATE");

TString waveFile[400];
	while(getline(inFile,line) && ifile<18000)		  //  evento 17 e' buono per Run4
	{
    	// define wafeforms input
    	waveFile[0] = (TString)line;
	
    	Ssiz_t s=waveFile[0].Sizeof();
    	waveFile[2] = (TString)waveFile[0].Replace(s-6,s,"3.txt"); // planacon
    	waveFile[1] = (TString)waveFile[0].Replace(s-6,s,"2.txt"); // B-MCP1 trigger
    	waveFile[0] = (TString)waveFile[0].Replace(s-6,s,"1.txt"); // B-MCP2 trigger
    	ifile++;
  
    	// read waveforms
		for (int iw=0;iw<3;iw++) 
		{
		// cout << "Reading... " << waveFile[iw] << endl;

      		float counts=0;
      		ampMax[iw] = 0;
	        tStamp[iw] = 0;

      		std::ifstream infile (waveFile[iw].Data(), std::ios::in);
      			for (int bin=1;bin<2001;bin++)		//1001 for all but run012 which uses 2001
      			{			
			infile >> counts;
			if (iw==2) counts = -counts; 	// planacon readout at the MCP output 
	                             			// plane and not at the anode 

				if (counts < ampMax[iw]) 
				{
				ampMax[iw] = counts;
				tStamp[iw] = bin; 
				}
			histo[iw] -> SetBinContent(bin, counts);
      			}
		//calculate baseline
		baseline[iw] = histo[iw] -> Integral(1,101)/100.;
		h_base[iw] -> Fill(baseline[iw]);
		
		//baseline correction
			for (int newbin = 1; newbin < 2001; newbin++)			//1001 for all but run012 which uses 2001
			histo[iw]->SetBinContent(newbin,histo[iw]->GetBinContent(newbin) - baseline[iw]);

		
		ampMax[iw] = ampMax[iw] - baseline[iw];

	        //histo[iw] -> Draw("same");

		int rif = 0;
			for (int scan = tStamp[iw]; scan > 0; --scan)
			{
			if(histo[iw]->GetBinContent(scan) > ampMax[iw]*frac/100) break;  
			tStamp[iw] = scan*0.05;
			rif = scan;
			}
		
		tStamp[iw] = trigger((rif-1)*0.05, histo[iw]->GetBinContent(rif-1),rif*0.05,histo[iw]->GetBinContent(rif), (rif+1)*0.05, histo[iw]->GetBinContent(rif+1), frac, ampMax[iw],ifile,iw);

			if (isnan(tStamp[iw]) == 1)
			cout<<"NAN!!!!!!!!!!!!!!!!!!!!!!!!   ==>   "<<ifile<<" , "<<iw<<endl;

		coeff[iw] = deriv((rif-1)*0.05, histo[iw]->GetBinContent(rif-1),rif*0.05,histo[iw]->GetBinContent(rif), (rif+1)*0.05, histo[iw]->GetBinContent(rif+1), frac, ampMax[iw],ifile,iw);
		h_coeff[iw]->Fill(coeff[iw]);
	
		}

		for (int iw=0;iw<3;iw++) 	
		{  
		sprintf (h1,"Ch%d_%03d",iw,isave);
	  	histo[iw]->SetName(h1);
	  	histo[iw]->SetTitle(h1);
	  	histo[iw]->Write();
		}
	
    	// fill ntuple
	nt->Fill();
 	}

//saving ntuple
for(int iw = 0; iw < 3; ++iw)
{
ampMax[iw] = 0.; 
tStamp[iw] = 0.; 
baseline[iw] = 0.; 
coeff[iw] = 0.; 
piedistallo[iw] = h_base[iw]->GetMean();
piedrms[iw] = h_base[iw]->GetRMS();
coeffmedio[iw] = h_coeff[iw]->GetMean();
}
nt->Fill();
nt->Write();
fw->Close();
}

//***************************************************************************************************************************************
//********************************************        READING NTUPLES      **************************************************************

int reading(TString runFolder, float HV, int point)
{
cout<<"num = "<<point<<endl;

TFile *Reading = new TFile(runFolder + ".root");
TTree *t1 = (TTree*)Reading->Get("nt");    	

float amplitude[3];
float time[3];
float offset[3];
float incl[3];
float inclmedia[3];
float base_mean[3];
float base_rms[3];

t1->SetBranchAddress("tStamp",&time);
t1->SetBranchAddress("ampMax",&amplitude);	
t1->SetBranchAddress("baseline",&offset);	
t1->SetBranchAddress("coeff",&incl);	
t1->SetBranchAddress("coeffmedio",&inclmedia);
t1->SetBranchAddress("piedistallo",&base_mean);
t1->SetBranchAddress("piedrms",&base_rms);

t1->GetEntry(t1->GetEntries()-1);
float cut0 = base_rms[0];
float cut1 = base_rms[1];
float cut2 = base_rms[2];

TH1F *h_diff = new TH1F ("h_diff","h_diff",1200,6.,18.);
TH1F *h_mean_mcp = new TH1F ("h_mean_mcp","h_mean_mcp",800,-20.,20.);
TH1F *h_1_mcp = new TH1F ("h_1_mcp","h_1_mcp",1200,-50.,50.);
TH1F *h_mcp_2 = new TH1F ("h_mcp_2","h_mcp_2",1200,-50.,50.);

TH1F *h_offset0 = new TH1F ("h_offset0","h_offset0",1000,-0.005,0.005);
TH1F *h_offset1 = new TH1F ("h_offset1","h_offset1",1000,-0.005,0.005);
TH1F *h_offset2 = new TH1F ("h_offset2","h_offset2",1000,-0.005,0.005);

TH1F *h_incl0 = new TH1F ("h_incl0","h_incl0",1000,-0.1,0.01);
TH1F *h_incl1 = new TH1F ("h_incl1","h_incl1",1000,-0.1,0.01);
TH1F *h_incl2 = new TH1F ("h_incl2","h_incl2",1000,-0.1,0.01);

TF1 *gauss = new TF1("gauss","[3]+[1]/sqrt(2*TMath::Pi()*[0]*[0])*exp(-(x-[2])*(x-[2])/(2*[0]*[0]))",-50,50);
gauss->SetParName(0,"sigma");
gauss->SetParName(1,"amplitude");
gauss->SetParName(2,"mean");
gauss->SetParName(3,"trasl");
gauss->SetNpx(1000000);

	for (int o = 0 ; o < t1->GetEntries() - 1; ++o)
	{
	t1->GetEntry(o);
        h_diff->Fill(time[1]-time[0]);
	h_mean_mcp->Fill((time[1]+time[0])*0.5-time[2]);
        h_1_mcp->Fill(time[1]-time[2]);
        h_mcp_2->Fill(time[2]-time[0]);

	h_offset0 -> Fill (offset[0]);
	h_offset1 -> Fill (offset[1]);
	h_offset2 -> Fill (offset[2]);
	
	h_incl0 -> Fill (incl[0]);
	h_incl1 -> Fill (incl[1]);
	h_incl2 -> Fill (incl[2]);
	}
	
//*************************************************************************************************************************************	
//***************************************************     DRAW CONCIDENCES     **********************************************************
  //0,1 -> external detectors
  //2 -> mcp detector
  
TCanvas *scatter = new TCanvas();
scatter -> cd();
dtvsamp->Draw("colz");
scatter->Print("scatter.root","root");

TCanvas *coincidenze = new TCanvas();
coincidenze -> Divide (4,1);

//h_diff => difference
float media = 0.;
float contenuto = 0.;
	for (int in = 300; in < 600; ++in)
		if (h_diff->GetBinContent(in) > contenuto)		
  		{
  		contenuto = h_diff->GetBinContent(in);
  		media = in*0.01 + 6.;
  		}
  
cout<<"il centro della gaussiana sta in "<<media<<endl;
 
h_diff->GetXaxis()->SetRangeUser(9.,12.);
h_diff->SetTitle("difference");	
coincidenze -> cd(1);
gauss->SetParameters(0.05,1.,media,1.);	
h_diff->Fit("gauss","","",9.,12.);		
h_diff->Draw(),
std::cout<<"sigma difference = "<<gauss->GetParameter(0)<<std::endl;
float sigma_diff = fabs(gauss->GetParameter(0));
float mean_diff = fabs(gauss->GetParameter(2));


//h_1_mcp => 1 - mcp
h_1_mcp->GetXaxis()->SetRangeUser(6.5,9.5);
h_1_mcp->SetTitle("1 - mcp");
coincidenze -> cd(3);
gauss->SetParameters(0.1,5.,7.6,1.);
h_1_mcp->Fit("gauss","","",-25.,25.);
h_1_mcp->Draw();
std::cout<<"sigma 1 - mcp = "<<gauss->GetParameter(0)<<std::endl;

//h_mcp_2 => mcp - 2
h_mcp_2->GetXaxis()->SetRangeUser(1.,5.);
h_mcp_2->SetTitle("mcp - 2");
coincidenze -> cd(4);
gauss->SetParameters(0.15,10.,3.,3.);	//Run 006
h_mcp_2->Fit("gauss","","",1.,5.);
h_mcp_2->Draw();
std::cout<<"sigma mcp - 2 = "<<gauss->GetParameter(0)<<std::endl;

//h_mean_mcp => mean - mcp
float mediamcp = 0.;
float contenutomcp = 0.;
	for (int inn = 430; inn < 480; ++inn)
		if (h_mean_mcp->GetBinContent(inn) > contenutomcp)		
		{
		contenutomcp = h_mean_mcp->GetBinContent(inn);
		mediamcp = inn*0.05 - 20.;
		}
  
cout<<"il centro della gaussiana mcp sta in "<<mediamcp<<endl;
  
h_mean_mcp->GetXaxis()->SetRangeUser(1.,4.);
h_mean_mcp->SetTitle("mean - mcp");
coincidenze -> cd(2);
gauss->SetParameters(0.05,10.,mediamcp,3.);	
h_mean_mcp->Fit("gauss","","",1.,4.);
h_mean_mcp->Draw();
std::cout<<"sigma mean - mcp = "<<gauss->GetParameter(0)<<std::endl;
float sigma_mcp = gauss->GetParameter(0);
float mean_mcp = gauss->GetParameter(2);

sigma_hv->SetPoint(point, HV,1000.*sqrt(sigma_mcp*sigma_mcp-sigma_diff*sigma_diff));

coincidenze->Print("coinc.root","root");
  	
chi->SetPoint(frac,frac,gauss->GetChisquare()/gauss->GetNDF());
sigma->SetPoint(frac,frac,fabs(gauss->GetParameter(0))); 

//***************************************************************************************************************************************
//*************************************************         EVALUATING EFFICIENCY        ************************************************

int bene = 0;
int benissimo = 0;
int fake = 0;
int trfake = 0;
int drake = 0;
	
TH1F *doppie = new TH1F ("doppie","doppie",1000,2.,12.);
TH1F *triple = new TH1F ("triple","triple",1000,2.,12.);
doppie -> SetFillColor (kGreen);
triple -> SetFillColor (kRed);

cout<<"cut0 vale "<<cut0<<endl;  
cout<<"cut1 vale "<<cut1<<endl;
cout<<"cut2 vale "<<cut2<<endl;

	for (int r = 0 ; r < t1->GetEntries() ; ++r)
	{
	t1->GetEntry(r);
//	dtvsamp->Fill(amplitude[1]+amplitude[0],time[1]-time[0]);	
	//double coincidence
	if (fabs((time[1]-time[0]) - mean_diff) < 3*sigma_diff && amplitude[1] < -5*cut1 && amplitude[0] < -5*cut0) 
		{					
		doppie -> Fill(time[1]-time[0]);		
		bene++;
//		p_dtvsamp->Fill(amplitude[2],(time[1]+time[0])*0.5-time[2]);	
		dtvsamp->Fill(amplitude[2],(time[1]+time[0])*0.5-time[2]);
	    	if (fabs((time[1]+time[0])*0.5-time[2] - 5.) < 1.5 && amplitude[2] < - fabs(5*cut2)) 	trfake++;	 
		if (fabs((time[1]+time[0])*0.5-time[2] - mean_mcp) < 3*sigma_mcp && amplitude[2] < -5*cut2)
			{
			benissimo++;
			p_dtvsamp->Fill(amplitude[2],(time[1]+time[0])*0.5-time[2]);	
			triple -> Fill((time[1]+time[0])*0.5-time[2]);		
			}
	 	}
	 if (fabs(time[1]-time[0] - 13.5) < 1.5 && amplitude[1] < -fabs(5*cut1) && amplitude[0] < -fabs(5*cut0))		fake++;	
         }

float rfake = (float)fake*3.*sigma_diff/1.5;
float rtrfake = (float)trfake*3.*sigma_mcp/1.5;

TCanvas *p_scatter = new TCanvas();
p_scatter -> cd();
TF1 *retta = new TF1 ("retta","pol1",-0.1,0.);
p_dtvsamp->Fit("retta","","",-0.02,-5*cut2);		//Run006
p_dtvsamp->Draw("colz");
cout<<"correggerei il problema dell'amplitude walk usando la retta y = "<<retta->GetParameter(1)<<" x + "<<retta->GetParameter(0)<<endl;
p_scatter->Print("scatter.root","root");

cout << "saved / triggered / double fake / triple fake: " << benissimo<< " / " << bene << " / " << rfake << " / "<< trfake<< endl;
cout << "fake = " <<fake<<endl;
float efficiency = (((float)benissimo-rtrfake)/((float)bene-rfake));
float error = sqrt (efficiency * (1. - efficiency) / ((float)bene-rfake));
cout.setf(std::ios::fixed);
cout.precision(5);
cout << "efficiency is (" <<efficiency*100.<< " +- " <<error*100.<<")\%" << endl;

eff_hv->SetPoint(point, HV, efficiency);
eff_hv->SetPointError(point, 5., error);

TCanvas *cfr = new TCanvas();
cfr->SetTitle("Doppie & Trilpe");
cfr->cd();
doppie -> Draw();
triple -> Draw("same");
cfr->Print("cfr.root","root");

Reading->Close();
//***************************************************************************************************************************************
//************************************************     DRAW SIGMA - CHI    **************************************************************
		
/*

  chi->SetMarkerStyle(20);
  sigma->SetMarkerStyle(20);
  chi->SetMarkerColor(kRed+1);
  sigma->SetMarkerColor(kGreen+1);
  TCanvas *tavola = new TCanvas();
  tavola->Divide(2,1);
  tavola->cd(1);	
  chi->Draw("AP");
  tavola->cd(2);	
  sigma->Draw("AP");
  TF1 *para = new TF1 ("para","pol2",0,60);	
  sigma->Fit("para","","",7.,55.);
  float minimum = para->GetMinimumX();
  cout<<"la frazione costante dell'ampiezza che sceglierei è "<<(int)minimum<<endl;
*/

//***************************************************************************************************************************************
//**********************************************      DRAW WAVEFORM	*****************************************************************
/*
	TCanvas *baseline = new TCanvas();
	baseline->Divide(3,1);

	TString printname;
	TFile *MyFile = new TFile("WFRun006_wf_save.root","READ");
	
	baseline->cd(1);	
	for (int i = 0 ; i < 85 ; ++i)
	{
		char count[10];
		sprintf(count, "%03d", i);
		TString help;
		help += count;
		TH1F* h_tmp = (TH1F*)MyFile->Get("Ch0_"+help);
		h_tmp->Draw("same");
	}
*/
//*********************************************************************************************************************			
return 0;
}

int main(int argc, char** argv)
{
TApplication *myApp = new TApplication("example",&argc, argv);
TString run;	
float voltage;		
int num;
/*
run = "WFRun005";	voltage = 2150;		num = 5;
//plottersimple(run, voltage, num);
reading(run, voltage, num);

run = "WFRun006";	voltage = 2251;		num = 6;
plottersimple(run, voltage, num);
reading(run, voltage, num);

run = "WFRun007";	voltage = 2050;		num = 7;
//plottersimple(run, voltage, num);
reading(run, voltage, num);

run = "WFRun008";	voltage = 1950;		num = 8;
//plottersimple(run, voltage, num);
reading(run, voltage, num);

run = "WFRun009";	voltage = 2350;		num = 9;
plottersimple(run, voltage, num);
reading(run, voltage, num);

run = "WFRun010";	voltage = 2250;		num = 10;
//plottersimple(run, voltage, num);
reading(run, voltage, num);

run = "WFRun012";	voltage = 1750;		num = 12;
//plottersimple(run, voltage, num);
reading(run, voltage, num);

run = "WFRun013";	voltage = 1750;		num = 13;
plottersimple(run, voltage, num);
reading(run, voltage, num);
*/

run = "WFRun015";	voltage = 1750;		num = 15;
plottersimple(run, voltage, num);
reading(run, voltage, num);

TCanvas *finale = new TCanvas();
finale->Divide(2,1);

finale-> cd(1);
eff_hv->GetXaxis()->SetTitle("High Voltage [V]");
eff_hv->GetYaxis()->SetTitle("Efficiency");
eff_hv->GetXaxis()->SetLimits(0.,3000.);
eff_hv->Draw("AP");

finale-> cd(2);
sigma_hv->SetMarkerStyle(34);
sigma_hv->SetMarkerSize(1);
sigma_hv->GetXaxis()->SetTitle("High Voltage [V]");
sigma_hv->GetYaxis()->SetTitle("Time Resolution [ps]");
sigma_hv->GetXaxis()->SetLimits(0.,3000.);
sigma_hv->Draw("AP");

finale->Print("eff12.root","root");

myApp->Run();
}
