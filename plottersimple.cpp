//c++ `root-config --cflags --glibs` -o plottersimple plottersimple.cpp 

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
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

//TH1F h_baseline0 = new TH1F("h_baseline0","h_baseline0",);
TH2F *dtvsamp = new TH2F("dtvsamp","",1000,-0.1,0.,1200,-20.,20.);
TGraph *chi = new TGraph();
TGraph *sigma = new TGraph();
std::vector<float> time0;
std::vector<float> time1;
std::vector<float> time2;

TProfile *p_dtvsamp = new TProfile("dtvsamp","",100,-0.1,0.);//,1200,-20.,20.);

int frac = 31;

float trigger (float x1, float y1, float x2, float y2, float x3, float y3, float frazione, float amp)
{
float ritorno;
float denominatore = (3*(x1*x1 + x2*x2 + x3*x3) - (x1+x2+x3)*(x1+x2+x3));
float m = (3*(x1*y1 + x2*y2 + x3*y3) - (x1+x2+x3)*(y1+y2+y3))/denominatore;
float q = ((y1+y2+y3)*(x1*x1 + x2*x2 + x3*x3) - (x1+x2+x3)*(x1*y1 + x2*y2 + x3*y3))/denominatore;
if (denominatore == 0.);
ritorno = x2;
ritorno = (amp*frazione/100. - q)/m; 
return ritorno;
}

void plottersimple(TString runFolder)
{
//for(int frac = 36 ; frac < 37 ; ++frac)
//{
	std::cout<<frac<<std::endl;
  // General reset
  gROOT->Reset();

  // Define histograms
  TH1F* histo[3];
  char h1[10];
    for (int iw=0;iw<3;iw++){
    sprintf (h1,"Ch%d",iw);
    histo[iw] = new TH1F(h1,h1,1000,0.,50.);
    histo[iw] -> SetXTitle ("Time (ns)");
    histo[iw] -> SetYTitle ("Amplitude (V)");
    histo[iw]->SetLineWidth(2);
    histo[iw]->SetMaximum(0.05);
    histo[iw]->SetMinimum(-0.05);
  }
  histo[0]->SetLineColor(kRed+1);   
  histo[1]->SetLineColor(kBlue+1);  
  histo[2]->SetLineColor(kGreen+1); 

  // Tree branches and tree strcuture
  float ampMax[3];
  float tStamp[3];
  float baseline[3]; 
  TTree *nt = new TTree("nt","nt");

  nt->Branch("ampMax",&ampMax,"ampMax[3]/F"); 
  nt->Branch("tStamp",&tStamp,"tStamp[3]/F"); 
  nt->Branch("baseline",&baseline,"baseline[3]/F"); 

  // define graphics
  TCanvas *c = new TCanvas("c","c");
  c ->cd(); 
  c->SetGrid();
  histo[0] -> Draw("");

  TLegend* leg = new TLegend(0.88,0.65,0.98,0.85);
  leg->SetFillColor(0);
  leg->SetTextFont(41);
  for (int iw=0;iw<3;iw++) {
    leg->AddEntry(histo[iw],histo[iw]->GetTitle(),"l");
  }
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
    for (int iw=0;iw<3;iw++) {
      // cout << "Reading... " << waveFile[iw] << endl;

      float counts=0;
      ampMax[iw] = 0;
      tStamp[iw] = 0;

      std::ifstream infile (waveFile[iw].Data(), std::ios::in);
      for (int bin=1;bin<1001;bin++){
	infile >> counts;
	if (iw==2) counts = -counts; // planacon readout at the MCP output 
	                             // plane and not at the anode 

	if (counts < ampMax[iw]) {
	  ampMax[iw] = counts;
	  tStamp[iw] = bin; 
	}
	histo[iw] -> SetBinContent(bin, counts);
      }
	//calculate baseline
	baseline[iw] = histo[iw] -> Integral(1,101)/100.;
	ampMax[iw] = ampMax[iw] - baseline[iw];

        histo[iw] -> Draw("same");

	int rif = 0;
	for (int scan = tStamp[iw]; scan > 0; --scan)
		{
		if(histo[iw]->GetBinContent(scan) > ampMax[iw]*frac/100) break;  
		tStamp[iw] = scan*0.05;
		rif = scan;
		}
	tStamp[iw] = trigger((rif-1)*0.05, histo[iw]->GetBinContent(rif-1),rif*0.05,histo[iw]->GetBinContent(rif), (rif+1)*0.05, histo[iw]->GetBinContent(rif+1), frac, ampMax[iw]);

	//amplitude walk correction for mcp
	if(runFolder == "WFRun012")
	if(iw == 2 && ampMax[iw] > -0.02 && ampMax[iw] < -5*8.11789e-05) tStamp[iw] = tStamp[iw] - 13.9*(-0.01 - ampMax[iw]);//(run012)
	if(runFolder == "WFRun010")
	if(iw == 2 && ampMax[iw] > -0.02 && ampMax[iw] < -0.004) tStamp[iw] = tStamp[iw] - 13.7*(-0.01 - ampMax[iw]);//(run010)
	if(runFolder == "WFRun008" || runFolder == "WFRun008")
	if(iw == 2 && ampMax[iw] > -0.02 && ampMax[iw] < -0.004) tStamp[iw] = tStamp[iw] - 15.34*(-0.01 - ampMax[iw]);// (run008,run009)
	if(runFolder == "WFRun007")
	if(iw == 2 && ampMax[iw] > -0.02 && ampMax[iw] < -0.004) tStamp[iw] = tStamp[iw] - 0.*26.85*(-0.01 - ampMax[iw]);// (run007)
	if(runFolder == "WFRun004")
	if(iw == 2 && ampMax[iw] > -0.02 && ampMax[iw] < -0.004) tStamp[iw] = tStamp[iw] - 13.9*(-0.01 - ampMax[iw]);// (run004)
	if(runFolder == "WFRun006")
	if(iw == 2 && ampMax[iw] > -0.02 && ampMax[iw] < -5*9.50036e-05) tStamp[iw] = tStamp[iw] - 19.5*(-0.01 - ampMax[iw]);//(run006)	

	}

	
    // save interesting wavefo (VEEEEERY ROUGH)
    if (fabs((tStamp[1]-tStamp[0])-11.15) < 0.25 ) 			//double coincidence
	 {
	 iwin++;
    	 if ((tStamp[2]-tStamp[0])>2.5 && (tStamp[2]-tStamp[0])<3.5 && ampMax[2]<-0.0026)
		{
		for (int iw=0;iw<3;iw++) 	
			{  
			sprintf (h1,"Ch%d_%03d",iw,isave);
	  		histo[iw]->SetName(h1);
	  		histo[iw]->SetTitle(h1);
	  		histo[iw]->Write();
			}
		isave++;
      		}	
    	}
	
    // fill ntuple
    nt->Fill();
 }
  //saving ntuple
  nt->Write();
  fw->Close();


//***************************************************************************************************************************************
//********************************************        READING NTUPLES      **************************************************************



	TFile *Reading = new TFile(runFolder + ".root");
	TTree *t1 = (TTree*)Reading->Get("nt");    	

	float amplitude[3];
	float time[3];
	float offset[3];
	
	t1->SetBranchAddress("tStamp",&time);
	t1->SetBranchAddress("ampMax",&amplitude);	
	t1->SetBranchAddress("baseline",&offset);	

	TH1F *h_diff = new TH1F ("h_diff","h_diff",1200,6.,18.);
	TH1F *h_mean_mcp = new TH1F ("h_mean_mcp","h_mean_mcp",600,-20.,20.);
	TH1F *h_1_mcp = new TH1F ("h_1_mcp","h_1_mcp",1200,-50.,50.);
	TH1F *h_mcp_2 = new TH1F ("h_mcp_2","h_mcp_2",1200,-50.,50.);

	TH1F *h_offset0 = new TH1F ("h_offset0","h_offset0",1000,-0.005,0.005);
	TH1F *h_offset1 = new TH1F ("h_offset1","h_offset1",1000,-0.005,0.005);
	TH1F *h_offset2 = new TH1F ("h_offset2","h_offset2",1000,-0.005,0.005);

	TF1 *gauss = new TF1("gauss","[3]+[1]/sqrt(2*TMath::Pi()*[0]*[0])*exp(-(x-[2])*(x-[2])/(2*[0]*[0]))",-50,50);
	gauss->SetParName(0,"sigma");
	gauss->SetParName(1,"amplitude");
	gauss->SetParName(2,"mean");
	gauss->SetParName(3,"trasl");
	gauss->SetNpx(1000000);

	for (int o = 0 ; o < t1->GetEntries() ; ++o)
	{
	t1->GetEntry(o);
        h_diff->Fill(time[1]-time[0]);
	h_mean_mcp->Fill((time[1]+time[0])*0.5-time[2]);
        h_1_mcp->Fill(time[1]-time[2]);
        h_mcp_2->Fill(time[2]-time[0]);

	h_offset0 -> Fill (offset[0]);
	h_offset1 -> Fill (offset[1]);
	h_offset2 -> Fill (offset[2]);

	}
	
//*************************************************************************************************************************************	
//*************************************************         EVALUATING BASELINE      ****************************************************

  TF1 *plain_gauss = new TF1("plain_gauss","gaus",0.,0.01);	

  TCanvas *base = new TCanvas();
  base -> Divide (3,1);

  base->cd(1);
  h_offset0->GetXaxis()->SetRangeUser(0.0008,0.0022);
  h_offset0->Fit("plain_gauss","","",0.0008,0.0022);
  float cut0 = plain_gauss->GetParameter(2);

  base->cd(2);
  h_offset1->GetXaxis()->SetRangeUser(0.0008,0.0022);
  h_offset1->Fit("plain_gauss","","",0.0008,0.0022);
  float cut1 = plain_gauss->GetParameter(2);
  
  base->cd(3);
  h_offset2->GetXaxis()->SetRangeUser(-0.0016,-0.0002);
  h_offset2->Fit("plain_gauss","","",-0.0016,-0.0002);
  float cut2 = plain_gauss->GetParameter(2);
  
  cout<<"cut0 vale "<<cut0<<endl;  
  cout<<"cut1 vale "<<cut1<<endl;
  cout<<"cut2 vale "<<cut2<<endl;
  
//***************************************************************************************************************************************
//***************************************************     DRAW CONCIDENCES     **********************************************************
  //0,1 -> external detectors
  //2 -> mcp detector
  
  TCanvas *scatter = new TCanvas();
  scatter -> cd();
  dtvsamp->Draw("colz");

  TCanvas *coincidenze = new TCanvas();
  coincidenze -> Divide (4,1);

  //h_diff => difference
  float media = 0.;
  float contenuto = 0.;
  cout<<h_diff->GetBinContent(300)<<endl;
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
//  gauss->SetParameters(0.04,2.,11.05,2.);	//Run006	
  gauss->SetParameters(0.015,5.,media,1.);	//Run 012
  h_diff->Fit("gauss","","",9.,12.);		//Run006
  std::cout<<"sigma difference = "<<gauss->GetParameter(0)<<std::endl;
  float sigma_diff = fabs(gauss->GetParameter(0));
  float mean_diff = fabs(gauss->GetParameter(2));

  //h_1_mcp => 1 - mcp
  h_1_mcp->GetXaxis()->SetRangeUser(6.5,9.5);
  h_1_mcp->SetTitle("1 - mcp");
  coincidenze -> cd(3);
  gauss->SetParameters(0.1,5.,7.6,1.);
  h_1_mcp->Fit("gauss","","",-25.,25.);
  std::cout<<"sigma 1 - mcp = "<<gauss->GetParameter(0)<<std::endl;

  //h_mcp_2 => mcp - 2
  h_mcp_2->GetXaxis()->SetRangeUser(1.,5.);
  h_mcp_2->SetTitle("mcp - 2");
  coincidenze -> cd(4);
  gauss->SetParameters(0.15,10.,3.,3.);	//Run 006
//  gauss->SetParameters(0.1,5.,3.1,1.);		//Run 012
  h_mcp_2->Fit("gauss","","",1.,5.);
  std::cout<<"sigma mcp - 2 = "<<gauss->GetParameter(0)<<std::endl;

  //h_mean_mcp => mean - mcp
  h_mean_mcp->GetXaxis()->SetRangeUser(1.,4.);
  h_mean_mcp->SetTitle("mean - mcp");
  coincidenze -> cd(2);
  gauss->SetParameters(0.05,10.,2.5,3.);	
  h_mean_mcp->Fit("gauss","","",-20.,20.);
  std::cout<<"sigma mean - mcp = "<<gauss->GetParameter(0)<<std::endl;
  float sigma_mcp = gauss->GetParameter(0);
  float mean_mcp = gauss->GetParameter(2);

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

	for (int r = 0 ; r < t1->GetEntries() ; ++r)
	{
	t1->GetEntry(r);
//	dtvsamp->Fill(amplitude[1]+amplitude[0],time[1]-time[0]);	
	//double coincidence
	if (fabs((time[1]-time[0]) - mean_diff) < 3*sigma_diff && amplitude[1] < -5*cut1 && amplitude[0] < -5*cut0) 
		{					
		doppie -> Fill(time[1]-time[0]);		
		bene++;
		p_dtvsamp->Fill(amplitude[2],(time[1]+time[0])*0.5-time[2]);	
		dtvsamp->Fill(amplitude[2],(time[1]+time[0])*0.5-time[2]);
	    	if (fabs((time[1]+time[0])*0.5-time[2] - 5.) < 1.5 && amplitude[2] < - fabs(5*cut2)) 	trfake++;	 
		if (fabs((time[1]+time[0])*0.5-time[2] - mean_mcp) < 3*sigma_mcp && amplitude[2] < -5*cut2)
			{
			benissimo++;
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
//  p_dtvsamp->Fit("retta","","",-0.035,-0.004);	//Run010
  p_dtvsamp->Fit("retta","","",-0.02,-5*cut2);		//Run006
  p_dtvsamp->Draw("colz");
  cout<<"correggerei il problema dell'amplitude walk usando la retta y = "<<retta->GetParameter(1)<<" x + "<<retta->GetParameter(0)<<endl;

cout << "saved / triggered / double fake / triple fake: " << benissimo<< " / " << bene << " / " << rfake << " / "<< trfake<< endl;
cout << "fake = " <<fake<<endl;
float efficiency = (((float)benissimo-rtrfake)/((float)bene-rfake))*100.;
cout << "efficiency is " <<efficiency<< "\%" << endl;

TCanvas *cfr = new TCanvas();
cfr->SetTitle("Doppie & Trilpe");
cfr->cd();
doppie -> Draw();
triple -> Draw("same");

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

}

int main(int argc, char** argv)
{
TApplication *myApp = new TApplication("example",&argc, argv);
TString run = "WFRun006";
plottersimple(run);
myApp->Run();
}

