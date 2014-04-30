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
#include "TProfile.h"
#include "TF1.h"
#include "TGraph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

TH2F *dtvsamp = new TH2F("dtvsamp","",100,-0.1,0.,1200,-20.,20.);
TGraph *chi = new TGraph();
TGraph *sigma = new TGraph();
std::vector<float> time0;
std::vector<float> time1;
std::vector<float> time2;

TProfile *p_dtvsamp = new TProfile("dtvsamp","",100,-0.1,0.);//,1200,-20.,20.);

int frac = 36;

float trigger (float x1, float y1, float x2, float y2, float x3, float y3, float frazione, float amp)
{
float ritorno;
float denominatore = (3*(x1*x1 + x2*x2 + x3*x3) - (x1+x2+x3)*(x1+x2+x3));
float m = (3*(x1*y1 + x2*y2 + x3*y3) - (x1+x2+x3)*(y1+y2+y3))/denominatore;
float q = ((y1+y2+y3)*(x1*x1 + x2*x2 + x3*x3) - (x1+x2+x3)*(x1*y1 + x2*y2 + x3*y3))/denominatore;
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
  TTree *nt = new TTree("nt","nt");

  nt->Branch("ampMax",&ampMax,"ampMax[3]/F"); 
  nt->Branch("tStamp",&tStamp,"tStamp[3]/F"); 


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
  string line;
  int ifile=0;
  int isave=0;
  int iwin =0;

    // define file for wf
    TString fileWF = runFolder + "_wf_save.root"; 
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
	if(iw == 2 && ampMax[iw] > -0.02 && ampMax[iw] < -0.004) tStamp[iw] = tStamp[iw] - 14.2*(-0.01 - ampMax[iw]);	
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
 /*   histos->Fill(tStamp[1]-tStamp[0]);
    h_signal->Fill((tStamp[1]+tStamp[0])*0.5-tStamp[2]);
    h_signal2->Fill(tStamp[1]-tStamp[2]);
    h_signal3->Fill(tStamp[2]-tStamp[0]);
	dtvsamp->Fill(ampMax[2],(tStamp[1]+tStamp[0])*0.5-tStamp[2]);	
	p_dtvsamp->Fill(ampMax[2],(tStamp[1]+tStamp[0])*0.5-tStamp[2]);	
	time0.push_back(tStamp[0]);  
	time1.push_back(tStamp[1]);
	time2.push_back(tStamp[2]);
*/

}
  fw->Close();

//  cout << "saved / triggered: " << isave << " / " << iwin << endl;
  
  TString filOut = runFolder + ".root"; 
  TFile *f = new TFile(filOut,"recreate"); 
  nt->Write();
  f->Close();
//}//fine del ciclo for
//***************************************************************************************************************************************
//********************************************        READING NTUPLES      **************************************************************

	TFile *Reading = new TFile("WFRun006.root");
	TTree *t1 = (TTree*)Reading->Get("nt");    	

	float amplitude[3];
	float time[3];
	
	t1->SetBranchAddress("tStamp",&time);
	t1->SetBranchAddress("ampMax",&amplitude);	

	TH1F *histos = new TH1F ("histos","histos",1200,6.,18.);
	TH1F *h_signal = new TH1F ("h_signal","h_signal",600,-20.,20.);
	TH1F *h_signal2 = new TH1F ("h_signal2","h_signal2",1200,-50.,50.);
	TH1F *h_signal3 = new TH1F ("h_signal3","h_signal3",1200,-50.,50.);
	/*
	TH1F* h_base[3];
	char hb[10];
	for (int i=0;i<3;i++)
		{
		sprintf (hb,"Baseline_%d",i);
		hbase[iw] = new TH1F(hb,hb,1000,0.,1.);
		}
	*/
	TF1 *gauss = new TF1("gauss","[3]+[1]/sqrt(2*TMath::Pi()*[0]*[0])*exp(-(x-[2])*(x-[2])/(2*[0]*[0]))",-50,50);
	gauss->SetParName(0,"sigma");
	gauss->SetParName(1,"amplitude");
	gauss->SetParName(2,"mean");
	gauss->SetParName(3,"trasl");
	gauss->SetNpx(1000000);

	for (int o = 0 ; o < t1->GetEntries() ; ++o)
	{
	t1->GetEntry(o);
        histos->Fill(time[1]-time[0]);
	h_signal->Fill((time[1]+time[0])*0.5-time[2]);
        h_signal2->Fill(time[1]-time[2]);
        h_signal3->Fill(time[2]-time[0]);
	//dtvsamp->Fill(amplitude[2],(time[1]+time[0])*0.5-time[2]);	
	p_dtvsamp->Fill(amplitude[2],(time[1]+time[0])*0.5-time[2]);	
//	time0.push_back(time[0]);  
//	time1.push_back(time[1]);
//	time2.push_back(time[2]);
	

	}

//*****************************************************************************************************************************************
//***************************************************     DRAW CONCIDENCES     ************************************************************

  //0,1 -> external detectors
  //2 -> mcp detector
  
  TCanvas *scatter = new TCanvas();
  scatter -> cd();
  dtvsamp->Draw("colz");

  TCanvas *p_scatter = new TCanvas();
  p_scatter -> cd();
  TF1 *retta = new TF1 ("retta","pol1",-0.02,0.);
  p_dtvsamp->Fit("retta","","",-0.02,-0.004);
  p_dtvsamp->Draw("colz");
  cout<<"correggerei il problema dell'amplitude walk usando la retta y = "<<retta->GetParameter(1)<<" x + "<<retta->GetParameter(0)<<endl;

  TCanvas *coincidenze = new TCanvas();
  coincidenze -> Divide (4,1);

  //histos => difference
  histos->GetXaxis()->SetRangeUser(10.8,11.3);
  histos->SetTitle("difference");	
  coincidenze -> cd(1);
  gauss->SetParameters(0.1,30.,11.2,10.);	
  histos->Fit("gauss","","",8.,16.);
  std::cout<<"sigma difference = "<<gauss->GetParameter(0)<<std::endl;
  float sigma_diff = gauss->GetParameter(0);
  float mean_diff = gauss->GetParameter(2);

  //h_signal2 => 1 - mcp
  h_signal2->GetXaxis()->SetRangeUser(6.5,9.5);
  h_signal2->SetTitle("1 - mcp");
  coincidenze -> cd(3);
  gauss->SetParameters(0.1,50.,8,3.);
  h_signal2->Fit("gauss","","",-25.,25.);
  std::cout<<"sigma 1 - mcp = "<<gauss->GetParameter(0)<<std::endl;

  //h_signal3 => mcp - 2
  h_signal3->GetXaxis()->SetRangeUser(1.,5.);
  h_signal3->SetTitle("mcp - 2");
  coincidenze -> cd(4);
  gauss->SetParameters(0.15,10.,3.,3.);	
  h_signal3->Fit("gauss","","",1.,5.);
  std::cout<<"sigma mcp - 2 = "<<gauss->GetParameter(0)<<std::endl;

  //h_signal => mean - mcp
  h_signal->GetXaxis()->SetRangeUser(1.,4.);
  h_signal->SetTitle("mean - mcp");
  coincidenze -> cd(2);
  gauss->SetParameters(0.05,10.,2.5,3.);	
  h_signal->Fit("gauss","","",-20.,20.);
  std::cout<<"sigma mean - mcp = "<<gauss->GetParameter(0)<<std::endl;
  float sigma_mcp = gauss->GetParameter(0);
  float mean_mcp = gauss->GetParameter(2);

  chi->SetPoint(frac,frac,gauss->GetChisquare()/gauss->GetNDF());
  sigma->SetPoint(frac,frac,fabs(gauss->GetParameter(0))); 

//*****************************************************************************************************************************************
//*************************************************         EVALUATING EFFICIENCY        *************************************************

int bene = 0;
int benissimo = 0;
int fake = 0;
int drake = 0;
	
	for (int r = 0 ; r < t1->GetEntries() ; ++r)
	{
	t1->GetEntry(r);
//	cout << time[1] <<endl;    
	//dtvsamp->Fill(amplitude[1]+amplitude[0],time[1]-time[0]);	
	if (fabs((time[1]-time[0]) - mean_diff) < 3*sigma_diff ) 			//double coincidence
		{					
		bene++;
		dtvsamp->Fill(amplitude[2],(time[1]+time[0])*0.5-time[2]);
	    	if (fabs((time[1]+time[0])*0.5-time[2] - 5.) < 1.5 && amplitude[2]<-0.002)	drake++;	 
		if (fabs((time[1]+time[0])*0.5-time[2] - mean_mcp) < 3*sigma_mcp && amplitude[2]<-0.002)
			{
			 	
			 benissimo++;

			}	      	 
		}
	// if (fabs(time[1]-time[0] - mean_diff + 12*sigma_diff) < 3*sigma_diff)// && amplitude[2]<-0.0026)
	 if (fabs(time[1]-time[0] - 13.5) < 1.5)// && amplitude[2]<-0.0026)
		fake++;	
			
	}

float rfake = (float)fake*3.*sigma_diff/1.5;
float rdrake = (float)drake*3.*sigma_mcp/1.5;

cout << "saved / triggered / fake: " << benissimo<< " / " << bene << " / " << rfake << " / "<< drake<< endl;
float efficiency = (((float)benissimo-rdrake)/((float)bene-(float)rfake))*100.;
cout << "efficiency is " <<efficiency<< "\%" << endl;

//*****************************************************************************************************************************************
//************************************************     DRAW SIGMA - CHI    ****************************************************************
		
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
//*************************************************************************************************************************************	
/*	sort(time0.begin(),time0.end());
	sort(time1.begin(),time1.end());
	sort(time2.begin(),time2.end());
	
	std::cout<<"min_0 = "<<time0.at(0)<<std::endl;
	std::cout<<"min_1 = "<<time1.at(0)<<std::endl;
	std::cout<<"min_2 = "<<time2.at(0)<<std::endl;
	std::cout<<"lunghezza = "<<time2.size()<<std::endl;
*/
//*****************************************************************************************************************************************
//**********************************************      DRAW WAVEFORM	*******************************************************************
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
