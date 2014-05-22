// g++ -o plottersimpleMacros_BTF `root-config --cflags --glibs` plottersimpleMacros_BTF.cpp
// ./plottersimpleMacros_BTF ../formatRoma/ run002 0  // analyze all events from ../formatRoma/run002/merged.root
// ./plottersimpleMacros_BTF ../formatRoma/ run002 1  X // analyze and plot event X from  ../formatRoma/run002/merged.root

#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"
#include "TApplication.h"
#include "TF1.h"
#include "TGraph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h> 

bool correctForAmplitudeWalk = false;
int frac = 36;

/*
float trigger (float x1, float y1, float x2, float y2, float x3, float y3, float AmpFraction, float amp)
{
  float den = ( 3 * (x1*x1 + x2*x2 + x3*x3) - (x1+x2+x3)*(x1+x2+x3) );
  float m = ( 3 * (x1*y1 + x2*y2 + x3*y3) - (x1+x2+x3) * (y1+y2+y3) ) / den;
  float q = ( (y1+y2+y3) * (x1*x1 + x2*x2 + x3*x3) - (x1+x2+x3) * (x1*y1 + x2*y2 + x3*y3) ) / den;
  float returnValue = (amp * AmpFraction/100. - q)/m; 
  return returnValue;
}
*/

double trigger(int x1, double y1, int x2, double y2, int x3, double y3, double AmpFraction, double amp, int Npar = 3)
{
  std::vector<double> xx;
  std::vector<double> xy;
  double Sx = 0.;
  double Sy = 0.;
  double Sxx = 0.;
  double Sxy = 0.;

  xx.clear();
  xy.clear();

  xx.push_back(x1*x1);  
  xx.push_back(x2*x2);
  xx.push_back(x3*x3);

  xy.push_back(x1*y1);  
  xy.push_back(x2*y2);
  xy.push_back(x3*y3);

  Sx = x1+x2+x3;
  Sy = y2+y2+y3;
  Sxx = xx.at(0) + xx.at(1) + xx.at(2);
  Sxy = xy.at(0) + xy.at(1) + xy.at(2);

  double Delta = 3*Sxx - Sx*Sx;
  double A = (Sxx*Sy - Sx*Sxy) / Delta;
  double B = (3*Sxy - Sx*Sy) / Delta;
 
  // A+Bx = AmpFraction * amp
  double interpolation = (amp * AmpFraction/100. - A) / B; 

  if(isnan(interpolation) == true) {
    std::cout << " >>> A = " << A << " B = " <<  B << std::endl;
    std::cout << " >>> Delta = " << Delta << std::endl;
    std::cout << " >>> 3*Sxx = " << 3.*Sxx << std::endl;
    std::cout << " >>> SxSx = " << Sx*Sx << std::endl;
    std::cout << " >>> Sx = " << Sx << " Sy = " <<  Sy << std::endl;
    std::cout << " >>> Sxx = " << Sxx << " Sxy = " <<  Sxy << std::endl;
    std::cout << " >>> xx.at(0) = " << xx.at(0) << " xx.at(1) = " << xx.at(1) << " xx.at(2) = " << xx.at(2) << std::endl;
    std::cout << " >>> x1, y1 = " << x1 << " " << y1 << " >>> x2, y2 = " << x2 << " " << y2
	      << " >>> x3, y3 = " << x3 << " " << y3 << std::endl;
  }
  return interpolation;
}

int main(int argc, char** argv)
{
  gROOT->Reset();


  // Settings
  std::string  Folder = argv[1];
  std::string  runFolder = argv[2];
  int drawEvent;
  if(argc > 3) drawEvent = atoi(argv[3]);
  int eventNumber;
  if(argc > 4) eventNumber = atoi(argv[4]);

  std::cout << "Folder = " << Folder << std::endl;
  std::cout << "runFolder = " << runFolder << std::endl;
  std::cout << "drawEvent = " << drawEvent << std::endl;
  std::cout << "eventNumber = " << eventNumber << std::endl;

  TApplication* myApp;
  if(drawEvent == 1) myApp = new TApplication("myApp", &argc, argv);

  const int nDigiCh = 9;
  const float timeToDigi = 0.2;
  int nSampling = 1024;

  int colorsExt[9] = {kRed+1, kBlue, kPink, kGreen+1, kCyan+2, kOrange+7, kAzure+7, kWhite, kBlack};
  std::vector<int> colors;
  for(int posVec = 0; posVec<nDigiCh; ++posVec){
    colors.push_back(colorsExt[posVec]);
  }


  // Define histograms
  TH1F* histo[nDigiCh];
  char h1[10];
  for (int iw=0;iw<nDigiCh;iw++){
    sprintf (h1,"Ch%d",iw);
    histo[iw] = new TH1F( h1, "", 1024, 0., 1024.);
    histo[iw] -> SetXTitle(Form("Time (channel * %f) ns", timeToDigi));
    histo[iw] -> SetYTitle("Amplitude (ADC)");
    histo[iw]->SetLineWidth(2);
    histo[iw]->SetLineColor(colors.at(iw));   
  }
  
  TH1F* dummyPed = new TH1F("dummyPed", "", 1024, 0., 1024);
  for(int ii=0; ii<nSampling; ++ii){
    dummyPed->SetBinContent(ii+1, 1);
  }

  TH1F* hBaseLine[nDigiCh];
  char h2[10];
  for (int iw=0;iw<nDigiCh;iw++){
    sprintf (h2,"Bs%d",iw);
    hBaseLine[iw] = new TH1F( h2, "", 400, 4000, -4000);
    hBaseLine[iw] -> SetXTitle ("BaseLine (V)");
    hBaseLine[iw]->SetLineColor(colors.at(iw));
  }

  
  // Define graphics
  TCanvas *cc = new TCanvas("cc","cc");
  cc->cd(); 
  cc->SetGrid();
  histo[0]->Draw("");

  std::vector<std::string> nameChannels;
  nameChannels.push_back("MIB1");
  nameChannels.push_back("MIB2");
  nameChannels.push_back("scintB");
  nameChannels.push_back("PLANA");
  nameChannels.push_back("MIB3");
  nameChannels.push_back("ROMAX");
  nameChannels.push_back("scintF");
  nameChannels.push_back("nan");
  nameChannels.push_back("trig0");


  TLegend* leg = new TLegend(0.88,0.65,0.98,0.85);
  leg->SetFillColor(0);
  leg->SetTextFont(41);
  for (int iw=0;iw<nDigiCh;iw++) {
    leg->AddEntry(histo[iw],nameChannels.at(iw).c_str(),"l");
  }
  leg->Draw();


  // Input tree
  std::string nameIn = std::string(Folder) + std::string(runFolder) + "/merged.root"; 
  TFile* inF = TFile::Open(nameIn.c_str());
  TTree* tntu = (TTree*)inF->Get("eventRawData");
  std::cout << " >> inFile->GetEntries() = " << tntu->GetEntries() << std::endl;

  unsigned int nDigiSamples;
  unsigned int digiGroup[9216];
  unsigned int digiChannel[9216];
  unsigned int digiSampleIndex[9216];
  float digiSampleValue[9216];
  unsigned int nAdcChannels;
  unsigned int adcData[40];
  unsigned int adcBoard[40];
  unsigned int adcChannel[40];
  
  tntu->SetBranchStatus("*",0);

  tntu->SetBranchStatus("nDigiSamples",1);     tntu->SetBranchAddress("nDigiSamples", &nDigiSamples);
  tntu->SetBranchStatus("digiGroup",1);        tntu->SetBranchAddress("digiGroup", digiGroup);
  tntu->SetBranchStatus("digiChannel",1);      tntu->SetBranchAddress("digiChannel",digiChannel);
  tntu->SetBranchStatus("digiSampleIndex",1);  tntu->SetBranchAddress("digiSampleIndex",digiSampleIndex);
  tntu->SetBranchStatus("digiSampleValue",1);  tntu->SetBranchAddress("digiSampleValue",digiSampleValue);

  tntu->SetBranchStatus("nAdcChannels",1);     tntu->SetBranchAddress("nAdcChannels",&nAdcChannels);
  tntu->SetBranchStatus("adcData",1);          tntu->SetBranchAddress("adcData",adcData);
  tntu->SetBranchStatus("adcBoard",1);         tntu->SetBranchAddress("adcBoard",adcBoard);
  tntu->SetBranchStatus("adcChannel",1);       tntu->SetBranchAddress("adcChannel",adcChannel);


  // Output tree
  TString filOut = "Ntuples/"+runFolder + "analyzed_btf.root";
  TFile* f = new TFile(filOut,"recreate");

  float ampMax[nDigiCh];
  float ped[nDigiCh];
  float integral[nDigiCh];
  float tStamp[nDigiCh]; 
  float tStampMax[nDigiCh]; 
  float bLine[nDigiCh];
  float bLineSig[nDigiCh];
  int nEle;

  TTree* nt = new TTree("nt", "nt");
  nt->Branch("ampMax",&ampMax,"ampMax[9]/F"); 
  nt->Branch("ped",&ped,"ped[9]/F"); 
  nt->Branch("integral",&integral,"integral[9]/F"); 
  nt->Branch("tStamp",&tStamp,"tStamp[9]/F"); 
  nt->Branch("tStampMax",&tStampMax,"tStampMax[9]/F"); 
  nt->Branch("bLine",&bLine,"bLine[9]/F");
  nt->Branch("bLineSig",&bLineSig,"bLineSig[9]/F"); 
  nt->Branch("nEle",&nEle,"nEle/I"); 


  // Loop over input tree
  for(int iEntry=0; iEntry<tntu->GetEntries();++iEntry){
    if(drawEvent == 1 && iEntry < eventNumber) continue;
    if(drawEvent == 1 && iEntry > eventNumber) break;

    if(iEntry%1000 == 0)    std::cout << " Reading entry " << iEntry << std::endl;
    tntu->GetEntry(iEntry);

    for(int iCount=0; iCount<nDigiCh; ++iCount){
      //clean before fill
      ampMax[iCount] = 40000;
      ped[iCount] = 0;
      integral[iCount] = 0;
      tStamp[iCount] = 0;
      tStampMax[iCount] = 0;
      bLine[iCount] = 0;
      bLineSig[iCount] = 0;
      nEle = -1;
    }

    for (int bin=0;bin<nDigiSamples;++bin){
       if(digiChannel[bin] == 7){
 	bin += 1024;
 	continue;      
       }
       //      std::cout << " bin, counts, time = " << bin << " " << digiSampleValue[bin] << " " << digiSampleIndex[bin] << std::endl;
      
       if(digiChannel[bin] == 3) digiSampleValue[bin] = -digiSampleValue[bin]; // planacon readout at the MCP output not anode 
            
      // compute pedestal over first 100 points = 100*0.2ns 
       if((bin+1)%1024 != 0 && digiSampleIndex[bin] >= 0. && digiSampleIndex[bin] <= 100.) {
	 ped[digiChannel[bin]] += digiSampleValue[bin];
       }
      
      if (digiSampleValue[bin] < ampMax[digiChannel[bin]]) {
	ampMax[digiChannel[bin]] = digiSampleValue[bin];
	tStamp[digiChannel[bin]] = float(digiSampleIndex[bin]*timeToDigi); 
      }
      histo[digiGroup[bin]+digiChannel[bin]]->SetBinContent(digiSampleIndex[bin]+1, digiSampleValue[bin]);

      //if last entry before moving to next channel
      if( (bin+1)%1024 == 0 ){
	//	std::cout << " >> moving to channel = " << digiChannel[bin] << std::endl;
	ped[digiChannel[bin]] = ped[digiChannel[bin]] / 100.;
	ampMax[digiChannel[bin]] = ampMax[digiChannel[bin]] - ped[digiChannel[bin]];
	tStampMax[digiChannel[bin]] = tStamp[digiChannel[bin]];

	if(ampMax[digiChannel[bin]] > 0. || tStampMax[digiChannel[bin]] == 1) continue;      
	if(ampMax[digiChannel[0]] > -100.) break;
	
	histo[digiChannel[bin]]->Add(dummyPed, (-1.*ped[digiChannel[bin]])); 
	hBaseLine[digiChannel[bin]]->Fill(ped[digiChannel[bin]]);
	
	cc->cd();
	if(digiChannel[bin] != 7)  histo[digiChannel[bin]]->Draw("sames");

	if((digiChannel[bin] == 0 || digiChannel[bin] == 1 || digiChannel[bin] == 3 || digiChannel[bin] == 4 || digiChannel[bin] == 5) &&
	   (ampMax[digiChannel[bin]] < -100.)){

	  int rif = 0;
	  for(int scan = int(tStamp[digiChannel[bin]]/timeToDigi); scan > 0; --scan)
	    {
	      if(histo[digiChannel[bin]]->GetBinContent(scan) > ampMax[digiChannel[bin]]*frac/100) break;  
	      tStamp[digiChannel[bin]] = scan*timeToDigi;
	      rif = scan;
	    }
	  
	  tStamp[digiChannel[bin]] = float(trigger((rif-1), histo[digiChannel[bin]]->GetBinContent(rif-1), 
						   rif, histo[digiChannel[bin]]->GetBinContent(rif), 
						   (rif+1), histo[digiChannel[bin]]->GetBinContent(rif+1), frac, ampMax[digiChannel[bin]]))*timeToDigi;
	}

	integral[digiChannel[bin]] = histo[digiChannel[bin]]->Integral(tStamp[digiChannel[bin]], tStampMax[digiChannel[bin]]);      

	if(drawEvent == 1){
	  sprintf (h1,"Ch%d",digiChannel[bin]);
	  histo[digiChannel[bin]]->SetName(h1);
	  histo[digiChannel[bin]]->SetTitle(h1);
	  histo[digiChannel[bin]]->Write();
	}
      }// ok next channel

      //       sprintf (h2,"Bs%d",digiChannel[bin]);
      //       hBaseLine[digiChannel[bin]]->SetName(h2);
      //       hBaseLine[digiChannel[bin]]->SetTitle(h2);
      //       hBaseLine[digiChannel[bin]]->Write();
    }// loop over bin

    for(int iCh=0; iCh<nAdcChannels; ++iCh){
      if(adcBoard[iCh] == 1 && adcChannel[iCh] == 0){
	if(adcData[iCh] < 500) nEle = 0;
	if(adcData[iCh] > 500 && adcData[iCh] < 1500) nEle = 1;
	if(adcData[iCh] > 1500 && adcData[iCh] < 2500) nEle = 2;
	if(adcData[iCh] > 2500 && adcData[iCh] < 3500) nEle = 3;
      }
    }

    nt->Fill(); 
  }// loop over entries	  

  // set pedestal RMS and mean for each channel computed over the all entries (only valid entry is the last one of output tree)
  for(int iw=0; iw<nDigiCh; ++iw){
    ampMax[iw] = 0;
    ped[iw] = 0;
    integral[iw] = 0;
    tStamp[iw] = 0;
    tStampMax[iw] = 0;
    nEle = -1;
    bLine[iw] = hBaseLine[iw]->GetMean();
    bLineSig[iw] = hBaseLine[iw]->GetRMS();
  }
  nt->Fill();

  std::cout << " END " << std::endl;
  nt->Write();
  f->Close();
  inF->Close();
  std::cout << f->GetName() << " created " << std::endl;
  if(drawEvent == 1)  myApp->Run(true);
  return 0;
}

