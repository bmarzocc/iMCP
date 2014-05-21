// g++ -o plottersimpleMacros_BTF `root-config --cflags --glibs` plottersimpleMacros_BTF.cpp
// ./plottersimpleMacros_BTF ../formatRoma/ run002 0  // analyze all events from ../formatRoma/run002/merged.root
// ./plottersimpleMacros_BTF ../formatRoma/ run002 1  X // analyze and plot event X from  ../formatRoma/run002/merged.root

/*
	Ricordati che la media sara' circa 10ns, cosi' recuperi la scala
sembra buono:	emailData147_140327184601
*/

// root
// .L plottersimpleMacros.C+
// plottersimpleMacros("path_from_local_to_WaveForms", "WFRunX")

// .L plottersimpleMacros.C+
// drawSimple("WFRunX")



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

TGraph* chi = new TGraph();
TGraph* sigma = new TGraph();
std::vector<float> time0;
std::vector<float> time1;
std::vector<float> time2;

bool drawSigmaChi2 = false;
bool correctForAmplitudeWalk = false;

//TProfile* p_dtvsamp = new TProfile("p_dtvsamp", "", 1000, -100., 100.);//,1200,-20.,20.);
//TH2F* dtvsamp = new TH2F("dtvsamp", "", 1000, -0.1, 0., 2000, -1000., 1000.);

TProfile* p_dtvsamp = new TProfile("p_dtvsamp", "", 100, -0.1, 0.);//,1200,-20.,20.);
//TProfile* p_dtvsamp = new TProfile("p_dtvsamp", "", 10000, -0.1, 0.);//,1200,-20.,20.);
TH2F* dtvsamp = new TH2F("dtvsamp", "", 100, -0.1, 0.0, 1200, -25., 25.);



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

 //double trigger(int& x1, double& y1, int& x2, double& y2, int& x3, double& y3, double& AmpFraction, double& amp, int Npar = 3)
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

  //std::cout << " >> interpolation = " << interpolation << std::endl;
  return interpolation;
}

//void plottersimpleMacros_BTF(TString Folder, TString runFolder)
int main(int argc, char** argv)
{
  gROOT->Reset();


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

  int colorsExt[9] = {kRed+1, kBlue, kPink, kGreen+1, kCyan+2, kOrange+7, kAzure+7, kWhite, kBlack};
  std::vector<int> colors;
  for(int posVec = 0; posVec<nDigiCh; ++posVec){
    colors.push_back(colorsExt[posVec]);
  }

  int nSampling = 1024;

  // Define histograms
  TH1F* histo[nDigiCh];
  char h1[10];
    for (int iw=0;iw<nDigiCh;iw++){
      sprintf (h1,"Ch%d",iw);
      histo[iw] = new TH1F( h1, "", 1024, 0., 1024.);
      histo[iw] -> SetXTitle ("Time (ns)");
      histo[iw] -> SetYTitle ("Amplitude (V)");
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



  // define graphics
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


  //input TTre
  std::string nameIn = std::string(Folder) + std::string(runFolder) + "/merged.root"; 
  TFile* inF = TFile::Open(nameIn.c_str());
  TTree* tntu = (TTree*)inF->Get("eventRawData");

  std::cout << " >> inFile->GetEntries() = " << tntu->GetEntries() << std::endl;
  //input tree
  unsigned int nDigiSamples;
  unsigned int digiGroup[9216];
  unsigned int digiChannel[9216];
  unsigned int digiSampleIndex[9216];
  float digiSampleValue[9216];
  
  tntu->SetBranchStatus("*",0);

  tntu->SetBranchStatus("nDigiSamples",1);     tntu->SetBranchAddress("nDigiSamples", &nDigiSamples);
  tntu->SetBranchStatus("digiGroup",1);        tntu->SetBranchAddress("digiGroup", digiGroup);
  tntu->SetBranchStatus("digiChannel",1);      tntu->SetBranchAddress("digiChannel",digiChannel);
  tntu->SetBranchStatus("digiSampleIndex",1);  tntu->SetBranchAddress("digiSampleIndex",digiSampleIndex);
  tntu->SetBranchStatus("digiSampleValue",1);  tntu->SetBranchAddress("digiSampleValue",digiSampleValue);

  ///OUTPUT TTREE
  // define file for wf
  TString filOut = "Ntuples/"+runFolder + "analyzed_btf.root";
  TFile* f = new TFile(filOut,"recreate");

  // Tree branches and tree strcuture
  float ampMax[nDigiCh];
  float ped[nDigiCh];
  float integral[nDigiCh];
  float tStamp[nDigiCh]; 
  float tStampMax[nDigiCh]; 
  float bLine[nDigiCh];
  float bLineSig[nDigiCh];

  TTree* nt = new TTree("nt", "nt");
  nt->Branch("ampMax",&ampMax,"ampMax[9]/F"); 
  nt->Branch("ped",&ped,"ped[9]/F"); 
  nt->Branch("integral",&integral,"integral[9]/F"); 
  nt->Branch("tStamp",&tStamp,"tStamp[9]/F"); 
  nt->Branch("tStampMax",&tStampMax,"tStampMax[9]/F"); 
  nt->Branch("bLine",&bLine,"bLine[9]/F"); 
  nt->Branch("bLineSig",&bLineSig,"bLineSig[9]/F"); 

  std::cout << " tntu->GetEntries() = " << tntu->GetEntries() << std::endl;
  //  for(int iEntry=4000; iEntry<4050;++iEntry){
  for(int iEntry=0; iEntry<tntu->GetEntries();++iEntry){
    if(drawEvent == 1 && iEntry < eventNumber) continue;
    if(drawEvent == 1 && iEntry > eventNumber) break;

    if(iEntry%100 == 0)    std::cout << " Reading entry " << iEntry << std::endl;
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
    }

    for (int bin=0;bin<nDigiSamples;++bin){

      //      std::cout << " bin, counts, time = " << bin << " " << digiSampleValue[bin] << " " << digiSampleIndex[bin] << std::endl;
      
      if(digiChannel[bin] == 3) digiSampleValue[bin] = -digiSampleValue[bin]; // planacon readout at the MCP output 
                                                                              // plane and not at the anode 
            
      //20ns /0.20 = 200/2 = 100 
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
	   (ampMax[digiChannel[digiChannel[bin]]] < -100.)){

	  int rif = 0;
	  for(int scan = int(tStamp[digiChannel[bin]]/timeToDigi); scan > 0; --scan)
	    {
	      //	    std::cout << "scan = " << scan << std::endl;
	      if(histo[digiChannel[bin]]->GetBinContent(scan) > ampMax[digiChannel[bin]]*frac/100) break;  
	      tStamp[digiChannel[bin]] = scan*timeToDigi;
	      rif = scan;
	    }
	  
	  tStamp[digiChannel[bin]] = float(trigger((rif-1), histo[digiChannel[bin]]->GetBinContent(rif-1), 
						   rif, histo[digiChannel[bin]]->GetBinContent(rif), 
						   (rif+1), histo[digiChannel[bin]]->GetBinContent(rif+1), frac, ampMax[digiChannel[bin]]))*timeToDigi;
	}

	if(tStamp[digiChannel[bin]] < 0 || isnan(tStamp[digiChannel[bin]]) == true) {
// 	  std::cout << " >>> iEntry = " << iEntry << std::endl;
// 	  std::cout << " >>> tStamp[digiChannel[bin]] = " << tStamp[digiChannel[bin]] << std::endl;
// 	  std::cout << " >>> tStampMax[digiChannel[bin]] = " << tStampMax[digiChannel[bin]] << std::endl;

	}
      
	integral[digiChannel[bin]] = histo[digiChannel[bin]]->Integral(tStamp[digiChannel[bin]], tStampMax[digiChannel[bin]]);      
      }// ok next channel
      

//       sprintf (h1,"Ch%d",digiChannel[bin]);
//       histo[digiChannel[bin]]->SetName(h1);
//       histo[digiChannel[bin]]->SetTitle(h1);
//       histo[digiChannel[bin]]->Write();
      
//       sprintf (h2,"Bs%d",digiChannel[bin]);
//       hBaseLine[digiChannel[bin]]->SetName(h2);
//       hBaseLine[digiChannel[bin]]->SetTitle(h2);
//       hBaseLine[digiChannel[bin]]->Write();


    }// loop over bin
    nt->Fill(); 
  }// loop over entries	  

  for(int iw=0; iw<nDigiCh; ++iw){
    ampMax[iw] = 0;
    ped[iw] = 0;
    integral[iw] = 0;
    tStamp[iw] = 0;
    tStampMax[iw] = 0;

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




/*
  TF1* bs0 = new TF1("bs0", "gaus", -5., 5.);
  bs0->SetParameters(20., -0.001, 0.0001);      
  TF1* bs1 = (TF1*)bs0->Clone("bs1");
  TF1* bs2 = (TF1*)bs0->Clone("bs2");

  for(int counter=0; counter<nDigiCh; ++counter){
    
  
  
  hBaseLine[0]->Fit("bs0");
  hBaseLine[1]->Fit("bs1");
  hBaseLine[2]->Fit("bs2");


       for(int iw=0; iw<3; ++iw){
	 ampMax[iw] = 0;
	 ped[iw] = 0;
	 integral[iw] = 0;
	 tStamp[iw] = 0;
	 tStampMax[iw] = 0;
       }

       bLine[0] = bs0->GetParameter(1);
       bLine[1] = bs1->GetParameter(1);
       bLine[2] = bs2->GetParameter(1);
       bLineSig[0] = bs0->GetParameter(2);
       bLineSig[1] = bs1->GetParameter(2);
       bLineSig[2] = bs2->GetParameter(2);
       nt->Fill();

       TCanvas* cBaseLine = new TCanvas("cBaseLine", "cBaseLine");
       cBaseLine->cd();
       hBaseLine[0]->Draw();
       hBaseLine[1]->Draw("same");
       hBaseLine[2]->Draw("same");
       

       bs0->Write();
       bs1->Write();
       bs2->Write();

       

}//fine del ciclo for

*/

//***************************************************************************************************************************************
//********************************************        READING NTUPLES      **************************************************************
void drawSimple(TString runFolder){

  TFile* Reading = new TFile(TString(runFolder+".root"), "read");
  TTree* t1 = (TTree*)Reading->Get("nt");    	

  float amplitude[3];
  float ped[3];
  float integral[3];
  float time[3];
  float timeMax[3];
  float bLine[3];
  float bLineSig[3];
  t1->SetBranchAddress("ampMax",amplitude);		
  t1->SetBranchAddress("ped",ped);		
  t1->SetBranchAddress("integral",integral);		
  t1->SetBranchAddress("tStamp",time);
  t1->SetBranchAddress("tStampMax",timeMax);
  t1->SetBranchAddress("bLine",bLine);
  t1->SetBranchAddress("bLineSig",bLineSig);
  float bLine_OK[3];
  float bLineSig_OK[3];

  TH1F* histos = new TH1F ("histos", "histos", 200, 6., 18.);
  //  TH1F* histos = new TH1F ("histos", "histos", 1200, -100., 100.);
  if(runFolder == "WFRun010") histos = new TH1F ("histos", "histos", 200, 6., 18.);
  //  if(runFolder == "WFRun008") histos = new TH1F ("histos", "histos", 1200, 6., 18.);
  //  if(runFolder == "WFRun007") histos = new TH1F ("histos", "histos", 1200, 6., 18.);

  TH1F* h_signal = new TH1F ("h_signal", "h_signal", 600, -20., 20.);
  if(runFolder == "WFRun012") h_signal = new TH1F ("h_signal", "h_signal", 1200, -40., 40.);
  if(runFolder == "WFRun010") h_signal = new TH1F ("h_signal", "h_signal", 1200, -40., 40.);
  TH1F* h_signal2 = new TH1F ("h_signal2", "h_signal2", 1200, -50., 50.);
  TH1F* h_signal3 = new TH1F ("h_signal3", "h_signal3", 1200, -50., 50.);

  TH1F* h_integral_Ch0 = new TH1F("h_integral_Ch0", "", 1000, -1., 50.);
  TH1F* h_integral_Ch1 = new TH1F("h_integral_Ch1", "", 1000, -1., 50.);
  TH1F* h_integral_Ch2 = new TH1F("h_integral_Ch2", "", 1000, -1., 50.);

  TH1F* h_amplitude_Ch0 = new TH1F("h_amplitude_Ch0", "", 1000, -0.1, 0.05);
  TH1F* h_amplitude_Ch1 = new TH1F("h_amplitude_Ch1", "", 1000, -0.1, 0.05);
  TH1F* h_amplitude_Ch2 = new TH1F("h_amplitude_Ch2", "", 1000, -0.1, 0.05);

  TH1F* h_triple = new TH1F("h_triple", "", 300, 0., 18);
  TH1F* h_double = new TH1F("h_double", "", 300, 0., 18);

  /*
    TH1F* h_base[3];
    char hb[10];
    for (int i=0;i<3;i++)
    {
    sprintf (hb,"Baseline_%d",i);
    hbase[iw] = new TH1F(hb,hb,1000,0.,1.);
    }
  */

  TF1* gGauss = new TF1("gGauss", "[3] + [1]/sqrt(2*TMath::Pi()*[0]*[0])*exp(-(x-[2])*(x-[2])/(2*[0]*[0]))", -50, 50);
  gGauss->SetParName(0, "sigma");
  gGauss->SetParName(1, "amplitude");
  gGauss->SetParName(2, "mean");
  gGauss->SetParName(3, "offset");
  gGauss->SetNpx(100000);
  
  t1->GetEntry(t1->GetEntries()-1);
  bLine_OK[0] = bLine[0];
  bLine_OK[1] = bLine[1];
  bLine_OK[2] = bLine[2];
  bLineSig_OK[0] = bLineSig[0];
  bLineSig_OK[1] = bLineSig[1];
  bLineSig_OK[2] = bLineSig[2];

  std::cout << " >>> bLineSig_OK[0] = " << bLineSig_OK[0] << std::endl;
  std::cout << " >>> bLineSig_OK[1] = " << bLineSig_OK[1] << std::endl;
  std::cout << " >>> bLineSig_OK[2] = " << bLineSig_OK[2] << std::endl;

  std::cout << " >>> number of WF = " << t1->GetEntries() << std::endl;
  for(int iCount=0; iCount<t1->GetEntries()-2; ++iCount){
    t1->GetEntry(iCount);

    if(amplitude[0] > 0. || timeMax[0] == 1 ||
       amplitude[1] > 0. || timeMax[1] == 1) continue;


    if(//integral[0] > -15. && integral[1] > -15. && 
	 amplitude[0] < -3.*bLineSig_OK[0]  && amplitude[1] < -3.*bLineSig_OK[1] &&
	 amplitude[0] < -0.002 && amplitude[1] < -0.002 && timeMax[0] > 1 && timeMax[1] > 1 ){

      histos->Fill(time[1]-time[0]);
      h_integral_Ch0->Fill(-integral[0]);
      h_integral_Ch1->Fill(-integral[1]);

      h_amplitude_Ch0->Fill(amplitude[0]);
      h_amplitude_Ch1->Fill(amplitude[1]);
      if(//integral[2] > -15. && 
	 amplitude[2] < -3.*bLineSig_OK[2] && 
	 amplitude[2] < -0.002 && timeMax[2] > 1){
	h_signal->Fill((time[1]+time[0])*0.5-time[2]);
	h_signal2->Fill(time[1]-time[2]);
	h_signal3->Fill(time[2]-time[0]);
	h_integral_Ch2->Fill(-integral[2]);
	h_amplitude_Ch2->Fill(amplitude[2]);
	dtvsamp->Fill(amplitude[2],(time[1]+time[0])*0.5-time[2]);	
	p_dtvsamp->Fill(amplitude[2],(time[1]+time[0])*0.5-time[2]);	
	if(isnan((time[1]+time[0])*0.5-time[2]) == true) std::cout << " >>>>>>>>>>>>>>>>>>>>>>>>>> " << iCount << std::endl;
      }
    }
  }

  h_integral_Ch0->SetLineColor(kBlue);
  h_integral_Ch1->SetLineColor(kRed);
  h_integral_Ch2->SetLineColor(kGreen);

  TCanvas* cIntegral = new TCanvas();
  h_integral_Ch0->GetXaxis()->SetTitle("integral 36%-Max");
  h_integral_Ch0->Draw();
  h_integral_Ch1->Draw("same");
  h_integral_Ch2->Draw("same");

  h_amplitude_Ch0->SetLineColor(kBlue);
  h_amplitude_Ch1->SetLineColor(kRed);
  h_amplitude_Ch2->SetLineColor(kGreen);

  TCanvas* cAmplitude = new TCanvas();
  h_amplitude_Ch0->GetXaxis()->SetTitle("amp-Max");
  h_amplitude_Ch0->Draw();
  h_amplitude_Ch1->Draw("same");
  h_amplitude_Ch2->Draw("same");


  std::cout << "Now draw concidences" << std::endl;

  //*****************************************************************************************************************************************
  //***************************************************     DRAW CONCIDENCES     ************************************************************

  //0,1 -> external detectors
  //2 -> mcp detector
  
  TCanvas *coincidenze = new TCanvas();
  coincidenze -> Divide (4,1);

  //histos => difference
  histos->GetXaxis()->SetTitle("t0 - t1 (ns)");
  histos->SetTitle("difference");	
  coincidenze -> cd(1);
  histos->GetXaxis()->SetRangeUser(10.,12.);
  gGauss->SetParameters(0.1, 1., 11.1, 3.);	
  histos->Fit("gGauss","","",9.5,12.5);
  if(runFolder == "WFRun004"){
    histos->GetXaxis()->SetRangeUser(10.,12.);
    gGauss->SetParameters(0.1, 1., 10.1, 3.);	
    histos->Fit("gGauss","","",9.5,12.5);
  }
  if(runFolder == "WFRun008"){
    histos->GetXaxis()->SetRangeUser(10.,12.);
    gGauss->SetParameters(0.5, 1., 11.2, 3.);	
    histos->Fit("gGauss","","",8.,16.);
  }
  if(runFolder == "WFRun010"){
    histos->GetXaxis()->SetRangeUser(10.,12.);
    gGauss->SetParameters(0.5, 500., 11.1, 10.);	
    histos->Fit("gGauss","","",9.5,12.5);
  }
  if(runFolder == "WFRun012"){
    histos->GetXaxis()->SetRangeUser(10.,12.);
    gGauss->SetParameters(0.5,1.,11.,3.);
    histos->Fit("gGauss","","",8.,16.);
  }

  std::cout<<"sigma difference = "<<gGauss->GetParameter(0)<<std::endl;
  float sigma_diff = fabs(gGauss->GetParameter(0));
  float sigma_diffE = fabs(gGauss->GetParError(0));
  float mean_diff = gGauss->GetParameter(2);

  //h_signal2 => 1 - mcp
  h_signal2->GetXaxis()->SetRangeUser(6.5,9.5);
  h_signal2->GetXaxis()->SetTitle("t1 - t2 (ns)");
  h_signal2->SetTitle("1 - mcp");
  coincidenze -> cd(3);
  gGauss->SetParameters(0.1,50.,8,3.);
  h_signal2->Fit("gGauss","","",-25.,25.);
  std::cout<<"sigma 1 - mcp = "<<gGauss->GetParameter(0)<<std::endl;
  float sigma_t2dat1 = fabs(gGauss->GetParameter(0));
  float sigma_t2dat1E = fabs(gGauss->GetParError(0));
  float mean_t2dat1 = gGauss->GetParameter(2);

  //h_signal3 => mcp - 2
  h_signal3->GetXaxis()->SetRangeUser(1.,5.);
  h_signal3->GetXaxis()->SetTitle("t0 - t2 (ns)");
  h_signal3->SetTitle("mcp - 2");
  coincidenze -> cd(4);
  gGauss->SetParameters(0.15,10.,3.,3.);	
  h_signal3->Fit("gGauss","","",1.,5.);
  std::cout<<"sigma mcp - 2 = "<<gGauss->GetParameter(0)<<std::endl;
  float sigma_t2dat0 = fabs(gGauss->GetParameter(0));
  float sigma_t2dat0E = fabs(gGauss->GetParError(0));
  float mean_t2dat0 = gGauss->GetParameter(2);

  //h_signal => mean - mcp
  h_signal->GetXaxis()->SetRangeUser(1.,4.);
  h_signal->GetXaxis()->SetTitle("t_{mean} - t_{pl} (ns)");
  h_signal->SetTitle("mean - mcp");
  coincidenze -> cd(2);
  gGauss->SetParameters(0.05,10.,2.5,3.);	
  h_signal->Fit("gGauss","","",-20.,20.);
  if(runFolder == "WFRun005"){
    gGauss->SetParameters(0.15,10.,2.,0.5);	
    h_signal->Fit("gGauss","","",-20.,20.);
  }
  if(runFolder == "WFRun010"){
    //    h_signal->GetXaxis()->SetRangeUser(0.,4.);
    gGauss->SetParameters(0.05,200.,2.2,1.);	
    h_signal->Fit("gGauss","","",-40.,40.);
  }
  if(runFolder == "WFRun012"){
    h_signal->GetXaxis()->SetRangeUser(0.,4.);
    gGauss->SetParameters(0.15,20.,2.,0.5);	
    h_signal->Fit("gGauss","","",-40.,40.);
  }
  std::cout<<"sigma mean - mcp = "<<gGauss->GetParameter(0)<<std::endl;
  float sigma_mcp = fabs(gGauss->GetParameter(0));
  float sigma_mcpE = fabs(gGauss->GetParError(0));
  float mean_mcp = gGauss->GetParameter(2);

  chi->SetPoint(frac,frac,gGauss->GetChisquare()/gGauss->GetNDF());
  sigma->SetPoint(frac,frac,fabs(gGauss->GetParameter(0))); 

//   TLegend *legCutsExample = new TLegend(0.70,0.80,0.99,1.,NULL,"brNDC");
//   legCutsExample->SetTextFont(42);
//   legCutsExample->SetFillColor(kWhite);
//   legCutsExample->SetLineColor(kWhite);
//   legCutsExample->SetShadowColor(kWhite);

  float sigmaMCP = fabs(sigma_diff/sqrt(2.));
  float sigmaMCPE = fabs(sigma_diffE/sqrt(2.));
  std::cout << " sigmaMCP = " << sigmaMCP << "+/-" << sigmaMCP << std::endl;

  float sigmaMEAN = fabs(sigmaMCP/sqrt(2.));
  float sigmaMEANE = fabs(sigmaMCPE/sqrt(2.));
  std::cout << " sigmaMEAN = " << sigmaMEAN << "+/-" << sigmaMEANE << std::endl;

  float sigmaPLANA = sqrt( pow(sigma_mcp,2.) - pow(sigmaMEAN,2.) );
  float sigmaPLANAE = (pow(2.*sigma_mcp*sigma_mcpE,2.)+pow(2.*sigmaMEAN*sigmaMEANE,2.))/2./sqrt(pow(sigma_mcp,2.)-pow(sigmaMEAN,2.));
  std::cout << " sigmaPLANA = " << sigmaPLANA << "+/-" << sigmaPLANAE << std::endl;

  TLatex* latex;

  TCanvas* Results = new TCanvas();
  Results->Divide(2,1);
  Results->cd(1);
  //  legCutsExample->AddEntry(histos, "t0 - t1"(typeSet.at(0)).c_str(), "pl");
  latex = new TLatex(0.1,0.90, Form("#sigma fit: %1.2e #pm %1.2e",sigma_diff, sigma_diffE));
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.04);
  latex -> SetTextColor(1);
  latex->Draw("same");
  latex = new TLatex(0.1,0.80, Form("#sigma mcp: %1.2e #pm %1.2e",sigmaMCP, sigmaMCPE));
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.04);
  latex -> SetTextColor(1);
  histos->SetLineWidth(2);
  histos->DrawCopy();
  Results->cd(2);
  latex = new TLatex(0.1,0.90, Form("#sigma fit: %1.2e #pm %1.2e",sigma_mcp, sigma_mcpE));
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.04);
  latex -> SetTextColor(1);
  latex->Draw("same");
  latex = new TLatex(0.1,0.80, Form("#sigma plan: %1.2e #pm %1.2e",sigmaPLANA, sigmaPLANAE));
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.04);
  latex -> SetTextColor(1);
  h_signal->SetLineWidth(2);
  h_signal->DrawCopy();

  std::cout << "t0-t1: #sigma = " << sigma_diff << "+/-" << sigma_diffE 
	    << " => t_{std}: #sigma = " << sigmaMCP << "+/-" << sigmaMCPE << std::endl;

  std::cout << "t2-t1: #sigma = " << sigma_t2dat1 << "+/-" << sigma_t2dat1E 
	    << " => t2: #sigma = " << sqrt( pow(sigma_t2dat1,2.) - pow(sigmaMCP,2.) ) 
	    << "+/-" << (pow(2.*sigma_t2dat1*sigma_t2dat1E,2.)+pow(2.*sigmaMCP*sigmaMCPE,2.))/2./sqrt(pow(sigma_t2dat1,2.)-pow(sigmaMCP,2.)) 
	    << std::endl;

  std::cout << "t2-t0: #sigma = " << sigma_t2dat0 << "+/-" << sigma_t2dat0E 
	    << " => t2: #sigma = " << sqrt( pow(sigma_t2dat0,2.) - pow(sigmaMCP,2.) )
	    << "+/-" << (pow(2.*sigma_t2dat0*sigma_t2dat0E,2.)+pow(2.*sigmaMCP*sigmaMCPE,2.))/2./sqrt(pow(sigma_t2dat0,2.)-pow(sigmaMCP,2.))
	    << std::endl;

  std::cout << "t_{mean}-t_{planacon}: #sigma = " << sigma_mcp << "+/-" << sigma_mcpE
	    << " => t_{planacon}: #sigma = " << sigmaPLANA
	    << "+/-" << sigmaPLANAE << std::endl;

  //*****************************************************************************************************************************************
  //*************************************************         EVALUATING EFFICIENCY        *************************************************

  int doubleC = 0;
  int tripleC = 0;
  int doubleFake = 0;
  int tripleFake = 0;
  
  for (int r=0 ; r<t1->GetEntries()-2 ; ++r)
    {
      t1->GetEntry(r);
      if(fabs((time[1]-time[0]) - mean_diff) < 3*sigma_diff && 
	 amplitude[0] < -3.*bLineSig_OK[0] && amplitude[1] < -3.*bLineSig_OK[1] &&
	 amplitude[0] < -0.002 && amplitude[1] < -0.002 && timeMax[0] > 1 && timeMax[1] > 1) //double coincidence
	{					
	  doubleC++;
	  h_double->Fill(time[1]-time[0]);
	  if (fabs((time[1]+time[0])*0.5-time[2] - 5.) < 1.5 && 
	      amplitude[2] < -3.*bLineSig_OK[2] && amplitude[2] < -0.002 && timeMax[2] > 1)   tripleFake++;	 
	  if (fabs((time[1]+time[0])*0.5-time[2] - mean_mcp) < 3*sigma_mcp && 
	      amplitude[2] < -3.*bLineSig_OK[2] && amplitude[2] < -0.002 && timeMax[2] > 1)
	    {
	      tripleC++;
	      h_triple->Fill( (time[1]+time[0])*0.5-time[2] - mean_mcp + mean_diff);
	      //	      std::cout << " fabs(sigma_diff) - fabs(sigma_mcp) = " << fabs(sigma_diff) - fabs(sigma_mcp) << std::endl;
	    }	      	 
	}
      if (fabs(time[1]-time[0] - 13.5) < 1.5 && 
	  amplitude[0] < -3.*bLineSig_OK[0] && amplitude[1] < -3.*bLineSig_OK[1] && 
	  amplitude[0] < -0.002 && amplitude[1] < -0.002 && timeMax[0] > 1 && timeMax[1] > 1)
	doubleFake++;	
    }

  float rateDoubleFake = (float)doubleFake * 3. * sigma_diff/1.5;
  float rateTripleFake = (float)tripleFake * 3. * sigma_mcp/1.5;
  
  std::cout << "saved / triggered / fakeTriple / fakeDouble: " 
	    << tripleC << " / " << doubleC << " / " << rateTripleFake << " / " << rateDoubleFake << std::endl;
  float efficiency = (((float)tripleC - (float)rateTripleFake)/ ((float)doubleC - (float)rateDoubleFake))*100.;
  std::cout << " Planacon efficiency = " << efficiency<< "\%" << std::endl;


  TCanvas *p_scatter = new TCanvas();
  p_scatter -> cd();
  TF1 *retta = new TF1 ("retta","pol1",-0.02,0.);
  retta->SetParameters(0, 1);
  //  p_dtvsamp->Fit("retta","","",-0.02,-0.004);
  p_dtvsamp->Fit("retta","","", -0.1 , 0.);
  p_dtvsamp->Draw();
  p_dtvsamp->GetYaxis()->SetTitle("t_{mean} - t_{planacon} (ns)");
  p_dtvsamp->GetXaxis()->SetTitle("Amplitude (mV)");
  std::cout << "Amplitude walk: parametrised as y = " << retta->GetParameter(1) << " x + " << retta->GetParameter(0) << std::endl;

  TCanvas* scatter = new TCanvas();
  scatter->cd();
  dtvsamp->GetYaxis()->SetTitle("t_{mean} - t_{planacon} (ns)");
  dtvsamp->GetXaxis()->SetTitle("Amplitude (mV)");
  dtvsamp->Draw("colz");

  h_double->SetLineColor(kBlue);
  h_triple->SetLineColor(kGreen+1);

  h_double->SetLineWidth(2);
  h_triple->SetLineWidth(2);

  TCanvas* cC = new TCanvas();
  h_double->Draw();
  h_triple->Draw("same");



//*****************************************************************************************************************************************
//************************************************     DRAW SIGMA - CHI    ****************************************************************
	
  if(drawSigmaChi2){
    chi->SetMarkerStyle(20);
    sigma->SetMarkerStyle(20);
    chi->SetMarkerColor(kRed+1);
    sigma->SetMarkerColor(kGreen+1);

    TCanvas* tavola = new TCanvas();
    tavola->Divide(2,1);
    tavola->cd(1);	
    chi->Draw("AP");
    tavola->cd(2);	
    sigma->Draw("AP");

    TF1 *para = new TF1 ("para","pol2",0,60);	
    sigma->Fit("para","","",7.,55.);
    float minimum = para->GetMinimumX();
    std::cout << "Best fit for constant fraction at " << (int)minimum << "% amp Max" << std::endl;

    //*************************************************************************************************************************************	
    /*
    sort(time0.begin(),time0.end());
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

}
