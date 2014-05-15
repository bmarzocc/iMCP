/*************************************************************

    simple program to read bin files
    compile with --->  c++ -o makeAnalysis `root-config --cflags --glibs` makeAnalysis.cpp
         or with --->  c++ -o makeAnalysis makeAnalysis.cpp `root-config --cflags --glibs`
    run with --->  ./makeAnalysis file.root 1000 0 // usage: ./makeAnalysis file.root nBins saveWFHistos thresAmp[1] thresAmp[2] thresAmp[3]
    
*************************************************************/

#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm> 
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>
#include <utility>

#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"

#include "InitTree.h"

using namespace std;

const float BinToTime = 0.2;

float computeMean(std::vector<float> sample);
float computeDispersion(std::vector<float> sample);
void subtractBkg(TH1F* h, TH1F* h_sub, float& bkg);

int main(int argc, char** argv)
{
    gROOT->ProcessLine("#include <vector>");
    
    char* inputFile = argv[1];

    int nBins = atoi(argv[2]);
    int saveWFHistos = atoi(argv[3]);

    float thresAmp[4];

    thresAmp[0] = -1.;
    thresAmp[1] = 100.;
    thresAmp[2] = 100.;
    thresAmp[3] = 35.;
    
    if(argc >= 5) thresAmp[1] = atof(argv[4]);
    if(argc >= 6) thresAmp[2] = atof(argv[5]);
    if(argc >= 7) thresAmp[3] = atof(argv[6]);
    
    std::cout << "inputFile    = " << inputFile << std::endl;
    std::cout << "nBins        = " << nBins << std::endl;
    std::cout << "saveWFHistos = " << saveWFHistos << std::endl;
    std::cout << "thresAmp[1]  = " << thresAmp[1] << std::endl;
    std::cout << "thresAmp[2]  = " << thresAmp[2] << std::endl;
    std::cout << "thresAmp[3]  = " << thresAmp[3] << std::endl;
     
    TFile* input = TFile::Open(inputFile);
    
    TTree* nt = (TTree*)input->Get("nt");

    InitTree(nt);

    int event_output;
    int run_output;
    int isDouble_output;
    int isTriple_output;
    float ampMax_Trigger_output;
    float timeAmpMax_Trigger_output;
    float timeConstFrac_Trigger_output;
    float slopeConstFrac_Trigger_output;
    float chi2ConstFrac_Trigger_output;
    float ampMax_channel1_output;
    float timeAmpMax_channel1_output;
    float timeConstFrac_channel1_output;
    float slopeConstFrac_channel1_output;
    float chi2ConstFrac_channel1_output;
    float ampMax_channel2_output;
    float timeAmpMax_channel2_output;
    float timeConstFrac_channel2_output;
    float slopeConstFrac_channel2_output;
    float chi2ConstFrac_channel2_output;
    float ampMax_channel3_output;
    float timeAmpMax_channel3_output;
    float timeConstFrac_channel3_output;
    float slopeConstFrac_channel3_output;
    float chi2ConstFrac_channel3_output;

    float bkg_double_output = 0.;
    float bkg_triple_output = 0.;
    float bkg_double13_output = 0.;
    float bkg_double23_output = 0.;
    float sigma_double = 0.;
    float sigmaError_double = 0.;
    float mean_double = 0.;
    float sigma_triple = 0.;
    float sigmaError_triple = 0.;
    float mean_triple = 0.;
    float sigma_double13 = 0.;
    float sigmaError_double13 = 0.;
    float mean_double13 = 0.;
    float sigma_double23 = 0.;
    float sigmaError_double23 = 0.;
    float mean_double23 = 0.;
    float nTriple_output = 0.;
    float nTriple_Fake_output = 0.;
    float nDouble_output = 0.;
    float nDouble_Fake_output = 0.;
    
    TTree *nt_output = new TTree("nt","nt");
    nt_output->SetDirectory(0);

    nt_output->Branch("event",&event_output,"event/I"); 
    nt_output->Branch("run",&run_output,"run/I"); 
    nt_output->Branch("isDouble",&isDouble_output,"isDouble/I"); 
    nt_output->Branch("isTriple",&isTriple_output,"isTriple/I"); 
    nt_output->Branch("ampMax_Trigger",&ampMax_Trigger_output,"ampMax_Trigger/F"); 
    nt_output->Branch("timeAmpMax_Trigger",&timeAmpMax_Trigger_output,"timeAmpMax_Trigger/F"); 
    nt_output->Branch("timeConstFrac_Trigger",&timeConstFrac_Trigger_output,"timeConstFrac_Trigger/F"); 
    nt_output->Branch("slopeConstFrac_Trigger",&slopeConstFrac_Trigger_output,"slopeConstFrac_Trigger/F"); 
    nt_output->Branch("chi2ConstFrac_Trigger",&chi2ConstFrac_Trigger_output,"chi2ConstFrac_Trigger/F"); 
    nt_output->Branch("ampMax_channel1",&ampMax_channel1_output,"ampMax_channel1/F"); 
    nt_output->Branch("timeAmpMax_channel1",&timeAmpMax_channel1_output,"timeAmpMax_channel1/F"); 
    nt_output->Branch("timeConstFrac_channel1",&timeConstFrac_channel1_output,"timeConstFrac_channel1/F"); 
    nt_output->Branch("slopeConstFrac_channel1",&slopeConstFrac_channel1_output,"slopeConstFrac_channel1/F"); 
    nt_output->Branch("chi2ConstFrac_channel1",&chi2ConstFrac_channel1_output,"chi2ConstFrac_channel1/F"); 
    nt_output->Branch("ampMax_channel2",&ampMax_channel2_output,"ampMax_channel2/F"); 
    nt_output->Branch("timeAmpMax_channel2",&timeAmpMax_channel2_output,"timeAmpMax_channel2/F"); 
    nt_output->Branch("timeConstFrac_channel2",&timeConstFrac_channel2_output,"timeConstFrac_channel2/F"); 
    nt_output->Branch("slopeConstFrac_channel2",&slopeConstFrac_channel2_output,"slopeConstFrac_channel2/F"); 
    nt_output->Branch("chi2ConstFrac_channel2",&chi2ConstFrac_channel2_output,"chi2ConstFrac_channel2/F"); 
    nt_output->Branch("ampMax_channel3",&ampMax_channel3_output,"ampMax_channel3/F"); 
    nt_output->Branch("timeAmpMax_channel3",&timeAmpMax_channel3_output,"timeAmpMax_channel3/F"); 
    nt_output->Branch("timeConstFrac_channel3",&timeConstFrac_channel3_output,"timeConstFrac_channel3/F"); 
    nt_output->Branch("slopeConstFrac_channel3",&slopeConstFrac_channel3_output,"slopeConstFrac_channel3/F"); 
    nt_output->Branch("chi2ConstFrac_channel3",&chi2ConstFrac_channel3_output,"chi2ConstFrac_channel3/F"); 

    nt_output->Branch("bkg_double_output",&bkg_double_output,"bkg_double_output/F"); 
    nt_output->Branch("bkg_triple_output",&bkg_triple_output,"bkg_triple_output/F"); 
    nt_output->Branch("bkg_double13_output",&bkg_double13_output,"bkg_double13_output/F"); 
    nt_output->Branch("bkg_double23_output",&bkg_double23_output,"bkg_double23_output/F"); 
    nt_output->Branch("sigma_double",&sigma_double,"sigma_double/F"); 
    nt_output->Branch("sigmaError_double",&sigmaError_double,"sigmaError_double/F"); 
    nt_output->Branch("mean_double",&mean_double,"mean_double/F"); 
    nt_output->Branch("sigma_triple",&sigma_triple,"sigma_triple/F"); 
    nt_output->Branch("sigmaError_triple",&sigmaError_triple,"sigmaError_triple/F"); 
    nt_output->Branch("mean_triple",&mean_triple,"mean_triple/F"); 
    nt_output->Branch("sigma_double13",&sigma_double13,"sigma_double13/F"); 
    nt_output->Branch("sigmaError_double13",&sigmaError_double13,"sigmaError_double13/F"); 
    nt_output->Branch("mean_double13",&mean_double13,"mean_double13/F"); 
    nt_output->Branch("sigma_double23",&sigma_double23,"sigma_double23/F"); 
    nt_output->Branch("sigmaError_double23",&sigmaError_double23,"sigmaError_double23/F"); 
    nt_output->Branch("mean_double23",&mean_double23,"mean_double23/F"); 
    nt_output->Branch("nTriple_output",&nTriple_output,"nTriple_output/F"); 
    nt_output->Branch("nTriple_Fake_output",&nTriple_Fake_output,"nTriple_Fake_output/F"); 
    nt_output->Branch("nDouble_output",&nDouble_output,"nDouble_output/F"); 
    nt_output->Branch("nDouble_Fake_output",&nDouble_Fake_output,"nDouble_Fake_output/F"); 


    TF1* gGauss = new TF1("gGauss", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))");
    
    TH1F* h_double = new TH1F("h_double","h_double",nBins,0.,100.);
    TH1F* h_double_sub = new TH1F("h_double_sub","h_double_sub",nBins,0.,100.);
    TH1F* h_triple = new TH1F("h_triple","h_triple",nBins,0.,100.);
    TH1F* h_triple_sub = new TH1F("h_triple_sub","h_triple_sub",nBins,0.,100.);
    TH1F* h_double13 = new TH1F("h_double13","h_double13",nBins,0.,100.);
    TH1F* h_double13_sub = new TH1F("h_double13_sub","h_double13_sub",nBins,0.,100.);
    TH1F* h_double23 = new TH1F("h_double23","h_double23",nBins,0.,100.);
    TH1F* h_double23_sub = new TH1F("h_double23_sub","h_double23_sub",nBins,0.,100.);

    std::map<int,std::map<int,TH1F*> > h_WF_Trigger;
    std::map<int,std::map<int,TH1F*> > h_WF_channel1; 
    std::map<int,std::map<int,TH1F*> > h_WF_channel2;
    std::map<int,std::map<int,TH1F*> > h_WF_channel3;

    std::vector<std::pair<int,int> > vec_run_event;
    std::pair<int,int> tmp_pair;

    char histoName[200];
   
    for(int ientry = 0; ientry < nt->GetEntries(); ientry++){
        if(ientry%100==0) std::cout<<"--- Reading entry = "<< ientry << " for filling the histos..." <<std::endl;
        nt->GetEntry(ientry);

        if(ampMax_channel1 <= thresAmp[1]) continue;
        if(ampMax_channel2 <= thresAmp[2]) continue;

        if(timeAmpMax_channel1/BinToTime < 2 || timeAmpMax_channel1/BinToTime > waveForm_channel1->size()-2) continue;
        if(timeAmpMax_channel2/BinToTime < 2 || timeAmpMax_channel2/BinToTime > waveForm_channel2->size()-2) continue;

        if((timeConstFrac_channel2-timeConstFrac_channel1) < 0.) continue;

        h_double->Fill(timeAmpMax_channel2-timeAmpMax_channel1);
        
        if(ampMax_channel3 > thresAmp[3] && timeAmpMax_channel3/BinToTime >= 2 && timeAmpMax_channel3/BinToTime <= waveForm_channel3->size()-2){
      
           h_triple->Fill((timeConstFrac_channel2+timeConstFrac_channel1)/2.-timeConstFrac_channel3);
           h_double13->Fill(timeAmpMax_channel1-timeAmpMax_channel3);
           h_double23->Fill(timeAmpMax_channel2-timeAmpMax_channel3);
        }
        
        int nSize = (int)waveForm_Trigger->size();

        if(saveWFHistos == 1){

            sprintf(histoName, "WF_Trigger_%d_%d",run,event);  
            h_WF_Trigger[run][event] = new TH1F(histoName,histoName,nSize,0.,nSize*BinToTime);

            for(unsigned int ii = 0; ii < waveForm_Trigger->size(); ii++)
                h_WF_Trigger[run][event]->SetBinContent(ii+1,waveForm_Trigger->at(ii));   

            sprintf(histoName, "WF_channel1_%d_%d",run,event);  
            h_WF_channel1[run][event] = new TH1F(histoName,histoName,nSize,0.,nSize*BinToTime);

            for(unsigned int ii = 0; ii < waveForm_channel1->size(); ii++)
                h_WF_channel1[run][event]->SetBinContent(ii+1,waveForm_channel1->at(ii));   
            
            sprintf(histoName, "WF_channel2_%d_%d",run,event);  
            h_WF_channel2[run][event] = new TH1F(histoName,histoName,nSize,0.,nSize*BinToTime);

            for(unsigned int ii = 0; ii < waveForm_channel2->size(); ii++)
                h_WF_channel2[run][event]->SetBinContent(ii+1,waveForm_channel2->at(ii));   

            sprintf(histoName, "WF_channel3_%d_%d",run,event);  
            h_WF_channel3[run][event] = new TH1F(histoName,histoName,nSize,0.,nSize*BinToTime);

            for(unsigned int ii = 0; ii < waveForm_channel3->size(); ii++)
                h_WF_channel3[run][event]->SetBinContent(ii+1,waveForm_channel3->at(ii));   
        }
        
        tmp_pair = std::make_pair(run,event);
        vec_run_event.push_back(tmp_pair);

    }

    float bkg_double = 0.;
    float bkg_triple = 0.;
    float bkg_double13 = 0.;
    float bkg_double23 = 0.;

    subtractBkg(h_double,h_double_sub,bkg_double);
    subtractBkg(h_triple,h_triple_sub,bkg_triple);
    subtractBkg(h_double13,h_double13_sub,bkg_double13);
    subtractBkg(h_double23,h_double23_sub,bkg_double23);

    bkg_double_output = bkg_double;
    bkg_triple_output = bkg_triple;
    bkg_double13_output = bkg_double13;
    bkg_double23_output = bkg_double23;
    
    float x_min = h_double_sub->GetBinCenter(h_double_sub->GetMaximumBin())-2*h_double_sub->GetRMS();
    float x_max = h_double_sub->GetBinCenter(h_double_sub->GetMaximumBin())+2*h_double_sub->GetRMS();
    gGauss->SetParameter(0,h_double_sub->Integral());
    gGauss->SetParameter(1,h_double_sub->GetMean());
    gGauss->SetParameter(2,h_double_sub->GetRMS());
    h_double_sub->Fit("gGauss","","",x_min*0.98,x_max*1.02);

    std::cout<<"sigma difference = "<<gGauss->GetParameter(2)<<std::endl;
    float sigma_diff = fabs(gGauss->GetParameter(2));
    float sigma_diffE = fabs(gGauss->GetParError(2));
    float mean_diff = gGauss->GetParameter(1);

    sigma_double = sigma_diff;
    sigmaError_double = sigma_diffE; 
    mean_double = mean_diff;

    x_min = h_triple_sub->GetBinCenter(h_triple_sub->GetMaximumBin())-2*h_triple_sub->GetRMS();
    x_max = h_triple_sub->GetBinCenter(h_triple_sub->GetMaximumBin())+2*h_triple_sub->GetRMS();
    gGauss->SetParameter(0,h_triple_sub->Integral());
    gGauss->SetParameter(1,h_triple_sub->GetMean());
    gGauss->SetParameter(2,h_triple_sub->GetRMS());
    h_triple_sub->Fit("gGauss","","",x_min*0.98,x_max*1.02);
    float sigma_mcp = fabs(gGauss->GetParameter(2));
    float sigma_mcpE = fabs(gGauss->GetParError(2));
    float mean_mcp = gGauss->GetParameter(1);

    sigma_triple = sigma_mcp;
    sigmaError_triple = sigma_mcpE; 
    mean_triple = mean_mcp;

    x_min = h_double13_sub->GetBinCenter(h_double13_sub->GetMaximumBin())-2*h_double13_sub->GetRMS();
    x_max = h_double13_sub->GetBinCenter(h_double13_sub->GetMaximumBin())+2*h_double13_sub->GetRMS();
    gGauss->SetParameter(0,h_double13_sub->Integral());
    gGauss->SetParameter(1,h_double13_sub->GetMean());
    gGauss->SetParameter(2,h_double13_sub->GetRMS());
    h_double13_sub->Fit("gGauss","","",x_min*0.98,x_max*1.02);
    float sigma_t2dat0 = fabs(gGauss->GetParameter(2));
    float sigma_t2dat0E = fabs(gGauss->GetParError(2));
    float mean_t2dat0 = gGauss->GetParameter(1);

    sigma_double13 = sigma_t2dat0;
    sigmaError_double13 = sigma_t2dat0E; 
    mean_double13 = mean_t2dat0;

    x_min = h_double23_sub->GetBinCenter(h_double23_sub->GetMaximumBin())-2*h_double23_sub->GetRMS();
    x_max = h_double23_sub->GetBinCenter(h_double23_sub->GetMaximumBin())+2*h_double23_sub->GetRMS();
    gGauss->SetParameter(0,h_double23_sub->Integral());
    gGauss->SetParameter(1,h_double23_sub->GetMean());
    gGauss->SetParameter(2,h_double23_sub->GetRMS());
    h_double23_sub->Fit("gGauss","","",x_min*0.98,x_max*1.02);
    float sigma_t2dat1 = fabs(gGauss->GetParameter(2));
    float sigma_t2dat1E = fabs(gGauss->GetParError(2));
    float mean_t2dat1 = gGauss->GetParameter(1);

    sigma_double23 = sigma_t2dat1;
    sigmaError_double23 = sigma_t2dat1E; 
    mean_double23 = mean_t2dat1;

    std::cout << std::endl;
    
    float sigmaMCP = fabs(sigma_diff/sqrt(2.));
    float sigmaMCPE = fabs(sigma_diffE/sqrt(2.));
    std::cout << " sigmaMCP = " << sigmaMCP*1000 << " +/- " << sigmaMCPE*1000 << " ps" << std::endl;

    float sigmaMEAN = fabs(sigmaMCP/sqrt(2.));
    float sigmaMEANE = fabs(sigmaMCPE/sqrt(2.));
    std::cout << " sigmaMEAN = " << sigmaMEAN*1000 << " +/- " << sigmaMEANE*1000 << " ps" << std::endl;

    float sigmaPLANA = sqrt( pow(sigma_mcp,2.) - pow(sigmaMEAN,2.) );
    float sigmaPLANAE = (pow(2.*sigma_mcp*sigma_mcpE,2.)+pow(2.*sigmaMEAN*sigmaMEANE,2.))/2./sqrt(pow(sigma_mcp,2.)-pow(sigmaMEAN,2.));
    std::cout << " sigmaPLANA = " << sigmaPLANA*1000 << " +/- " << sigmaPLANAE*1000 << " ps" << std::endl;

    std::cout << std::endl;

    int nDouble = 0;
    float nDouble_fake = 0;
    int nTriple = 0;
    float nTriple_fake = 0;
   
    for(int ientry = 0; ientry < nt->GetEntries(); ientry++){
        if(ientry%100==0) std::cout<<"--- Reading entry = "<< ientry  << " for computing efficiencies..."<<std::endl;
        nt->GetEntry(ientry);
        
        bool isDouble = false;
        bool isTriple = false;

        if(ampMax_channel1 <= thresAmp[1]) continue;
        if(ampMax_channel2 <= thresAmp[2]) continue;

        if(timeAmpMax_channel1/BinToTime < 2 || timeAmpMax_channel1/BinToTime > waveForm_channel1->size()-2) continue;
        if(timeAmpMax_channel2/BinToTime < 2 || timeAmpMax_channel2/BinToTime > waveForm_channel2->size()-2) continue;

        if(fabs(timeConstFrac_channel2-timeConstFrac_channel1 - mean_diff) >= 3*sigma_diff) continue;

        nDouble++;
       
        if(ampMax_channel3 <= thresAmp[3]) continue;
        if(timeAmpMax_channel3/BinToTime < 2 || timeAmpMax_channel3/BinToTime > waveForm_channel3->size()-2) continue;

        if(fabs((timeConstFrac_channel2+timeConstFrac_channel1)/2.- timeConstFrac_channel3 - mean_mcp) >= 3*sigma_mcp) continue;

        nTriple++;
        
    }

    std::cout << std::endl;

    nTriple_fake = 6*sigma_mcp*bkg_triple/h_triple->GetBinWidth(1);
    nDouble_fake = 6*sigma_diff*bkg_double/h_double->GetBinWidth(1);
    float eff = float(nTriple-nTriple_fake)/float(nDouble-nDouble_fake);
    float eff_err = sqrt(eff*(1-eff)/float(nDouble-nDouble_fake));

    std::cout << "nTriple = " << float(nTriple) << " - nTriple_fake = " << nTriple_fake << std::endl;
    std::cout << "nDouble = " << float(nDouble) << " - nDouble_fake = " << nDouble_fake << std::endl;
    std::cout << "Efficiencies = " << eff*100 << " +/- "<< eff_err*100 << " %" << std::endl;

    nTriple_output = nTriple;
    nTriple_Fake_output = nTriple_fake;
    nDouble_output = nDouble;
    nDouble_Fake_output = nDouble_fake;

    std::cout << std::endl;

    for(int ientry = 0; ientry < nt->GetEntries(); ientry++){
        if(ientry%100==0) std::cout<<"--- Reading entry = "<< ientry << " for filling the tree..." <<std::endl;
        nt->GetEntry(ientry);

        event_output = event;
        run_output = run; 

        ampMax_Trigger_output = ampMax_Trigger;
        ampMax_channel1_output = ampMax_channel1;
        ampMax_channel2_output = ampMax_channel2;
        ampMax_channel3_output = ampMax_channel3;
        
        timeAmpMax_Trigger_output = timeAmpMax_Trigger;
        timeAmpMax_channel1_output = timeAmpMax_channel1;
        timeAmpMax_channel2_output = timeAmpMax_channel2;
        timeAmpMax_channel3_output = timeAmpMax_channel3;
        
        timeConstFrac_Trigger_output = timeConstFrac_Trigger;          
        timeConstFrac_channel1_output = timeConstFrac_channel1;
        timeConstFrac_channel2_output = timeConstFrac_channel2;
        timeConstFrac_channel3_output = timeConstFrac_channel3;
        
        slopeConstFrac_Trigger_output = slopeConstFrac_Trigger;          
        slopeConstFrac_channel1_output = slopeConstFrac_channel1;
        slopeConstFrac_channel2_output = slopeConstFrac_channel2;
        slopeConstFrac_channel3_output = slopeConstFrac_channel3;
        
        chi2ConstFrac_Trigger_output = chi2ConstFrac_Trigger;          
        chi2ConstFrac_channel1_output = chi2ConstFrac_channel1;
        chi2ConstFrac_channel2_output = chi2ConstFrac_channel2;
        chi2ConstFrac_channel3_output = chi2ConstFrac_channel3;

        isDouble_output = 0;
        isTriple_output = 0;
        
        int isDouble = 0;
        int isTriple = 0;

        if(ampMax_channel1 <= thresAmp[1]) continue;
        if(ampMax_channel2 <= thresAmp[2]) continue;

        if(timeAmpMax_channel1/BinToTime < 2 || timeAmpMax_channel1/BinToTime > waveForm_channel1->size()-2) continue;
        if(timeAmpMax_channel2/BinToTime < 2 || timeAmpMax_channel2/BinToTime > waveForm_channel2->size()-2) continue;

        if((timeConstFrac_channel2-timeConstFrac_channel1) < 0.) continue;

        isDouble = 1;
        h_double->Fill(timeAmpMax_channel2-timeAmpMax_channel1);
        
        if(ampMax_channel3 > thresAmp[3] && timeAmpMax_channel3/BinToTime >= 2 && timeAmpMax_channel3/BinToTime <= waveForm_channel3->size()-2){
      
           isTriple = 1;
           h_triple->Fill((timeConstFrac_channel2+timeConstFrac_channel1)/2.-timeConstFrac_channel3);
           h_double13->Fill(timeAmpMax_channel1-timeAmpMax_channel3);
           h_double23->Fill(timeAmpMax_channel2-timeAmpMax_channel3);
        }

        isDouble_output = isDouble;
        isTriple_output = isTriple;

        nt_output->Fill();
        
    }

    std::cout << std::endl;
    
    TFile* f1 = new TFile(("analyzed_"+std::string(inputFile)).c_str(),"RECREATE");
    f1->cd();

    nt_output->Write();
    h_double->Write();
    h_double_sub->Write();
    h_triple->Write();
    h_triple_sub->Write();
    h_double13->Write();
    h_double13_sub->Write();
    h_double23->Write();
    h_double23_sub->Write();

    if(saveWFHistos == 1){

      for(unsigned int ii = 0; ii < vec_run_event.size(); ii++)
        {

          h_WF_Trigger[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          h_WF_channel1[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          h_WF_channel2[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          h_WF_channel3[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
            
        }
    }

    f1->Close();
}


float computeMean(std::vector<float> sample)
{
  float mean = 0.;
  
  for(unsigned int ii = 0; ii < sample.size(); ii++)
      mean = mean + sample.at(ii);

  return mean/sample.size();
}

float computeDispersion(std::vector<float> sample)
{
  float error = 0.; 
  float mean = computeMean(sample);

  for(unsigned int ii = 0; ii < sample.size(); ii++)
      error = error + (sample.at(ii)-mean)*(sample.at(ii)-mean);  

  //return sqrt(error/(sample.size()-1))/sqrt(sample.size()); 
  return sqrt(error/(sample.size()-1)); 
}

void subtractBkg(TH1F* h, TH1F* h_sub, float& bkg)
{

    std::vector<float> tmp;
    for(int ii = 1; ii <= h->GetNbinsX(); ii++)
        if((h->GetBinCenter(ii) < (h->GetBinCenter(h->GetMaximumBin())-2*h->GetRMS()) || h->GetBinCenter(ii) > (h->GetBinCenter(h->GetMaximumBin())+2*h->GetRMS())) && h->GetBinContent(ii) != 0.) tmp.push_back(h->GetBinContent(ii));

    bkg = computeMean(tmp);
    for(int ii = 1; ii <= h->GetNbinsX(); ii++)
        if(h->GetBinContent(ii)-bkg >= 0.) h_sub->SetBinContent(ii,h->GetBinContent(ii)-bkg);
        else h_sub->SetBinContent(ii,0.);
}

