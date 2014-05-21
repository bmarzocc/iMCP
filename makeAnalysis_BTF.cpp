/*************************************************************

    simple program to read bin files
    compile with --->  c++ -o makeAnalysis_BTF `root-config --cflags --glibs` makeAnalysis_BTF.cpp
         or with --->  c++ -o makeAnalysis_BTF makeAnalysis_BTF.cpp `root-config --cflags --glibs`
    run with --->  ./makeAnalysis file.root 1000 0 // usage: ./makeAnalysis file.root nBins saveWFHistos thresAmp[0] thresAmp[1] thresAmp[2]
    
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
const int nSize = 1024;

const int timingChannels[3] = {0,5,3};

float computeMean(std::vector<float> sample);
float computeDispersion(std::vector<float> sample);
void subtractBkg(TH1F* h, TH1F* h_sub, float& bkg);
void ChooseChannels(std::map<int,int> map_channels,float& ampMax1,float& ampMax2,float& ampMax3,float& timeAmpMax1,float& timeAmpMax2,float& timeAmpMax3,float& timeConstFrac1,float& timeConstFrac2,float& timeConstFrac3);

int main(int argc, char** argv)
{
    gROOT->ProcessLine("#include <vector>");
    
    char* inputFile = argv[1];

    int nBins = atoi(argv[2]);
    int saveWFHistos = atoi(argv[3]);

    float thresAmp[3];

    thresAmp[0] = 150.;
    thresAmp[1] = 150.;
    thresAmp[2] = 0.;
    
    if(argc >= 5) thresAmp[0] = atof(argv[4]);
    if(argc >= 6) thresAmp[1] = atof(argv[5]);
    if(argc >= 7) thresAmp[2] = atof(argv[6]);
    
    std::cout << "inputFile    = " << inputFile << std::endl;
    std::cout << "nBins        = " << nBins << std::endl;
    std::cout << "saveWFHistos = " << saveWFHistos << std::endl;
    std::cout << "thresAmp[0]  = " << thresAmp[0] << std::endl;
    std::cout << "thresAmp[1]  = " << thresAmp[1] << std::endl;
    std::cout << "thresAmp[2]  = " << thresAmp[2] << std::endl;
     
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
    float sigmaMCP_output = 0.;
    float sigmaMEAN_output = 0.;
    float sigmaPLANA_output = 0.;
    float nTriple_output = 0.;
    float nTriple_Fake_output = 0.;
    float nDouble_output = 0.;
    float nDouble_Fake_output = 0.;
    float Efficiency = 0.;
    float EfficiencyError = 0.;
    
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
    nt_output->Branch("sigmaMCP",&sigmaMCP_output,"sigmaMCP/F"); 
    nt_output->Branch("sigmaMEAN",&sigmaMEAN_output,"sigmaMEAN/F"); 
    nt_output->Branch("sigmaPLANA",&sigmaPLANA_output,"sigmaPLANA/F"); 
    nt_output->Branch("nTriple",&nTriple_output,"nTriple/F"); 
    nt_output->Branch("nTriple_Fake",&nTriple_Fake_output,"nTriple_Fake/F"); 
    nt_output->Branch("nDouble",&nDouble_output,"nDouble/F"); 
    nt_output->Branch("nDouble_Fake",&nDouble_Fake_output,"nDouble_Fake/F"); 
    nt_output->Branch("Efficiency",&Efficiency,"Efficiency/F"); 
    nt_output->Branch("EfficiencyError",&EfficiencyError,"EfficiencyError/F"); 


    TF1* gGauss = new TF1("gGauss", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))");
    
    TH1F* h_double = new TH1F("h_double","h_double",nBins,-500.,500.);
    TH1F* h_triple = new TH1F("h_triple","h_triple",nBins,-500.,500.);
    TH1F* h_double13 = new TH1F("h_double13","h_double13",nBins,-500.,500.);
    TH1F* h_double23 = new TH1F("h_double23","h_double23",nBins,-500.,500.);

    std::map<int,std::map<int,TH1F*> > h_WF_Trigger;
    std::map<int,std::map<int,TH1F*> > h_WF_channel1; 
    std::map<int,std::map<int,TH1F*> > h_WF_channel2;
    std::map<int,std::map<int,TH1F*> > h_WF_channel3;

    std::vector<std::pair<int,int> > vec_run_event;
    std::pair<int,int> tmp_pair;

    char histoName[200];

    float ampMax1 = 0.;
    float ampMax2 = 0.;
    float ampMax3 = 0.;
    float timeAmpMax1 = 0.;
    float timeAmpMax2 = 0.;
    float timeAmpMax3 = 0.;
    float timeConstFrac1 = 0.;
    float timeConstFrac2 = 0.;
    float timeConstFrac3 = 0.;

    std::map<int,int> map_channels;
    map_channels[0] = timingChannels[0];
    map_channels[1] = timingChannels[1];
    map_channels[2] = timingChannels[2];

    float nDouble_simple = 0.;
    float nTriple_simple = 0.;
   
    for(int ientry = 0; ientry < nt->GetEntries(); ientry++){
        if(ientry%100==0) std::cout<<"--- Reading entry = "<< ientry << " for filling the histos..." <<std::endl;
        nt->GetEntry(ientry);
        
        ChooseChannels(map_channels,ampMax1,ampMax2,ampMax3,timeAmpMax1,timeAmpMax2,timeAmpMax3,timeConstFrac1,timeConstFrac2,timeConstFrac3);
        
        if(ampMax1 <= thresAmp[0]) continue;
        if(ampMax2 <= thresAmp[1]) continue;

        h_double->Fill(timeConstFrac2-timeConstFrac1);
        nDouble_simple++;
        
        if(ampMax3 > thresAmp[2]){
      
           h_triple->Fill((timeConstFrac2+timeConstFrac1)/2.-timeConstFrac3);
           h_double13->Fill(timeConstFrac1-timeConstFrac3);
           h_double23->Fill(timeConstFrac2-timeConstFrac3);
           nTriple_simple++;
        }
        
        if(saveWFHistos == 1){

            sprintf(histoName, "WF_Trigger_%d_%d",run,event);  
            h_WF_Trigger[run][event] = new TH1F(histoName,histoName,nSize,0.,nSize*BinToTime);

            for(unsigned int ii = 0; ii < waveForm_Trigger->size(); ii++)
                h_WF_Trigger[run][event]->SetBinContent(ii+1,waveForm_Trigger->at(ii));   

            sprintf(histoName, "WF_channel1_%d_%d",run,event);  
            h_WF_channel1[run][event] = new TH1F(histoName,histoName,nSize,0.,nSize*BinToTime);

            if(timingChannels[0] == 0){
             for(unsigned int ii = 0; ii < waveForm_channel0->size(); ii++)
                 h_WF_channel1[run][event]->SetBinContent(ii+1,waveForm_channel0->at(ii));   
            }
            else if(timingChannels[0] == 1){
             for(unsigned int ii = 0; ii < waveForm_channel1->size(); ii++)
                 h_WF_channel1[run][event]->SetBinContent(ii+1,waveForm_channel1->at(ii));   
            }
            else if(timingChannels[0] == 2){
             for(unsigned int ii = 0; ii < waveForm_channel2->size(); ii++)
                 h_WF_channel1[run][event]->SetBinContent(ii+1,waveForm_channel2->at(ii));   
            }
            else if(timingChannels[0] == 3){
             for(unsigned int ii = 0; ii < waveForm_channel3->size(); ii++)
                 h_WF_channel1[run][event]->SetBinContent(ii+1,waveForm_channel3->at(ii));   
            } 
            else if(timingChannels[0] == 4){
             for(unsigned int ii = 0; ii < waveForm_channel4->size(); ii++)
                 h_WF_channel1[run][event]->SetBinContent(ii+1,waveForm_channel4->at(ii));   
            }
            else if(timingChannels[0] == 5){
             for(unsigned int ii = 0; ii < waveForm_channel5->size(); ii++)
                 h_WF_channel1[run][event]->SetBinContent(ii+1,waveForm_channel5->at(ii));   
            }
            else if(timingChannels[0] == 6){
             for(unsigned int ii = 0; ii < waveForm_channel6->size(); ii++)
                 h_WF_channel1[run][event]->SetBinContent(ii+1,waveForm_channel6->at(ii));   
            }
            else if(timingChannels[0] == 7){
             for(unsigned int ii = 0; ii < waveForm_channel7->size(); ii++)
                 h_WF_channel1[run][event]->SetBinContent(ii+1,waveForm_channel7->at(ii));   
            }
            else if(timingChannels[0] == 8){
             for(unsigned int ii = 0; ii < waveForm_channel8->size(); ii++)
                 h_WF_channel1[run][event]->SetBinContent(ii+1,waveForm_channel8->at(ii));   
            }
             
            sprintf(histoName, "WF_channel2_%d_%d",run,event);  
            h_WF_channel2[run][event] = new TH1F(histoName,histoName,nSize,0.,nSize*BinToTime);

            if(timingChannels[1] == 0){
             for(unsigned int ii = 0; ii < waveForm_channel0->size(); ii++)
                 h_WF_channel2[run][event]->SetBinContent(ii+1,waveForm_channel0->at(ii));   
            }
            else if(timingChannels[1] == 1){
             for(unsigned int ii = 0; ii < waveForm_channel1->size(); ii++)
                 h_WF_channel2[run][event]->SetBinContent(ii+1,waveForm_channel1->at(ii));   
            }
            else if(timingChannels[1] == 2){
             for(unsigned int ii = 0; ii < waveForm_channel2->size(); ii++)
                 h_WF_channel2[run][event]->SetBinContent(ii+1,waveForm_channel2->at(ii));   
            }
            else if(timingChannels[1] == 3){
             for(unsigned int ii = 0; ii < waveForm_channel3->size(); ii++)
                 h_WF_channel2[run][event]->SetBinContent(ii+1,waveForm_channel3->at(ii));   
            } 
            else if(timingChannels[1] == 4){
             for(unsigned int ii = 0; ii < waveForm_channel4->size(); ii++)
                 h_WF_channel2[run][event]->SetBinContent(ii+1,waveForm_channel4->at(ii));   
            }
            else if(timingChannels[1] == 5){
             for(unsigned int ii = 0; ii < waveForm_channel5->size(); ii++)
                 h_WF_channel2[run][event]->SetBinContent(ii+1,waveForm_channel5->at(ii));   
            }
            else if(timingChannels[1] == 6){
             for(unsigned int ii = 0; ii < waveForm_channel6->size(); ii++)
                 h_WF_channel2[run][event]->SetBinContent(ii+1,waveForm_channel6->at(ii));   
            }
            else if(timingChannels[1] == 7){
             for(unsigned int ii = 0; ii < waveForm_channel7->size(); ii++)
                 h_WF_channel2[run][event]->SetBinContent(ii+1,waveForm_channel7->at(ii));   
            }
            else if(timingChannels[1] == 8){
             for(unsigned int ii = 0; ii < waveForm_channel8->size(); ii++)
                 h_WF_channel2[run][event]->SetBinContent(ii+1,waveForm_channel8->at(ii));   
            } 

            sprintf(histoName, "WF_channel3_%d_%d",run,event);  
            h_WF_channel3[run][event] = new TH1F(histoName,histoName,nSize,0.,nSize*BinToTime);

            if(timingChannels[2] == 0){
             for(unsigned int ii = 0; ii < waveForm_channel0->size(); ii++)
                 h_WF_channel3[run][event]->SetBinContent(ii+1,waveForm_channel0->at(ii));   
            }
            else if(timingChannels[2] == 1){
             for(unsigned int ii = 0; ii < waveForm_channel1->size(); ii++)
                 h_WF_channel3[run][event]->SetBinContent(ii+1,waveForm_channel1->at(ii));   
            }
            else if(timingChannels[2] == 2){
             for(unsigned int ii = 0; ii < waveForm_channel2->size(); ii++)
                 h_WF_channel3[run][event]->SetBinContent(ii+1,waveForm_channel2->at(ii));   
            }
            else if(timingChannels[2] == 3){
             for(unsigned int ii = 0; ii < waveForm_channel3->size(); ii++)
                 h_WF_channel3[run][event]->SetBinContent(ii+1,waveForm_channel3->at(ii));   
            } 
            else if(timingChannels[2] == 4){
             for(unsigned int ii = 0; ii < waveForm_channel4->size(); ii++)
                 h_WF_channel3[run][event]->SetBinContent(ii+1,waveForm_channel4->at(ii));   
            }
            else if(timingChannels[2] == 5){
             for(unsigned int ii = 0; ii < waveForm_channel5->size(); ii++)
                 h_WF_channel3[run][event]->SetBinContent(ii+1,waveForm_channel5->at(ii));   
            }
            else if(timingChannels[2] == 6){
             for(unsigned int ii = 0; ii < waveForm_channel6->size(); ii++)
                 h_WF_channel3[run][event]->SetBinContent(ii+1,waveForm_channel6->at(ii));   
            }
            else if(timingChannels[2] == 7){
             for(unsigned int ii = 0; ii < waveForm_channel7->size(); ii++)
                 h_WF_channel3[run][event]->SetBinContent(ii+1,waveForm_channel7->at(ii));   
            }
            else if(timingChannels[2] == 8){
             for(unsigned int ii = 0; ii < waveForm_channel8->size(); ii++)
                 h_WF_channel3[run][event]->SetBinContent(ii+1,waveForm_channel8->at(ii));   
            }  
        }
        
        tmp_pair = std::make_pair(run,event);
        vec_run_event.push_back(tmp_pair);

    }

    float x_min = h_double->GetBinCenter(h_double->GetMaximumBin())-h_double->GetRMS();
    float x_max = h_double->GetBinCenter(h_double->GetMaximumBin())+h_double->GetRMS();
    gGauss->SetParameter(0,h_double->Integral(h_double->FindBin(x_min),h_double->FindBin(x_max)));
    gGauss->SetParameter(1,h_double->GetMean(h_double->GetBinCenter(h_double->GetMaximumBin())));
    gGauss->SetParameter(2,h_double->GetRMS());
    h_double->Fit("gGauss","","",x_min/1.3,x_max/1.3);

    std::cout<<"sigma difference = "<<gGauss->GetParameter(2)<<std::endl;
    float sigma_diff = fabs(gGauss->GetParameter(2));
    float sigma_diffE = fabs(gGauss->GetParError(2));
    float mean_diff = gGauss->GetParameter(1);

    sigma_double = sigma_diff;
    sigmaError_double = sigma_diffE; 
    mean_double = mean_diff;

    x_min = h_triple->GetBinCenter(h_triple->GetMaximumBin())-h_triple->GetRMS();
    x_max = h_triple->GetBinCenter(h_triple->GetMaximumBin())+h_triple->GetRMS();
    gGauss->SetParameter(0,h_triple->Integral(h_triple->FindBin(x_min),h_triple->FindBin(x_max)));
    gGauss->SetParameter(1,h_triple->GetMean(h_triple->GetBinCenter(h_triple->GetMaximumBin())));
    gGauss->SetParameter(2,h_triple->GetRMS());
    h_triple->Fit("gGauss","","",x_min/1.3,x_max/1.3);
    float sigma_mcp = fabs(gGauss->GetParameter(2));
    float sigma_mcpE = fabs(gGauss->GetParError(2));
    float mean_mcp = gGauss->GetParameter(1);

    sigma_triple = sigma_mcp;
    sigmaError_triple = sigma_mcpE; 
    mean_triple = mean_mcp;

    x_min = h_double13->GetBinCenter(h_double13->GetMaximumBin())-h_double13->GetRMS();
    x_max = h_double13->GetBinCenter(h_double13->GetMaximumBin())+h_double13->GetRMS();
    gGauss->SetParameter(0,h_double13->Integral(h_double13->FindBin(x_min),h_double13->FindBin(x_max)));
    gGauss->SetParameter(1,h_double13->GetMean(h_double13->GetBinCenter(h_double13->GetMaximumBin())));
    gGauss->SetParameter(2,h_double13->GetRMS());
    h_double13->Fit("gGauss","","",x_min/1.3,x_max/1.3);
    float sigma_t2dat0 = fabs(gGauss->GetParameter(2));
    float sigma_t2dat0E = fabs(gGauss->GetParError(2));
    float mean_t2dat0 = gGauss->GetParameter(1);

    sigma_double13 = sigma_t2dat0;
    sigmaError_double13 = sigma_t2dat0E; 
    mean_double13 = mean_t2dat0;

    x_min = h_double23->GetBinCenter(h_double23->GetMaximumBin())-h_double23->GetRMS();
    x_max = h_double23->GetBinCenter(h_double23->GetMaximumBin())+h_double23->GetRMS();
    gGauss->SetParameter(0,h_double23->Integral(h_double23->FindBin(x_min),h_double23->FindBin(x_max)));
    gGauss->SetParameter(1,h_double23->GetMean(h_double23->GetBinCenter(h_double23->GetMaximumBin())));
    gGauss->SetParameter(2,h_double23->GetRMS());
    h_double23->Fit("gGauss","","",x_min/1.3,x_max/1.3);
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
    sigmaMCP_output = sigmaMCP;

    float sigmaMEAN = fabs(sigmaMCP/sqrt(2.));
    float sigmaMEANE = fabs(sigmaMCPE/sqrt(2.));
    std::cout << " sigmaMEAN = " << sigmaMEAN*1000 << " +/- " << sigmaMEANE*1000 << " ps" << std::endl;
    sigmaMEAN_output = sigmaMEAN;

    float sigmaPLANA = sqrt( pow(sigma_mcp,2.) - pow(sigmaMEAN,2.) );
    float sigmaPLANAE = (pow(2.*sigma_mcp*sigma_mcpE,2.)+pow(2.*sigmaMEAN*sigmaMEANE,2.))/2./sqrt(pow(sigma_mcp,2.)-pow(sigmaMEAN,2.));
    std::cout << " sigmaPLANA = " << sigmaPLANA*1000 << " +/- " << sigmaPLANAE*1000 << " ps" << std::endl;
    sigmaPLANA_output = sigmaPLANA;

    std::cout << std::endl;

    int nDouble = 0;
    float nDouble_fake = 0;
    int nTriple = 0;
    float nTriple_fake = 0;

    ampMax1 = 0.;
    ampMax2 = 0.;
    ampMax3 = 0.;
    timeAmpMax1 = 0.;
    timeAmpMax2 = 0.;
    timeAmpMax3 = 0.;
    timeConstFrac1 = 0.;
    timeConstFrac2 = 0.;
    timeConstFrac3 = 0.;

    std::vector<float> vec_double_fake;
    std::vector<float> vec_triple_fake;
   
    for(int ientry = 0; ientry < nt->GetEntries(); ientry++){
        if(ientry%100==0) std::cout<<"--- Reading entry = "<< ientry  << " for computing efficiencies..."<<std::endl;
        nt->GetEntry(ientry);
        
        bool isDouble = false;
        bool isTriple = false;

        ChooseChannels(map_channels,ampMax1,ampMax2,ampMax3,timeAmpMax1,timeAmpMax2,timeAmpMax3,timeConstFrac1,timeConstFrac2,timeConstFrac3);
          
        if(ampMax1 <= thresAmp[0]) continue;
        if(ampMax2 <= thresAmp[1]) continue;

        if(fabs(timeConstFrac2-timeConstFrac1 - mean_diff) < 5*sigma_diff) nDouble++;
        else vec_double_fake.push_back(timeConstFrac2-timeConstFrac1);
       
        if(ampMax3 <= thresAmp[2]) continue;
        
        if(fabs((timeConstFrac2+timeConstFrac1)/2.- timeConstFrac3 - mean_mcp) < 5*sigma_mcp) nTriple++;
        else vec_triple_fake.push_back((timeConstFrac2+timeConstFrac2)/2.- timeConstFrac_channel3);
        
    }

    

    std::cout << std::endl;

    std::sort(vec_double_fake.begin(),vec_double_fake.end());
    std::sort(vec_triple_fake.begin(),vec_triple_fake.end());
    
    
    if(vec_triple_fake.size() > 1) nTriple_fake = 10*sigma_mcp*vec_triple_fake.size()/fabs(vec_triple_fake.at(0)-vec_triple_fake.at(vec_triple_fake.size()-1));
    else if (vec_triple_fake.size() == 1) nTriple_fake = 10*sigma_mcp;
    else if (vec_triple_fake.size() == 0) nTriple_fake = 0.;

    if(vec_double_fake.size() > 1) nDouble_fake = 10*sigma_diff*vec_double_fake.size()/fabs(vec_double_fake.at(0)-vec_double_fake.at(vec_double_fake.size()-1));
    else if (vec_double_fake.size() == 1) nDouble_fake = 10*sigma_diff;
    else if (vec_double_fake.size() == 0) nDouble_fake = 0.;

    float eff = float(nTriple-nTriple_fake)/float(nDouble-nDouble_fake);
    float eff_err = sqrt(eff*(1-eff)/float(nDouble-nDouble_fake));

    std::cout << "nTriple = " << float(nTriple) << " - nTriple_fake = " << nTriple_fake << std::endl;
    std::cout << "nDouble = " << float(nDouble) << " - nDouble_fake = " << nDouble_fake << std::endl;
    std::cout << "Efficiency = " << eff*100 << " +/- "<< eff_err*100 << " %" << std::endl;
    std::cout << "nTriple_simple = " << float(nTriple_simple) << std::endl;
    std::cout << "nDouble_simple = " << float(nDouble_simple) << std::endl;
    std::cout << "Efficiency_simple = " << nTriple_simple/nDouble_simple*100 << " %" << std::endl;

    nTriple_output = (float)nTriple;
    nTriple_Fake_output = (float)nTriple_fake;
    nDouble_output = (float)nDouble;
    nDouble_Fake_output = (float)nDouble_fake;
    Efficiency = eff;
    EfficiencyError = eff_err;

    std::cout << std::endl;

    ampMax1 = 0.;
    ampMax2 = 0.;
    ampMax3 = 0;
    timeAmpMax1 = 0.;
    timeAmpMax2 = 0.;
    timeAmpMax3 = 0.;
    timeConstFrac1 = 0.;
    timeConstFrac2 = 0.;
    timeConstFrac3 = 0.;

    for(int ientry = 0; ientry < nt->GetEntries(); ientry++){
        if(ientry%100==0) std::cout<<"--- Reading entry = "<< ientry << " for filling the tree..." <<std::endl;
        nt->GetEntry(ientry);

        event_output = event;
        run_output = run; 

        ChooseChannels(map_channels,ampMax1,ampMax2,ampMax3,timeAmpMax1,timeAmpMax2,timeAmpMax3,timeConstFrac1,timeConstFrac2,timeConstFrac3);

        ampMax_Trigger_output = ampMax_Trigger;
        ampMax_channel1_output = ampMax1;
        ampMax_channel2_output = ampMax2;
        ampMax_channel3_output = ampMax3;
        
        timeAmpMax_Trigger_output = timeAmpMax_Trigger;
        timeAmpMax_channel1_output = timeAmpMax1;
        timeAmpMax_channel2_output = timeAmpMax2;
        timeAmpMax_channel3_output = timeAmpMax3;
        
        timeConstFrac_Trigger_output = timeConstFrac_Trigger;          
        timeConstFrac_channel1_output = timeConstFrac1;
        timeConstFrac_channel2_output = timeConstFrac2;
        timeConstFrac_channel3_output = timeConstFrac3;
       
        isDouble_output = 0;
        isTriple_output = 0;
        
        int isDouble = 0;
        int isTriple = 0;

        if(ampMax1 <= thresAmp[0]) continue;
        if(ampMax2 <= thresAmp[1]) continue;

        isDouble = 1;
        
        if(ampMax3 <= thresAmp[2]) continue;
        
        isTriple = 1;
        
        isDouble_output = isDouble;
        isTriple_output = isTriple;

        nt_output->Fill();
        
    }

    std::cout << std::endl;
    
    TFile* f1 = new TFile(("analyzed_"+std::string(inputFile)).c_str(),"RECREATE");
    f1->cd();

    nt_output->Write();
    h_double->Write();
    h_triple->Write();
    h_double13->Write();
    h_double23->Write();
    
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


void ChooseChannels(std::map<int,int> map_channels,float& ampMax1,float& ampMax2,float& ampMax3,float& timeAmpMax1,float& timeAmpMax2,float& timeAmpMax3,float& timeConstFrac1,float& timeConstFrac2,float& timeConstFrac3){
         
        if(map_channels[0] == 0){
           ampMax1 = ampMax_channel0;
           timeAmpMax1 = timeAmpMax_channel0;
           timeConstFrac1 = timeConstFrac_channel0;
        }
        else if(map_channels[0] == 1){
           ampMax1 = ampMax_channel1;
           timeAmpMax1 = timeAmpMax_channel1;
           timeConstFrac1 = timeConstFrac_channel1;
        }
        else if(map_channels[0] == 2){
           ampMax1 = ampMax_channel2;
           timeAmpMax1 = timeAmpMax_channel2;
           timeConstFrac1 = timeConstFrac_channel2;
        }
        else if(map_channels[0] == 3){
           ampMax1 = ampMax_channel3;
           timeAmpMax1 = timeAmpMax_channel3;
           timeConstFrac1 = timeConstFrac_channel3;
        }
        else if(map_channels[0] == 4){
           ampMax1 = ampMax_channel4;
           timeAmpMax1 = timeAmpMax_channel4;
           timeConstFrac1 = timeConstFrac_channel4;
        }
        else if(map_channels[0] == 5){
           ampMax1 = ampMax_channel5;
           timeAmpMax1 = timeAmpMax_channel5;
           timeConstFrac1 = timeConstFrac_channel5;
        }
        else if(map_channels[0] == 6){
           ampMax1 = ampMax_channel6;
           timeAmpMax1 = timeAmpMax_channel6;
           timeConstFrac1 = timeConstFrac_channel6;
        }
        else if(map_channels[0] == 7){
           ampMax1 = ampMax_channel7;
           timeAmpMax1 = timeAmpMax_channel7;
           timeConstFrac1 = timeConstFrac_channel7;
        }
        else if(map_channels[0] == 8){
           ampMax1 = ampMax_channel8;
           timeAmpMax1 = timeAmpMax_channel8;
           timeConstFrac1 = timeConstFrac_channel8;
        }

        if(map_channels[1] == 0){
           ampMax2 = ampMax_channel0;
           timeAmpMax2 = timeAmpMax_channel0;
           timeConstFrac2 = timeConstFrac_channel0;
        }
        else if(map_channels[1] == 1){
           ampMax2 = ampMax_channel1;
           timeAmpMax2 = timeAmpMax_channel1;
           timeConstFrac2 = timeConstFrac_channel1;
        }
        else if(map_channels[1] == 2){
           ampMax2 = ampMax_channel2;
           timeAmpMax2 = timeAmpMax_channel2;
           timeConstFrac2= timeConstFrac_channel2;
        }
        else if(map_channels[1] == 3){
           ampMax2 = ampMax_channel3;
           timeAmpMax2 = timeAmpMax_channel3;
           timeConstFrac2 = timeConstFrac_channel3;
        }
        else if(map_channels[1] == 4){
           ampMax2 = ampMax_channel4;
           timeAmpMax2 = timeAmpMax_channel4;
           timeConstFrac2 = timeConstFrac_channel4;
        }
        else if(map_channels[1] == 5){
           ampMax2 = ampMax_channel5;
           timeAmpMax2 = timeAmpMax_channel5;
           timeConstFrac2 = timeConstFrac_channel5;
        }
        else if(map_channels[1] == 6){
           ampMax2 = ampMax_channel6;
           timeAmpMax2 = timeAmpMax_channel6;
           timeConstFrac2 = timeConstFrac_channel6;
        }
        else if(map_channels[1] == 7){
           ampMax2 = ampMax_channel7;
           timeAmpMax2 = timeAmpMax_channel7;
           timeConstFrac2 = timeConstFrac_channel7;
        }
        else if(map_channels[1] == 8){
           ampMax2 = ampMax_channel8;
           timeAmpMax2 = timeAmpMax_channel8;
           timeConstFrac2 = timeConstFrac_channel8;
        }

        if(map_channels[2] == 0){
           ampMax3 = ampMax_channel0;
           timeAmpMax3 = timeAmpMax_channel0;
           timeConstFrac3 = timeConstFrac_channel0;
        }
        else if(map_channels[2] == 1){
           ampMax3 = ampMax_channel1;
           timeAmpMax3 = timeAmpMax_channel1;
           timeConstFrac3 = timeConstFrac_channel1;
        }
        else if(map_channels[2] == 2){
           ampMax3 = ampMax_channel2;
           timeAmpMax3 = timeAmpMax_channel2;
           timeConstFrac3 = timeConstFrac_channel2;
        }
        else if(map_channels[2] == 3){
           ampMax3 = ampMax_channel3;
           timeAmpMax3 = timeAmpMax_channel3;
           timeConstFrac3 = timeConstFrac_channel3;
        }
        else if(map_channels[2] == 4){
           ampMax3 = ampMax_channel4;
           timeAmpMax3 = timeAmpMax_channel4;
           timeConstFrac3 = timeConstFrac_channel4;
        }
        else if(map_channels[2] == 5){
           ampMax3 = ampMax_channel5;
           timeAmpMax3 = timeAmpMax_channel5;
           timeConstFrac3 = timeConstFrac_channel5;
        }
        else if(map_channels[2] == 6){
           ampMax3 = ampMax_channel6;
           timeAmpMax3 = timeAmpMax_channel6;
           timeConstFrac3 = timeConstFrac_channel6;
        }
        else if(map_channels[2] == 7){
           ampMax3 = ampMax_channel7;
           timeAmpMax3 = timeAmpMax_channel7;
           timeConstFrac3 = timeConstFrac_channel7;
        }
        else if(map_channels[2] == 8){
           ampMax3 = ampMax_channel8;
           timeAmpMax3 = timeAmpMax_channel8;
           timeConstFrac3 = timeConstFrac_channel8;
        }
}

