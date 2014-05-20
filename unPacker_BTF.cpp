/*************************************************************

    simple program to read bin files
    compile with --->  c++ -o unPacker_BTF `root-config --cflags --glibs` unPacker_BTF.cpp
         or with --->  c++ -o unPacker_BTF unPacker_BTF.cpp `root-config --cflags --glibs`
    run with --->  ./unPacker_BTF WaveForms_BTF/run_000.root 0 5 10 1 0 // usage: ./unPacker_BTF directory/file.bin TriggerChannel nPoints-interpolation baseline-time saveWaveForm saveWFHistos
    
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

#include "TimeConstFraction.h"
#include "InitTree_BTF.h"
#include "SetOutputTree.h"

using namespace std;

bool computeBestFrac = false;

const float BinToTime = 0.2;
const float thresAmp[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
const float isNegative[9] = {1,1,1,0,1,1,1,1,1};

float computeMean(std::vector<float> sample);
float computeDispersion(std::vector<float> sample);

int main(int argc, char** argv)
{

    gROOT->ProcessLine("#include <vector>");

    char* runFolder = argv[1];
    int triggerChannel = atoi(argv[2]) ;
    int nPointsInterpolation = atoi(argv[3]);
    float timeBaseLine = atof(argv[4]);
    int saveWF = atoi(argv[5]);
    int saveWFHistos = atoi(argv[6]);
  
    std::cout << "runFolder             = " << runFolder << std::endl;
    std::cout << "triggerChannel        = " << triggerChannel << std::endl;
    std::cout << "nPointsInterpolation  = " << nPointsInterpolation << std::endl;
    std::cout << "timeBaseLine          = " << timeBaseLine << std::endl;
    std::cout << "saveWF                = " << saveWF << std::endl;
    std::cout << "saveWFHistos          = " << saveWFHistos << std::endl;

    int binBaseLine = int(timeBaseLine/BinToTime);
    
    TString command;
    std::string tmpRunFolder = std::string(runFolder);

    if(tmpRunFolder.find(".root") != std::string::npos) command = "ls "+ tmpRunFolder + " > input.tmp"; 
    else command = "ls "+ tmpRunFolder + "/*.root > input.tmp"; 

    gSystem -> Exec(command); 
 
    TTree *nt = new TTree("nt","nt");
    nt->SetDirectory(0);
   
    outTree(nt,triggerChannel);

    // loop over the file and read data
    ifstream inFile("input.tmp");
    string line;
    int ifile=0;
    int nCh = 9;
    int nSize = 1024;

    std::vector<int> vec_run;
    std::vector<float> vec_slope;
    std::vector<float> vec_slope_error;
    std::vector<float> vec_chi2;
    std::vector<float> vec_chi2_error;
    std::vector<float> fitted_fraction;
    std::map<int,std::vector<float> > channels;
    std::map<int,std::vector<float> > tmp;
    std::map<int,std::vector<float> > tmp_baseline;
    std::map<int,std::map<float,int> > map_AmpTm;

    std::map<int,std::vector<float> > baseline_global;
    float baselineMean_global[9], baselineDispersion_global[9];

    std::map<int,std::map<int,TH1F*> > h_WF_Trigger;
    std::map<int,std::map<int,TH1F*> > h_WF_channel0; 
    std::map<int,std::map<int,TH1F*> > h_WF_channel1; 
    std::map<int,std::map<int,TH1F*> > h_WF_channel2;
    std::map<int,std::map<int,TH1F*> > h_WF_channel3;
    std::map<int,std::map<int,TH1F*> > h_WF_channel4;
    std::map<int,std::map<int,TH1F*> > h_WF_channel5;
    std::map<int,std::map<int,TH1F*> > h_WF_channel6;
    std::map<int,std::map<int,TH1F*> > h_WF_channel7;
    std::map<int,std::map<int,TH1F*> > h_WF_channel8;

    std::vector<std::pair<int,int> > vec_run_event;
    std::pair<int,int> tmp_pair;

    char histoName[200];

    std::map<int,std::map<int,TH1F*> > h_Slope;
    for(int jj = 0; jj < 9; jj++)
        for(int ii = 0; ii <= 100; ii++){
            sprintf(histoName, "h_Slope_%d_%d",jj,ii);  
            h_Slope[jj][ii] = new TH1F(histoName,histoName,400000,-2000.,2000.);
        }

    std::map<int,std::map<int,TH1F*> > h_Chi2;
    for(int jj = 0; jj < 9; jj++)
        for(int ii = 0; ii <= 100; ii++){
            sprintf(histoName, "h_Chi2_%d_%d",jj,ii);  
            h_Chi2[jj][ii] = new TH1F(histoName,histoName,200000,0.,2000.);
        }

    std::map<int,TGraphErrors*> gr_slope;
    for(int jj = 0; jj < 9; jj++)
        gr_slope[jj] = new TGraphErrors();

    std::map<int,TGraphErrors*> gr_chi2;
    for(int jj = 0; jj < 9; jj++)
        gr_chi2[jj] = new TGraphErrors();

    TFile* inputFile;
    TTree* inputTree;
    
    nPointsInterpolation_output = nPointsInterpolation;
     
  if(computeBestFrac == true){
    while(getline(inFile,line)){

      ifile++;
      std::cout << "Reading File: " << ifile << " - " << line << " for computing best fraction..." <<std::endl;

      //read data
      inputFile = TFile::Open(line.c_str());
      inputTree = (TTree*)inputFile->Get("eventRawData");
      InitTree(inputTree);

      for(int ievent=0; ievent<inputTree->GetEntries(); ievent++)
      {

        event = evtNumber;

        if(ievent%100==0) cout << "--- Reading entry = " << ievent << endl;
        inputTree->GetEntry(ievent);

        for(int iSample=0; iSample<nDigiSamples; iSample++){
            channels[digiChannel[iSample]].push_back(digiSampleValue[iSample]);
            tmp[digiChannel[iSample]].push_back(digiSampleValue[iSample]);
            map_AmpTm[digiChannel[iSample]][digiSampleValue[iSample]] = digiSampleIndex[iSample];
            if(digiSampleIndex[iSample] <= binBaseLine) tmp_baseline[digiChannel[iSample]].push_back(digiSampleValue[iSample]);
        }

        for(int iCh=0; iCh<nCh; iCh++)
        {   
            float baseline = computeMean(tmp_baseline[iCh]);
            
            baseline_global[iCh].push_back(baseline);

            std::sort(tmp[iCh].begin(),tmp[iCh].end());
            
            float tmpAmp = 0.;
            int tmpTime = 0;
            float tmpTimeConstFrac = 0.;

            //std::cout << "isNegative[" << iCh << "] = " << isNegative[iCh] << std::endl;
            if(isNegative[iCh] == 1){
               tmpAmp = tmp[iCh].at(0);
               tmpTime = map_AmpTm[iCh][tmp[iCh].at(0)];
            }else{
               tmpAmp = tmp[iCh].at(tmp.size()-1);
               tmpTime = map_AmpTm[iCh][tmp[iCh].at(tmp[iCh].size()-1)];
            }
            
            for(int ii = 0; ii <= 100; ii++){

                int ref = 0;
                float fraction = ii/100.;
                
	        for(int iSample = tmpTime; iSample >= 0; --iSample)
                {
                    float absAmp = fabs(channels[iCh].at(iSample)-baseline);
                    float absFrac = fabs(tmpAmp-baseline)*fraction;
		    ref = iSample;
                    if(absAmp <= absFrac) break;
                }
                if(fabs(channels[iCh].at(tmpTime)-baseline) > thresAmp[iCh]){
                   TimeConstFraction timeCF(ref,tmpTime,channels[iCh],fraction,baseline,nPointsInterpolation,BinToTime);
                   h_Slope[iCh][ii]->Fill(timeCF.getSlope());
                   h_Chi2[iCh][ii]->Fill(timeCF.getChi2());
                }
            }

        }
        
        for(int iCh = 0; iCh < nCh; iCh++){
            tmp[iCh].clear();
            tmp_baseline[iCh].clear();
            channels[iCh].clear();
        }
        
      }
    }

    for(int iCh=0; iCh<nCh; iCh++)
    {       
          for(unsigned int ii = 0; ii < h_Slope[iCh].size(); ii++){
              vec_slope.push_back(h_Slope[iCh][ii]->GetMean());
              vec_slope_error.push_back(h_Slope[iCh][ii]->GetRMS()/sqrt(h_Slope[iCh][ii]->GetEntries()));
          } 

          for(unsigned int ii = 0; ii < h_Chi2[iCh].size(); ii++){
              vec_chi2.push_back(h_Chi2[iCh][ii]->GetMean());
              vec_chi2_error.push_back(h_Chi2[iCh][ii]->GetRMS()/sqrt(h_Chi2[iCh][ii]->GetEntries()));
          }    

          for(unsigned int ii = 0; ii < vec_slope.size(); ii++){
              gr_slope[iCh]->SetPoint(ii,ii/100.,vec_slope.at(ii));
              gr_slope[iCh]->SetPointError(ii,0,vec_slope_error.at(ii));
          }

          for(unsigned int ii = 0; ii < vec_slope.size(); ii++){
              gr_chi2[iCh]->SetPoint(ii,ii/100.,vec_chi2.at(ii));
              gr_chi2[iCh]->SetPointError(ii,0,vec_chi2_error.at(ii));
          }

          baselineMean_global[iCh] = computeMean(baseline_global[iCh]);
          baselineDispersion_global[iCh] = computeDispersion(baseline_global[iCh]);
          
          vec_slope.clear();
          vec_slope_error.clear();
          vec_chi2.clear();
          vec_chi2_error.clear();
    }

    for(int iCh=0; iCh<nCh; iCh++)
    {
     
        TF1* func = new TF1("func","pol2",0.,1.);
        gr_slope[iCh]->Fit(func,"","",0.2,0.9);
          
        fitted_fraction.push_back(-1*func->GetParameter(1)/(2*func->GetParameter(2)));

    }

  }else{
    for(int iCh=0; iCh<nCh; iCh++)
        fitted_fraction.push_back(0.55);
  }

    ifile=0;
    line = "";
    ifstream inFile_fill("input.tmp");
 
     while(getline(inFile_fill,line)){

      ifile++;
      std::cout << "Reading File: " << ifile << " - " << line << " for filling the ntuple..." <<std::endl;

      // get run number
      char split_char = '/';
      std::vector<std::string> tokens;
      std::istringstream split(line);
      for(std::string each; getline(split, each, split_char); tokens.push_back(each));

      int icut_run = -1;
      for(unsigned int ii = 0; ii < tokens.size(); ii++)
        if(tokens.at(ii).find("run") != std::string::npos) icut_run = ii;

      split_char = '_';
      std::vector<std::string> tokens_run;
      std::istringstream split_run(tokens.at(icut_run));
      for(std::string each; getline(split_run, each, split_char); tokens_run.push_back(each));
    
      run = ::atoi(tokens_run.at(1).c_str());
      vec_run.push_back(run);

      //read data
      inputFile = TFile::Open(line.c_str());
      inputTree = (TTree*)inputFile->Get("eventRawData");
      InitTree(inputTree);

      int nCh = 9;
      vector<float> channels[nCh];
      int isNegative[9];
      for(int ii = 0; ii < 9; ii++)
          isNegative[ii] = -1;

      for(int ievent=0; ievent<inputTree->GetEntries(); ievent++)
      {

        event = evtNumber;

        if(ievent%100==0) cout << "--- Reading entry = " << ievent << endl;
        inputTree->GetEntry(ievent);

        tmp_pair = std::make_pair(run,event);
        vec_run_event.push_back(tmp_pair);

        baseLine_Trigger = 0.;
        baseLine_channel0 = 0.;
        baseLine_channel1 = 0.;
        baseLine_channel2 = 0.;
        baseLine_channel3 = 0.;
        baseLine_channel4 = 0.;
        baseLine_channel5 = 0.;
        baseLine_channel6 = 0.;
        baseLine_channel7 = 0.;
        baseLine_channel8 = 0.;

        baselineDispersion_Trigger = 0.;
        baselineDispersion_channel0 = 0.;
        baselineDispersion_channel1 = 0.;
        baselineDispersion_channel2 = 0.;
        baselineDispersion_channel3 = 0.;
        baselineDispersion_channel4 = 0.;
        baselineDispersion_channel5 = 0.;
        baselineDispersion_channel6 = 0.;
        baselineDispersion_channel7 = 0.;
        baselineDispersion_channel8 = 0.;

        baseLineGlobal_Trigger = 0.;
        baseLineGlobal_channel0 = 0.;
        baseLineGlobal_channel1 = 0.;
        baseLineGlobal_channel2 = 0.;
        baseLineGlobal_channel3 = 0.;
        baseLineGlobal_channel4 = 0.;
        baseLineGlobal_channel5 = 0.;
        baseLineGlobal_channel6 = 0.;
        baseLineGlobal_channel7 = 0.;
        baseLineGlobal_channel8 = 0.;

        baselineDispersionGlobal_Trigger = 0.; 
        baselineDispersionGlobal_channel0 = 0.;
        baselineDispersionGlobal_channel1 = 0.;
        baselineDispersionGlobal_channel2 = 0.;
        baselineDispersionGlobal_channel3 = 0.;
        baselineDispersionGlobal_channel4 = 0.;
        baselineDispersionGlobal_channel5 = 0.;
        baselineDispersionGlobal_channel6 = 0.;
        baselineDispersionGlobal_channel7 = 0.;
        baselineDispersionGlobal_channel8 = 0.;

        ampMax_Trigger = 0.;
        ampMax_channel0 = 0.;
        ampMax_channel1 = 0.;
        ampMax_channel2 = 0.;
        ampMax_channel3 = 0.;
        ampMax_channel4 = 0.;
        ampMax_channel5 = 0.;
        ampMax_channel6 = 0.;
        ampMax_channel7 = 0.;
        ampMax_channel8 = 0.;

        timeAmpMax_Trigger = 0.;
        timeAmpMax_channel0 = 0.;
        timeAmpMax_channel1 = 0.;
        timeAmpMax_channel2 = 0.;
        timeAmpMax_channel3 = 0.;
        timeAmpMax_channel4 = 0.;
        timeAmpMax_channel5 = 0.;
        timeAmpMax_channel6 = 0.;
        timeAmpMax_channel7 = 0.;
        timeAmpMax_channel8 = 0.;

        constFrac_Trigger = 0.;     
        constFrac_channel0 = 0.;     
        constFrac_channel1 = 0.;
        constFrac_channel2 = 0.;
        constFrac_channel3 = 0.;
        constFrac_channel4 = 0.;
        constFrac_channel5 = 0.;
        constFrac_channel6 = 0.;
        constFrac_channel7 = 0.;
        constFrac_channel8 = 0.;

        timeConstFrac_Trigger = 0.;     
        timeConstFrac_channel0 = 0.;     
        timeConstFrac_channel1 = 0.;
        timeConstFrac_channel2 = 0.;
        timeConstFrac_channel3 = 0.;
        timeConstFrac_channel4 = 0.;
        timeConstFrac_channel5 = 0.;
        timeConstFrac_channel6 = 0.;
        timeConstFrac_channel7 = 0.;
        timeConstFrac_channel8 = 0.;

        slopeConstFrac_Trigger = 0.;  
        slopeConstFrac_channel0 = 0.;        
        slopeConstFrac_channel1 = 0.;
        slopeConstFrac_channel2 = 0.;
        slopeConstFrac_channel3 = 0.;
        slopeConstFrac_channel4 = 0.;
        slopeConstFrac_channel5 = 0.;
        slopeConstFrac_channel6 = 0.;
        slopeConstFrac_channel7 = 0.;
        slopeConstFrac_channel8 = 0.;

        chi2ConstFrac_Trigger = 0.;         
        chi2ConstFrac_channel0 = 0.; 
        chi2ConstFrac_channel1 = 0.;
        chi2ConstFrac_channel2 = 0.;
        chi2ConstFrac_channel3 = 0.;
        chi2ConstFrac_channel4 = 0.;
        chi2ConstFrac_channel5 = 0.;
        chi2ConstFrac_channel6 = 0.;
        chi2ConstFrac_channel7 = 0.;
        chi2ConstFrac_channel8 = 0.;
        
        for(int iSample=0; iSample<nDigiSamples; iSample++){
            channels[digiChannel[iSample]].push_back(digiSampleValue[iSample]);
            tmp[digiChannel[iSample]].push_back(digiSampleValue[iSample]);
            map_AmpTm[digiChannel[iSample]][digiSampleValue[iSample]] = digiSampleIndex[iSample];
            if(digiSampleIndex[iSample] <= binBaseLine) tmp_baseline[digiChannel[iSample]].push_back(digiSampleValue[iSample]);
        }
        
        for(int iCh=0; iCh<nCh; iCh++)
        {   
            float baseline = computeMean(tmp_baseline[iCh]);
            float baseline_dispersion = computeDispersion(tmp_baseline[iCh]);
            
            baseline_global[iCh].push_back(baseline);

            std::sort(tmp[iCh].begin(),tmp[iCh].end());
            
            float tmpAmp = 0.;
            int tmpTime = 0;
            float tmpFrac = 0.;
            float tmpTimeConstFrac = 0.;
            float tmpSlopeConstFrac = 0.;
            float tmpChi2ConstFrac = 0.;

            //std::cout << "isNegative[" << iCh << "] = " << isNegative[iCh] << std::endl;
            if(isNegative[iCh] == -1){
               tmpAmp = tmp[iCh].at(0);
               tmpTime = map_AmpTm[iCh][tmp[iCh].at(0)];
            }else{
               tmpAmp = tmp[iCh].at(tmp.size()-1);
               tmpTime = map_AmpTm[iCh][tmp[iCh].at(tmp[iCh].size()-1)];
            }
            
            int ref = 0;
            
	    for(int iSample = tmpTime; iSample >= 0; --iSample)
            {
                float absAmp = fabs(channels[iCh].at(iSample)-baseline);
                float absFrac = fabs(tmpAmp-baseline)*fitted_fraction.at(iCh);
		ref = iSample;
                if(absAmp <= absFrac) break;
            }
            
            TimeConstFraction timeCF(ref,tmpTime,channels[iCh],fitted_fraction.at(iCh),baseline,nPointsInterpolation,BinToTime);
            tmpFrac = fitted_fraction.at(iCh);
            tmpTimeConstFrac = timeCF.getTime();
            tmpSlopeConstFrac = timeCF.getSlope();
            tmpChi2ConstFrac = timeCF.getChi2();

            if(iCh == triggerChannel)
            {
               for(unsigned int iSample=0; iSample<channels[iCh].size(); iSample++)
                   if(saveWF == 1) waveForm_Trigger.push_back(channels[iCh].at(iSample));

               baseLine_Trigger = baseline;
               baselineDispersion_Trigger = baseline_dispersion;
               baseLineGlobal_Trigger = baselineMean_global[iCh];
               baselineDispersionGlobal_Trigger = baselineDispersion_global[iCh];
               ampMax_Trigger = fabs(tmpAmp-baseline);
               timeAmpMax_Trigger = BinToTime*tmpTime;   
               constFrac_Trigger = tmpFrac;     
               timeConstFrac_Trigger = tmpTimeConstFrac;  
               slopeConstFrac_Trigger = tmpSlopeConstFrac;     
               chi2ConstFrac_Trigger = tmpChi2ConstFrac;    
                
               if(saveWFHistos == 1){

                  sprintf(histoName, "WF_Trigger_%d_%d",run,event);  
                  h_WF_Trigger[run][event] = new TH1F(histoName,histoName,nSize,0.,float(nSize)*BinToTime);

                  for(unsigned int ii = 0; ii < waveForm_Trigger.size(); ii++)
                      h_WF_Trigger[run][event]->SetBinContent(ii+1,waveForm_Trigger.at(ii));   
               }
            }
            else if(iCh == 0)
            {
               for(unsigned int iSample=0; iSample<channels[iCh].size(); iSample++)
                   if(saveWF == 1) waveForm_channel0.push_back(channels[iCh].at(iSample));

               baseLine_channel0 = baseline;
               baselineDispersion_channel0 = baseline_dispersion;
               baseLineGlobal_channel0 = baselineMean_global[iCh];
               baselineDispersionGlobal_channel0 = baselineDispersion_global[iCh];
               ampMax_channel0 = fabs(tmpAmp-baseline);
               timeAmpMax_channel0 = BinToTime*tmpTime;   
               constFrac_channel0 = tmpFrac;     
               timeConstFrac_channel0 = tmpTimeConstFrac;  
               slopeConstFrac_channel0 = tmpSlopeConstFrac;     
               chi2ConstFrac_channel0 = tmpChi2ConstFrac;    
                
               if(saveWFHistos == 1){

                  sprintf(histoName, "WF_channel0_%d_%d",run,event);  
                  h_WF_channel0[run][event] = new TH1F(histoName,histoName,nSize,0.,float(nSize)*BinToTime);

                  for(unsigned int ii = 0; ii < waveForm_channel0.size(); ii++)
                      h_WF_channel0[run][event]->SetBinContent(ii+1,waveForm_channel0.at(ii));   
               }
            }
            else if(iCh == 1)
            {
               for(unsigned int iSample=0; iSample<channels[iCh].size(); iSample++)
                   if(saveWF == 1) waveForm_channel1.push_back(channels[iCh].at(iSample));

               baseLine_channel1 = baseline;
               baselineDispersion_channel1 = baseline_dispersion;
               baseLineGlobal_channel1 = baselineMean_global[iCh];
               baselineDispersionGlobal_channel1 = baselineDispersion_global[iCh];
               ampMax_channel1 = fabs(tmpAmp-baseline);
               timeAmpMax_channel1 = BinToTime*tmpTime; 
               constFrac_channel1 = tmpFrac;     
               timeConstFrac_channel1 = tmpTimeConstFrac;  
               slopeConstFrac_channel1 = tmpSlopeConstFrac;     
               chi2ConstFrac_channel1 = tmpChi2ConstFrac;     

               if(saveWFHistos == 1){

                  sprintf(histoName, "WF_channel1_%d_%d",run,event);  
                  h_WF_channel1[run][event] = new TH1F(histoName,histoName,nSize,0.,float(nSize)*BinToTime);

                  for(unsigned int ii = 0; ii < waveForm_channel1.size(); ii++)
                      h_WF_channel1[run][event]->SetBinContent(ii+1,waveForm_channel1.at(ii));   
               }  
 
            }
            else if(iCh == 2)
            {
               for(unsigned int iSample=0; iSample<channels[iCh].size(); iSample++)
                   if(saveWF == 1) waveForm_channel2.push_back(channels[iCh].at(iSample));

               baseLine_channel2 = baseline;
               baselineDispersion_channel2 = baseline_dispersion;
               baseLineGlobal_channel2 = baselineMean_global[iCh];
               baselineDispersionGlobal_channel2 = baselineDispersion_global[iCh];
               ampMax_channel2 = fabs(tmpAmp-baseline);
               timeAmpMax_channel2 = BinToTime*tmpTime;  
               constFrac_channel2 = tmpFrac;     
               timeConstFrac_channel2 = tmpTimeConstFrac;  
               slopeConstFrac_channel2 = tmpSlopeConstFrac;     
               chi2ConstFrac_channel2 = tmpChi2ConstFrac;     

               if(saveWFHistos == 1){

                  sprintf(histoName, "WF_channel2_%d_%d",run,event);  
                  h_WF_channel2[run][event] = new TH1F(histoName,histoName,nSize,0.,float(nSize)*BinToTime);

                  for(unsigned int ii = 0; ii < waveForm_channel2.size(); ii++)
                      h_WF_channel2[run][event]->SetBinContent(ii+1,waveForm_channel2.at(ii));   
               } 

            }
            else if(iCh == 3)
            {
               for(unsigned int iSample=0; iSample<channels[iCh].size(); iSample++)
                   if(saveWF == 1) waveForm_channel3.push_back(channels[iCh].at(iSample));

               baseLine_channel3 = baseline;
               baselineDispersion_channel3 = baseline_dispersion;
               baseLineGlobal_channel3 = baselineMean_global[iCh];
               baselineDispersionGlobal_channel3 = baselineDispersion_global[iCh];
               ampMax_channel3 = fabs(tmpAmp-baseline);
               timeAmpMax_channel3 = BinToTime*tmpTime;  
               constFrac_channel3 = tmpFrac;     
               timeConstFrac_channel3 = tmpTimeConstFrac;  
               slopeConstFrac_channel3 = tmpSlopeConstFrac;     
               chi2ConstFrac_channel3 = tmpChi2ConstFrac;   

               if(saveWFHistos == 1){

                  sprintf(histoName, "WF_channel3_%d_%d",run,event);  
                  h_WF_channel3[run][event] = new TH1F(histoName,histoName,nSize,0.,float(nSize)*BinToTime);

                  for(unsigned int ii = 0; ii < waveForm_channel3.size(); ii++)
                      h_WF_channel3[run][event]->SetBinContent(ii+1,waveForm_channel3.at(ii));   
               } 
               
            }
            else if(iCh == 4)
            {
               for(unsigned int iSample=0; iSample<channels[iCh].size(); iSample++)
                   if(saveWF == 1) waveForm_channel4.push_back(channels[iCh].at(iSample));

               baseLine_channel4 = baseline;
               baselineDispersion_channel4 = baseline_dispersion;
               baseLineGlobal_channel4 = baselineMean_global[iCh];
               baselineDispersionGlobal_channel4 = baselineDispersion_global[iCh];
               ampMax_channel4 = fabs(tmpAmp-baseline);
               timeAmpMax_channel4 = BinToTime*tmpTime;
               constFrac_channel4 = tmpFrac;     
               timeConstFrac_channel4 = tmpTimeConstFrac;  
               slopeConstFrac_channel4 = tmpSlopeConstFrac;     
               chi2ConstFrac_channel4 = tmpChi2ConstFrac;   

               if(saveWFHistos == 1){

                  sprintf(histoName, "WF_channel4_%d_%d",run,event);  
                  h_WF_channel4[run][event] = new TH1F(histoName,histoName,nSize,0.,float(nSize)*BinToTime);

                  for(unsigned int ii = 0; ii < waveForm_channel4.size(); ii++)
                      h_WF_channel4[run][event]->SetBinContent(ii+1,waveForm_channel4.at(ii));   
               }   

            }
            else if(iCh == 5)
            {
               for(unsigned int iSample=0; iSample<channels[iCh].size(); iSample++)
                   if(saveWF == 1) waveForm_channel5.push_back(channels[iCh].at(iSample));

               baseLine_channel5 = baseline;
               baselineDispersion_channel5 = baseline_dispersion;
               baseLineGlobal_channel5 = baselineMean_global[iCh];
               baselineDispersionGlobal_channel5 = baselineDispersion_global[iCh];
               ampMax_channel5 = fabs(tmpAmp-baseline);
               timeAmpMax_channel5 = BinToTime*tmpTime;  
               constFrac_channel1 = tmpFrac;     
               timeConstFrac_channel5 = tmpTimeConstFrac;  
               slopeConstFrac_channel5 = tmpSlopeConstFrac;     
               chi2ConstFrac_channel5 = tmpChi2ConstFrac;   

               if(saveWFHistos == 1){

                  sprintf(histoName, "WF_channel5_%d_%d",run,event);  
                  h_WF_channel5[run][event] = new TH1F(histoName,histoName,nSize,0.,float(nSize)*BinToTime);

                  for(unsigned int ii = 0; ii < waveForm_channel5.size(); ii++)
                      h_WF_channel5[run][event]->SetBinContent(ii+1,waveForm_channel5.at(ii));   
               }  

            }
            else if(iCh == 6)
            {
               for(unsigned int iSample=0; iSample<channels[iCh].size(); iSample++)
                   if(saveWF == 1) waveForm_channel6.push_back(channels[iCh].at(iSample));

               baseLine_channel6 = baseline;
               baselineDispersion_channel6 = baseline_dispersion;
               baseLineGlobal_channel6 = baselineMean_global[iCh];
               baselineDispersionGlobal_channel6 = baselineDispersion_global[iCh];
               ampMax_channel6 = fabs(tmpAmp-baseline);
               timeAmpMax_channel6 = BinToTime*tmpTime;   
               constFrac_channel6 = tmpFrac;     
               timeConstFrac_channel6 = tmpTimeConstFrac;  
               slopeConstFrac_channel6 = tmpSlopeConstFrac;     
               chi2ConstFrac_channel6 = tmpChi2ConstFrac;   

               if(saveWFHistos == 1){

                  sprintf(histoName, "WF_channel6_%d_%d",run,event);  
                  h_WF_channel6[run][event] = new TH1F(histoName,histoName,nSize,0.,float(nSize)*BinToTime);

                  for(unsigned int ii = 0; ii < waveForm_channel6.size(); ii++)
                      h_WF_channel6[run][event]->SetBinContent(ii+1,waveForm_channel6.at(ii));   
               }  

            }
            else if(iCh == 7)
            {
               for(unsigned int iSample=0; iSample<channels[iCh].size(); iSample++)
                   if(saveWF == 1) waveForm_channel7.push_back(channels[iCh].at(iSample));

               baseLine_channel7 = baseline;
               baselineDispersion_channel7 = baseline_dispersion;
               baseLineGlobal_channel7 = baselineMean_global[iCh];
               baselineDispersionGlobal_channel7 = baselineDispersion_global[iCh];
               ampMax_channel7 = fabs(tmpAmp-baseline);
               timeAmpMax_channel7 = BinToTime*tmpTime;
               constFrac_channel7 = tmpFrac;     
               timeConstFrac_channel7 = tmpTimeConstFrac;  
               slopeConstFrac_channel7 = tmpSlopeConstFrac;     
               chi2ConstFrac_channel7 = tmpChi2ConstFrac;   

               if(saveWFHistos == 1){

                  sprintf(histoName, "WF_channel7_%d_%d",run,event);  
                  h_WF_channel7[run][event] = new TH1F(histoName,histoName,nSize,0.,float(nSize)*BinToTime);

                  for(unsigned int ii = 0; ii < waveForm_channel7.size(); ii++)
                      h_WF_channel7[run][event]->SetBinContent(ii+1,waveForm_channel7.at(ii));   
               } 

            }
            else if(iCh == 8)
            {
               for(unsigned int iSample=0; iSample<channels[iCh].size(); iSample++)
                   if(saveWF == 1) waveForm_channel8.push_back(channels[iCh].at(iSample));

               baseLine_channel8 = baseline;
               baselineDispersion_channel8 = baseline_dispersion;
               baseLineGlobal_channel8 = baselineMean_global[iCh];
               baselineDispersionGlobal_channel8 = baselineDispersion_global[iCh];
               ampMax_channel8 = fabs(tmpAmp-baseline);
               timeAmpMax_channel8 = BinToTime*tmpTime;  
               constFrac_channel8 = tmpFrac;     
               timeConstFrac_channel8 = tmpTimeConstFrac;  
               slopeConstFrac_channel8 = tmpSlopeConstFrac;     
               chi2ConstFrac_channel8 = tmpChi2ConstFrac;    

               if(saveWFHistos == 1){

                  sprintf(histoName, "WF_channel8_%d_%d",run,event);  
                  h_WF_channel8[run][event] = new TH1F(histoName,histoName,nSize,0.,float(nSize)*BinToTime);

                  for(unsigned int ii = 0; ii < waveForm_channel8.size(); ii++)
                      h_WF_channel8[run][event]->SetBinContent(ii+1,waveForm_channel8.at(ii));   
               }  
          
            }
            
        }

        if(saveWF == 0){

            waveForm_Trigger.push_back(0.);
            waveForm_channel0.push_back(0.);
            waveForm_channel1.push_back(0.);
            waveForm_channel2.push_back(0.);
            waveForm_channel3.push_back(0.);
            waveForm_channel4.push_back(0.);
            waveForm_channel5.push_back(0.);
            waveForm_channel6.push_back(0.);
            waveForm_channel7.push_back(0.);
            waveForm_channel8.push_back(0.);

        }
        else{
            if(triggerChannel == 0) waveForm_channel0.push_back(0.);
            if(triggerChannel == 1) waveForm_channel1.push_back(0.);
            if(triggerChannel == 2) waveForm_channel2.push_back(0.);
            if(triggerChannel == 3) waveForm_channel3.push_back(0.);
            if(triggerChannel == 4) waveForm_channel4.push_back(0.);
            if(triggerChannel == 5) waveForm_channel5.push_back(0.);
            if(triggerChannel == 6) waveForm_channel6.push_back(0.);
            if(triggerChannel == 7) waveForm_channel7.push_back(0.);
            if(triggerChannel == 8) waveForm_channel8.push_back(0.);
        }
   
        nt->Fill();  

        waveForm_Trigger.clear();
        waveForm_channel0.clear(); 
        waveForm_channel1.clear(); 
        waveForm_channel2.clear(); 
        waveForm_channel3.clear(); 
        waveForm_channel4.clear(); 
        waveForm_channel5.clear(); 
        waveForm_channel6.clear(); 
        waveForm_channel7.clear(); 
        waveForm_channel8.clear(); 


        for(int ii = 0; ii < nCh; ii++){
            tmp[ii].clear();
            tmp_baseline[ii].clear();
            channels[ii].clear();
        }
   
      }
  
    }

    TFile *f1;
    std::string tmpRunFolder_cp = tmpRunFolder;
    if(tmpRunFolder.find(".root") != std::string::npos){
       tmpRunFolder_cp.erase(tmpRunFolder_cp.size()-5,tmpRunFolder_cp.size()-1);
       replace(tmpRunFolder_cp.begin(),tmpRunFolder_cp.end(),'/','_');
       f1 = new TFile((tmpRunFolder_cp+"_tree.root").c_str(),"RECREATE");
    }
    else f1 = new TFile((std::string(runFolder)+"_tree.root").c_str(),"RECREATE"); 

    f1->cd();
    nt->Write("nt");

    for(int ii = 0; ii < nCh; ii++){
        if(ii == 0) sprintf(histoName, "Slope_graph_Trigger"); 
        else sprintf(histoName, "Slope_graph_channel%d",ii); 
        gr_slope[ii]->Write(histoName);
    }
    for(int ii = 0; ii < nCh; ii++){
        if(ii == 0) sprintf(histoName, "Chi2_graph_Trigger"); 
        else sprintf(histoName, "Chi2_graph_channel%d",ii); 
        gr_chi2[ii]->Write(histoName);
    }

    f1->Close();
    
    if(saveWFHistos == 1){
       
      tmpRunFolder_cp = tmpRunFolder;
      TFile *f2;
      if(tmpRunFolder.find(".root") != std::string::npos){
         tmpRunFolder_cp.erase(tmpRunFolder_cp.size()-5,tmpRunFolder_cp.size()-1);
         replace(tmpRunFolder_cp.begin(),tmpRunFolder_cp.end(),'/','_');
         f2 = new TFile(("histos_"+tmpRunFolder_cp+"_tree.root").c_str(),"RECREATE");
      }
      else f2 = new TFile(("histos_"+std::string(runFolder)+".root").c_str(),"RECREATE"); 
      f2->cd();
      for(unsigned int ii = 0; ii < vec_run_event.size(); ii++){
          h_WF_Trigger[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          if(triggerChannel != 0) h_WF_channel0[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          if(triggerChannel != 1) h_WF_channel1[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          if(triggerChannel != 2) h_WF_channel2[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          if(triggerChannel != 3) h_WF_channel3[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          if(triggerChannel != 4) h_WF_channel4[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          if(triggerChannel != 5) h_WF_channel5[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          if(triggerChannel != 6) h_WF_channel6[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          if(triggerChannel != 7) h_WF_channel7[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
          if(triggerChannel != 8) h_WF_channel8[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Write();
      }

      f2->Close();
    }

    if(saveWFHistos == 1){
       
       gStyle->SetOptStat(0000); 

       std::sort(vec_run.begin(),vec_run.end());

       std::vector<float> max, min;
       
       if(h_WF_Trigger.size() != 0){

          TCanvas* c0 = new TCanvas("c0","c0");
          c0 -> cd();
            
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++){
               max.push_back(h_WF_Trigger[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMaximum()); 
               min.push_back(h_WF_Trigger[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_Trigger[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*1.05);
          else h_WF_Trigger[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*0.95);
          if(h_min > 0.) h_WF_Trigger[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*0.95); 
          else h_WF_Trigger[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*1.05); 
          
          h_WF_Trigger[vec_run_event.at(0).first][vec_run_event.at(0).second]->Draw();
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++)
               h_WF_Trigger[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Draw("same"); 
               
          c0 -> Print("WF_trigger_total.png","png");
          c0 -> Print("WF_trigger_total.pdf","pdf");

          max.clear();
          min.clear();

          delete c0;
       } 
 
       if(h_WF_channel0.size() != 0){

          TCanvas* c0 = new TCanvas("c0","c0");
          c0 -> cd();
            
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++){
               max.push_back(h_WF_channel0[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMaximum()); 
               min.push_back(h_WF_channel0[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel0[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*1.05);
          else h_WF_channel0[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*0.95);
          if(h_min > 0.) h_WF_channel0[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*0.95); 
          else h_WF_channel0[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*1.05); 
          
          h_WF_channel0[vec_run_event.at(0).first][vec_run_event.at(0).second]->Draw();
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++)
               h_WF_channel0[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Draw("same"); 
               
          c0 -> Print("WF_channel0_total.png","png");
          c0 -> Print("WF_channel0_total.pdf","pdf");

          max.clear();
          min.clear();

          delete c0;
       } 

       if(h_WF_channel1.size() != 0){

          TCanvas* c0 = new TCanvas("c0","c0");
          c0 -> cd();
            
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++){
               max.push_back(h_WF_channel1[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMaximum()); 
               min.push_back(h_WF_channel1[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel1[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*1.05);
          else h_WF_channel1[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*0.95);
          if(h_min > 0.) h_WF_channel1[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*0.95); 
          else h_WF_channel1[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*1.05); 
          
          h_WF_channel1[vec_run_event.at(0).first][vec_run_event.at(0).second]->Draw();
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++)
               h_WF_channel1[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Draw("same"); 
               
          c0 -> Print("WF_channel1_total.png","png");
          c0 -> Print("WF_channel1_total.pdf","pdf");

          max.clear();
          min.clear();

          delete c0;
       } 

       if(h_WF_channel2.size() != 0){

          TCanvas* c0 = new TCanvas("c0","c0");
          c0 -> cd();
            
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++){
               max.push_back(h_WF_channel2[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMaximum()); 
               min.push_back(h_WF_channel2[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel2[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*1.05);
          else h_WF_channel2[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*0.95);
          if(h_min > 0.) h_WF_channel2[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*0.95); 
          else h_WF_channel2[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*1.05); 
          
          h_WF_channel2[vec_run_event.at(0).first][vec_run_event.at(0).second]->Draw();
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++)
               h_WF_channel2[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Draw("same"); 
               
          c0 -> Print("WF_channel2_total.png","png");
          c0 -> Print("WF_channel2_total.pdf","pdf");

          max.clear();
          min.clear();

          delete c0;
       } 

       if(h_WF_channel3.size() != 0){

          TCanvas* c0 = new TCanvas("c0","c0");
          c0 -> cd();
            
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++){
               max.push_back(h_WF_channel3[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMaximum()); 
               min.push_back(h_WF_channel3[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel3[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*1.05);
          else h_WF_channel3[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*0.95);
          if(h_min > 0.) h_WF_channel3[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*0.95); 
          else h_WF_channel3[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*1.05); 
          
          h_WF_channel3[vec_run_event.at(0).first][vec_run_event.at(0).second]->Draw();
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++)
               h_WF_channel3[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Draw("same"); 
               
          c0 -> Print("WF_channel3_total.png","png");
          c0 -> Print("WF_channel3_total.pdf","pdf");

          max.clear();
          min.clear();

          delete c0;
       } 
       
       if(h_WF_channel4.size() != 0){

          TCanvas* c0 = new TCanvas("c0","c0");
          c0 -> cd();
            
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++){
               max.push_back(h_WF_channel4[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMaximum()); 
               min.push_back(h_WF_channel4[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel4[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*1.05);
          else h_WF_channel4[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*0.95);
          if(h_min > 0.) h_WF_channel4[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*0.95); 
          else h_WF_channel4[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*1.05); 
          
          h_WF_channel4[vec_run_event.at(0).first][vec_run_event.at(0).second]->Draw();
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++)
               h_WF_channel4[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Draw("same"); 
               
          c0 -> Print("WF_channel4_total.png","png");
          c0 -> Print("WF_channel4_total.pdf","pdf");

          max.clear();
          min.clear();

          delete c0;
       } 
       
       if(h_WF_channel5.size() != 0){

          TCanvas* c0 = new TCanvas("c0","c0");
          c0 -> cd();
            
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++){
               max.push_back(h_WF_channel5[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMaximum()); 
               min.push_back(h_WF_channel5[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel5[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*1.05);
          else h_WF_channel5[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*0.95);
          if(h_min > 0.) h_WF_channel5[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*0.95); 
          else h_WF_channel5[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*1.05); 
          
          h_WF_channel5[vec_run_event.at(0).first][vec_run_event.at(0).second]->Draw();
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++)
               h_WF_channel5[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Draw("same"); 
               
          c0 -> Print("WF_channel5_total.png","png");
          c0 -> Print("WF_channel5_total.pdf","pdf");

          max.clear();
          min.clear();

          delete c0;
       } 
       
       if(h_WF_channel6.size() != 0){

          TCanvas* c0 = new TCanvas("c0","c0");
          c0 -> cd();
            
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++){
               max.push_back(h_WF_channel6[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMaximum()); 
               min.push_back(h_WF_channel6[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel6[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*1.05);
          else h_WF_channel6[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*0.95);
          if(h_min > 0.) h_WF_channel6[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*0.95); 
          else h_WF_channel6[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*1.05); 
          
          h_WF_channel6[vec_run_event.at(0).first][vec_run_event.at(0).second]->Draw();
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++)
               h_WF_channel6[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Draw("same"); 
               
          c0 -> Print("WF_channel6_total.png","png");
          c0 -> Print("WF_channel6_total.pdf","pdf");

          max.clear();
          min.clear();

          delete c0;
       } 
       
       if(h_WF_channel7.size() != 0){

          TCanvas* c0 = new TCanvas("c0","c0");
          c0 -> cd();
            
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++){
               max.push_back(h_WF_channel7[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMaximum()); 
               min.push_back(h_WF_channel7[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel7[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*1.05);
          else h_WF_channel7[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*0.95);
          if(h_min > 0.) h_WF_channel7[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*0.95); 
          else h_WF_channel7[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*1.05); 
          
          h_WF_channel7[vec_run_event.at(0).first][vec_run_event.at(0).second]->Draw();
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++)
               h_WF_channel7[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Draw("same"); 
               
          c0 -> Print("WF_channel7_total.png","png");
          c0 -> Print("WF_channel7_total.pdf","pdf");

          max.clear();
          min.clear();

          delete c0;
       } 
       
       if(h_WF_channel8.size() != 0){

          TCanvas* c0 = new TCanvas("c0","c0");
          c0 -> cd();
            
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++){
               max.push_back(h_WF_channel8[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMaximum()); 
               min.push_back(h_WF_channel8[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel8[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*1.05);
          else h_WF_channel8[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMaximum(h_max*0.95);
          if(h_min > 0.) h_WF_channel8[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*0.95); 
          else h_WF_channel8[vec_run_event.at(0).first][vec_run_event.at(0).second]->SetMinimum(h_min*1.05); 
          
          h_WF_channel8[vec_run_event.at(0).first][vec_run_event.at(0).second]->Draw();
          for(unsigned int ii = 0; ii < vec_run_event.size(); ii++)
               h_WF_channel8[vec_run_event.at(ii).first][vec_run_event.at(ii).second]->Draw("same"); 
               
          c0 -> Print("WF_channel8_total.png","png");
          c0 -> Print("WF_channel8_total.pdf","pdf");

          max.clear();
          min.clear();

          delete c0;
       } 
 
    }
    
    gSystem -> Exec("rm input.tmp"); 
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
