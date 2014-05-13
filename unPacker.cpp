/*************************************************************

    simple program to read bin files
    compile with --->  c++ -o unPacker `root-config --cflags --glibs` unPacker.cpp
         or with --->  c++ -o unPacker unPacker.cpp `root-config --cflags --glibs`
    run with --->  ./unPacker WaveForms 5 10 1 0 // usage: ./unPacker directory nPoints-interpolation baseline-time saveWaveForm saveWFHistos
    
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

using namespace std;

const float BinToTime = 0.2;
const float thresAmp = 200.;
const float thresAmpCh3 = 100.;

float computeMean(std::vector<float> sample);
float computeDispersion(std::vector<float> sample);

int main(int argc, char** argv)
{

    gROOT->ProcessLine("#include <vector>");

    char* runFolder = argv[1];
    int nPointsInterpolation = atoi(argv[2]);
    float timeBaseLine = atof(argv[3]);
    int saveWF = atoi(argv[4]);
    int saveWFHistos = atoi(argv[5]);
  
    std::cout << "runFolder    = " << runFolder << std::endl;
    std::cout << "nPointsInterpolation  = " << nPointsInterpolation << std::endl;
    std::cout << "timeBaseLine = " << timeBaseLine << std::endl;
    std::cout << "saveWF       = " << saveWF << std::endl;
    std::cout << "saveWFHistos = " << saveWFHistos << std::endl;

    int binBaseLine = int(timeBaseLine/BinToTime);

    TString command = "ls "+ std::string(runFolder) + "/*.bin > input.tmp"; 
    gSystem -> Exec(command); 

    int nChannels = 0;

    // Tree branches and tree structure
    int event;
    int run;
    int nPointsInterpolation_output;
    float timeBaseLine_output;

    vector<float> waveForm_Trigger;
    float baseLine_Trigger;
    float baselineDispersion_Trigger;
    float baseLineGlobal_Trigger;
    float baselineDispersionGlobal_Trigger;
    float ampMax_Trigger;
    float timeAmpMax_Trigger;
    float constFrac_Trigger;
    float timeConstFrac_Trigger;
    float slopeConstFrac_Trigger;
    float chi2ConstFrac_Trigger;
    vector<float> waveForm_channel1;
    float baseLine_channel1;
    float baselineDispersion_channel1;
    float baseLineGlobal_channel1;
    float baselineDispersionGlobal_channel1;
    float ampMax_channel1;
    float timeAmpMax_channel1;
    float constFrac_channel1;
    float timeConstFrac_channel1;
    float slopeConstFrac_channel1;
    float chi2ConstFrac_channel1;
    vector<float> waveForm_channel2;
    float baseLine_channel2;
    float baselineDispersion_channel2;
    float baseLineGlobal_channel2;
    float baselineDispersionGlobal_channel2;
    float ampMax_channel2;
    float timeAmpMax_channel2;
    float constFrac_channel2;
    float timeConstFrac_channel2;
    float slopeConstFrac_channel2;
    float chi2ConstFrac_channel2;
    vector<float> waveForm_channel3;
    float baseLine_channel3;
    float baselineDispersion_channel3;
    float baseLineGlobal_channel3;
    float baselineDispersionGlobal_channel3;
    float ampMax_channel3;
    float timeAmpMax_channel3;
    float constFrac_channel3;
    float timeConstFrac_channel3;
    float slopeConstFrac_channel3;
    float chi2ConstFrac_channel3;
    vector<float> waveForm_channel4;
    float baseLine_channel4;
    float baselineDispersion_channel4;
    float baseLineGlobal_channel4;
    float baselineDispersionGlobal_channel4;
    float ampMax_channel4;
    float timeAmpMax_channel4;
    float constFrac_channel4;
    float timeConstFrac_channel4;
    float slopeConstFrac_channel4;
    float chi2ConstFrac_channel4;
    vector<float> waveForm_channel5;
    float baseLine_channel5;
    float baselineDispersion_channel5;
    float baseLineGlobal_channel5;
    float baselineDispersionGlobal_channel5;
    float ampMax_channel5;
    float timeAmpMax_channel5;
    float constFrac_channel5;
    float timeConstFrac_channel5;
    float slopeConstFrac_channel5;
    float chi2ConstFrac_channel5;
    vector<float> waveForm_channel6;
    float baseLine_channel6;
    float baselineDispersion_channel6;
    float baseLineGlobal_channel6;
    float baselineDispersionGlobal_channel6;
    float ampMax_channel6;
    float timeAmpMax_channel6;
    float constFrac_channel6;
    float timeConstFrac_channel6;
    float slopeConstFrac_channel6;
    float chi2ConstFrac_channel6;
    vector<float> waveForm_channel7;
    float baseLine_channel7;
    float baselineDispersion_channel7;
    float baseLineGlobal_channel7;
    float baselineDispersionGlobal_channel7;
    float ampMax_channel7;
    float timeAmpMax_channel7;
    float constFrac_channel7;
    float timeConstFrac_channel7;
    float slopeConstFrac_channel7;
    float chi2ConstFrac_channel7;
    vector<float> waveForm_channel8;
    float baseLine_channel8; 
    float baselineDispersion_channel8;
    float baseLineGlobal_channel8;
    float baselineDispersionGlobal_channel8;
    float ampMax_channel8;
    float timeAmpMax_channel8;
    float constFrac_channel8;
    float timeConstFrac_channel8;
    float slopeConstFrac_channel8;
    float chi2ConstFrac_channel8;

    TTree *nt = new TTree("nt","nt");
    nt->SetDirectory(0);

    nt->Branch("event",&event,"event/I"); 
    nt->Branch("run",&run,"run/I"); 
    nt->Branch("nPointsInterpolation",&nPointsInterpolation,"nPointsInterpolation/I"); 
    nt->Branch("timeBaseLine",&timeBaseLine_output,"timeBaseLine/F"); 
    nt->Branch("waveForm_Trigger","std::vector<float>",&waveForm_Trigger); 
    nt->Branch("baseLine_Trigger",&baseLine_Trigger,"baseLine_Trigger/F");  
    nt->Branch("baselineDispersion_Trigger",&baselineDispersion_Trigger,"baselineDispersion_Trigger/F"); 
    nt->Branch("baseLineGlobal_Trigger",&baseLineGlobal_Trigger,"baseLineGlobal_Trigger/F");  
    nt->Branch("baselineDispersionGlobal_Trigger",&baselineDispersionGlobal_Trigger,"baselineDispersionGlobal_Trigger/F"); 
    nt->Branch("ampMax_Trigger",&ampMax_Trigger,"ampMax_Trigger/F"); 
    nt->Branch("timeAmpMax_Trigger",&timeAmpMax_Trigger,"timeAmpMax_Trigger/F"); 
    nt->Branch("constFrac_Trigger",&constFrac_Trigger,"constFrac_Trigger/F"); 
    nt->Branch("timeConstFrac_Trigger",&timeConstFrac_Trigger,"timeConstFrac_Trigger/F"); 
    nt->Branch("slopeConstFrac_Trigger",&slopeConstFrac_Trigger,"slopeConstFrac_Trigger/F"); 
    nt->Branch("chi2ConstFrac_Trigger",&chi2ConstFrac_Trigger,"chi2ConstFrac_Trigger/F"); 
    nt->Branch("waveForm_channel1","std::vector<float>",&waveForm_channel1);  
    nt->Branch("baseLine_channel1",&baseLine_channel1,"baseLine_channel1/F"); 
    nt->Branch("baselineDispersion_channel1",&baselineDispersion_channel1,"baselineDispersion_channel1/F"); 
    nt->Branch("baseLineGlobal_channel1",&baseLineGlobal_channel1,"baseLineGlobal_channel1/F");  
    nt->Branch("baselineDispersionGlobal_channel1",&baselineDispersionGlobal_channel1,"baselineDispersionGlobal_channel1/F"); 
    nt->Branch("ampMax_channel1",&ampMax_channel1,"ampMax_channel1/F"); 
    nt->Branch("timeAmpMax_channel1",&timeAmpMax_channel1,"timeAmpMax_channel1/F"); 
    nt->Branch("constFrac_channel1",&constFrac_channel1,"constFrac_channel1/F"); 
    nt->Branch("timeConstFrac_channel1",&timeConstFrac_channel1,"timeConstFrac_channel1/F"); 
    nt->Branch("slopeConstFrac_channel1",&slopeConstFrac_channel1,"slopeConstFrac_channel1/F"); 
    nt->Branch("chi2ConstFrac_channel1",&chi2ConstFrac_channel1,"chi2ConstFrac_channel1/F"); 
    nt->Branch("waveForm_channel2","std::vector<float>",&waveForm_channel2); 
    nt->Branch("baseLine_channel2",&baseLine_channel2,"baseLine_channel2/F"); 
    nt->Branch("baselineDispersion_channel2",&baselineDispersion_channel2,"baselineDispersion_channel2/F"); 
    nt->Branch("baseLineGlobal_channel2",&baseLineGlobal_channel2,"baseLineGlobal_channel2/F");  
    nt->Branch("baselineDispersionGlobal_channel2",&baselineDispersionGlobal_channel2,"baselineDispersionGlobal_channel2/F"); 
    nt->Branch("ampMax_channel2",&ampMax_channel2,"ampMax_channel2/F"); 
    nt->Branch("timeAmpMax_channel2",&timeAmpMax_channel2,"timeAmpMax_channel2/F"); 
    nt->Branch("constFrac_channel2",&constFrac_channel2,"constFrac_channel2/F"); 
    nt->Branch("timeConstFrac_channel2",&timeConstFrac_channel2,"timeConstFrac_channel2/F"); 
    nt->Branch("slopeConstFrac_channel2",&slopeConstFrac_channel2,"slopeConstFrac_channel2/F"); 
    nt->Branch("chi2ConstFrac_channel2",&chi2ConstFrac_channel2,"chi2ConstFrac_channel2/F"); 
    nt->Branch("waveForm_channel3","std::vector<float>",&waveForm_channel3); 
    nt->Branch("baseLine_channel3",&baseLine_channel3,"baseLine_channel3/F"); 
    nt->Branch("baselineDispersion_channel3",&baselineDispersion_channel3,"baselineDispersion_channel3/F"); 
    nt->Branch("baseLineGlobal_channel3",&baseLineGlobal_channel3,"baseLineGlobal_channel3/F");  
    nt->Branch("baselineDispersionGlobal_channel3",&baselineDispersionGlobal_channel3,"baselineDispersionGlobal_channel3/F"); 
    nt->Branch("ampMax_channel3",&ampMax_channel3,"ampMax_channel3/F"); 
    nt->Branch("timeAmpMax_channel3",&timeAmpMax_channel3,"timeAmpMax_channel3/F"); 
    nt->Branch("constFrac_channel3",&constFrac_channel3,"constFrac_channel3/F"); 
    nt->Branch("timeConstFrac_channel3",&timeConstFrac_channel3,"timeConstFrac_channel3/F"); 
    nt->Branch("slopeConstFrac_channel3",&slopeConstFrac_channel3,"slopeConstFrac_channel3/F"); 
    nt->Branch("chi2ConstFrac_channel3",&chi2ConstFrac_channel3,"chi2ConstFrac_channel3/F"); 
    nt->Branch("waveForm_channel4","std::vector<float>",&waveForm_channel4); 
    nt->Branch("baseLine_channel4",&baseLine_channel4,"baseLine_channel4/F"); 
    nt->Branch("baselineDispersion_channel4",&baselineDispersion_channel4,"baselineDispersion_channel4/F"); 
    nt->Branch("baseLineGlobal_channel4",&baseLineGlobal_channel4,"baseLineGlobal_channel4/F");  
    nt->Branch("baselineDispersionGlobal_channel4",&baselineDispersionGlobal_channel4,"baselineDispersionGlobal_channel4/F"); 
    nt->Branch("ampMax_channel4",&ampMax_channel4,"ampMax_channel4/F"); 
    nt->Branch("timeAmpMax_channel4",&timeAmpMax_channel4,"timeAmpMax_channel4/F"); 
    nt->Branch("constFrac_channel4",&constFrac_channel4,"constFrac_channel4/F"); 
    nt->Branch("timeConstFrac_channel4",&timeConstFrac_channel4,"timeConstFrac_channel4/F"); 
    nt->Branch("slopeConstFrac_channel4",&slopeConstFrac_channel4,"slopeConstFrac_channel4/F"); 
    nt->Branch("chi2ConstFrac_channel4",&chi2ConstFrac_channel4,"chi2ConstFrac_channel4/F"); 
    nt->Branch("waveForm_channel5","std::vector<float>",&waveForm_channel5);
    nt->Branch("baseLine_channel5",&baseLine_channel5,"baseLine_channel5/F"); 
    nt->Branch("baselineDispersion_channel5",&baselineDispersion_channel5,"baselineDispersion_channel5/F"); 
    nt->Branch("baseLineGlobal_channel5",&baseLineGlobal_channel5,"baseLineGlobal_channel5/F");  
    nt->Branch("baselineDispersionGlobal_channel5",&baselineDispersionGlobal_channel5,"baselineDispersionGlobal_channel5/F"); 
    nt->Branch("ampMax_channel5",&ampMax_channel5,"ampMax_channel5/F");  
    nt->Branch("timeAmpMax_channel5",&timeAmpMax_channel5,"timeAmpMax_channel5/F"); 
    nt->Branch("constFrac_channel5",&constFrac_channel5,"constFrac_channel5/F"); 
    nt->Branch("timeConstFrac_channel5",&timeConstFrac_channel5,"timeConstFrac_channel5/F"); 
    nt->Branch("slopeConstFrac_channel5",&slopeConstFrac_channel5,"slopeConstFrac_channel5/F"); 
    nt->Branch("chi2ConstFrac_channel5",&chi2ConstFrac_channel5,"chi2ConstFrac_channel5/F"); 
    nt->Branch("waveForm_channel6","std::vector<float>",&waveForm_channel6); 
    nt->Branch("baseLine_channel6",&baseLine_channel6,"baseLine_channel6/F"); 
    nt->Branch("baselineDispersion_channel6",&baselineDispersion_channel6,"baselineDispersion_channel6/F"); 
    nt->Branch("baseLineGlobal_channel6",&baseLineGlobal_channel6,"baseLineGlobal_channel6/F");  
    nt->Branch("baselineDispersionGlobal_channel6",&baselineDispersionGlobal_channel6,"baselineDispersionGlobal_channel6/F"); 
    nt->Branch("ampMax_channel6",&ampMax_channel6,"ampMax_channel6/F"); 
    nt->Branch("timeAmpMax_channel6",&timeAmpMax_channel6,"timeAmpMax_channel6/F"); 
    nt->Branch("constFrac_channel6",&constFrac_channel6,"constFrac_channel6/F"); 
    nt->Branch("timeConstFrac_channel6",&timeConstFrac_channel6,"timeConstFrac_channel6/F"); 
    nt->Branch("slopeConstFrac_channel6",&slopeConstFrac_channel6,"slopeConstFrac_channel6/F"); 
    nt->Branch("chi2ConstFrac_channel6",&chi2ConstFrac_channel6,"chi2ConstFrac_channel6/F"); 
    nt->Branch("waveForm_channel7","std::vector<float>",&waveForm_channel7); 
    nt->Branch("baseLine_channel7",&baseLine_channel7,"baseLine_channel7/F"); 
    nt->Branch("baselineDispersion_channel7",&baselineDispersion_channel7,"baselineDispersion_channel7/F"); 
    nt->Branch("baseLineGlobal_channel7",&baseLineGlobal_channel7,"baseLineGlobal_channel7/F");  
    nt->Branch("baselineDispersionGlobal_channel7",&baselineDispersionGlobal_channel7,"baselineDispersionGlobal_channel7/F"); 
    nt->Branch("ampMax_channel7",&ampMax_channel7,"ampMax_channel7/F"); 
    nt->Branch("timeAmpMax_channel7",&timeAmpMax_channel7,"timeAmpMax_channel7/F"); 
    nt->Branch("constFrac_channel7",&constFrac_channel7,"constFrac_channel7/F"); 
    nt->Branch("timeConstFrac_channel7",&timeConstFrac_channel7,"timeConstFrac_channel7/F"); 
    nt->Branch("slopeConstFrac_channel7",&slopeConstFrac_channel7,"slopeConstFrac_channel7/F"); 
    nt->Branch("chi2ConstFrac_channel7",&chi2ConstFrac_channel7,"chi2ConstFrac_channel7/F"); 
    nt->Branch("waveForm_channel8","std::vector<float>",&waveForm_channel8);
    nt->Branch("baseLine_channel8",&baseLine_channel8,"baseLine_channel8/F"); 
    nt->Branch("baselineDispersion_channel8",&baselineDispersion_channel8,"baselineDispersion_channel8/F"); 
    nt->Branch("baseLineGlobal_channel8",&baseLineGlobal_channel8,"baseLineGlobal_channel8/F");  
    nt->Branch("baselineDispersionGlobal_channel8",&baselineDispersionGlobal_channel8,"baselineDispersionGlobal_channel8/F"); 
    nt->Branch("ampMax_channel8",&ampMax_channel8,"ampMax_channel8/F"); 
    nt->Branch("timeAmpMax_channel8",&timeAmpMax_channel8,"timeAmpMax_channel8/F"); 
    nt->Branch("constFrac_channel8",&constFrac_channel8,"constFrac_channel8/F"); 
    nt->Branch("timeConstFrac_channel8",&timeConstFrac_channel8,"timeConstFrac_channel8/F"); 
    nt->Branch("slopeConstFrac_channel8",&slopeConstFrac_channel8,"slopeConstFrac_channel8/F"); 
    nt->Branch("chi2ConstFrac_channel8",&chi2ConstFrac_channel8,"chi2ConstFrac_channel8/F"); 
    
    waveForm_Trigger.clear();
    waveForm_channel1.clear();
    waveForm_channel2.clear();
    waveForm_channel3.clear();
    waveForm_channel4.clear();
    waveForm_channel5.clear();
    waveForm_channel6.clear();
    waveForm_channel7.clear();
    waveForm_channel8.clear();


    // loop over the file and read data
    ifstream inFile("input.tmp");
    string line;
    int ifile=0;

    std::vector<int> vec_run;
    std::vector<float> vec_slope;
    std::vector<float> vec_slope_error;
    std::vector<float> vec_chi2;
    std::vector<float> vec_chi2_error;
    std::vector<float> fitted_fraction;
    std::vector<float> tmp, tmp_WF;
    std::vector<float> tmp_baseline;
    std::map<float,int> map_AmpTm;

    std::map<int,std::vector<float> > baseline_global;
    float baselineMean_global[9], baselineDispersion_global[9];

    std::map<int,std::map<int,TH1F*> > h_WF_Trigger;
    std::map<int,std::map<int,TH1F*> > h_WF_channel1; 
    std::map<int,std::map<int,TH1F*> > h_WF_channel2;
    std::map<int,std::map<int,TH1F*> > h_WF_channel3;
    std::map<int,std::map<int,TH1F*> > h_WF_channel4;
    std::map<int,std::map<int,TH1F*> > h_WF_channel5;
    std::map<int,std::map<int,TH1F*> > h_WF_channel6;
    std::map<int,std::map<int,TH1F*> > h_WF_channel7;
    std::map<int,std::map<int,TH1F*> > h_WF_channel8;

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

    while(getline(inFile,line)){

      ifile++;
      std::cout << "Reading File: " << ifile << " - " << line << " for computing best fraction..." <<std::endl;

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
      FILE* input;
      vector<float> channels[9];
      vector<int> eventNumber;
      int *BinHeader, nCh=0, nSize=0, iEvent=0, nEvent=0;

      if((input = fopen(line.c_str(),"rb")) == NULL)
      {
        cout << "input file not found" << endl;
        return -1;
      }   

      BinHeader = (int*) malloc (sizeof(int)*10);
      fread(BinHeader, sizeof(int), 10, input);
      nSize = BinHeader[0];
      nCh = BinHeader[1];

      nChannels = nCh;

      int isNegative[9];
      isNegative[0] = -1;
      for(int ii = 1; ii < 9; ii++)
          isNegative[ii] = BinHeader[ii+1];


      //for(int ii = 0; ii < 10; ii++)
      //std::cout << "BinHeader[" << ii << "] = " << BinHeader[ii] << std::endl;

      fread(&iEvent, sizeof(int),1,input);
      
      float* buffer = (float*) malloc (sizeof(float)*nSize*(nCh+1));
    
      while (fread(buffer, sizeof(float), nSize*(nCh+1), input) == nSize*(nCh+1))
      {
        eventNumber.push_back(iEvent);
        for(int iCh=0; iCh<=nCh; iCh++)
        {
            for(int iSample=0; iSample<nSize; iSample++)
            {
                int pos = iCh*nSize + iSample;
                channels[iCh].push_back(buffer[pos]);
            }
        }
        nEvent++;
        fread(&iEvent, sizeof(int), 1, input);
      }
    
      for(int i=0; i<nEvent; i++)
      {

        event = eventNumber.at(i);

        if(event%100==0) cout << "--- Reading entry = " << event << endl;

        for(int iCh=0; iCh<=nCh; iCh++)
        {   
            for(int iSample=0; iSample<nSize; iSample++){

                if(iSample <= binBaseLine) tmp_baseline.push_back(channels[iCh].at(iSample+nSize*i));

                tmp.push_back(channels[iCh].at(iSample+nSize*i));
                tmp_WF.push_back(channels[iCh].at(iSample+nSize*i));
                map_AmpTm[channels[iCh].at(iSample+nSize*i)] = iSample;
            }

            float baseline = computeMean(tmp_baseline);
            
            baseline_global[iCh].push_back(baseline);

            std::sort(tmp.begin(),tmp.end());
            
            float tmpAmp = 0.;
            int tmpTime = 0;
            float tmpTimeConstFrac = 0.;

            //std::cout << "isNegative[" << iCh << "] = " << isNegative[iCh] << std::endl;
            if(isNegative[iCh] == -1){
               tmpAmp = tmp.at(0);
               tmpTime = map_AmpTm[tmp.at(0)];
            }else{
               tmpAmp = tmp.at(tmp.size()-1);
               tmpTime = map_AmpTm[tmp.at(tmp.size()-1)];
            }
            
            for(int ii = 0; ii <= 100; ii++){

                int ref = 0;
                float fraction = ii/100.;
                
	        for(int iSample = tmpTime; iSample >= 0; --iSample)
                {
                    float absAmp = fabs(channels[iCh].at(iSample+nSize*i)-baseline);
                    float absFrac = fabs(tmpAmp-baseline)*fraction;
		    ref = iSample;
                    if(absAmp <= absFrac) break;
                }
                if(iCh != 3){
                   if(fabs(tmp_WF.at(tmpTime)-baseline) > thresAmp){
                      TimeConstFraction timeCF(ref,tmpTime,tmp_WF,fraction,baseline,nPointsInterpolation,BinToTime);
                      h_Slope[iCh][ii]->Fill(timeCF.getSlope());
                      h_Chi2[iCh][ii]->Fill(timeCF.getChi2());
                   }
                }else{
                   if(fabs(tmp_WF.at(tmpTime)-baseline) > thresAmpCh3){
                      TimeConstFraction timeCF(ref,tmpTime,tmp_WF,fraction,baseline,nPointsInterpolation,BinToTime);
                      h_Slope[iCh][ii]->Fill(timeCF.getSlope());
                      h_Chi2[iCh][ii]->Fill(timeCF.getChi2());
                   }
                }
            }

            tmp.clear();
            tmp_WF.clear();
            tmp_baseline.clear();
        }
      }
    }

    for(int iCh=0; iCh<=8; iCh++)
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

    for(int iCh=0; iCh<=nChannels; iCh++)
    {
     
        TF1* func = new TF1("func","pol2",0.,1.);
        gr_slope[iCh]->Fit(func,"","",0.2,0.9);
          
        fitted_fraction.push_back(-1*func->GetParameter(1)/(2*func->GetParameter(2)));

    }

    vec_run.clear();

    ifile=0;
    line = "";
    ifstream inFile_fill("input.tmp");
    
    while(getline(inFile_fill,line)){

      ifile++;
      std::cout << "Reading File: " << ifile << " - " << line << " for filling the ntuple..." <<std::endl;

      nPointsInterpolation_output = nPointsInterpolation;
      timeBaseLine_output = timeBaseLine;

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
      FILE* input;
      vector<float> channels[9];
      vector<int> eventNumber;
      int *BinHeader, nCh=0, nSize=0, iEvent=0, nEvent=0;

      if((input = fopen(line.c_str(),"rb")) == NULL)
      {
        cout << "input file not found" << endl;
        return -1;
      }   

      BinHeader = (int*) malloc (sizeof(int)*10);
      fread(BinHeader, sizeof(int), 10, input);
      nSize = BinHeader[0];
      nCh = BinHeader[1];

      int isNegative[9];
      isNegative[0] = -1;
      for(int ii = 1; ii < 9; ii++)
          isNegative[ii] = BinHeader[ii+1];


      //for(int ii = 0; ii < 10; ii++)
      //std::cout << "BinHeader[" << ii << "] = " << BinHeader[ii] << std::endl;

      fread(&iEvent, sizeof(int),1,input);
      
      float* buffer = (float*) malloc (sizeof(float)*nSize*(nCh+1));
    
      while (fread(buffer, sizeof(float), nSize*(nCh+1), input) == nSize*(nCh+1))
      {
        eventNumber.push_back(iEvent);
        for(int iCh=0; iCh<=nCh; iCh++)
        {
            for(int iSample=0; iSample<nSize; iSample++)
            {
                int pos = iCh*nSize + iSample;
                channels[iCh].push_back(buffer[pos]);
            }
        }
        nEvent++;
        fread(&iEvent, sizeof(int), 1, input);
      }
    
      //Fill tree
      for(int i=0; i<nEvent; i++)
      {

        event = eventNumber.at(i);

        if(event%100==0) cout << "--- Reading entry = " << event << endl;

        baseLine_Trigger = 0.;
        baseLine_channel1 = 0.;
        baseLine_channel2 = 0.;
        baseLine_channel3 = 0.;
        baseLine_channel4 = 0.;
        baseLine_channel5 = 0.;
        baseLine_channel6 = 0.;
        baseLine_channel7 = 0.;
        baseLine_channel8 = 0.;

        baselineDispersion_Trigger = 0.;
        baselineDispersion_channel1 = 0.;
        baselineDispersion_channel2 = 0.;
        baselineDispersion_channel3 = 0.;
        baselineDispersion_channel4 = 0.;
        baselineDispersion_channel5 = 0.;
        baselineDispersion_channel6 = 0.;
        baselineDispersion_channel7 = 0.;
        baselineDispersion_channel8 = 0.;

        baseLineGlobal_Trigger = 0.;
        baseLineGlobal_channel1 = 0.;
        baseLineGlobal_channel2 = 0.;
        baseLineGlobal_channel3 = 0.;
        baseLineGlobal_channel4 = 0.;
        baseLineGlobal_channel5 = 0.;
        baseLineGlobal_channel6 = 0.;
        baseLineGlobal_channel7 = 0.;
        baseLineGlobal_channel8 = 0.;

        baselineDispersionGlobal_Trigger = 0.;
        baselineDispersionGlobal_channel1 = 0.;
        baselineDispersionGlobal_channel2 = 0.;
        baselineDispersionGlobal_channel3 = 0.;
        baselineDispersionGlobal_channel4 = 0.;
        baselineDispersionGlobal_channel5 = 0.;
        baselineDispersionGlobal_channel6 = 0.;
        baselineDispersionGlobal_channel7 = 0.;
        baselineDispersionGlobal_channel8 = 0.;

        ampMax_Trigger = 0.;
        ampMax_channel1 = 0.;
        ampMax_channel2 = 0.;
        ampMax_channel3 = 0.;
        ampMax_channel4 = 0.;
        ampMax_channel5 = 0.;
        ampMax_channel6 = 0.;
        ampMax_channel7 = 0.;
        ampMax_channel8 = 0.;

        timeAmpMax_Trigger = 0.;
        timeAmpMax_channel1 = 0.;
        timeAmpMax_channel2 = 0.;
        timeAmpMax_channel3 = 0.;
        timeAmpMax_channel4 = 0.;
        timeAmpMax_channel5 = 0.;
        timeAmpMax_channel6 = 0.;
        timeAmpMax_channel7 = 0.;
        timeAmpMax_channel8 = 0.;

        constFrac_Trigger = 0.;          
        constFrac_channel1 = 0.;
        constFrac_channel2 = 0.;
        constFrac_channel3 = 0.;
        constFrac_channel4 = 0.;
        constFrac_channel5 = 0.;
        constFrac_channel6 = 0.;
        constFrac_channel7 = 0.;
        constFrac_channel8 = 0.;

        timeConstFrac_Trigger = 0.;          
        timeConstFrac_channel1 = 0.;
        timeConstFrac_channel2 = 0.;
        timeConstFrac_channel3 = 0.;
        timeConstFrac_channel4 = 0.;
        timeConstFrac_channel5 = 0.;
        timeConstFrac_channel6 = 0.;
        timeConstFrac_channel7 = 0.;
        timeConstFrac_channel8 = 0.;

        slopeConstFrac_Trigger = 0.;          
        slopeConstFrac_channel1 = 0.;
        slopeConstFrac_channel2 = 0.;
        slopeConstFrac_channel3 = 0.;
        slopeConstFrac_channel4 = 0.;
        slopeConstFrac_channel5 = 0.;
        slopeConstFrac_channel6 = 0.;
        slopeConstFrac_channel7 = 0.;
        slopeConstFrac_channel8 = 0.;

        chi2ConstFrac_Trigger = 0.;          
        chi2ConstFrac_channel1 = 0.;
        chi2ConstFrac_channel2 = 0.;
        chi2ConstFrac_channel3 = 0.;
        chi2ConstFrac_channel4 = 0.;
        chi2ConstFrac_channel5 = 0.;
        chi2ConstFrac_channel6 = 0.;
        chi2ConstFrac_channel7 = 0.;
        chi2ConstFrac_channel8 = 0.;

        for(int iCh=0; iCh<=nCh; iCh++)
        {   
            for(int iSample=0; iSample<nSize; iSample++){

                if(iSample <= binBaseLine) tmp_baseline.push_back(channels[iCh].at(iSample+nSize*i));

                tmp.push_back(channels[iCh].at(iSample+nSize*i));
                tmp_WF.push_back(channels[iCh].at(iSample+nSize*i));
                map_AmpTm[channels[iCh].at(iSample+nSize*i)] = iSample;
            }
            
            float baseline = computeMean(tmp_baseline);
            float baseline_dispersion = computeDispersion(tmp_baseline);

            std::sort(tmp.begin(),tmp.end());
            
            float tmpAmp = 0.;
            int tmpTime = 0;
            float tmpFrac = 0.;
            float tmpTimeConstFrac = 0.;
            float tmpSlopeConstFrac = 0.;
            float tmpChi2ConstFrac = 0.;
            
            //std::cout << "isNegative[" << iCh << "] = " << isNegative[iCh] << std::endl;
            if(isNegative[iCh] == -1){
               tmpAmp = tmp.at(0);
               tmpTime = map_AmpTm[tmp.at(0)];
            }else{
               tmpAmp = tmp.at(tmp.size()-1);
               tmpTime = map_AmpTm[tmp.at(tmp.size()-1)];
            }   
             

            int ref = 0;
              
	    for(int iSample = tmpTime; iSample >= 0; --iSample)
            {
                float absAmp = fabs(channels[iCh].at(iSample+nSize*i)-baseline);
                float absFrac = fabs(tmpAmp-baseline)*fitted_fraction.at(iCh);
		ref = iSample;
                if(absAmp <= absFrac) break;
            }
           
            TimeConstFraction timeCF(ref,tmpTime,tmp_WF,fitted_fraction.at(iCh),baseline,nPointsInterpolation,BinToTime);
            tmpFrac = fitted_fraction.at(iCh);
            tmpTimeConstFrac = timeCF.getTime();
            tmpSlopeConstFrac = timeCF.getSlope();
            tmpChi2ConstFrac = timeCF.getChi2();
            
            if(iCh == 0)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_Trigger.push_back(channels[iCh].at(iSample+nSize*i));

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

            if(iCh == 1)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel1.push_back(channels[iCh].at(iSample+nSize*i));

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

            if(iCh == 2)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel2.push_back(channels[iCh].at(iSample+nSize*i));

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
             
            if(iCh == 3)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel3.push_back(channels[iCh].at(iSample+nSize*i));

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

            if(iCh == 4)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel4.push_back(channels[iCh].at(iSample+nSize*i));

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

            if(iCh == 5)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel5.push_back(channels[iCh].at(iSample+nSize*i));

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

            if(iCh == 6)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel6.push_back(channels[iCh].at(iSample+nSize*i));

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

            if(iCh == 7)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel7.push_back(channels[iCh].at(iSample+nSize*i));

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

            if(iCh == 8)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel8.push_back(channels[iCh].at(iSample+nSize*i));

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
 
            tmp.clear();
            tmp_WF.clear();
            tmp_baseline.clear();

            
        }

        for(int iCh=nCh+1; iCh<=8; iCh++){

            if(iCh == 0) waveForm_Trigger.push_back(0.);
            if(iCh == 1) waveForm_channel1.push_back(0.);
            if(iCh == 2) waveForm_channel2.push_back(0.);
            if(iCh == 3) waveForm_channel3.push_back(0.);
            if(iCh == 4) waveForm_channel4.push_back(0.);
            if(iCh == 5) waveForm_channel5.push_back(0.);
            if(iCh == 6) waveForm_channel6.push_back(0.);
            if(iCh == 7) waveForm_channel7.push_back(0.);
            if(iCh == 8) waveForm_channel8.push_back(0.);
        }
        
        if(saveWF == 0){

            waveForm_Trigger.push_back(0.);
            waveForm_channel1.push_back(0.);
            waveForm_channel2.push_back(0.);
            waveForm_channel3.push_back(0.);
            waveForm_channel4.push_back(0.);
            waveForm_channel5.push_back(0.);
            waveForm_channel6.push_back(0.);
            waveForm_channel7.push_back(0.);
            waveForm_channel8.push_back(0.);

        }
   
        nt->Fill();  

        waveForm_Trigger.clear();
        waveForm_channel1.clear(); 
        waveForm_channel2.clear(); 
        waveForm_channel3.clear(); 
        waveForm_channel4.clear(); 
        waveForm_channel5.clear(); 
        waveForm_channel6.clear(); 
        waveForm_channel7.clear(); 
        waveForm_channel8.clear(); 
      }
  
    }


    TFile *f1 = new TFile((std::string(runFolder)+"_tree.root").c_str(),"RECREATE"); 
    f1->cd();
    nt->Write("nt");
    f1->Close();

    if(saveWFHistos == 1){
      TFile *f2 = new TFile(("histos_"+std::string(runFolder)+".root").c_str(),"RECREATE"); 
      f2->cd();
      for(unsigned int ii = 0; ii < h_WF_Trigger.size(); ii++)
        for(unsigned int jj = 0; jj < h_WF_Trigger[ii].size(); jj++)
            h_WF_Trigger[ii][jj]->Write();
      for(unsigned int ii = 0; ii < h_WF_channel1.size(); ii++)
        for(unsigned int jj = 0; jj < h_WF_channel1[ii].size(); jj++)
            h_WF_channel1[ii][jj]->Write();
      for(unsigned int ii = 0; ii < h_WF_channel2.size(); ii++)
        for(unsigned int jj = 0; jj < h_WF_channel2[ii].size(); jj++)
            h_WF_channel2[ii][jj]->Write();
      for(unsigned int ii = 0; ii < h_WF_channel3.size(); ii++)
        for(unsigned int jj = 0; jj < h_WF_channel3[ii].size(); jj++)
            h_WF_channel3[ii][jj]->Write();
      for(unsigned int ii = 0; ii < h_WF_channel4.size(); ii++)
        for(unsigned int jj = 0; jj < h_WF_channel4[ii].size(); jj++)
            h_WF_channel4[ii][jj]->Write();
      for(unsigned int ii = 0; ii < h_WF_channel5.size(); ii++)
        for(unsigned int jj = 0; jj < h_WF_channel5[ii].size(); jj++)
            h_WF_channel5[ii][jj]->Write();
      for(unsigned int ii = 0; ii < h_WF_channel6.size(); ii++)
        for(unsigned int jj = 0; jj < h_WF_channel6[ii].size(); jj++)
            h_WF_channel6[ii][jj]->Write();
      for(unsigned int ii = 0; ii < h_WF_channel7.size(); ii++)
        for(unsigned int jj = 0; jj < h_WF_channel7[ii].size(); jj++)
            h_WF_channel7[ii][jj]->Write();
      for(unsigned int ii = 0; ii < h_WF_channel8.size(); ii++)
        for(unsigned int jj = 0; jj < h_WF_channel8[ii].size(); jj++)
            h_WF_channel8[ii][jj]->Write();
 
      for(int ii = 0; ii < 9; ii++){
          if(ii == 0) sprintf(histoName, "Slope_graph_Trigger"); 
          else sprintf(histoName, "Slope_graph_channel%d",ii); 
          gr_slope[ii]->Write(histoName);
      }
      for(int ii = 0; ii < 9; ii++){
          if(ii == 0) sprintf(histoName, "Chi2_graph_Trigger"); 
          else sprintf(histoName, "Chi2_graph_channel%d",ii); 
          gr_chi2[ii]->Write(histoName);
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
            
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_Trigger[ii].size(); jj++){
               max.push_back(h_WF_Trigger[ii][jj]->GetMaximum()); 
               min.push_back(h_WF_Trigger[ii][jj]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_Trigger[vec_run.at(0)][0]->SetMaximum(h_max*1.05);
          else h_WF_Trigger[vec_run.at(0)][0]->SetMaximum(h_max*0.95);
          if(h_min > 0.)h_WF_Trigger[vec_run.at(0)][0]->SetMinimum(h_min*0.95); 
          else h_WF_Trigger[vec_run.at(0)][0]->SetMinimum(h_min*1.05); 
          
          h_WF_Trigger[vec_run.at(0)][0]->Draw();
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_Trigger[ii].size(); jj++)
               h_WF_Trigger[ii][jj]->Draw("same"); 
               
          c0 -> Print("WF_trigger_total.png","png");
          c0 -> Print("WF_trigger_total.pdf","pdf");

          max.clear();
          min.clear();
       }
      
       if(h_WF_channel1.size() != 0){

          TCanvas* c1 = new TCanvas("c1","c1");
          c1 -> cd();

          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel1[ii].size(); jj++){
               max.push_back(h_WF_channel1[ii][jj]->GetMaximum()); 
               min.push_back(h_WF_channel1[ii][jj]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel1[vec_run.at(0)][0]->SetMaximum(h_max*1.05);
          else h_WF_channel1[vec_run.at(0)][0]->SetMaximum(h_max*0.95);
          if(h_min > 0.)h_WF_channel1[vec_run.at(0)][0]->SetMinimum(h_min*0.95); 
          else h_WF_channel1[vec_run.at(0)][0]->SetMinimum(h_min*1.05); 
          
          h_WF_channel1[vec_run.at(0)][0]->Draw();
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel1[ii].size(); jj++)
               h_WF_channel1[ii][jj]->Draw("same");

          c1 -> Print("WF_channel1_total.png","png");
          c1 -> Print("WF_channel1_total.pdf","pdf");

          max.clear();
          min.clear(); 
       }
       
       if(h_WF_channel2.size() != 0){

          TCanvas* c2 = new TCanvas("c2","c2");
          c2 -> cd();
          
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel2[ii].size(); jj++){
               max.push_back(h_WF_channel2[ii][jj]->GetMaximum()); 
               min.push_back(h_WF_channel2[ii][jj]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel2[vec_run.at(0)][0]->SetMaximum(h_max*1.05);
          else h_WF_channel2[vec_run.at(0)][0]->SetMaximum(h_max*0.95);
          if(h_min > 0.)h_WF_channel2[vec_run.at(0)][0]->SetMinimum(h_min*0.95); 
          else h_WF_channel2[vec_run.at(0)][0]->SetMinimum(h_min*1.05); 
          
          h_WF_channel2[vec_run.at(0)][0]->Draw();
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel2[ii].size(); jj++)
               h_WF_channel2[ii][jj]->Draw("same");

          c2 -> Print("WF_channel2_total.png","png");
          c2 -> Print("WF_channel2_total.pdf","pdf");

          max.clear();
          min.clear(); 
       }
        
       if(h_WF_channel3.size() != 0){   

          TCanvas* c3 = new TCanvas("c3","c3");
          c3 -> cd();

          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel3[ii].size(); jj++){
               max.push_back(h_WF_channel3[ii][jj]->GetMaximum()); 
               min.push_back(h_WF_channel3[ii][jj]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel3[vec_run.at(0)][0]->SetMaximum(h_max*1.05);
          else h_WF_channel3[vec_run.at(0)][0]->SetMaximum(h_max*0.95);
          if(h_min > 0.)h_WF_channel3[vec_run.at(0)][0]->SetMinimum(h_min*0.95); 
          else h_WF_channel3[vec_run.at(0)][0]->SetMinimum(h_min*1.05); 
          
          h_WF_channel3[vec_run.at(0)][0]->Draw();
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel3[ii].size(); jj++)
               h_WF_channel3[ii][jj]->Draw("same");

          c3 -> Print("WF_channel3_total.png","png");
          c3 -> Print("WF_channel3_total.pdf","pdf");

          max.clear();
          min.clear(); 
       }
       

       if(h_WF_channel4.size() != 0){
       
          TCanvas* c4 = new TCanvas("c4","c4");
          c4 -> cd();

          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel4[ii].size(); jj++){
               max.push_back(h_WF_channel4[ii][jj]->GetMaximum()); 
               min.push_back(h_WF_channel4[ii][jj]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel4[vec_run.at(0)][0]->SetMaximum(h_max*1.05);
          else h_WF_channel4[vec_run.at(0)][0]->SetMaximum(h_max*0.95);
          if(h_min > 0.)h_WF_channel4[vec_run.at(0)][0]->SetMinimum(h_min*0.95); 
          else h_WF_channel4[vec_run.at(0)][0]->SetMinimum(h_min*1.05); 
          
          h_WF_channel4[vec_run.at(0)][0]->Draw();
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel4[ii].size(); jj++)
               h_WF_channel4[ii][jj]->Draw("same");

          c4 -> Print("WF_channel4_total.png","png");
          c4 -> Print("WF_channel4_total.pdf","pdf");

          max.clear();
          min.clear(); 
       }
  
       if(h_WF_channel5.size() != 0){

          TCanvas* c5 = new TCanvas("c5","c5");
          c5 -> cd();
 
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel5[ii].size(); jj++){
               max.push_back(h_WF_channel5[ii][jj]->GetMaximum()); 
               min.push_back(h_WF_channel5[ii][jj]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel5[vec_run.at(0)][0]->SetMaximum(h_max*1.05);
          else h_WF_channel5[vec_run.at(0)][0]->SetMaximum(h_max*0.95);
          if(h_min > 0.)h_WF_channel5[vec_run.at(0)][0]->SetMinimum(h_min*0.95); 
          else h_WF_channel5[vec_run.at(0)][0]->SetMinimum(h_min*1.05); 
          
          h_WF_channel5[vec_run.at(0)][0]->Draw();
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel5[ii].size(); jj++)
               h_WF_channel5[ii][jj]->Draw("same");

          c5 -> Print("WF_channel5_total.png","png");
          c5 -> Print("WF_channel5_total.pdf","pdf");

          max.clear();
          min.clear(); 
       }
      
       if(h_WF_channel6.size() != 0){

          TCanvas* c6 = new TCanvas("c6","c6");
          c6 -> cd();

          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel6[ii].size(); jj++){
               max.push_back(h_WF_channel6[ii][jj]->GetMaximum()); 
               min.push_back(h_WF_channel6[ii][jj]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel6[vec_run.at(0)][0]->SetMaximum(h_max*1.05);
          else h_WF_channel6[vec_run.at(0)][0]->SetMaximum(h_max*0.95);
          if(h_min > 0.)h_WF_channel6[vec_run.at(0)][0]->SetMinimum(h_min*0.95); 
          else h_WF_channel6[vec_run.at(0)][0]->SetMinimum(h_min*1.05); 
          
          h_WF_channel6[vec_run.at(0)][0]->Draw();
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel6[ii].size(); jj++)
               h_WF_channel6[ii][jj]->Draw("same");

          c6 -> Print("WF_channel6_total.png","png");
          c6 -> Print("WF_channel6_total.pdf","pdf");

          max.clear();
          min.clear(); 
       }
 
       if(h_WF_channel7.size() != 0){
       
          TCanvas* c7 = new TCanvas("c7","c7");
          c7 -> cd();

          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel7[ii].size(); jj++){
               max.push_back(h_WF_channel7[ii][jj]->GetMaximum()); 
               min.push_back(h_WF_channel7[ii][jj]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel7[vec_run.at(0)][0]->SetMaximum(h_max*1.05);
          else h_WF_channel7[vec_run.at(0)][0]->SetMaximum(h_max*0.95);
          if(h_min > 0.)h_WF_channel7[vec_run.at(0)][0]->SetMinimum(h_min*0.95); 
          else h_WF_channel7[vec_run.at(0)][0]->SetMinimum(h_min*1.05); 
          
          h_WF_channel7[vec_run.at(0)][0]->Draw();
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel7[ii].size(); jj++)
               h_WF_channel7[ii][jj]->Draw("same");

          c7 -> Print("WF_channel7_total.png","png");
          c7 -> Print("WF_channel7_total.pdf","pdf");

          max.clear();
          min.clear(); 
       }
       
       if(h_WF_channel8.size() != 0){
       
          TCanvas* c8 = new TCanvas("c8","c8");
          c8 -> cd();

          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel8[ii].size(); jj++){
               max.push_back(h_WF_channel8[ii][jj]->GetMaximum()); 
               min.push_back(h_WF_channel8[ii][jj]->GetMinimum()); 
           }
           
          std::sort(max.begin(),max.end());
          std::sort(min.begin(),min.end());
 
          float h_min = min.at(0);
          float h_max = max.at(max.size()-1);
          if(h_max > 0.) h_WF_channel8[vec_run.at(0)][0]->SetMaximum(h_max*1.05);
          else h_WF_channel8[vec_run.at(0)][0]->SetMaximum(h_max*0.95);
          if(h_min > 0.)h_WF_channel8[vec_run.at(0)][0]->SetMinimum(h_min*0.95); 
          else h_WF_channel8[vec_run.at(0)][0]->SetMinimum(h_min*1.05); 
          
          h_WF_channel8[vec_run.at(0)][0]->Draw();
          for(unsigned int ii = vec_run.at(0); ii <= vec_run.at(vec_run.size()-1); ii++)
           for(unsigned int jj = 0; jj < h_WF_channel8[ii].size(); jj++)
               h_WF_channel8[ii][jj]->Draw("same");

          c8 -> Print("WF_channel8_total.png","png");
          c8 -> Print("WF_channel8_total.pdf","pdf");

          max.clear();
          min.clear(); 
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
