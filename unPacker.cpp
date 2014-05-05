/*************************************************************

    simple program to read bin files
    compile with --->  c++ -o unPacker `root-config --cflags --glibs` unPacker.cpp
         or with --->  c++ -o unPacker unPacker.cpp `root-config --cflags --glibs`
    run with --->  ./unPacker WaveForms 36 5 1 // usage: ./unPacker directory fraction baseline-time saveWaveForm
    
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
#include "TTree.h"
#include "TFile.h"

using namespace std;

const float BinToTime = 0.2;

float findTimeConstFrac(float x1, float y1, float x2, float y2, float x3, float y3, float frac, float amp);
float computeMean(std::vector<float> sample);

int main(int argc, char** argv)
{

    gROOT->ProcessLine("#include <vector>");

    char* runFolder = argv[1];
    float ampFraction = atof(argv[2]);
    float timeBaseLine = atof(argv[3]);
    int saveWF = atoi(argv[4]);

    int binBaseLine = int(timeBaseLine/BinToTime);

    TString command = "ls "+ std::string(runFolder) + "/*.bin > input.tmp"; 
    gSystem -> Exec(command); 

    // Tree branches and tree structure
    int event;
    int run;
    float ampFraction_output;
    float timeBaseLine_output;

    vector<float> waveForm_Trigger;
    float baseLine_Trigger;
    float ampMax_Trigger;
    float timeAmpMax_Trigger;
    float timeConstFrac_Trigger;
    vector<float> waveForm_channel1;
    float baseLine_channel1;
    float ampMax_channel1;
    float timeAmpMax_channel1;
    float timeConstFrac_channel1;
    vector<float> waveForm_channel2;
    float baseLine_channel2;
    float ampMax_channel2;
    float timeAmpMax_channel2;
    float timeConstFrac_channel2;
    vector<float> waveForm_channel3;
    float baseLine_channel3;
    float ampMax_channel3;
    float timeAmpMax_channel3;
    float timeConstFrac_channel3;
    vector<float> waveForm_channel4;
    float baseLine_channel4;
    float ampMax_channel4;
    float timeAmpMax_channel4;
    float timeConstFrac_channel4;
    vector<float> waveForm_channel5;
    float baseLine_channel5;
    float ampMax_channel5;
    float timeAmpMax_channel5;
    float timeConstFrac_channel5;
    vector<float> waveForm_channel6;
    float baseLine_channel6;
    float ampMax_channel6;
    float timeAmpMax_channel6;
    float timeConstFrac_channel6;
    vector<float> waveForm_channel7;
    float baseLine_channel7;
    float ampMax_channel7;
    float timeAmpMax_channel7;
    float timeConstFrac_channel7;
    vector<float> waveForm_channel8;
    float baseLine_channel8;
    float ampMax_channel8;
    float timeAmpMax_channel8;
    float timeConstFrac_channel8;

    TTree *nt = new TTree("nt","nt");
    nt->SetDirectory(0);

    nt->Branch("event",&event,"event/I"); 
    nt->Branch("run",&run,"run/I"); 
    nt->Branch("ampFraction",&ampFraction_output,"ampFraction/F"); 
    nt->Branch("timeBaseLine",&timeBaseLine_output,"timeBaseLine/F"); 
    nt->Branch("waveForm_Trigger","std::vector<float>",&waveForm_Trigger); 
    nt->Branch("baseLine_Trigger",&baseLine_Trigger,"baseLine_Trigger/F"); 
    nt->Branch("ampMax_Trigger",&ampMax_Trigger,"ampMax_Trigger/F"); 
    nt->Branch("timeAmpMax_Trigger",&timeAmpMax_Trigger,"timeAmpMax_Trigger/F"); 
    nt->Branch("timeConstFrac_Trigger",&timeConstFrac_Trigger,"timeConstFrac_Trigger/F"); 
    nt->Branch("waveForm_channel1","std::vector<float>",&waveForm_channel1);  
    nt->Branch("baseLine_channel1",&baseLine_channel1,"baseLine_channel1/F"); 
    nt->Branch("ampMax_channel1",&ampMax_channel1,"ampMax_channel1/F"); 
    nt->Branch("timeAmpMax_channel1",&timeAmpMax_channel1,"timeAmpMax_channel1/F"); 
    nt->Branch("timeConstFrac_channel1",&timeConstFrac_channel1,"timeConstFrac_channel1/F"); 
    nt->Branch("waveForm_channel2","std::vector<float>",&waveForm_channel2); 
    nt->Branch("baseLine_channel2",&baseLine_channel2,"baseLine_channel2/F"); 
    nt->Branch("ampMax_channel2",&ampMax_channel2,"ampMax_channel2/F"); 
    nt->Branch("timeAmpMax_channel2",&timeAmpMax_channel2,"timeAmpMax_channel2/F"); 
    nt->Branch("timeConstFrac_channel2",&timeConstFrac_channel2,"timeConstFrac_channel2/F");
    nt->Branch("waveForm_channel3","std::vector<float>",&waveForm_channel3); 
    nt->Branch("baseLine_channel3",&baseLine_channel3,"baseLine_channel3/F"); 
    nt->Branch("ampMax_channel3",&ampMax_channel3,"ampMax_channel3/F"); 
    nt->Branch("timeAmpMax_channel3",&timeAmpMax_channel3,"timeAmpMax_channel3/F"); 
    nt->Branch("timeConstFrac_channel3",&timeConstFrac_channel3,"timeConstFrac_channel3/F");
    nt->Branch("waveForm_channel4","std::vector<float>",&waveForm_channel4); 
    nt->Branch("baseLine_channel4",&baseLine_channel4,"baseLine_channel4/F"); 
    nt->Branch("ampMax_channel4",&ampMax_channel4,"ampMax_channel4/F"); 
    nt->Branch("timeAmpMax_channel4",&timeAmpMax_channel4,"timeAmpMax_channel4/F"); 
    nt->Branch("timeConstFrac_channel4",&timeConstFrac_channel4,"timeConstFrac_channel4/F");
    nt->Branch("waveForm_channel5","std::vector<float>",&waveForm_channel5);
    nt->Branch("baseLine_channel5",&baseLine_channel5,"baseLine_channel5/F"); 
    nt->Branch("ampMax_channel5",&ampMax_channel5,"ampMax_channel5/F");  
    nt->Branch("timeAmpMax_channel5",&timeAmpMax_channel5,"timeAmpMax_channel5/F"); 
    nt->Branch("timeConstFrac_channel5",&timeConstFrac_channel5,"timeConstFrac_channel5/F");
    nt->Branch("waveForm_channel6","std::vector<float>",&waveForm_channel6); 
    nt->Branch("baseLine_channel6",&baseLine_channel6,"baseLine_channel6/F"); 
    nt->Branch("ampMax_channel6",&ampMax_channel6,"ampMax_channel6/F"); 
    nt->Branch("timeAmpMax_channel6",&timeAmpMax_channel6,"timeAmpMax_channel6/F"); 
    nt->Branch("timeConstFrac_channel6",&timeConstFrac_channel6,"timeConstFrac_channel6/F");
    nt->Branch("waveForm_channel7","std::vector<float>",&waveForm_channel7); 
    nt->Branch("baseLine_channel7",&baseLine_channel7,"baseLine_channel7/F"); 
    nt->Branch("ampMax_channel7",&ampMax_channel7,"ampMax_channel7/F"); 
    nt->Branch("timeAmpMax_channel7",&timeAmpMax_channel7,"timeAmpMax_channel7/F"); 
    nt->Branch("timeConstFrac_channel7",&timeConstFrac_channel7,"timeConstFrac_channel7/F");
    nt->Branch("waveForm_channel8","std::vector<float>",&waveForm_channel8);
    nt->Branch("baseLine_channel8",&baseLine_channel8,"baseLine_channel8/F"); 
    nt->Branch("ampMax_channel8",&ampMax_channel8,"ampMax_channel8/F"); 
    nt->Branch("timeAmpMax_channel8",&timeAmpMax_channel8,"timeAmpMax_channel8/F"); 
    nt->Branch("timeConstFrac_channel8",&timeConstFrac_channel8,"timeConstFrac_channel8/F");
    
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

    std::vector<float> tmp;
    std::vector<float> tmp_baseline;
    std::map<float,int> map_AmpTm;
    
    bool isNegative[9];
    for(int ii = 0; ii < 9; ii++)
        isNegative[ii] = true;
    isNegative[3] = false;

    while(getline(inFile,line)){

      ifile++;
      std::cout << "Reading File: " << ifile << " - " << line << std::endl;

      ampFraction_output = ampFraction;
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

      BinHeader = (int*) malloc (sizeof(int)*3);
      fread(BinHeader, sizeof(int), 3, input);
      iEvent = BinHeader[0];
      nSize = BinHeader[1];
      nCh = BinHeader[2];
      
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
        fread(BinHeader, sizeof(int), 3, input);
        iEvent = BinHeader[0];
      }

      //Fill tree
      for(int i=0; i<nEvent; i++)
      {

        event = eventNumber.at(i);

        cout << "Event: " << event << endl;

        baseLine_Trigger = 0.;
        baseLine_channel1 = 0.;
        baseLine_channel2 = 0.;
        baseLine_channel3 = 0.;
        baseLine_channel4 = 0.;
        baseLine_channel5 = 0.;
        baseLine_channel6 = 0.;
        baseLine_channel7 = 0.;
        baseLine_channel8 = 0.;

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

        timeConstFrac_Trigger = 0.;          
        timeConstFrac_channel1 = 0.;
        timeConstFrac_channel2 = 0.;
        timeConstFrac_channel3 = 0.;
        timeConstFrac_channel4 = 0.;
        timeConstFrac_channel5 = 0.;
        timeConstFrac_channel6 = 0.;
        timeConstFrac_channel7 = 0.;
        timeConstFrac_channel8 = 0.;

        for(int iCh=0; iCh<=nCh; iCh++)
        {   
            for(int iSample=0; iSample<nSize; iSample++){

                if(iSample <= binBaseLine) tmp_baseline.push_back(channels[iCh].at(iSample+nSize*i));

                tmp.push_back(channels[iCh].at(iSample+nSize*i));
                map_AmpTm[channels[iCh].at(iSample+nSize*i)] = iSample;
            }

            float baseline = computeMean(tmp_baseline);

            std::sort(tmp.begin(),tmp.end());
            
            float tmpAmp = 0.;
            int tmpTime = 0;
            float tmpTimeConstFrac = 0.;

            if(isNegative[iCh] == true){
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
                float absFrac = fabs(tmpAmp-baseline)*ampFraction/100.;
		if(absAmp < absFrac) break;
		ref = iSample;
            }

	    tmpTimeConstFrac = findTimeConstFrac((ref-1)*BinToTime,channels[iCh].at((ref-1)+nSize*i),ref*BinToTime,channels[iCh].at(ref+nSize*i),(ref+1)*BinToTime,channels[iCh].at((ref+1)+nSize*i),ampFraction,tmpAmp);

            if(iCh == 0)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_Trigger.push_back(channels[iCh].at(iSample+nSize*i));

               baseLine_Trigger = baseline;
               ampMax_Trigger = fabs(tmpAmp-baseline);
               timeAmpMax_Trigger = BinToTime*tmpTime;   
               timeConstFrac_Trigger = tmpTimeConstFrac;     
            }

            if(iCh == 1)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel1.push_back(channels[iCh].at(iSample+nSize*i));

               baseLine_channel1 = baseline;
               ampMax_channel1 = fabs(tmpAmp-baseline);
               timeAmpMax_channel1 = BinToTime*tmpTime; 
               timeConstFrac_channel1 = tmpTimeConstFrac;    
 
            }

            if(iCh == 2)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel2.push_back(channels[iCh].at(iSample+nSize*i));

               baseLine_channel2 = baseline;
               ampMax_channel2 = fabs(tmpAmp-baseline);
               timeAmpMax_channel2 = BinToTime*tmpTime;  
               timeConstFrac_channel2 = tmpTimeConstFrac;  

            }
             
            if(iCh == 3)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel3.push_back(channels[iCh].at(iSample+nSize*i));

               baseLine_channel3 = baseline;
               ampMax_channel3 = fabs(tmpAmp-baseline);
               timeAmpMax_channel3 = BinToTime*tmpTime;  
               timeConstFrac_channel3 = tmpTimeConstFrac;  
               
            }

            if(iCh == 4)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel4.push_back(channels[iCh].at(iSample+nSize*i));

               baseLine_channel4 = baseline;
               ampMax_channel4 = fabs(tmpAmp-baseline);
               timeAmpMax_channel4 = BinToTime*tmpTime;
               timeConstFrac_channel4 = tmpTimeConstFrac;     

            }

            if(iCh == 5)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel5.push_back(channels[iCh].at(iSample+nSize*i));

               baseLine_channel5 = baseline;
               ampMax_channel5 = fabs(tmpAmp-baseline);
               timeAmpMax_channel5 = BinToTime*tmpTime;  
               timeConstFrac_channel5 = tmpTimeConstFrac;  

            }

            if(iCh == 6)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel6.push_back(channels[iCh].at(iSample+nSize*i));

               baseLine_channel6 = baseline;
               ampMax_channel6 = fabs(tmpAmp-baseline);
               timeAmpMax_channel6 = BinToTime*tmpTime;   
               timeConstFrac_channel6 = tmpTimeConstFrac;  

            }

            if(iCh == 7)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel7.push_back(channels[iCh].at(iSample+nSize*i));

               baseLine_channel7 = baseline;
               ampMax_channel7 = fabs(tmpAmp-baseline);
               timeAmpMax_channel7 = BinToTime*tmpTime;
               timeConstFrac_channel7 = tmpTimeConstFrac;     

            }

            if(iCh == 8)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel8.push_back(channels[iCh].at(iSample+nSize*i));

               baseLine_channel8 = baseline;
               ampMax_channel8 = fabs(tmpAmp-baseline);
               timeAmpMax_channel8 = BinToTime*tmpTime;  
               timeConstFrac_channel8 = tmpTimeConstFrac;  
          
            }
 
            tmp.clear();
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

    TFile *f = new TFile((std::string(runFolder)+"_tree.root").c_str(),"RECREATE"); 
    f->cd();
    nt->Write("nt");
    f->Close();

    gSystem -> Exec("rm input.tmp"); 
}

float findTimeConstFrac(float x1, float y1, float x2, float y2, float x3, float y3, float frac, float amp)
{
  float denom = (3*(x1*x1 + x2*x2 + x3*x3) - (x1+x2+x3)*(x1+x2+x3));
  float m = (3*(x1*y1 + x2*y2 + x3*y3) - (x1+x2+x3)*(y1+y2+y3))/denom;
  float q = ((y1+y2+y3)*(x1*x1 + x2*x2 + x3*x3) - (x1+x2+x3)*(x1*y1 + x2*y2 + x3*y3))/denom;
  float time = (amp*frac/100. - q)/m; 
  return time;
}

float computeMean(std::vector<float> sample)
{
  float mean = 0.;
  
  for(unsigned int ii = 0; ii < sample.size(); ii++)
      mean = mean + sample.at(ii);

  return mean/sample.size();
}
    
