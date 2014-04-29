/*************************************************************

    simple program to read bin files
    compile with --->  c++ -o unPacker `root-config --cflags --glibs` unPacker.cpp
    run with --->  ./unPacker WaveForms
    
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

#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

int main(int argc, char** argv)
{

    gROOT->ProcessLine("#include <vector>");

    char* runFolder = argv[1];
    TString command = "ls "+ std::string(runFolder) + "/*.bin > input.tmp"; 
    //std::cout  << command << std::endl;
    gSystem -> Exec(command); 

    // Tree branches and tree structure
    int event;
    int run;
    vector<float> waveForm_channel0;
    float ampMax_channel0;
    vector<float> waveForm_channel1;
    float ampMax_channel1;
    vector<float> waveForm_channel2;
    float ampMax_channel2;
    vector<float> waveForm_channel3;
    float ampMax_channel3;
    vector<float> waveForm_channel4;
    float ampMax_channel4;
    vector<float> waveForm_channel5;
    float ampMax_channel5;
    vector<float> waveForm_channel6;
    float ampMax_channel6;
    vector<float> waveForm_channel7;
    float ampMax_channel7;
    vector<float> waveForm_channel8;
    float ampMax_channel8;

    TTree *nt = new TTree("nt","nt");
    nt->SetDirectory(0);

    nt->Branch("event",&event,"event/I"); 
    nt->Branch("run",&run,"run/I"); 
    nt->Branch("waveForm_channel0","std::vector<float>",&waveForm_channel0); 
    nt->Branch("ampMax_channel0",&ampMax_channel0,"ampMax_channel0/F"); 
    nt->Branch("waveForm_channel1","std::vector<float>",&waveForm_channel1);  
    nt->Branch("ampMax_channel1",&ampMax_channel1,"ampMax_channel1/F"); 
    nt->Branch("waveForm_channel2","std::vector<float>",&waveForm_channel2); 
    nt->Branch("ampMax_channel2",&ampMax_channel2,"ampMax_channel2/F"); 
    nt->Branch("waveForm_channel3","std::vector<float>",&waveForm_channel3); 
    nt->Branch("ampMax_channel3",&ampMax_channel3,"ampMax_channel3/F"); 
    nt->Branch("waveForm_channel4","std::vector<float>",&waveForm_channel4); 
    nt->Branch("ampMax_channel4",&ampMax_channel4,"ampMax_channel4/F"); 
    nt->Branch("waveForm_channel5","std::vector<float>",&waveForm_channel5);
    nt->Branch("ampMax_channel5",&ampMax_channel5,"ampMax_channel5/F");  
    nt->Branch("waveForm_channel6","std::vector<float>",&waveForm_channel6); 
    nt->Branch("ampMax_channel6",&ampMax_channel6,"ampMax_channel6/F"); 
    nt->Branch("waveForm_channel7","std::vector<float>",&waveForm_channel7); 
    nt->Branch("ampMax_channel7",&ampMax_channel7,"ampMax_channel7/F"); 
    nt->Branch("waveForm_channel8","std::vector<float>",&waveForm_channel8);
    nt->Branch("ampMax_channel8",&ampMax_channel8,"ampMax_channel8/F"); 
    
    // loop over the file and read data
    ifstream inFile("input.tmp");
    string line;
    int ifile=0;

    while(getline(inFile,line)){

      ifile++;
      std::cout << "Reading File: " << ifile << " - " << line << std::endl;

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

        ampMax_channel0 = 0.;
        ampMax_channel1 = 0.;
        ampMax_channel2 = 0.;
        ampMax_channel3 = 0.;
        ampMax_channel4 = 0.;
        ampMax_channel5 = 0.;
        ampMax_channel6 = 0.;
        ampMax_channel7 = 0.;
        ampMax_channel8 = 0.;

        for(int iCh=0; iCh<=nCh; iCh++)
        {   
            std::sort(channels[iCh].begin(),channels[iCh].end());

            if(iCh == 0)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   waveForm_channel0.push_back(channels[iCh].at(iSample+nSize*i));

               ampMax_channel0 = channels[iCh].at(0);
               
            }

            if(iCh == 1)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   waveForm_channel1.push_back(channels[iCh].at(iSample+nSize*i));

               ampMax_channel1 = channels[iCh].at(0);
 
            }

            if(iCh == 2)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   waveForm_channel2.push_back(channels[iCh].at(iSample+nSize*i));

               ampMax_channel2 = channels[iCh].at(0);

            }
             
            if(iCh == 3)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   waveForm_channel3.push_back(channels[iCh].at(iSample+nSize*i));

               ampMax_channel3 = channels[iCh].at(0);  
               
            }

            if(iCh == 4)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   waveForm_channel4.push_back(channels[iCh].at(iSample+nSize*i));

               ampMax_channel4 = channels[iCh].at(0);

            }

            if(iCh == 5)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   waveForm_channel5.push_back(channels[iCh].at(iSample+nSize*i));

               ampMax_channel5 = channels[iCh].at(0);

            }

            if(iCh == 6)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   waveForm_channel6.push_back(channels[iCh].at(iSample+nSize*i));

               ampMax_channel6 = channels[iCh].at(0);

            }

            if(iCh == 7)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   waveForm_channel7.push_back(channels[iCh].at(iSample+nSize*i));

               ampMax_channel7 = channels[iCh].at(0);

            }

            if(iCh == 8)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   waveForm_channel8.push_back(channels[iCh].at(iSample+nSize*i));

               ampMax_channel8 = channels[iCh].at(0);
          
            }
            
        }

        for(int iCh=nCh+1; iCh<=8; iCh++){

            if(iCh == 0) waveForm_channel0.push_back(0.);
            if(iCh == 1) waveForm_channel1.push_back(0.);
            if(iCh == 2) waveForm_channel2.push_back(0.);
            if(iCh == 3) waveForm_channel3.push_back(0.);
            if(iCh == 4) waveForm_channel4.push_back(0.);
            if(iCh == 5) waveForm_channel5.push_back(0.);
            if(iCh == 6) waveForm_channel6.push_back(0.);
            if(iCh == 7) waveForm_channel7.push_back(0.);
            if(iCh == 8) waveForm_channel8.push_back(0.);
        }
        
        nt->Fill();   
      }

    }

    TFile *f = new TFile((std::string(runFolder)+"_tree.root").c_str(),"RECREATE"); 
    f->cd();
    nt->Write("nt");
    f->Close();

    gSystem -> Exec("rm input.tmp"); 
}
    
