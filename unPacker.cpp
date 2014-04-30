/*************************************************************

    simple program to read bin files
    compile with --->  c++ -o unPacker `root-config --cflags --glibs` unPacker.cpp
         or with --->  c++ -o unPacker unPacker.cpp `root-config --cflags --glibs`
    run with --->  ./unPacker WaveForms 0 // usage: ./unPacker directory saveWaveForm
    
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
    int saveWF = atoi(argv[2]);

    TString command = "ls "+ std::string(runFolder) + "/*.bin > input.tmp"; 
    gSystem -> Exec(command); 

    // Tree branches and tree structure
    int event;
    int run;
    vector<float> waveForm_Trigger;
    float ampMax_Trigger;
    int timeAmpMax_Trigger;
    vector<float> waveForm_channel1;
    float ampMax_channel1;
    int timeAmpMax_channel1;
    vector<float> waveForm_channel2;
    float ampMax_channel2;
    int timeAmpMax_channel2;
    vector<float> waveForm_channel3;
    float ampMax_channel3;
    int timeAmpMax_channel3;
    vector<float> waveForm_channel4;
    float ampMax_channel4;
    int timeAmpMax_channel4;
    vector<float> waveForm_channel5;
    float ampMax_channel5;
    int timeAmpMax_channel5;
    vector<float> waveForm_channel6;
    float ampMax_channel6;
    int timeAmpMax_channel6;
    vector<float> waveForm_channel7;
    float ampMax_channel7;
    int timeAmpMax_channel7;
    vector<float> waveForm_channel8;
    float ampMax_channel8;
    int timeAmpMax_channel8;

    TTree *nt = new TTree("nt","nt");
    nt->SetDirectory(0);

    nt->Branch("event",&event,"event/I"); 
    nt->Branch("run",&run,"run/I"); 
    nt->Branch("waveForm_Trigger","std::vector<float>",&waveForm_Trigger); 
    nt->Branch("ampMax_Trigger",&ampMax_Trigger,"ampMax_Trigger/F"); 
    nt->Branch("timeAmpMax_Trigger",&timeAmpMax_Trigger,"timeAmpMax_Trigger/I"); 
    nt->Branch("waveForm_channel1","std::vector<float>",&waveForm_channel1);  
    nt->Branch("ampMax_channel1",&ampMax_channel1,"ampMax_channel1/F"); 
    nt->Branch("timeAmpMax_channel1",&timeAmpMax_channel1,"timeAmpMax_channel1/I"); 
    nt->Branch("waveForm_channel2","std::vector<float>",&waveForm_channel2); 
    nt->Branch("ampMax_channel2",&ampMax_channel2,"ampMax_channel2/F"); 
    nt->Branch("timeAmpMax_channel2",&timeAmpMax_channel2,"timeAmpMax_channel2/I"); 
    nt->Branch("waveForm_channel3","std::vector<float>",&waveForm_channel3); 
    nt->Branch("ampMax_channel3",&ampMax_channel3,"ampMax_channel3/F"); 
    nt->Branch("timeAmpMax_channel3",&timeAmpMax_channel3,"timeAmpMax_channel3/I"); 
    nt->Branch("waveForm_channel4","std::vector<float>",&waveForm_channel4); 
    nt->Branch("ampMax_channel4",&ampMax_channel4,"ampMax_channel4/F"); 
    nt->Branch("timeAmpMax_channel4",&timeAmpMax_channel4,"timeAmpMax_channel4/I"); 
    nt->Branch("waveForm_channel5","std::vector<float>",&waveForm_channel5);
    nt->Branch("ampMax_channel5",&ampMax_channel5,"ampMax_channel5/F");  
    nt->Branch("timeAmpMax_channel5",&timeAmpMax_channel5,"timeAmpMax_channel5/I"); 
    nt->Branch("waveForm_channel6","std::vector<float>",&waveForm_channel6); 
    nt->Branch("ampMax_channel6",&ampMax_channel6,"ampMax_channel6/F"); 
    nt->Branch("timeAmpMax_channel6",&timeAmpMax_channel6,"timeAmpMax_channel6/I"); 
    nt->Branch("waveForm_channel7","std::vector<float>",&waveForm_channel7); 
    nt->Branch("ampMax_channel7",&ampMax_channel7,"ampMax_channel7/F"); 
    nt->Branch("timeAmpMax_channel7",&timeAmpMax_channel7,"timeAmpMax_channel7/I"); 
    nt->Branch("waveForm_channel8","std::vector<float>",&waveForm_channel8);
    nt->Branch("ampMax_channel8",&ampMax_channel8,"ampMax_channel8/F"); 
    nt->Branch("timeAmpMax_channel8",&timeAmpMax_channel8,"timeAmpMax_channel8/I"); 
    
    // loop over the file and read data
    ifstream inFile("input.tmp");
    string line;
    int ifile=0;

    std::vector<float> tmp;
    std::map<float,int> map_AmpTm;
    
    bool isNegative[9];
    for(int ii = 0; ii < 9; ii++)
        isNegative[ii] = true;
    isNegative[3] = false;

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

        cout << "Event: " << event << endl;

        ampMax_Trigger = 0.;
        ampMax_channel1 = 0.;
        ampMax_channel2 = 0.;
        ampMax_channel3 = 0.;
        ampMax_channel4 = 0.;
        ampMax_channel5 = 0.;
        ampMax_channel6 = 0.;
        ampMax_channel7 = 0.;
        ampMax_channel8 = 0.;

        timeAmpMax_Trigger = 0;
        timeAmpMax_channel1 = 0;
        timeAmpMax_channel2 = 0;
        timeAmpMax_channel3 = 0;
        timeAmpMax_channel4 = 0;
        timeAmpMax_channel5 = 0;
        timeAmpMax_channel6 = 0;
        timeAmpMax_channel7 = 0;
        timeAmpMax_channel8 = 0;

        for(int iCh=0; iCh<=nCh; iCh++)
        {   
            for(int iSample=0; iSample<nSize; iSample++){
                tmp.push_back(channels[iCh].at(iSample+nSize*i));
                map_AmpTm[channels[iCh].at(iSample+nSize*i)] = iSample;
            }

            std::sort(tmp.begin(),tmp.end());

            if(iCh == 0)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_Trigger.push_back(channels[iCh].at(iSample+nSize*i));

               if(isNegative[iCh] == true){
                  ampMax_Trigger = tmp.at(0);
                  timeAmpMax_Trigger = 0.2*map_AmpTm[tmp.at(0)]; 
               }   
               else{
                  ampMax_Trigger = tmp.at(tmp.size()-1);
                  timeAmpMax_Trigger = 0.2*map_AmpTm[tmp.at(tmp.size()-1)]; 
               }            
            }

            if(iCh == 1)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel1.push_back(channels[iCh].at(iSample+nSize*i));

               if(isNegative[iCh] == true){
                  ampMax_channel1 = tmp.at(0);
                  timeAmpMax_channel1 = 0.2*map_AmpTm[tmp.at(0)]; 
               }   
               else{
                  ampMax_channel1 = tmp.at(tmp.size()-1);
                  timeAmpMax_channel1 = 0.2*map_AmpTm[tmp.at(tmp.size()-1)]; 
               }    
 
            }

            if(iCh == 2)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel2.push_back(channels[iCh].at(iSample+nSize*i));

               if(isNegative[iCh] == true){
                  ampMax_channel2 = tmp.at(0);
                  timeAmpMax_channel2 = 0.2*map_AmpTm[tmp.at(0)]; 
               }   
               else{
                  ampMax_channel2 = tmp.at(tmp.size()-1);
                  timeAmpMax_channel2 = 0.2*map_AmpTm[tmp.at(tmp.size()-1)]; 
               }   

            }
             
            if(iCh == 3)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel3.push_back(channels[iCh].at(iSample+nSize*i));

               if(isNegative[iCh] == true){
                  ampMax_channel3 = tmp.at(0);
                  timeAmpMax_channel3 = 0.2*map_AmpTm[tmp.at(0)]; 
               }   
               else{
                  ampMax_channel3 = tmp.at(tmp.size()-1);
                  timeAmpMax_channel3 = 0.2*map_AmpTm[tmp.at(tmp.size()-1)]; 
               }   
               
            }

            if(iCh == 4)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel4.push_back(channels[iCh].at(iSample+nSize*i));

               if(isNegative[iCh] == true){
                  ampMax_channel4 = tmp.at(0);
                  timeAmpMax_channel4 = 0.2*map_AmpTm[tmp.at(0)]; 
               }   
               else{
                  ampMax_channel4 = tmp.at(tmp.size()-1);
                  timeAmpMax_channel4 = 0.2*map_AmpTm[tmp.at(tmp.size()-1)]; 
               }   

            }

            if(iCh == 5)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel5.push_back(channels[iCh].at(iSample+nSize*i));

               if(isNegative[iCh] == true){
                  ampMax_channel5 = tmp.at(0);
                  timeAmpMax_channel5 = 0.2*map_AmpTm[tmp.at(0)]; 
               }   
               else{
                  ampMax_channel5 = tmp.at(tmp.size()-1);
                  timeAmpMax_channel5 = 0.2*map_AmpTm[tmp.at(tmp.size()-1)]; 
               } 

            }

            if(iCh == 6)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel6.push_back(channels[iCh].at(iSample+nSize*i));

               if(isNegative[iCh] == true){
                  ampMax_channel6 = tmp.at(0);
                  timeAmpMax_channel6 = 0.2*map_AmpTm[tmp.at(0)]; 
               }   
               else{
                  ampMax_channel6 = tmp.at(tmp.size()-1);
                  timeAmpMax_channel6 = 0.2*map_AmpTm[tmp.at(tmp.size()-1)]; 
               }    

            }

            if(iCh == 7)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel7.push_back(channels[iCh].at(iSample+nSize*i));

               if(isNegative[iCh] == true){
                  ampMax_channel7 = tmp.at(0);
                  timeAmpMax_channel7 = 0.2*map_AmpTm[tmp.at(0)]; 
               }   
               else{
                  ampMax_channel7 = tmp.at(tmp.size()-1);
                  timeAmpMax_channel7 = 0.2*map_AmpTm[tmp.at(tmp.size()-1)]; 
               }   

            }

            if(iCh == 8)
            {
               for(int iSample=0; iSample<nSize; iSample++)
                   if(saveWF == 1) waveForm_channel8.push_back(channels[iCh].at(iSample+nSize*i));

               if(isNegative[iCh] == true){
                  ampMax_channel8 = tmp.at(0);
                  timeAmpMax_channel8 = 0.2*map_AmpTm[tmp.at(0)]; 
               }   
               else{
                  ampMax_channel8 = tmp.at(tmp.size()-1);
                  timeAmpMax_channel8 = 0.2*map_AmpTm[tmp.at(tmp.size()-1)]; 
               }   
          
            }
 
            tmp.clear();
            
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
    
