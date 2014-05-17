/*************************************************************

    simple program to read bin files
    compile with --->  g++ -Wall -o makeAnalysis_unbinned `root-config --cflags --glibs` -lRooFit -lRooFitCore -lMinuit2 makeAnalysis_unbinned.cpp
         or with --->  g++ -Wall -o makeAnalysis_unbinned makeAnalysis_unbinned.cpp `root-config --cflags --glibs` -lRooFit -lRooFitCore -lMinuit2
    run with --->  ./makeAnalysis file.root 0 // usage: ./makeAnalysis file.root saveWFHistos thresAmp[1] thresAmp[2] thresAmp[3]
    
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
#include "TApplication.h"

#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooExponential.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooPolynomial.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif 

#include "InitTree.h"

using namespace std;
using namespace RooFit ;

const float BinToTime = 0.2;

float computeMean(std::vector<float> sample);
float computeDispersion(std::vector<float> sample);
void MakeFit(RooRealVar* deltaT,RooDataSet* data, float& sigma, float& sigmaErr, float& Mean, std::string Name);

int main(int argc, char** argv)
{
    gROOT->ProcessLine("#include <vector>");

    char* inputFile = argv[1];

    int saveWFHistos = atoi(argv[2]);

    float thresAmp[4];

    thresAmp[0] = -1.;
    thresAmp[1] = 100.;
    thresAmp[2] = 100.;
    thresAmp[3] = 35.;
    
    if(argc >= 4) thresAmp[1] = atof(argv[3]);
    if(argc >= 5) thresAmp[2] = atof(argv[4]);
    if(argc >= 6) thresAmp[3] = atof(argv[5]);
    
    std::cout << "inputFile    = " << inputFile << std::endl;
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

    
    std::vector<float> vec_double;
    std::vector<float> vec_triple;
    std::vector<float> vec_double13;
    std::vector<float> vec_double23;

    std::map<int,std::map<int,TH1F*> > h_WF_Trigger;
    std::map<int,std::map<int,TH1F*> > h_WF_channel1; 
    std::map<int,std::map<int,TH1F*> > h_WF_channel2;
    std::map<int,std::map<int,TH1F*> > h_WF_channel3;

    std::vector<std::pair<int,int> > vec_run_event;
    std::pair<int,int> tmp_pair;

    RooRealVar* deltaTDouble = new RooRealVar("deltaTDouble","deltaTDouble",0.,30.);  
    RooDataSet* dataDouble = new RooDataSet("dataDouble","dataset with deltaTDouble",RooArgSet(*deltaTDouble));
 
    RooRealVar* deltaTTriple = new RooRealVar("deltaTTriple","deltaTTriple",0.,30.);  
    RooDataSet* dataTriple = new RooDataSet("dataTriple","dataset with deltaTTriple",RooArgSet(*deltaTTriple));

    RooRealVar* deltaTDouble13 = new RooRealVar("deltaTDouble13","deltaTDouble13",0.,30.);  
    RooDataSet* dataDouble13 = new RooDataSet("dataDouble13","dataset with deltaTDouble13",RooArgSet(*deltaTDouble13));

   RooRealVar* deltaTDouble23 = new RooRealVar("deltaTDouble23","deltaTDouble23",0.,30.);  
    RooDataSet* dataDouble23 = new RooDataSet("dataDouble23","dataset with deltaTDouble23",RooArgSet(*deltaTDouble23));
    
    char histoName[200];
   
    for(int ientry = 0; ientry < nt->GetEntries(); ientry++){
        if(ientry%100==0) std::cout<<"--- Reading entry = "<< ientry << " for filling the histos..." <<std::endl;
        nt->GetEntry(ientry);

        if(ampMax_channel1 <= thresAmp[1]) continue;
        if(ampMax_channel2 <= thresAmp[2]) continue;

        if(timeAmpMax_channel1/BinToTime < 2 || timeAmpMax_channel1/BinToTime > waveForm_channel1->size()-2) continue;
        if(timeAmpMax_channel2/BinToTime < 2 || timeAmpMax_channel2/BinToTime > waveForm_channel2->size()-2) continue;

        if((timeConstFrac_channel2-timeConstFrac_channel1) < 0.) continue;

        vec_double.push_back(timeConstFrac_channel2-timeConstFrac_channel1);
        deltaTDouble->setVal(timeConstFrac_channel2-timeConstFrac_channel1);
        dataDouble->add(RooArgSet(*deltaTDouble));
        
        if(ampMax_channel3 > thresAmp[3] && timeAmpMax_channel3/BinToTime >= 2 && timeAmpMax_channel3/BinToTime <= waveForm_channel3->size()-2){
      
           vec_triple.push_back((timeConstFrac_channel2+timeConstFrac_channel1)/2.-timeConstFrac_channel3);
           deltaTTriple->setVal((timeConstFrac_channel2+timeConstFrac_channel1)/2.-timeConstFrac_channel3);
           dataTriple->add(RooArgSet(*deltaTTriple));
        
           vec_double13.push_back(timeConstFrac_channel1-timeConstFrac_channel3);
           deltaTDouble13->setVal(timeConstFrac_channel1-timeConstFrac_channel3);
           dataDouble13->add(RooArgSet(*deltaTDouble13));
        
           vec_double23.push_back(timeConstFrac_channel2-timeConstFrac_channel3);
           deltaTDouble23->setVal(timeConstFrac_channel2-timeConstFrac_channel3);
           dataDouble23->add(RooArgSet(*deltaTDouble23));
        
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

    float sigma_diff = 0.;
    float sigma_diffE = 0.;
    float mean_diff = 0.;
    MakeFit(deltaTDouble,dataDouble,sigma_diff,sigma_diffE,mean_diff,std::string("Double"));

    float sigma_mcp = 0.;
    float sigma_mcpE = 0.;
    float mean_mcp = 0.;
    MakeFit(deltaTTriple,dataTriple,sigma_mcp,sigma_mcpE,mean_mcp,std::string("Triple"));

    float sigma_t2dat0 = 0.;
    float sigma_t2dat0E = 0.;
    float mean_t2dat0 = 0.;
    MakeFit(deltaTDouble13,dataDouble13,sigma_t2dat0,sigma_t2dat0E,mean_t2dat0,std::string("Double13"));

    float sigma_t2dat1 = 0.;
    float sigma_t2dat1E = 0.;
    float mean_t2dat1 = 0.;
    MakeFit(deltaTDouble23,dataDouble23,sigma_t2dat1,sigma_t2dat1E,mean_t2dat1,std::string("Double23"));

    sigma_double = sigma_diff;
    sigmaError_double = sigma_diffE; 
    mean_double = mean_diff;

    sigma_triple = sigma_mcp;
    sigmaError_triple = sigma_mcpE; 
    mean_triple = mean_mcp;

    sigma_double13 = sigma_t2dat0;
    sigmaError_double13 = sigma_t2dat0E; 
    mean_double13 = mean_t2dat0;

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
    
    std::vector<float> vec_double_fake;
    std::vector<float> vec_triple_fake;
   
    for(int ientry = 0; ientry < nt->GetEntries(); ientry++){
        if(ientry%100==0) std::cout<<"--- Reading entry = "<< ientry  << " for computing efficiencies..."<<std::endl;
        nt->GetEntry(ientry);
        
        if(ampMax_channel1 <= thresAmp[1]) continue;
        if(ampMax_channel2 <= thresAmp[2]) continue;

        if(timeAmpMax_channel1/BinToTime < 2 || timeAmpMax_channel1/BinToTime > waveForm_channel1->size()-2) continue;
        if(timeAmpMax_channel2/BinToTime < 2 || timeAmpMax_channel2/BinToTime > waveForm_channel2->size()-2) continue;

        if(fabs(timeConstFrac_channel2-timeConstFrac_channel1 - mean_diff) < 5*sigma_diff) nDouble++;
        else vec_double_fake.push_back(timeConstFrac_channel2-timeConstFrac_channel1);
       
        if(ampMax_channel3 <= thresAmp[3]) continue;
        if(timeAmpMax_channel3/BinToTime < 2 || timeAmpMax_channel3/BinToTime > waveForm_channel3->size()-2) continue;

        if(fabs((timeConstFrac_channel2+timeConstFrac_channel1)/2.- timeConstFrac_channel3 - mean_mcp) < 5*sigma_mcp) nTriple++;
        else vec_triple_fake.push_back((timeConstFrac_channel2+timeConstFrac_channel1)/2.- timeConstFrac_channel3);
        
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

    nTriple_output = (float)nTriple;
    nTriple_Fake_output = (float)nTriple_fake;
    nDouble_output = (float)nDouble;
    nDouble_Fake_output = (float)nDouble_fake;
    Efficiency = eff;
    EfficiencyError = eff_err;

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
       
        if(ampMax_channel3 > thresAmp[3] && timeAmpMax_channel3/BinToTime >= 2 && timeAmpMax_channel3/BinToTime <= waveForm_channel3->size()-2) isTriple = 1;
          

        isDouble_output = isDouble;
        isTriple_output = isTriple;

        nt_output->Fill();
        
    }

    std::cout << std::endl;
    
    TFile* f1 = new TFile(("analyzed_"+std::string(inputFile)).c_str(),"RECREATE");
    f1->cd();

    nt_output->Write();

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

void MakeFit(RooRealVar* deltaT,RooDataSet* data, float& sigma, float& sigmaErr, float& Mean, std::string Name)
{
   RooRealVar mean("mean","mean",10.,0.,100.);
   RooRealVar width("width","width",4.,0.,100.);
   RooGaussian Gauss("Gauss", "Gaussian PDF",*deltaT,mean,width);
  
   RooPolynomial pol0("pol0", "pol0 PDF",*deltaT,RooArgList());

   RooRealVar nsig("nsig","#signal events",0.5,0.,1.) ;
   RooAddPdf sum("sum","g+a",RooArgList(Gauss,pol0),RooArgList(nsig)) ; 
   
   RooFitResult* res;
   res = sum.fitTo(*data,Save(kTRUE),SumW2Error(kTRUE),Minimizer("Minuit2"));

   std::cout << std::endl;
   std::cout << std::endl;
   std::cout << "MakeFit: " << Name << std::endl;
   
   float MEAN = mean.getVal();
   float MEANError = mean.getError();
   std::cout << "MEAN = " << MEAN << " +/- " << MEANError << std::endl;

   Mean = MEAN;

   float WIDTH = width.getVal();
   float WIDTHError = width.getError();
   std::cout << "WIDTH = " << WIDTH << " +/- " << WIDTHError << std::endl;

   sigma = WIDTH;
   sigmaErr = WIDTHError;

   float frac = nsig.getVal();
   float fracError = nsig.getError();
   std::cout << "NSIG = " << frac << " +/- " << fracError << std::endl;

   std::cout << std::endl;
   std::cout << std::endl;
   
   RooPlot* xframe = deltaT->frame(Bins(200));
   data->plotOn(xframe); 
   sum.plotOn(xframe,VisualizeError(*res,1,kTRUE),DrawOption("F"),FillColor(kOrange),VLines());
   sum.plotOn(xframe);
   
   xframe->GetYaxis()->SetRangeUser(0.0000001,100.);
   
   TCanvas* c = new TCanvas("c","c", 800,600);
   xframe->Draw();

   c->SaveAs((Name+".root").c_str());
   delete c;

}

