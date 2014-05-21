#ifndef __SetOutputTree__
#define __SetOutputTree__

#include "TTree.h"
#include "TROOT.h"

using namespace std;
  
    // Tree branches and tree structure
    int event;
    int run;
    int nPointsInterpolation_output;
    int adcSciFront;
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
    vector<float> waveForm_channel0;
    float baseLine_channel0;
    float baselineDispersion_channel0;
    float baseLineGlobal_channel0;
    float baselineDispersionGlobal_channel0;
    float ampMax_channel0;
    float timeAmpMax_channel0;
    float constFrac_channel0;
    float timeConstFrac_channel0;
    float slopeConstFrac_channel0;
    float chi2ConstFrac_channel0;
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

    void outTree(TTree* nt, const int& triggerChannel){
      nt->Branch("event",&event,"event/I"); 
      nt->Branch("run",&run,"run/I"); 
      nt->Branch("nPointsInterpolation",&nPointsInterpolation_output,"nPointsInterpolation/I"); 
      nt->Branch("adcSciFront",&adcSciFront,"adcSciFront/F"); 
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
      if(triggerChannel != 0){
         nt->Branch("waveForm_channel0","std::vector<float>",&waveForm_channel0);  
         nt->Branch("baseLine_channel0",&baseLine_channel0,"baseLine_channel0/F"); 
         nt->Branch("baselineDispersion_channel0",&baselineDispersion_channel0,"baselineDispersion_channel0/F"); 
         nt->Branch("baseLineGlobal_channel0",&baseLineGlobal_channel0,"baseLineGlobal_channel0/F");  
         nt->Branch("baselineDispersionGlobal_channel0",&baselineDispersionGlobal_channel0,"baselineDispersionGlobal_channel0/F"); 
         nt->Branch("ampMax_channel0",&ampMax_channel0,"ampMax_channel0/F"); 
         nt->Branch("timeAmpMax_channel0",&timeAmpMax_channel0,"timeAmpMax_channel0/F"); 
         nt->Branch("constFrac_channel0",&constFrac_channel0,"constFrac_channel0/F"); 
         nt->Branch("timeConstFrac_channel0",&timeConstFrac_channel0,"timeConstFrac_channel0/F"); 
         nt->Branch("slopeConstFrac_channel0",&slopeConstFrac_channel0,"slopeConstFrac_channel0/F"); 
         nt->Branch("chi2ConstFrac_channel0",&chi2ConstFrac_channel0,"chi2ConstFrac_channel0/F"); 
      }
      if(triggerChannel != 1){
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
      }
      if(triggerChannel != 2){
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
      }
      if(triggerChannel != 3){
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
      }
      if(triggerChannel != 4){
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
      }
      if(triggerChannel != 5){
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
      }
      if(triggerChannel != 6){
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
      }
      if(triggerChannel != 7){
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
      }
      if(triggerChannel != 8){
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
      }
    
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
    }

#endif
