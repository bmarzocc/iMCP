#ifndef __InitTree__
#define __InitTree__

#include "TTree.h"

    // Declaration of leaf types
   int           event;
   int           run;
   int           nPointsInterpolation;
   float         timeBaseLine;
   std::vector<float>   *waveForm_Trigger;
   float         baseLine_Trigger;
   float         baselineDispersion_Trigger;
   float         baseLineGlobal_Trigger;
   float         baselineDispersionGlobal_Trigger;
   float         ampMax_Trigger;
   float         timeAmpMax_Trigger;
   float         constFrac_Trigger;
   float         timeConstFrac_Trigger;
   float         slopeConstFrac_Trigger;
   float         chi2ConstFrac_Trigger;
   std::vector<float>   *waveForm_channel1;
   float         baseLine_channel1;
   float         baselineDispersion_channel1;
   float         baseLineGlobal_channel1;
   float         baselineDispersionGlobal_channel1;
   float         ampMax_channel1;
   float         timeAmpMax_channel1;
   float         constFrac_channel1;
   float         timeConstFrac_channel1;
   float         slopeConstFrac_channel1;
   float         chi2ConstFrac_channel1;
   std::vector<float>   *waveForm_channel2;
   float         baseLine_channel2;
   float         baselineDispersion_channel2;
   float         baseLineGlobal_channel2;
   float         baselineDispersionGlobal_channel2;
   float         ampMax_channel2;
   float         timeAmpMax_channel2;
   float         constFrac_channel2;
   float         timeConstFrac_channel2;
   float         slopeConstFrac_channel2;
   float         chi2ConstFrac_channel2;
   std::vector<float>   *waveForm_channel3;
   float         baseLine_channel3;
   float         baselineDispersion_channel3;
   float         baseLineGlobal_channel3;
   float         baselineDispersionGlobal_channel3;
   float         ampMax_channel3;
   float         timeAmpMax_channel3;
   float         constFrac_channel3;
   float         timeConstFrac_channel3;
   float         slopeConstFrac_channel3;
   float         chi2ConstFrac_channel3;
   std::vector<float>   *waveForm_channel4;
   float         baseLine_channel4;
   float         baselineDispersion_channel4;
   float         baseLineGlobal_channel4;
   float         baselineDispersionGlobal_channel4;
   float         ampMax_channel4;
   float         timeAmpMax_channel4;
   float         constFrac_channel4;
   float         timeConstFrac_channel4;
   float         slopeConstFrac_channel4;
   float         chi2ConstFrac_channel4;
   std::vector<float>   *waveForm_channel5;
   float         baseLine_channel5;
   float         baselineDispersion_channel5;
   float         baseLineGlobal_channel5;
   float         baselineDispersionGlobal_channel5;
   float         ampMax_channel5;
   float         timeAmpMax_channel5;
   float         constFrac_channel5;
   float         timeConstFrac_channel5;
   float         slopeConstFrac_channel5;
   float         chi2ConstFrac_channel5;
   std::vector<float>   *waveForm_channel6;
   float         baseLine_channel6;
   float         baselineDispersion_channel6;
   float         baseLineGlobal_channel6;
   float         baselineDispersionGlobal_channel6;
   float         ampMax_channel6;
   float         timeAmpMax_channel6;
   float         constFrac_channel6;
   float         timeConstFrac_channel6;
   float         slopeConstFrac_channel6;
   float         chi2ConstFrac_channel6;
   std::vector<float>   *waveForm_channel7;
   float         baseLine_channel7;
   float         baselineDispersion_channel7;
   float         baseLineGlobal_channel7;
   float         baselineDispersionGlobal_channel7;
   float         ampMax_channel7;
   float         timeAmpMax_channel7;
   float         constFrac_channel7;
   float         timeConstFrac_channel7;
   float         slopeConstFrac_channel7;
   float         chi2ConstFrac_channel7;
   std::vector<float>   *waveForm_channel8;
   float         baseLine_channel8;
   float         baselineDispersion_channel8;
   float         baseLineGlobal_channel8;
   float         baselineDispersionGlobal_channel8;
   float         ampMax_channel8;
   float         timeAmpMax_channel8;
   float         constFrac_channel8;
   float         timeConstFrac_channel8;
   float         slopeConstFrac_channel8;
   float         chi2ConstFrac_channel8;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_nPointsInterpolation;   //!
   TBranch        *b_timeBaseLine;   //!
   TBranch        *b_waveForm_Trigger;   //!
   TBranch        *b_baseLine_Trigger;   //!
   TBranch        *b_baselineDispersion_Trigger;   //!
   TBranch        *b_baseLineGlobal_Trigger;   //!
   TBranch        *b_baselineDispersionGlobal_Trigger;   //!
   TBranch        *b_ampMax_Trigger;   //!
   TBranch        *b_timeAmpMax_Trigger;   //!
   TBranch        *b_constFrac_Trigger;   //!
   TBranch        *b_timeConstFrac_Trigger;   //!
   TBranch        *b_slopeConstFrac_Trigger;   //!
   TBranch        *b_chi2ConstFrac_Trigger;   //!
   TBranch        *b_waveForm_channel1;   //!
   TBranch        *b_baseLine_channel1;   //!
   TBranch        *b_baselineDispersion_channel1;   //!
   TBranch        *b_baseLineGlobal_channel1;   //!
   TBranch        *b_baselineDispersionGlobal_channel1;   //!
   TBranch        *b_ampMax_channel1;   //!
   TBranch        *b_timeAmpMax_channel1;   //!
   TBranch        *b_constFrac_channel1;   //!
   TBranch        *b_timeConstFrac_channel1;   //!
   TBranch        *b_slopeConstFrac_channel1;   //!
   TBranch        *b_chi2ConstFrac_channel1;   //!
   TBranch        *b_waveForm_channel2;   //!
   TBranch        *b_baseLine_channel2;   //!
   TBranch        *b_baselineDispersion_channel2;   //!
   TBranch        *b_baseLineGlobal_channel2;   //!
   TBranch        *b_baselineDispersionGlobal_channel2;   //!
   TBranch        *b_ampMax_channel2;   //!
   TBranch        *b_timeAmpMax_channel2;   //!
   TBranch        *b_constFrac_channel2;   //!
   TBranch        *b_timeConstFrac_channel2;   //!
   TBranch        *b_slopeConstFrac_channel2;   //!
   TBranch        *b_chi2ConstFrac_channel2;   //!
   TBranch        *b_waveForm_channel3;   //!
   TBranch        *b_baseLine_channel3;   //!
   TBranch        *b_baselineDispersion_channel3;   //!
   TBranch        *b_baseLineGlobal_channel3;   //!
   TBranch        *b_baselineDispersionGlobal_channel3;   //!
   TBranch        *b_ampMax_channel3;   //!
   TBranch        *b_timeAmpMax_channel3;   //!
   TBranch        *b_constFrac_channel3;   //!
   TBranch        *b_timeConstFrac_channel3;   //!
   TBranch        *b_slopeConstFrac_channel3;   //!
   TBranch        *b_chi2ConstFrac_channel3;   //!
   TBranch        *b_waveForm_channel4;   //!
   TBranch        *b_baseLine_channel4;   //!
   TBranch        *b_baselineDispersion_channel4;   //!
   TBranch        *b_baseLineGlobal_channel4;   //!
   TBranch        *b_baselineDispersionGlobal_channel4;   //!
   TBranch        *b_ampMax_channel4;   //!
   TBranch        *b_timeAmpMax_channel4;   //!
   TBranch        *b_constFrac_channel4;   //!
   TBranch        *b_timeConstFrac_channel4;   //!
   TBranch        *b_slopeConstFrac_channel4;   //!
   TBranch        *b_chi2ConstFrac_channel4;   //!
   TBranch        *b_waveForm_channel5;   //!
   TBranch        *b_baseLine_channel5;   //!
   TBranch        *b_baselineDispersion_channel5;   //!
   TBranch        *b_baseLineGlobal_channel5;   //!
   TBranch        *b_baselineDispersionGlobal_channel5;   //!
   TBranch        *b_ampMax_channel5;   //!
   TBranch        *b_timeAmpMax_channel5;   //!
   TBranch        *b_constFrac_channel5;   //!
   TBranch        *b_timeConstFrac_channel5;   //!
   TBranch        *b_slopeConstFrac_channel5;   //!
   TBranch        *b_chi2ConstFrac_channel5;   //!
   TBranch        *b_waveForm_channel6;   //!
   TBranch        *b_baseLine_channel6;   //!
   TBranch        *b_baselineDispersion_channel6;   //!
   TBranch        *b_baseLineGlobal_channel6;   //!
   TBranch        *b_baselineDispersionGlobal_channel6;   //!
   TBranch        *b_ampMax_channel6;   //!
   TBranch        *b_timeAmpMax_channel6;   //!
   TBranch        *b_constFrac_channel6;   //!
   TBranch        *b_timeConstFrac_channel6;   //!
   TBranch        *b_slopeConstFrac_channel6;   //!
   TBranch        *b_chi2ConstFrac_channel6;   //!
   TBranch        *b_waveForm_channel7;   //!
   TBranch        *b_baseLine_channel7;   //!
   TBranch        *b_baselineDispersion_channel7;   //!
   TBranch        *b_baseLineGlobal_channel7;   //!
   TBranch        *b_baselineDispersionGlobal_channel7;   //!
   TBranch        *b_ampMax_channel7;   //!
   TBranch        *b_timeAmpMax_channel7;   //!
   TBranch        *b_constFrac_channel7;   //!
   TBranch        *b_timeConstFrac_channel7;   //!
   TBranch        *b_slopeConstFrac_channel7;   //!
   TBranch        *b_chi2ConstFrac_channel7;   //!
   TBranch        *b_waveForm_channel8;   //!
   TBranch        *b_baseLine_channel8;   //!
   TBranch        *b_baselineDispersion_channel8;   //!
   TBranch        *b_baseLineGlobal_channel8;   //!
   TBranch        *b_baselineDispersionGlobal_channel8;   //!
   TBranch        *b_ampMax_channel8;   //!
   TBranch        *b_timeAmpMax_channel8;   //!
   TBranch        *b_constFrac_channel8;   //!
   TBranch        *b_timeConstFrac_channel8;   //!
   TBranch        *b_slopeConstFrac_channel8;   //!
   TBranch        *b_chi2ConstFrac_channel8;   //!

   void InitTree(TTree* nt){

     // Set object pointer
     waveForm_Trigger = 0;
     waveForm_channel1 = 0;
     waveForm_channel2 = 0;
     waveForm_channel3 = 0;
     waveForm_channel4 = 0;
     waveForm_channel5 = 0;
     waveForm_channel6 = 0;
     waveForm_channel7 = 0;
     waveForm_channel8 = 0;

     // Set branch addresses and branch pointers
     nt->SetBranchAddress("event", &event, &b_event);
     nt->SetBranchAddress("run", &run, &b_run);
     nt->SetBranchAddress("nPointsInterpolation", &nPointsInterpolation, &b_nPointsInterpolation);
     nt->SetBranchAddress("timeBaseLine", &timeBaseLine, &b_timeBaseLine);
     nt->SetBranchAddress("waveForm_Trigger", &waveForm_Trigger, &b_waveForm_Trigger);
     nt->SetBranchAddress("baseLine_Trigger", &baseLine_Trigger, &b_baseLine_Trigger);
     nt->SetBranchAddress("baselineDispersion_Trigger", &baselineDispersion_Trigger, &b_baselineDispersion_Trigger);
     nt->SetBranchAddress("baseLineGlobal_Trigger", &baseLineGlobal_Trigger, &b_baseLineGlobal_Trigger);
     nt->SetBranchAddress("baselineDispersionGlobal_Trigger", &baselineDispersionGlobal_Trigger, &b_baselineDispersionGlobal_Trigger);
     nt->SetBranchAddress("ampMax_Trigger", &ampMax_Trigger, &b_ampMax_Trigger);
     nt->SetBranchAddress("timeAmpMax_Trigger", &timeAmpMax_Trigger, &b_timeAmpMax_Trigger);
     nt->SetBranchAddress("constFrac_Trigger", &constFrac_Trigger, &b_constFrac_Trigger);
     nt->SetBranchAddress("timeConstFrac_Trigger", &timeConstFrac_Trigger, &b_timeConstFrac_Trigger);
     nt->SetBranchAddress("slopeConstFrac_Trigger", &slopeConstFrac_Trigger, &b_slopeConstFrac_Trigger);
     nt->SetBranchAddress("chi2ConstFrac_Trigger", &chi2ConstFrac_Trigger, &b_chi2ConstFrac_Trigger);
     nt->SetBranchAddress("waveForm_channel1", &waveForm_channel1, &b_waveForm_channel1);
     nt->SetBranchAddress("baseLine_channel1", &baseLine_channel1, &b_baseLine_channel1);
     nt->SetBranchAddress("baselineDispersion_channel1", &baselineDispersion_channel1, &b_baselineDispersion_channel1);
     nt->SetBranchAddress("baseLineGlobal_channel1", &baseLineGlobal_channel1, &b_baseLineGlobal_channel1);
     nt->SetBranchAddress("baselineDispersionGlobal_channel1", &baselineDispersionGlobal_channel1, &b_baselineDispersionGlobal_channel1);
     nt->SetBranchAddress("ampMax_channel1", &ampMax_channel1, &b_ampMax_channel1);
     nt->SetBranchAddress("timeAmpMax_channel1", &timeAmpMax_channel1, &b_timeAmpMax_channel1);
     nt->SetBranchAddress("constFrac_channel1", &constFrac_channel1, &b_constFrac_channel1);
     nt->SetBranchAddress("timeConstFrac_channel1", &timeConstFrac_channel1, &b_timeConstFrac_channel1);
     nt->SetBranchAddress("slopeConstFrac_channel1", &slopeConstFrac_channel1, &b_slopeConstFrac_channel1);
     nt->SetBranchAddress("chi2ConstFrac_channel1", &chi2ConstFrac_channel1, &b_chi2ConstFrac_channel1);
     nt->SetBranchAddress("waveForm_channel2", &waveForm_channel2, &b_waveForm_channel2);
     nt->SetBranchAddress("baseLine_channel2", &baseLine_channel2, &b_baseLine_channel2);
     nt->SetBranchAddress("baselineDispersion_channel2", &baselineDispersion_channel2, &b_baselineDispersion_channel2);
     nt->SetBranchAddress("baseLineGlobal_channel2", &baseLineGlobal_channel2, &b_baseLineGlobal_channel2);
     nt->SetBranchAddress("baselineDispersionGlobal_channel2", &baselineDispersionGlobal_channel2, &b_baselineDispersionGlobal_channel2);
     nt->SetBranchAddress("ampMax_channel2", &ampMax_channel2, &b_ampMax_channel2);
     nt->SetBranchAddress("timeAmpMax_channel2", &timeAmpMax_channel2, &b_timeAmpMax_channel2);
     nt->SetBranchAddress("constFrac_channel2", &constFrac_channel2, &b_constFrac_channel2);
     nt->SetBranchAddress("timeConstFrac_channel2", &timeConstFrac_channel2, &b_timeConstFrac_channel2);
     nt->SetBranchAddress("slopeConstFrac_channel2", &slopeConstFrac_channel2, &b_slopeConstFrac_channel2);
     nt->SetBranchAddress("chi2ConstFrac_channel2", &chi2ConstFrac_channel2, &b_chi2ConstFrac_channel2);
     nt->SetBranchAddress("waveForm_channel3", &waveForm_channel3, &b_waveForm_channel3);
     nt->SetBranchAddress("baseLine_channel3", &baseLine_channel3, &b_baseLine_channel3);
     nt->SetBranchAddress("baselineDispersion_channel3", &baselineDispersion_channel3, &b_baselineDispersion_channel3);
     nt->SetBranchAddress("baseLineGlobal_channel3", &baseLineGlobal_channel3, &b_baseLineGlobal_channel3);
     nt->SetBranchAddress("baselineDispersionGlobal_channel3", &baselineDispersionGlobal_channel3, &b_baselineDispersionGlobal_channel3);
     nt->SetBranchAddress("ampMax_channel3", &ampMax_channel3, &b_ampMax_channel3);
     nt->SetBranchAddress("timeAmpMax_channel3", &timeAmpMax_channel3, &b_timeAmpMax_channel3);
     nt->SetBranchAddress("constFrac_channel3", &constFrac_channel3, &b_constFrac_channel3);
     nt->SetBranchAddress("timeConstFrac_channel3", &timeConstFrac_channel3, &b_timeConstFrac_channel3);
     nt->SetBranchAddress("slopeConstFrac_channel3", &slopeConstFrac_channel3, &b_slopeConstFrac_channel3);
     nt->SetBranchAddress("chi2ConstFrac_channel3", &chi2ConstFrac_channel3, &b_chi2ConstFrac_channel3);
     nt->SetBranchAddress("waveForm_channel4", &waveForm_channel4, &b_waveForm_channel4);
     nt->SetBranchAddress("baseLine_channel4", &baseLine_channel4, &b_baseLine_channel4);
     nt->SetBranchAddress("baselineDispersion_channel4", &baselineDispersion_channel4, &b_baselineDispersion_channel4);
     nt->SetBranchAddress("baseLineGlobal_channel4", &baseLineGlobal_channel4, &b_baseLineGlobal_channel4);
     nt->SetBranchAddress("baselineDispersionGlobal_channel4", &baselineDispersionGlobal_channel4, &b_baselineDispersionGlobal_channel4);
     nt->SetBranchAddress("ampMax_channel4", &ampMax_channel4, &b_ampMax_channel4);
     nt->SetBranchAddress("timeAmpMax_channel4", &timeAmpMax_channel4, &b_timeAmpMax_channel4);
     nt->SetBranchAddress("constFrac_channel4", &constFrac_channel4, &b_constFrac_channel4);
     nt->SetBranchAddress("timeConstFrac_channel4", &timeConstFrac_channel4, &b_timeConstFrac_channel4);
     nt->SetBranchAddress("slopeConstFrac_channel4", &slopeConstFrac_channel4, &b_slopeConstFrac_channel4);
     nt->SetBranchAddress("chi2ConstFrac_channel4", &chi2ConstFrac_channel4, &b_chi2ConstFrac_channel4);
     nt->SetBranchAddress("waveForm_channel5", &waveForm_channel5, &b_waveForm_channel5);
     nt->SetBranchAddress("baseLine_channel5", &baseLine_channel5, &b_baseLine_channel5);
     nt->SetBranchAddress("baselineDispersion_channel5", &baselineDispersion_channel5, &b_baselineDispersion_channel5);
     nt->SetBranchAddress("baseLineGlobal_channel5", &baseLineGlobal_channel5, &b_baseLineGlobal_channel5);
     nt->SetBranchAddress("baselineDispersionGlobal_channel5", &baselineDispersionGlobal_channel5, &b_baselineDispersionGlobal_channel5);
     nt->SetBranchAddress("ampMax_channel5", &ampMax_channel5, &b_ampMax_channel5);
     nt->SetBranchAddress("timeAmpMax_channel5", &timeAmpMax_channel5, &b_timeAmpMax_channel5);
     nt->SetBranchAddress("constFrac_channel5", &constFrac_channel5, &b_constFrac_channel5);
     nt->SetBranchAddress("timeConstFrac_channel5", &timeConstFrac_channel5, &b_timeConstFrac_channel5);
     nt->SetBranchAddress("slopeConstFrac_channel5", &slopeConstFrac_channel5, &b_slopeConstFrac_channel5);
     nt->SetBranchAddress("chi2ConstFrac_channel5", &chi2ConstFrac_channel5, &b_chi2ConstFrac_channel5);
     nt->SetBranchAddress("waveForm_channel6", &waveForm_channel6, &b_waveForm_channel6);
     nt->SetBranchAddress("baseLine_channel6", &baseLine_channel6, &b_baseLine_channel6);
     nt->SetBranchAddress("baselineDispersion_channel6", &baselineDispersion_channel6, &b_baselineDispersion_channel6);
     nt->SetBranchAddress("baseLineGlobal_channel6", &baseLineGlobal_channel6, &b_baseLineGlobal_channel6);
     nt->SetBranchAddress("baselineDispersionGlobal_channel6", &baselineDispersionGlobal_channel6, &b_baselineDispersionGlobal_channel6);
     nt->SetBranchAddress("ampMax_channel6", &ampMax_channel6, &b_ampMax_channel6);
     nt->SetBranchAddress("timeAmpMax_channel6", &timeAmpMax_channel6, &b_timeAmpMax_channel6);
     nt->SetBranchAddress("constFrac_channel6", &constFrac_channel6, &b_constFrac_channel6);
     nt->SetBranchAddress("timeConstFrac_channel6", &timeConstFrac_channel6, &b_timeConstFrac_channel6);
     nt->SetBranchAddress("slopeConstFrac_channel6", &slopeConstFrac_channel6, &b_slopeConstFrac_channel6);
     nt->SetBranchAddress("chi2ConstFrac_channel6", &chi2ConstFrac_channel6, &b_chi2ConstFrac_channel6);
     nt->SetBranchAddress("waveForm_channel7", &waveForm_channel7, &b_waveForm_channel7);
     nt->SetBranchAddress("baseLine_channel7", &baseLine_channel7, &b_baseLine_channel7);
     nt->SetBranchAddress("baselineDispersion_channel7", &baselineDispersion_channel7, &b_baselineDispersion_channel7);
     nt->SetBranchAddress("baseLineGlobal_channel7", &baseLineGlobal_channel7, &b_baseLineGlobal_channel7);
     nt->SetBranchAddress("baselineDispersionGlobal_channel7", &baselineDispersionGlobal_channel7, &b_baselineDispersionGlobal_channel7);
     nt->SetBranchAddress("ampMax_channel7", &ampMax_channel7, &b_ampMax_channel7);
     nt->SetBranchAddress("timeAmpMax_channel7", &timeAmpMax_channel7, &b_timeAmpMax_channel7);
     nt->SetBranchAddress("constFrac_channel7", &constFrac_channel7, &b_constFrac_channel7);
     nt->SetBranchAddress("timeConstFrac_channel7", &timeConstFrac_channel7, &b_timeConstFrac_channel7);
     nt->SetBranchAddress("slopeConstFrac_channel7", &slopeConstFrac_channel7, &b_slopeConstFrac_channel7);
     nt->SetBranchAddress("chi2ConstFrac_channel7", &chi2ConstFrac_channel7, &b_chi2ConstFrac_channel7);
     nt->SetBranchAddress("waveForm_channel8", &waveForm_channel8, &b_waveForm_channel8);
     nt->SetBranchAddress("baseLine_channel8", &baseLine_channel8, &b_baseLine_channel8);
     nt->SetBranchAddress("baselineDispersion_channel8", &baselineDispersion_channel8, &b_baselineDispersion_channel8);
     nt->SetBranchAddress("baseLineGlobal_channel8", &baseLineGlobal_channel8, &b_baseLineGlobal_channel8);
     nt->SetBranchAddress("baselineDispersionGlobal_channel8", &baselineDispersionGlobal_channel8, &b_baselineDispersionGlobal_channel8);
     nt->SetBranchAddress("ampMax_channel8", &ampMax_channel8, &b_ampMax_channel8);
     nt->SetBranchAddress("timeAmpMax_channel8", &timeAmpMax_channel8, &b_timeAmpMax_channel8);
     nt->SetBranchAddress("constFrac_channel8", &constFrac_channel8, &b_constFrac_channel8);
     nt->SetBranchAddress("timeConstFrac_channel8", &timeConstFrac_channel8, &b_timeConstFrac_channel8);
     nt->SetBranchAddress("slopeConstFrac_channel8", &slopeConstFrac_channel8, &b_slopeConstFrac_channel8);
     nt->SetBranchAddress("chi2ConstFrac_channel8", &chi2ConstFrac_channel8, &b_chi2ConstFrac_channel8);
   }

#endif
