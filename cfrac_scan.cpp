/*************************************************************

        simple program to read bin files
        compile with ---> c++ -o cfrac_scan cfrac_scan.cpp `root-config --cflags --glibs`

*************************************************************/

#include <cstdio>
#include <cstdlib>
#include <stdint.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TString.h"
#include "TApplication.h"
#include "TGraphErrors.h"

using namespace std;

#define DIGITIZER_SAMPLING_UNIT 0.2 //digitizer samples width (ns)

//---estimate time (ns) with CFD, samples must be a negative signal and baseline subtract
float TimeConstFrac(int x1, int x2, const vector<float>* samples, float AmpFraction, int Nsamples = 5)
{

    float xx= 0.;
    float xy= 0.;
    float Sx = 0.;
    float Sy = 0.;
    float Sxx = 0.;
    float Sxy = 0.;
    float Chi2 = 0.;
    int minSample=0;
    int cfSample=0; // first sample over AmpMax*CF 
    float minValue=0;
    
    for(int iSample=x1; iSample<x2; iSample++)
    {
        if(samples->at(iSample) < samples->at(maxSample)) minSample = iSample;
    }
    minValue = samples->at(minSample);
    for(int iSample=x1; iSample<x2; iSample++)
    {
        if(samples->at(iSample) > minValue*AmpFraction) cfSample = iSample;
    }
    for(int n=-(Nsample-1)/2; n<=(Nsample-1)/2; n++)
    {
        if(cfSample+n<0) continue;
        xx = (cfSample+n)*(cfSample+n)*0.2*0.2;
        xy = (cfSample+n)*0.2*(samples->at(cfSample+n));
        Sx = Sx + (cfSample+n)*0.2;
        Sy = Sy + samples->at(cfSample+n);
        Sxx = Sxx + xx;
        Sxy = Sxy + xy;
    }

    float Delta = Nsample*Sxx - Sx*Sx;
    float A = (Sxx*Sy - Sx*Sxy) / Delta;
    float B = (Nsample*Sxy - Sx*Sy) / Delta;

    float sigma2 = pow(0.2/sqrt(12)*B,2);
 
    for(int n=-(Nsample-1)/2; n<=(Nsample-1)/2; n++)
    {
        if(cfSample+n<0) continue;
        Chi2 = Chi2 + pow(samples->at(cfSample+n) - A - B*((cfSample+n)*0.2),2)/sigma2;
    } 
    // A+Bx = AmpFraction * amp
    float interpolation = (samples->at(xm) * AmpFraction - A) / B;
    return interpolation;
}


