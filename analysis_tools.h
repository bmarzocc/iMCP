/*************************************************************

*************************************************************/

#ifndef __analysis_tools__
#define __analysis_tools__

#include <cstdio>
#include <cstdlib>
#include <stdint.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TTree.h"

#include "histoFunc.h"

using namespace std;

#define DIGITIZER_SAMPLING_UNIT 0.2 //digitizer samples width (ns)
#define MS_SAMPLING_UNIT 0.1 //mean signal samples width (ns)
#define MS_LOW_TIME -15
#define MS_HIGH_TIME 30

//---Integral treshold th[iScan][Ch_n] Roma1 is on Ch=2
float _th[4][6];

//------------------------------------------------------------------------------
void init()
{
    _th[0][0] = -110;
    _th[1][0] = -110;
    _th[2][0] = -110;
    _th[3][0] = -110;

    _th[0][1] = -130;
    _th[1][1] = -130;
    _th[2][1] = -130;
    _th[3][1] = -130;

    _th[0][2] = -0;
    _th[1][2] = -0;
    _th[2][2] = -0;
    _th[3][2] = -0;
    
    _th[0][3] = -55;
    _th[1][3] = -55;
    _th[2][3] = -55;
    _th[3][3] = -55;

    _th[0][4] = -140;
    _th[1][4] = -140;
    _th[2][4] = -140;
    _th[3][4] = -140;

    _th[0][5] = -70;
    _th[1][5] = -70;
    _th[2][5] = -70;
    _th[3][5] = -45;
}

//------------------------------------------------------------------------------
void DFT_lowCut(vector<float>* samples, float f_cut)
{
    float* a = (float*)malloc(sizeof(float)*samples->size());
    float* b = (float*)malloc(sizeof(float)*samples->size());
    for (int k=0; k<samples->size(); k++) 
    { 
        a[k] = 0;
        b[k] = 0; 
        for (int t=0; t<samples->size(); t++) 
        { 
            float angle = 2 * M_PI * t * k / samples->size(); 
            a[k] += samples->at(t) * cos(angle); 
            b[k] += samples->at(t) * sin(angle); 
        } 
    }
}

//---estimate the baseline in a given range and then subtract it from the signal 
float SubtractBaseline(int tb1, int tb2, vector<float>* samples)
{
    float baseline=0;
    int minSample=tb1;
    //---compute baseline
    for(int iSample=tb1; iSample<tb2; iSample++)
    {
        baseline += samples->at(iSample);
        if(samples->at(iSample) < samples->at(minSample)) 
            minSample = iSample;
    }
    baseline = baseline/(float)(tb2-tb1);
    //---subtract baseline
    for(int iSample=0; iSample<samples->size(); iSample++)
    {
        samples->at(iSample) = samples->at(iSample) - baseline;
    }
    return samples->at(minSample);
}

//---estimate time (ns) with CFD, samples must be a negative signal and baseline subtract
float TimeConstFrac(int t1, int t2, const vector<float>* samples, float AmpFraction, 
                    float step=0.2, int Nsamples = 5)
{
    float xx= 0.;
    float xy= 0.;
    float Sx = 0.;
    float Sy = 0.;
    float Sxx = 0.;
    float Sxy = 0.;
    float Chi2 = 0.;
    int minSample=t1;
    int cfSample=t1; // first sample over AmpMax*CF 
    float minValue=0;

    for(int iSample=t1; iSample<t2; iSample++)
    {
        if(samples->at(iSample) < samples->at(minSample)) minSample = iSample;
    }
    minValue = samples->at(minSample);
    if(AmpFraction == 1) 
        return minSample*step;
    for(int iSample=minSample; iSample>t1; iSample--)
    {
        if(samples->at(iSample) > minValue*AmpFraction) 
        {
            cfSample = iSample;
            break;
        }
    }
    for(int n=-(Nsamples-1)/2; n<=(Nsamples-1)/2; n++)
    {
        if(cfSample+n<0) continue;
        xx = (cfSample+n)*(cfSample+n)*step*step;
        xy = (cfSample+n)*step*(samples->at(cfSample+n));
        Sx = Sx + (cfSample+n)*step;
        Sy = Sy + samples->at(cfSample+n);
        Sxx = Sxx + xx;
        Sxy = Sxy + xy;
    }

    float Delta = Nsamples*Sxx - Sx*Sx;
    float A = (Sxx*Sy - Sx*Sxy) / Delta;
    float B = (Nsamples*Sxy - Sx*Sy) / Delta;

    float sigma2 = pow(step/sqrt(12)*B,2);
 
    for(int n=-(Nsamples-1)/2; n<=(Nsamples-1)/2; n++)
    {
        if(cfSample+n<0) continue;
        Chi2 = Chi2 + pow(samples->at(cfSample+n) - A - B*((cfSample+n)*step),2)/sigma2;
    } 
    // A+Bx = AmpFraction * amp
    float interpolation = (samples->at(minSample) * AmpFraction - A) / B;
    return interpolation;
}

float AmpMax(int t1, int t2, const vector<float>* samples)
{
    int minSample=t1;
    for(int iSample=t1; iSample<t2; iSample++)
    {
        if(samples->at(iSample) < samples->at(minSample)) minSample = iSample;
    }
    return samples->at(minSample);
}

float ComputeIntegral(int t1, int t2, const vector<float>* samples)
{
    float integral=0;
    for(int bin=t1; bin<t2; bin++)
        integral += samples->at(bin);

    return integral;
}

#endif
