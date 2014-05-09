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

double Trigger(int x1, int xm, vector<float> samples, float AmpFraction, int Npar = 3, float sampling = 0.2)
{
    double xx= 0.;
    double xy= 0.;
    double Sx = 0.;
    double Sy = 0.;
    double Sxx = 0.;
    double Sxy = 0.;

    for(int n=-(Npar-1)/2; n<=(Npar-1)/2; n++)
    {
        if(x1+n<0) break;
        xx = (x1+n)*(x1+n)*0.2*0.2;
        xy = (x1+n)*0.2*(samples.at(x1+n));
        Sx = Sx + (x1+n)*0.2;
        Sy = Sy + samples.at(x1+n);
        Sxx = Sxx + xx;
        Sxy = Sxy + xy;
    }

    double Delta = Npar*Sxx - Sx*Sx;
    double A = (Sxx*Sy - Sx*Sxy) / Delta;
    double B = (Npar*Sxy - Sx*Sy) / Delta;

    // A+Bx = AmpFraction * amp
    double interpolation = (samples.at(xm) * AmpFraction - A) / B;
    return B;
}

//******************************************************************************
// main

int main(int argc, char** argv)
{
    FILE* input;
    vector<float> channels[9];
    vector<int> eventNumber;
    vector<float> trig[9], B, sigmaB;
    float *buffer, frac=0;
    int *BinHeader, isNegative[9], nCh=0, nSize=0, iEvent=0, nEvent=0, c=0;
    int x1=0;
    
    TH1F* double_coinc = new TH1F("difference","difference",2000,-50,50);
    TGraphErrors* scan = new TGraphErrors();
    scan->SetMarkerStyle(7);
    scan->SetMarkerSize(20);
    scan->SetMarkerColor(kBlue);
    scan->SetTitle("scan; costant fraction [% ampMax]; Slope [ADC ch/ns]");
    TApplication* app = new TApplication("app",0,0);
    
    //-----scan loop
    for(c=20; c<80; c++)
    {
        system("ls WaveForms/ > tmp");
        ifstream infile ("tmp");
        string file;
        string path = "WaveForms/";
        
        TH1F* B_histo = new TH1F("B_histo"+c,"B_histo",2000,-2000,10);
        //-----loop on all the Runs
        while(getline(infile, file))
        {
            path = "WaveForms/";
            path.append(file);
            if((input = fopen(path.c_str(),"rb")) == NULL)
            {
                cout << "input file not found" << endl;
                return -1;
            }   
    
            frac = c/100.;
    
            BinHeader = (int*) malloc (sizeof(int)*10);
            fread(BinHeader, sizeof(int), 10, input);
            nSize = BinHeader[0];
            nCh = BinHeader[1];
    
            isNegative[0] = -1;
            for(int ii = 1; ii < 9; ii++)
                isNegative[ii] = BinHeader[ii+1];

            fread(&iEvent, sizeof(int),1,input);
    
            buffer = (float*) malloc (sizeof(float)*nSize*(nCh+1));
    
	        while(fread(buffer, sizeof(float), nSize*(nCh+1), input) == nSize*(nCh+1))
            {
    	       eventNumber.push_back(iEvent);
               for(int iCh=0; iCh<=nCh; iCh++)
               {
                   channels[iCh].clear();
                   for(int iSample=0; iSample<nSize; iSample++)
                   {
                       int pos = iCh*nSize + iSample;
                       channels[iCh].push_back(buffer[pos]);
              	   }
                }
    
                //-----channels loop (only Russian MCP)
                for(int iCh=1; iCh<nCh; iCh++)
                {
                    //---compute baseline (5ns)
                    float baseline = 0.;
            	    for(int iBase=0; iBase<25; iBase++)
                    	baseline = baseline + channels[iCh].at(iBase);
        
            	    baseline = baseline/25;

            	    for(int iSample=0; iSample<nSize; iSample++)
            	    	channels[iCh].at(iSample) = channels[iCh].at(iSample) - baseline;
	        	
	        	    int maxSample = 0;
	        	    if (isNegative[iCh]<0) 
	        	    {
                        for(int iSample=0; iSample<nSize; iSample++)
                        {
                            if(channels[iCh].at(iSample)<channels[iCh].at(maxSample)) maxSample=iSample;
                        }    
	        	        for (int iSample=maxSample; iSample>0; iSample--)
	        	        {    
	        	            //---compute costant fraction time 
	        	            if(channels[iCh].at(iSample) > frac*channels[iCh].at(maxSample) &&
	        	               channels[iCh].at(iSample-1) > channels[iCh].at(iSample) &&
	        	               channels[iCh].at(iSample) > channels[iCh].at(iSample+1))
	        	            {
	        	                x1 = iSample;
	        	                break;
	        	            }
	        	        }
	        	        if(channels[iCh].at(maxSample) < -200)
                            B_histo->Fill(Trigger(x1, maxSample, channels[iCh], frac, 5));
      	            }
                }
                nEvent++;
                fread(&iEvent, sizeof(int), 1, input);
            }
            fclose(input);
        }            
        B_histo->SetAxisRange(-500,0,"X");
        B.push_back(B_histo->GetMean());
        sigmaB.push_back(B_histo->GetRMS()/sqrt(B_histo->GetEntries()));
        cout << B_histo->GetEntries() << endl; 
        delete B_histo;
    }

    for(int j=0; j<B.size(); j++)
    {
        scan->SetPoint(j,j+20,B.at(j)); 
        scan->SetPointError(j,0,sigmaB.at(j));
    }
    TCanvas* c1 = new TCanvas();
    scan->Draw("AP");    
    app->Run();    
}


