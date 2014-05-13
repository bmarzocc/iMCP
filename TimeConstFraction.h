#ifndef __TIMECONSTFRACTION__
#define __TIMECONSTFRACTION__

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

class TimeConstFraction
{

public:

    TimeConstFraction(const int& x1, const int& xm, std::vector<float> samples, const float& AmpFraction, const float& baseline, const int& Npoint = 3, const float&  sampling = 0.2);
    virtual ~TimeConstFraction();
    float getTime();
    float getSlope();
    float getIntercept();
    float getChi2();
    

protected:
    
    float xx_;
    float xy_;
    float Sx_;
    float Sy_;
    float Sxx_;
    float Sxy_;
    float Delta_;
    float q_;
    float m_;
    float time_;
    float sigma2_;
    float Chi2_;
    
};

//---------------------------------------------------------------------------------------------------------------------------

//ctor
TimeConstFraction::TimeConstFraction(const int& x1, const int& xm, std::vector<float> samples, const float& AmpFraction, const float& baseline, const int& Npoint, const float&  sampling)
{
  xx_ = 0.;
  xy_ = 0.;
  Sx_ = 0.;
  Sy_ = 0.;
  Sxx_ = 0.;
  Sxy_ = 0.;
  time_ = 0.;
  sigma2_ = 0.;
  Chi2_ = 0.;
  
  for(int n=-(Npoint-1)/2; n<=(Npoint-1)/2; n++)
  {

        if(x1+n<0) continue;
        if(x1+n>=samples.size()) continue;
        xx_ = (x1+n)*(x1+n)*sampling*sampling;
        xy_ = (x1+n)*sampling*(samples.at(x1+n)-baseline);
        Sx_ = Sx_ + (x1+n)*sampling;
        Sy_ = Sy_ + (samples.at(x1+n)-baseline);
        Sxx_ = Sxx_ + xx_;
        Sxy_ = Sxy_ + xy_;

  }

  Delta_ = Npoint*Sxx_ - Sx_*Sx_;
  q_ = (Sxx_*Sy_ - Sx_*Sxy_) / Delta_;
  m_ = (Npoint*Sxy_ - Sx_*Sy_) / Delta_;
  time_ = ((samples.at(xm)-baseline) * AmpFraction - q_) / m_;

  sigma2_ = pow(sampling/sqrt(12)*m_,2);
  for(int n=-(Npoint-1)/2; n<=(Npoint-1)/2; n++)
  {
      if(x1+n<0) continue;
      if(x1+n>=samples.size()) continue;
      Chi2_ = Chi2_ + pow((samples.at(x1+n)-baseline) - q_ - m_*((x1+n)*sampling),2)/sigma2_;
  } 
}

//---------------------------------------------------------------------------------------------------------------------------

//dtor
TimeConstFraction::~TimeConstFraction()
{

}

//---------------------------------------------------------------------------------------------------------------------------

float TimeConstFraction::getTime()
{
   float time = time_;
   return time;
}

//---------------------------------------------------------------------------------------------------------------------------

float TimeConstFraction::getSlope()
{
   float slope = m_;
   return slope;
}

//---------------------------------------------------------------------------------------------------------------------------

float TimeConstFraction::getIntercept()
{
   float intercept = q_;
   return intercept;
}

//---------------------------------------------------------------------------------------------------------------------------

float TimeConstFraction::getChi2()
{
  float Chi2 = Chi2_;
  return Chi2;
}

#endif
