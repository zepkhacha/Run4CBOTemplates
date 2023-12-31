#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include "gm2util/blinders/Blinders.hh"

static blinding::Blinders myblinders(blinding::Blinders::kOmega_a, "reunificationisblind");

//histogram binning
int          startBin = 203;
double       startTime = 0.0;
const int    nBins  = 4357 ;
double 	     n	    = 0.108;

// histograms
static TH1D* wiggle = new TH1D("wiggle", "wiggle", nBins, startTime, 650.0644);
static TH1D* lambda = new TH1D("lambda", "lambda", nBins, startTime, 650.0644);
static TH1D* varcorr = new TH1D("varcorr", "varcorr", nBins, startTime, 650.0644);
static TH1D* wiggleraw = new TH1D("wiggleraw", "wiggleraw", nBins, startTime, 650.0644);
static TH1D* cboIsolate = new TH1D("cboIsolate", "cboIsolate", nBins, startTime, 650.0644);

// kevin-style switches
//static double expcorr[nBins];
bool   fixtau = false;
bool   constraintau = true;
bool   correcttau   = false;
bool   asymmetry    = false;
int    fitrangelow   = startBin;
int    fitrangehigh  = nBins;
int windowBins;
double frlifetime = 0;
// zep-style switches
int desiredSeed = 0;
int desiredCalo = 0;
int windowNo    = 0;
bool includeMopTerm = false;
bool constrainMop = false;
bool pBinTau = false  ; //true ; 
bool pBinPhi = false  ; //true ; 
// bin-dependent variables
std::pair <double, double> magicTauConstraints (0.0,0.0);
std::pair <int,int> binRange (0,0); 
int magicBin = -10;
std::vector<double> binWeights;
std::vector<double> dp_p0;
std::vector<double> gammas;
std::vector<double> dPhi;
double weightedAvgTau;
// e-field 
double c_e_scale = 1.0; // set from 0 - 1
std::vector<double> c_e; //in ppm
std::vector<double> x_e;

double invert(double c, double s){
  //double i = atan(s/c);
  //if(c < 0){
  //  i += M_PI;
  //  return i;
  //}
  //if(s < 0){
  //  i += 2*M_PI;
  //  return i;
  //}
  //return i;
  return atan2(s,c);
}

double cbo(double time, double amp){
// For the endgame
//return 1.0 + 2.927*exp(-time/79.05)/time + 2.239*exp(-time/6.94)/time;
// For Run 2/3
return 1.0 + amp*exp(-time/24.4)/time;
}

double calcnuCboIsolate(double par[], double i){
  double time = (i-0.5) * 0.1492 + startTime;
  double nu = (1.0 + par[0]*cos((par[1])*time - par[2]));
  return nu; 
}

void fitCboIsolate( __attribute__((unused)) int& nDim,
                     __attribute__((unused)) double* gout,
                     double& result, double par[],
                     __attribute__((unused)) int flg){
  result = 0;
  for(int i=fitrangelow; i<fitrangehigh; i++){
    double n = cboIsolate->GetBinContent(i);
    double nu = calcnuCboIsolate(par, i);
    //double chisq = pow(n-nu, 2)/n;
    result += pow(n-nu, 2)/n;
    //printf("fitCboIsolate: i %i time %f n %f nu %f varcorr %f wiggleraw %f chisq %f \n", 
    //     i, (i-0.5)*0.1492, n, nu, varcorr->GetBinContent(i), wiggleraw->GetBinContent(i), result);
  }
}

double calcnu(double par[], double i){

  double time = (i-0.5) * 0.1492 + startTime;
  double blindr = myblinders.paramToFreq(par[4]);
  double phi = par[3];

  double nu = 0.0;

  nu = par[0]*exp(-1.*time/par[1]) * (1 + par[2]*cos((blindr)*time - phi));
  // now add x-terms
  nu *= 1.0 + par[5]*cos(par[6]*(time) - par[7]) + par[8]*cos(2*par[7]*(time)- par[9]);
  // now add y-terms
  nu *= (1.0 + par[10] * cos(par[12]*time - par[11])
             + par[13] * cos(par[15]*time - par[14]));
  // now add LM 
  nu *= (1.0 - par[16]*lambda->GetBinContent(i));
  return nu; 

}

