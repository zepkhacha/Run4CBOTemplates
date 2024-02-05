#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
//#include "TPaveText.h"
#include "TPaveStats.h"
#include "TGraphErrors.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include "gm2util/blinders/Blinders.hh"
#include "slidingFitFunc.hh"

void printhelp(){
    printf(
            "Analysis Program\n"
            "Supply arguments as follows:\n"
            "-h        : Prints this help text\n"
            "-i [s]    : Input filename (required)\n"
            "-o [s]    : Output filename (required)\n"
            "-b [s]    : Boost filename (required)\n"
            "-f [s]    : Format filename (required)\n"
            "-q [s]    : CBO reference filename (required)\n"
            "-a [i]    : Bool for asymm (1) or threshold (0)\n"
            "-c [d]    : Double for scale 0-1 to apply to c_e (default 1.0)\n"
            "-s [i]    : Int for seed (default 2)\n"
            "-w [i]    : Window number (required)\n"
            "-n [i]    : Int for calo (default 0)\n"
          );
}

void parse_cmdline(int argc, char** argv, 
        char* &formatFilename, 
        char* &inputFilename, char* &outputFilename, char* &boostFilename, 
        char* &fullFitFilename, 
        bool* asymmetry, double* c_e_scale, int* seed, int* windowNo,
        int* desiredCalo){
    const char* const opts = "hf:o:i:b:a:c:s:q:w:n:";
    bool done = false;
    while(!done){
        const char ch = getopt(argc, argv, opts);
        switch(ch){
            case -1: done = true; break;
            case 'h': printhelp(); exit(0);
            case 'f': formatFilename = optarg; break;
            case 'i': inputFilename = optarg; break;
            case 'o': outputFilename = optarg; break;
            case 'q': fullFitFilename = optarg; break;
            case 'b': boostFilename = optarg; break;
            case 'a': *asymmetry = bool(atoi(optarg)); break;
            case 'c': *c_e_scale = std::stod(optarg); break;
            case 's': *seed   = atoi(optarg); break;
            case 'w': *windowNo = atoi(optarg); break;
            case 'n': *desiredCalo = atoi(optarg); break;
            default: printhelp(); exit(1);
        }
    }
}

void minuitFunction( __attribute__((unused)) int& nDim,
        __attribute__((unused)) double* gout,
        double& result, double par[],
        __attribute__((unused)) int flg){
    result = 0;
    for(int i=fitrangelow; i<fitrangehigh; i++){
        double n = wiggle->GetBinContent(i);
        double nu = calcnu(par, i);
        if(asymmetry){
            result += pow(n-nu, 2)/varcorr->GetBinContent(i);
        }
        else{
            result += pow(n-nu, 2)/(wiggleraw->GetBinContent(i) + varcorr->GetBinContent(i));
        }
    }
}

int main(int argc, char* argv[]){
    printf("In main()...\n");
    // get input arguments
    char* formatFilename;
    char* outputFilename;
    char* inputFilename ;
    char* boostFilename ;
    char* fullFitFilename;

    parse_cmdline(argc,argv,
            formatFilename, 
            inputFilename, outputFilename, boostFilename, fullFitFilename,
            &asymmetry, &c_e_scale, &desiredSeed, &windowNo, &desiredCalo);

    // read in parameters for fit type
    printf("reading format from %s\n", formatFilename);
    std::fstream format;
    format.open(formatFilename, std::ios::in);
    if (format.is_open()) {
        std::string line;
        while (getline(format, line)){
            printf("%s\n",line.c_str());
            int datIndex=line.find("=")+1;
            if 	(line.find("pBinTau") != std::string::npos) {
                pBinTau = bool(atoi(line.substr(datIndex,100).c_str()));
            }else if 	(line.find("pBinPhi") != std::string::npos) {
                pBinPhi = bool(atoi(line.substr(datIndex,100).c_str()));
            }else if 	(line.find("constraintau") != std::string::npos) {
                constraintau = bool(atoi(line.substr(datIndex,100).c_str()));
            }else if 	(line.find("includeMopTerm") != std::string::npos) {
                includeMopTerm = bool(atoi(line.substr(datIndex,100).c_str()));
            }else if 	(line.find("constrainMop") != std::string::npos) {
                constrainMop = bool(atoi(line.substr(datIndex,100).c_str()));
            }else{}
        }
    } 
    format.close();

    printf("Parameters after format file: pBinTau %d pBinPhi %d constraintau %d includeMopTerm %d constrainMop %d c_e_scale %f \n", pBinTau, pBinPhi, constraintau, includeMopTerm, constrainMop, c_e_scale);
    printf("Opening file %s.\n", boostFilename);
    TFile* boostdata  = new TFile(boostFilename, "READ");
    TH1D*  boostgamma = (TH1D*) boostdata->Get("transform_gamma");
    frlifetime = boostgamma->GetMean();
    /*
       double* boost = boostgamma->GetX();
       double* freq = boostgamma->GetY();
       double sum = 0;
       double den = 0;
       for(int i=0; i<boostgamma->GetN(); i++){
       sum += boost[i]*freq[i];
       den += freq[i];
       }
       frlifetime = sum/den;
       */
    for(int i=0; i<nBins; i++){
        double tt = startTime + (i+0.5)*0.1492;
        double mdl = exp(-tt/(2.1969811*frlifetime));
        double xact = 0;
        for(int j=0; j<boostgamma->GetNbinsX(); j++){
            xact += boostgamma->GetBinContent(j+1)/boostgamma->Integral()*exp(-tt/(2.1969811*boostgamma->GetBinCenter(j+1)));
        }
        //expcorr[i] = xact/mdl;
    }
    boostdata->Close();

    TFile* output = new TFile(outputFilename, "RECREATE");
    TH1D* bestfit = new TH1D("bestfit", "bestfit", nBins, startTime, 650.0644);
    TTree* fitresults = new TTree("fitresults", "fitresults");

    int sidx;
    double chisq, rchisq;

    double N0, tau, A0, phi, R;
    double N0err, tauerr, A0err, phierr, Rerr;

    double alpha_CBO, beta_CBO, w_CBO, A_CBO, phi_CBO;
    double alpha_CBOerr, beta_CBOerr, w_CBOerr, A_CBOerr, phi_CBOerr;

    double alpha_2CBO, beta_2CBO, A_2CBO, phi_2CBO;
    double alpha_2CBOerr, beta_2CBOerr, A_2CBOerr, phi_2CBOerr;

    double alpha_y, beta_y, w_y, A_y, phi_y;
    double alpha_yerr, beta_yerr, w_yerr, A_yerr, phi_yerr;

    double alpha_vw, beta_vw, w_vw, A_vw, phi_vw;
    double alpha_vwerr, beta_vwerr, w_vwerr, A_vwerr, phi_vwerr;

    double LM, LMerr;

    double fitStart = fitrangelow  * 0.1492;
    double fitStop  = fitrangehigh * 0.1492;


    // read in full-fit results to fix 5-param + LM
    printf("reading full fit data from %s\n", fullFitFilename);
    TFile* cboFullFit = new TFile(fullFitFilename, "READ");

    double fullFit_Ax1, fullFit_wx1, fullFit_phix1;
    TTree* cboFullFitResults = (TTree*) cboFullFit->Get("fitresults");

    cboFullFitResults->SetBranchAddress("N0"   , &N0     );
    cboFullFitResults->SetBranchAddress("A0"   , &A0     );
    cboFullFitResults->SetBranchAddress("R"    , &R      );
    cboFullFitResults->SetBranchAddress("phi"  , &phi    );
    cboFullFitResults->SetBranchAddress("tau"  , &tau    );

    cboFullFitResults->SetBranchAddress("w_CBO", &fullFit_wx1);
    cboFullFitResults->SetBranchAddress("ANx1" , &fullFit_Ax1);
    cboFullFitResults->SetBranchAddress("pNx1" , &fullFit_phix1);
    cboFullFitResults->SetBranchAddress("LM"   , &LM     );

    cboFullFitResults->GetEntry(0);
    cboFullFit->Close();

    // declare tree to hold new fit results
    fitresults->Branch("fitStart", &fitStart);
    fitresults->Branch("fitStop",  &fitStop );
    fitresults->Branch("windowBins", &windowBins);
    fitresults->Branch("caloNum", &desiredCalo);

    fitresults->Branch("pBinTau", &pBinTau);
    fitresults->Branch("pBinPhi", &pBinPhi);
    fitresults->Branch("includeMopTerm", &includeMopTerm);
    fitresults->Branch("constrainMop", &constrainMop);
    fitresults->Branch("constrainTau", &constraintau);

    fitresults->Branch("seed", &sidx);
    fitresults->Branch("windowNo", &windowNo);

    fitresults->Branch("chisq", &chisq);
    fitresults->Branch("rchisq", &rchisq);

    //5-param parameters
    fitresults->Branch("N0", &N0);
    fitresults->Branch("A0", &A0);
    fitresults->Branch("tau", &tau);
    fitresults->Branch("R", &R);
    fitresults->Branch("phi", &phi);

    fitresults->Branch("N0err", &N0err);
    fitresults->Branch("A0err", &A0err);
    fitresults->Branch("tauerr", &tauerr);
    fitresults->Branch("Rerr", &Rerr);
    fitresults->Branch("phierr", &phierr);

    // 1-CBO parameters
    fitresults->Branch("w_CBO", &w_CBO);
    fitresults->Branch("alpha_CBO", &alpha_CBO);
    fitresults->Branch("beta_CBO", &beta_CBO);
    fitresults->Branch("A_CBO", A_CBO);
    fitresults->Branch("phi_CBO", &phi_CBO);

    fitresults->Branch("w_CBOerr", &w_CBOerr);
    fitresults->Branch("alpha_CBOerr", &alpha_CBOerr);
    fitresults->Branch("beta_CBOerr", &beta_CBOerr);
    fitresults->Branch("phi_CBOerr", &phi_CBOerr);

    // 2-CBO parameters
    fitresults->Branch("alpha_2CBO", &alpha_2CBO);
    fitresults->Branch("beta_2CBO", &beta_2CBO);
    fitresults->Branch("A_2CBO", A_2CBO);
    fitresults->Branch("phi_2CBO", &phi_2CBO);

    fitresults->Branch("alpha_2CBOerr", &alpha_2CBOerr);
    fitresults->Branch("beta_2CBOerr", &beta_2CBOerr);
    fitresults->Branch("phi_2CBOerr", &phi_2CBOerr);

    // 1-Y parameters
    fitresults->Branch("w_y", &w_y);
    fitresults->Branch("alpha_y", &alpha_y);
    fitresults->Branch("beta_y", &beta_y);
    fitresults->Branch("A_y", A_y);
    fitresults->Branch("phi_y", &phi_y);

    fitresults->Branch("w_yerr", &w_yerr);
    fitresults->Branch("alpha_yerr", &alpha_yerr);
    fitresults->Branch("beta_yerr", &beta_yerr);
    fitresults->Branch("phi_yerr", &phi_yerr);

    // 2-Y parameters
    fitresults->Branch("w_vw", &w_vw);
    fitresults->Branch("alpha_vw", &alpha_vw);
    fitresults->Branch("beta_vw", &beta_vw);
    fitresults->Branch("A_vw", A_vw);
    fitresults->Branch("phi_vw", &phi_vw);

    fitresults->Branch("w_vwerr", &w_vwerr);
    fitresults->Branch("alpha_vwerr", &alpha_vwerr);
    fitresults->Branch("beta_vwerr", &beta_vwerr);
    fitresults->Branch("phi_vwerr", &phi_vwerr);

    // LM parameters
    fitresults->Branch("LM", &LM);
    fitresults->Branch("LMerr", &LMerr);

    // assign number of bins to this window
    // use longer windows at later times
    // maintain equal statistics in each window

    fitrangelow  = startBin + (windowNo);

    if (windowNo == 0){
        windowBins = 240;
    }else{
        double c = 0.55;
        double insideLog = 1.0 - c*c*exp(-((startBin-fitrangelow)*0.1492)/tau) 
            + c*c*exp(-((startBin+117-fitrangelow)*0.1492)/tau);
        // if insideLog term is below threshold, set deltaT to some max window size
        double deltaT;
        if (insideLog<0.08){
            deltaT = 1100*0.1492;
        }else{
            deltaT = -tau*log( 1.0 
                    - c*c*exp(-((startBin-fitrangelow)*0.1492)/tau) 
                    + c*c*exp(-((startBin+117-fitrangelow)*0.1492)/tau) );
        }

        printf("deltaT %f\n", deltaT);
        windowBins = int(deltaT/0.1492);
        // now round to nearest 10.0
        windowBins = int (windowBins/10.0) * 10;
        // if deltaT < 17us, set it back to 17us
        if (windowBins < 240){
            windowBins = 240;
        }

    }
    windowBins = 240; // hard-coding 34us bins for now

    printf("windowBins = %i\n", windowBins);

    fitrangehigh = fitrangelow+windowBins;
    fitStart     = startTime + fitrangelow*0.1492;
    fitStop      = startTime + fitrangehigh*0.1492;

    // if we are extending too far out, stop fitting
    if (fitrangehigh >= nBins){
        exit(0);
    }

    printf("fitting sliding window %i from bin %i at %f to bin %i at %f (%i bins)\n\n", 
            windowNo, fitrangelow, fitStart, fitrangehigh, fitStop, windowBins);

    int seed = desiredSeed;
    sidx = desiredSeed;

    // Read in the wiggle histograms from a file
    printf("Reading data from %s\n", inputFilename);

    // read in histograms for positrons, pileup, muonloss, etc.
    TFile* input = new TFile(inputFilename,"READ");
    TH1I* wiggletemp   = (TH1I*) input->Get(Form("fitHist_seed%d_calo%i" , seed, desiredCalo));
    TH1F* pileup2temp  = (TH1F*) input->Get(Form("pileup2_seed%d_calo%i" , seed, desiredCalo));
    TH1F* pileup3temp  = (TH1F*) input->Get(Form("pileup3_seed%d_calo%i" , seed, desiredCalo));
    TH1D* muonlosstemp = (TH1D*) input->Get(Form("muonLoss_seed%d"       , seed              ));
    TH1F* corrtemp     = (TH1F*) input->Get(Form("varcorr_seed%d_calo%i" , seed, desiredCalo));

    // create the corrected wiggle with positron + pileup data
    for(int j=0; j<nBins; j++){
        double counts = 0;
        counts += wiggletemp->GetBinContent(j+1);
        counts -= pileup2temp->GetBinContent(j+1);
        counts -= pileup3temp->GetBinContent(j+1);
        wiggle->SetBinContent(j+1, counts);
        wiggleraw->SetBinContent(j+1, wiggletemp->GetBinContent(j+1));
        if(corrtemp){
            varcorr->SetBinContent(j+1, corrtemp->GetBinContent(j+1));
        }
    }

    // populate the lambda histogram
    double integral = 0;
    for(int j=0; j<startBin; j++){
        double ml = muonlosstemp->GetBinContent(j+1)/muonlosstemp->Integral();
        double wt = exp((j+0.5)*0.1492/64.44);
        integral += ml*wt;
    }

    for(int j=startBin; j<nBins; j++){
        double ml = muonlosstemp->GetBinContent(j+1)/muonlosstemp->Integral();
        double wt = exp((j+0.5)*0.1492/64.44);
        integral += ml*wt;
        lambda->SetBinContent(j+1, integral);
    }
    input->Close();
    fprintf(stderr, "Done Reading data\n");

    // Now do a chisq fit
    int noParams = 17; 
    double par[noParams];
    double errorplus[noParams];
    double errorminus[noParams];

    TMinuit minimizer(noParams);
    minimizer.Command("SET PRINTOUT 1"); // change to level 1
    minimizer.Command("SET NOWARNINGS");
    minimizer.SetFCN(minuitFunction);
    minimizer.fGraphicsMode = false;

    // allow basic 5-param fit to float except R and tau
    minimizer.DefineParameter(0, "N0", N0, 100, 0, 0);
    minimizer.DefineParameter(1, "tau", tau, 0.0, 60, 70); // FIX
    minimizer.DefineParameter(2, "A0", A0, 0.0001, 0, 0); 
    minimizer.DefineParameter(3, "phi", phi, 0.01, 0, 0); 
    minimizer.DefineParameter(4, "R", R, 0, -1000, 1000); // FIX

    // 1-CBO terms
    minimizer.DefineParameter(5, "alpha_CBO", 0.0, 0.001, 0, 0);
    minimizer.DefineParameter(6, "w_CBO", fullFit_wx1, 0.001, 0, 0);
    minimizer.DefineParameter(7, "beta_CBO", 0.0, 0.001, 0, 0);

    // 2-CBO terms
    minimizer.DefineParameter(8, "alpha_2CBO", 0.0, 0.001, 0, 0);
    minimizer.DefineParameter(9, "beta_2CBO", 0.0, 0.001, 0, 0);

    // 1-Y terms
    minimizer.DefineParameter(10, "alpha_y", 0.0, 0.001, 0, 0);
    minimizer.DefineParameter(11, "w_y", 0.0, 0.001, 0, 0);
    minimizer.DefineParameter(12, "beta_y", 0.0, 0.001, 0, 0);

    // 2-Y terms
    minimizer.DefineParameter(13, "alpha_vw", 0.0, 0.001, 0, 0);
    minimizer.DefineParameter(14, "w_vw", 0.0, 0.001, 0, 0);
    minimizer.DefineParameter(15, "beta_vw", 0.0, 0.001, 0, 0);
 
    // LM term
    minimizer.DefineParameter(16, "LM", LM, 0.0, -0.1, 0.1); // FIX

    // now go through various stages of fitting 
    printf("MINUIT - FIT ONLY WIGGLE\n");
    minimizer.Command("FIX 5");
    minimizer.Command("FIX 6");
    minimizer.Command("FIX 7");
    minimizer.Command("FIX 8");
    minimizer.Command("FIX 9");
    minimizer.Command("FIX 10");
    minimizer.Command("FIX 11");
    minimizer.Command("FIX 12");
    minimizer.Command("FIX 13");
    minimizer.Command("FIX 14");
    minimizer.Command("FIX 15");
    minimizer.Command("FIX 16");
    minimizer.Migrad();

    printf("MINUIT - FIT ONLY 1-CBO + wiggle\n");
    minimizer.Command("RES");
    minimizer.Command("FIX 8");
    minimizer.Command("FIX 9");
    minimizer.Command("FIX 11");
    minimizer.Command("FIX 12");
    minimizer.Command("FIX 14");
    minimizer.Command("FIX 15");
    minimizer.Command("FIX 16");
    minimizer.Migrad();

    printf("MINUIT - FIT ONLY 2-CBO\n");
    // fix all but 2-CBO
    minimizer.Command("RES");
    minimizer.Command("FIX 0");
    minimizer.Command("FIX 1");
    minimizer.Command("FIX 2");
    minimizer.Command("FIX 3");
    minimizer.Command("FIX 4");
    minimizer.Command("FIX 5");
    minimizer.Command("FIX 6");
    minimizer.Command("FIX 7");
    minimizer.Command("FIX 10");
    minimizer.Command("FIX 11");
    minimizer.Command("FIX 12");
    minimizer.Command("FIX 13");
    minimizer.Command("FIX 14");
    minimizer.Command("FIX 15");
    minimizer.Command("FIX 16");
    minimizer.Migrad();

    printf("MINUIT - FIT wiggle + 1-CBO + 2-CBO\n");
    minimizer.Command("RES");
    minimizer.Command("FIX 10");
    minimizer.Command("FIX 11");
    minimizer.Command("FIX 12");
    minimizer.Command("FIX 13");
    minimizer.Command("FIX 14");
    minimizer.Command("FIX 15");
    minimizer.Command("FIX 16");
    minimizer.Migrad();

    printf("MINUIT - FIT ONLY 1-y\n");
    // fix all but 1-y
    minimizer.Command("RES");
    minimizer.Command("FIX 0");
    minimizer.Command("FIX 1");
    minimizer.Command("FIX 2");
    minimizer.Command("FIX 3");
    minimizer.Command("FIX 4");
    minimizer.Command("FIX 5");
    minimizer.Command("FIX 6");
    minimizer.Command("FIX 7");
    minimizer.Command("FIX 8");
    minimizer.Command("FIX 9");
    minimizer.Command("FIX 13");
    minimizer.Command("FIX 14");
    minimizer.Command("FIX 15");
    minimizer.Command("FIX 16");
    minimizer.Migrad();

    printf("MINUIT - FIT ONLY 2-y\n");
    // fix all but 2-y
    minimizer.Command("RES");
    minimizer.Command("FIX 0");
    minimizer.Command("FIX 1");
    minimizer.Command("FIX 2");
    minimizer.Command("FIX 3");
    minimizer.Command("FIX 4");
    minimizer.Command("FIX 5");
    minimizer.Command("FIX 6");
    minimizer.Command("FIX 7");
    minimizer.Command("FIX 8");
    minimizer.Command("FIX 9");
    minimizer.Command("FIX 10");
    minimizer.Command("FIX 11");
    minimizer.Command("FIX 12");
    minimizer.Command("FIX 16");
    minimizer.Migrad();

    printf("MINUIT - FIT ALTOGETHER\n");
    // fit all together
    minimizer.Command("RES");
    minimizer.Migrad();

    printf("MINUIT - MINOS\n");
    minimizer.Command("RES");
    minimizer.Command("MIG 25000 0.01");
    minimizer.Command("MINO 25000");

    TString varname[noParams];
    TString chnam;
    double val, err, xlolim, xuplim;
    int iuint;
    double fmin, fedm, errdef;
    int npari, nparx, istat;
    double eplus, eminus, eparab, globc;

    for(int k=0; k<noParams; k++){
        minimizer.mnpout(k, chnam, val, err, xlolim, xuplim, iuint);
        par[k] = val;
        varname[k] = chnam;
        minimizer.mnerrs(k, eplus, eminus, eparab, globc);
        errorplus[k] = eplus;
        errorminus[k] = eminus;
    }

    for(int k=0; k<noParams; k++){
        fprintf(stderr, "%s \t %12.5e + %9.3e %9.3e\n",
                varname[k].Data(), par[k], errorplus[k], errorminus[k]);
    }
    minimizer.mnstat(fmin, fedm, errdef, npari, nparx, istat);
    int noFreeParams = minimizer.GetNumFreePars();
    fprintf(stderr, "chisq/ndf : %f / %d\n", fmin, fitrangehigh-fitrangelow - noFreeParams);
    fprintf(stderr, "reduced chisq: %f +- %f\n", fmin/(fitrangehigh-fitrangelow - noFreeParams), sqrt(2./(fitrangehigh-fitrangelow-noFreeParams)));
    minimizer.Command("SHO COR");
    double covariance[noParams][noParams];
    minimizer.mnemat(&covariance[0][0],noParams);

    chisq  = fmin;
    rchisq = fmin/(fitrangehigh-fitrangelow-noFreeParams);

    // set 5-param variables and their errors
    N0 = par[0];
    tau = par[1];
    A0 = par[2];
    phi = par[3];
    R = par[4];

    N0err = sqrt(-errorplus[0]*errorminus[0]);
    tauerr = sqrt(-errorplus[1]*errorminus[1]);
    A0err = sqrt(-errorplus[2]*errorminus[2]);
    phierr = sqrt(-errorplus[3]*errorminus[3]);
    Rerr = sqrt(-errorplus[4]*errorminus[4]);

    // set 1-CBO variables and their errors
    alpha_CBO = par[5];
    w_CBO = par[6];
    beta_CBO = par[7];
    A_CBO = sqrt(alpha_CBO*alpha_CBO + beta_CBO*beta_CBO);
    phi_CBO = invert(alpha_CBO, beta_CBO);

    alpha_CBOerr = sqrt(-errorplus[5]*errorminus[5]);
    w_CBOerr = sqrt(-errorplus[6]*errorminus[6]);
    beta_CBOerr = sqrt(-errorplus[7]*errorminus[7]);
    A_CBOerr = 0.0;
    phi_CBOerr = 0.0;

    // set 2-CBO variables and their errors
    alpha_2CBO = par[8];
    beta_2CBO = par[9];
    A_2CBO = sqrt(alpha_2CBO*alpha_2CBO + beta_2CBO*beta_2CBO);
    phi_2CBO - invert(alpha_2CBO, beta_2CBO);

    alpha_2CBOerr = sqrt(-errorplus[8]*errorminus[8]);
    beta_2CBOerr = sqrt(-errorplus[9]*errorminus[9]);
    A_2CBOerr = 0.0;
    phi_2CBOerr = 0.0;

    // set 1-Y variables and their errors
    alpha_y = par[10];
    w_y = par[11];
    beta_y = par[12];
    A_y = sqrt(alpha_y*alpha_y + beta_y*beta_y);
    phi_y = invert(alpha_y, beta_y);

    alpha_yerr = sqrt(-errorplus[10]*errorminus[10]);
    w_yerr  = sqrt(-errorplus[11]*errorminus[11]);
    beta_yerr = sqrt(-errorplus[12]*errorminus[12]);
    A_yerr = 0.0;
    phi_yerr = 0.0;

    // set 2-Y variables and their errors
    alpha_vw = par[13];
    w_vw = par[14];
    beta_vw = par[15];
    A_vw = sqrt(alpha_vw*alpha_vw + beta_vw*beta_vw);
    phi_vw = invert(alpha_vw, beta_vw);

    alpha_vwerr = sqrt(-errorplus[13]*errorminus[13]);
    w_vwerr = sqrt(-errorplus[14]*errorminus[14]);
    beta_vwerr = sqrt(-errorplus[15]*errorminus[15]);
    A_vwerr = 0.0;
    phi_vwerr = 0.0;

    // set LM term
    LM = par[16];
    LMerr = sqrt(-errorplus[16]*errorminus[16]);

    fitresults->Fill();

    double sum = 0;
    for(int i=fitrangelow; i<fitrangehigh; i++){
        bestfit->Fill((i-0.5)*0.1492 + startTime, calcnu(par, i));
        double pucorr = calcnu(par, i)*calcnu(par, i)/(calcnu(par, i-0.5)*calcnu(par, i+0.5));
        sum += pucorr;
    }

    fprintf(stderr, "Mean corr %f\n", sum/nBins);
    //}; // end loop over seeds

    output->cd();
    TH1D* residua = new TH1D("residua", "residua", nBins, startTime, 650.0644);
    TTree* pulls = new TTree("pulls", "pulls");
    double pull;
    pulls->Branch("pull", &pull);
    for(int i=fitrangelow; i<fitrangehigh; i++){
        double residual = wiggle->GetBinContent(i) - bestfit->GetBinContent(i);
        residua->SetBinContent(i, residual);
        pull = residual/sqrt(wiggleraw->GetBinContent(i) + varcorr->GetBinContent(i));
        pulls->Fill();
    }

    output->cd();
    residua->Write();
    pulls->Write();
    fitresults->Write();  
    lambda->Write();
    bestfit->Write();
    wiggle->Write();
    output->Close();

}
