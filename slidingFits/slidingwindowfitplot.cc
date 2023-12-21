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
            "-p [s]    : CBO isolate filename (required)\n"
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
        char* &inputFilename, char* &outputFilename, char* &boostFilename, char* &cboIsolateFilename,
        char* &cboRefFilename, 
        bool* asymmetry, double* c_e_scale, int* seed, int* windowNo,
        char* &windowFilename){
    const char* const opts = "hf:o:i:b:a:c:s:p:q:w:n:W:";
    bool done = false;
    while(!done){
        const char ch = getopt(argc, argv, opts);
        switch(ch){
            case -1: done = true; break;
            case 'h': printhelp(); exit(0);
            case 'f': formatFilename = optarg; break;
            case 'i': inputFilename = optarg; break;
            case 'o': outputFilename = optarg; break;
            case 'p': cboIsolateFilename = optarg; break;
            case 'q': cboRefFilename = optarg; break;
            case 'b': boostFilename = optarg; break;
            case 'a': *asymmetry = bool(atoi(optarg)); break;
            case 'c': *c_e_scale = std::stod(optarg); break;
            case 's': *seed   = atoi(optarg); break;
            case 'w': *windowNo = atoi(optarg); break;
            case 'n': desiredCalo = atoi(optarg); break;
            case 'W': windowFilename = optarg; break;
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
    char* cboIsolateFilename;
    char* cboRefFilename;
    char* windowFilename;

    parse_cmdline(argc,argv,
            formatFilename, 
            inputFilename, outputFilename, boostFilename, cboIsolateFilename,cboRefFilename,
            &asymmetry, &c_e_scale, &desiredSeed, &windowNo, windowFilename);

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
    // printf("Opening file %s.\n", boostFilename);
    // TFile* boostdata  = new TFile(boostFilename, "READ");
    // TH1D*  boostgamma = (TH1D*) boostdata->Get("transform_gamma");
    // frlifetime = boostgamma->GetMean();
    // /*
    // double* boost = boostgamma->GetX();
    // double* freq = boostgamma->GetY();
    // double sum = 0;
    // double den = 0;
    // for(int i=0; i<boostgamma->GetN(); i++){
    //   sum += boost[i]*freq[i];
    //   den += freq[i];
    // }
    // frlifetime = sum/den;
    // */
    // for(int i=0; i<nBins; i++){
    //   double tt = startTime + (i+0.5)*0.1492;
    //   double mdl = exp(-tt/(2.1969811*frlifetime));
    //   double xact = 0;
    //   for(int j=0; j<boostgamma->GetNbinsX(); j++){
    //     xact += boostgamma->GetBinContent(j+1)/boostgamma->Integral()*exp(-tt/(2.1969811*boostgamma->GetBinCenter(j+1)));
    //   }
    //   expcorr[i] = xact/mdl;
    // }

    //*************
    // get momentum-binned weights
    //*************

    // TH1F* transform_dp_p0 = (TH1F*) boostdata->Get("transform_dp_p0");
    // TH1F* transform_gamma = (TH1F*) boostdata->Get("transform_gamma");
    // TH1F* transform_x     = (TH1F*) boostdata->Get("transform_x");

    // int pBins = int(transform_dp_p0->GetNbinsX());
    // std::pair <double,double> physicalRange (-0.0056799,0.0056085);  // in dp/p0

    // // first isolate the physical range of the spectra
    // // get the integral and bin numbers 
    // std::vector<double> pWeights(pBins, 0.0); 

    // for (int i=0; i<pBins; i++){
    //    double lowEdge = transform_dp_p0->GetBinLowEdge(i);
    //    if (physicalRange.first <= lowEdge) {
    // binRange.first = i;
    //       break;
    //    }
    // }
    // for (int i=pBins-1; i>=0; i--){
    //    double upEdge = transform_dp_p0->GetBinLowEdge(i+1);
    //    if (physicalRange.second >= upEdge) {
    // binRange.second = i;
    //       break;
    //    }
    // }

    // beta_0 = sqrt(1.0 - (1.0 / pow(29.3,2.0))); // this uses design gamma
    // R_0    = 7112; // design R_0 in mm
    // // now go through bins and store values for dp_p0, dPhi, gammas, tau, x
    // for (int i=0; i<pBins; i++){
    //     dp_p0 .push_back(transform_dp_p0       ->GetBinCenter(i));
    //     gammas.push_back(transform_gamma       ->GetBinCenter(i));
    //     dPhi  .push_back(-1.0 * transform_dp_p0 ->GetBinCenter(i)); // negative if fitting cos(wt - phi), pos if wt + phi
    //     x_e   .push_back(transform_x->GetBinCenter(i)); //convert from radial offset to radial position
    //     c_e   .push_back( c_e_scale * (2.0 * ((beta_0*beta_0)/R_0) * n 
    // 			* transform_x->GetBinCenter(i) 
    // 			* transform_dp_p0->GetBinCenter(i)) );
    //     printf("bin %i c_e %f\n", i, c_e.at(i));
    // }

    // printf("Physical range of dp/p0 is (%f,%f) from bins (%i,%i).\n",
    //   transform_dp_p0->GetBinLowEdge(binRange.first), transform_dp_p0->GetBinLowEdge(binRange.second+1),
    //   binRange.first, binRange.second);

    // double totalContent = transform_dp_p0->Integral(binRange.first, binRange.second);

    // // now find weights in only the physical range
    // totalContent = transform_dp_p0->Integral(binRange.first, binRange.second);
    // double sumWeights = 0.0;
    // for (int i=0; i<=pBins; i++){
    //    // find magicBin 
    //    double lowEdge = transform_dp_p0->GetBinLowEdge(i);
    //    double uppEdge = transform_dp_p0->GetBinLowEdge(i+1);

    //    if ( lowEdge <= 0.0  and 0.0 < uppEdge){
    //       magicBin = i;
    //       printf("magic bin is %i with tau %f\n", magicBin, 2.1969811*gammas[magicBin]);
    //    }

    //    // find weights according to momentum
    //    double thisWeight = transform_dp_p0->GetBinContent(i) / totalContent;

    //    // is this bin in a physical momentum bin?
    //    bool isPhysical = (i>=binRange.first and i<= binRange.second);

    //    if (pBinTau or pBinPhi){
    //       binWeights.push_back( isPhysical? thisWeight : 0.0 );
    //    }else{ // if we're not binning, set weight=1.0 in magic bin and 0.0 otherwise
    //       binWeights.push_back( (i==magicBin) ? 1.0 : 0.0 );
    //    }
    //    sumWeights += binWeights[i];
    //    weightedAvgTau += binWeights[i]*transform_gamma->GetBinCenter(i)*2.1969811;
    //    printf("bin %i momentum %f isPhysical? %d weight %f\n",i,dp_p0[i],isPhysical,binWeights[i]);
    // }

    // boostdata->Close();

    //*************
    // end get momentum-binned weights + taus
    //*************

    TH1D* bestfit = new TH1D("bestfit", "bestfit", nBins, startTime, 650.0644);
    TTree* fitresults = new TTree("fitresults", "fitresults");

    int sidx;
    double chisq, rchisq;
    double N0, tau, A0, phi, R;
    double w_CBO, A_CBO, phi_CBO;
    double w_CBOerr, A_CBOerr, phi_CBOerr;
    double LM;
    double N0err, tauerr, A0err, phierr, Rerr;
    double LMerr;
    double wCBO_phiCBO_cov;
    double ANy1, pNy1, Ky, A_ct;
    double ANx2, pNx2; 
    double ANx2err, pNx2err;
    double ANy2, pNy2;
    double ANy1err, pNy1err, Kyerr, A_cterr;
    double ANy2err, pNy2err;
    double w_y, w_yerr;
    double w_vw, w_vwerr;

    double fitStart = fitrangelow  * 0.1492;
    double fitStop  = fitrangehigh * 0.1492;

    double wCBO_expCoeff, wCBO_expCoefferr;
    double wCBO_expOffset, wCBO_expOffseterr;
    double wCBO_expT, wCBO_expTerr;
    double wCBO_linCoeff, wCBO_linCoefferr;
    double wCBO_const, wCBO_consterr;
    

    //*************
    // get CBO terms from CBO Isolate
    //*************

    printf("reading cbo reference file (full fit) from %s\n", cboRefFilename);
    TFile* cboFullFit = new TFile(cboRefFilename, "READ");

    double guess_ANx1, guess_w, guess_p;
    double guess_AVW, guess_pVW, guess_wVW;

    TTree* cboFullFitResults = (TTree*) cboFullFit->Get("fitresults");

    cboFullFitResults->SetBranchAddress("N0"   , &N0     );
    cboFullFitResults->SetBranchAddress("A0"   , &A0     );
    cboFullFitResults->SetBranchAddress("R"    , &R      );
    cboFullFitResults->SetBranchAddress("phi"  , &phi    );
    cboFullFitResults->SetBranchAddress("tau"  , &tau    );

    cboFullFitResults->SetBranchAddress("w_CBO", &guess_w);
    cboFullFitResults->SetBranchAddress("ANx1" , &guess_ANx1);
    cboFullFitResults->SetBranchAddress("ANx2" , &ANx2);
    cboFullFitResults->SetBranchAddress("pNx1" , &guess_p);
    cboFullFitResults->SetBranchAddress("pNx2" , &pNx2);
    cboFullFitResults->SetBranchAddress("ANy1" , &ANy1);
    cboFullFitResults->SetBranchAddress("pNy1" , &pNy1);
    cboFullFitResults->SetBranchAddress("ANy2" , &ANy2);
    cboFullFitResults->SetBranchAddress("pNy2" , &pNy2);
    cboFullFitResults->SetBranchAddress("Ky"   , &Ky  );
    cboFullFitResults->SetBranchAddress("A_ct" , &A_ct);
    cboFullFitResults->SetBranchAddress("LM"   , &LM     );

    cboFullFitResults->GetEntry(0);

    guess_AVW = ANy2;
    guess_pVW = pNy2;
    guess_wVW = wy(Ky, guess_w);

    printf("A_CBO %f w_CBO %f phi_CBO %f N0 %f R %f tau %f \n", guess_ANx1, guess_w, guess_p, N0, R, tau);
    // we are setting these results as starting guesses for fitCboIsolate
    cboFullFit->Close();

    printf("Reading cboRef file %s\n", cboIsolateFilename);
    TFile* cboRef = new TFile(cboIsolateFilename, "READ");
    TH1D* cboHisto = (TH1D*) cboRef->Get("dataDivFnc_wo_CBO");

    for(int j=1; j<nBins; j++){
        //printf("j %i cboHisto->GetBinContent(%i) = %f \n", j, j, cboHisto->GetBinContent(j));
        cboIsolate->SetBinContent(j, cboHisto->GetBinContent(j));
    }
    cboRef->Close();

    //*************
    // end: get CBO terms from CBO Isolate
    //*************

    fitresults->Branch("fitStart", &fitStart);
    fitresults->Branch("fitStop",  &fitStop );
    fitresults->Branch("windowBins", &windowBins);

    fitresults->Branch("pBinTau", &pBinTau);
    fitresults->Branch("pBinPhi", &pBinPhi);
    fitresults->Branch("includeMopTerm", &includeMopTerm);
    fitresults->Branch("constrainMop", &constrainMop);
    fitresults->Branch("constrainTau", &constraintau);

    fitresults->Branch("seed", &sidx);
    fitresults->Branch("windowNo", &windowNo);

    fitresults->Branch("chisq", &chisq);
    fitresults->Branch("rchisq", &rchisq);

    fitresults->Branch("N0", &N0);
    fitresults->Branch("A0", &A0);
    fitresults->Branch("tau", &tau);
    fitresults->Branch("R", &R);
    fitresults->Branch("phi", &phi);

    fitresults->Branch("w_CBO", &w_CBO);
    fitresults->Branch("A_CBO", &A_CBO);
    fitresults->Branch("phi_CBO", &phi_CBO);

    fitresults->Branch("ANx2", &ANx2);
    fitresults->Branch("pNx2", &pNx2);

    fitresults->Branch("ANy1", &ANy1);
    fitresults->Branch("pNy1", &pNy1);
    fitresults->Branch("w_y" , &w_y );

    fitresults->Branch("ANy2", &ANy2);
    fitresults->Branch("pNy2", &pNy2);
    fitresults->Branch("w_vw", &w_vw);

    fitresults->Branch("LM", &LM);

    fitresults->Branch("wCBO_expCoeff", &wCBO_expCoeff);
    fitresults->Branch("wCBO_expOffset", &wCBO_expOffset);
    fitresults->Branch("wCBO_expT", &wCBO_expT);
    fitresults->Branch("wCBO_linCoeff", &wCBO_linCoeff);
    fitresults->Branch("wCBO_const", &wCBO_const);

    fitresults->Branch("N0err", &N0err);
    fitresults->Branch("A0err", &A0err);
    fitresults->Branch("tauerr", &tauerr);
    fitresults->Branch("Rerr", &Rerr);
    fitresults->Branch("phierr", &phierr);

    fitresults->Branch("w_CBOerr", &w_CBOerr);
    fitresults->Branch("A_CBOerr", &A_CBOerr);
    fitresults->Branch("phi_CBOerr", &phi_CBOerr);

    fitresults->Branch("LMerr", &LMerr);

    fitresults->Branch("wCBO_phiCBO_cov", &wCBO_phiCBO_cov);

    fitresults->Branch("ANx2err", &ANx2err);
    fitresults->Branch("pNx2err", &pNx2err);

    fitresults->Branch("ANy1err", &ANy1err);
    fitresults->Branch("pNy1err", &pNy1err);
    fitresults->Branch("w_yerr" , &w_yerr );

    fitresults->Branch("ANy2err", &ANy2err);
    fitresults->Branch("pNy2err", &pNy2err);
    fitresults->Branch("w_vwerr", &w_vwerr);

    fitresults->Branch("wCBO_expCoefferr", &wCBO_expCoefferr);
    fitresults->Branch("wCBO_expOffseterr", &wCBO_expOffseterr);
    fitresults->Branch("wCBO_expTerr", &wCBO_expTerr);
    fitresults->Branch("wCBO_linCoefferr", &wCBO_linCoefferr);
    fitresults->Branch("wCBO_consterr", &wCBO_consterr);

    // assign number of bins to this window
    // use longer windows at later times
    // maintain equal statistics in each window

    fitrangelow  = startBin + (windowNo);

    if (windowNo == 0){
        windowBins = 120;
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
        if (windowBins < 120){
            windowBins = 120;
        }

    }
    printf("windowBins = %i\n", windowBins);

    fitrangehigh = fitrangelow+windowBins;
    fitStart     = startTime + fitrangelow*0.1492; //TODO: add startTime
    fitStop      = startTime + fitrangehigh*0.1492;

    // if we are extending too far out, stop fitting
    if (fitrangehigh >= nBins){
        exit(0);
    }

    printf("fitting sliding window %i from bin %i at %f to bin %i at %f (%i bins)\n\n", 
            windowNo, fitrangelow, fitStart, fitrangehigh, fitStop, windowBins);

    //for(int seed=desiredSeed; seed<desiredSeed+1; seed++){
    int seed = desiredSeed;
    sidx = desiredSeed; //seed;
    fprintf(stderr, "Reading data\n");
    printf("Reading data from %s\n", inputFilename);
    // Read in the wiggle histograms from a file
    TFile* input = new TFile(inputFilename,"READ");
    TH1I* wiggletemp   = (TH1I*) input->Get(Form("fitHist_seed%d_calo%i" , seed, desiredCalo));
    TH1F* pileup2temp  = (TH1F*) input->Get(Form("pileup2_seed%d_calo%i" , seed, desiredCalo));
    TH1F* pileup3temp  = (TH1F*) input->Get(Form("pileup3_seed%d_calo%i" , seed, desiredCalo));
    TH1D* muonlosstemp = (TH1D*) input->Get(Form("muonLoss_seed%d"       , seed              ));
    TH1F* corrtemp     = (TH1F*) input->Get(Form("varcorr_seed%d_calo%i" , seed, desiredCalo));

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

    double integral = 0;
    // go from bin 0 to start time bin, remember that muonloss goes from 0 - 650us
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
    printf("fitting cbo isolate\n");
    //*************
    // get CBO terms from CBO Isolate
    //*************
    //

    TH1D* cboIsolateData = new TH1D(Form("cboIsolateWindowData_window%i",windowNo),
            Form("cboIsolateWindowData_window%i; time [#mus]; arb units",windowNo),windowBins,fitStart,fitStop);
    for(int iso_i = 0; iso_i < windowBins; iso_i ++ ){
        cboIsolateData->SetBinContent(iso_i+1, cboIsolate->GetBinContent(fitrangelow+iso_i));
    }
    TF1* isoFit = new TF1("isoFit", "(1.0 + [0]*cos(([1])*x - [2]))", fitrangelow*0.1492, fitrangehigh*0.1492);

    TMinuit cboMinimizer(3);
    cboMinimizer.Command("SET PRINTOUT 1"); // change to level 1
    cboMinimizer.Command("SET NOWARNINGS");
    cboMinimizer.SetFCN(fitCboIsolate);
    cboMinimizer.fGraphicsMode = false;

    // make a guess for phi based on previous window...
    // but use CBO isolate for A and w
    cboMinimizer.DefineParameter(0, "A_CBO"  , guess_ANx1 , 1E-5, 0  ,   0);//-1E-5, 0.002);
    cboMinimizer.DefineParameter(1, "w_CBO"  , guess_w , 0.01, 2.0, 3.0);
    cboMinimizer.DefineParameter(2, "phi_CBO", guess_p , 0.01, 0  ,2*M_PI);
    cboMinimizer.Migrad();

    double cbo_err, cbo_xlolim, cbo_xuplim;
    int cbo_iuint;
    TString cbo_chnam;

    cboMinimizer.mnpout(0, cbo_chnam, A_CBO, cbo_err, cbo_xlolim, cbo_xuplim, cbo_iuint);
    cboMinimizer.mnpout(1, cbo_chnam, w_CBO, cbo_err, cbo_xlolim, cbo_xuplim, cbo_iuint);
    cboMinimizer.mnpout(2, cbo_chnam, phi_CBO, cbo_err, cbo_xlolim, cbo_xuplim, cbo_iuint);

    printf("setting A_CBO=%f w_CBO=%f, phi_CBO=%f\n", A_CBO, w_CBO, phi_CBO);
    // set same results into TF1 for plotting
    isoFit->SetParameter(0,A_CBO);
    isoFit->SetParameter(1,w_CBO);
    isoFit->SetParameter(2,phi_CBO);

    //*************
    // end: get CBO terms from CBO Isolate
    //*************
    //
    // read window 1300 to get wCBO
    double window_wCBO;
    if (windowFilename){
        TFile *windowFile = new TFile(windowFilename, "READ");
        TTree *windowTree = (TTree*)windowFile->Get("fitresults");
        windowTree->SetBranchAddress("w_CBO", &window_wCBO);
        windowTree->GetEntry(0);
    }else{
        window_wCBO = w_CBO;
    }

    printf("done fitting cbo isolate\n");

    // Now do a 28 parameter fit
    double par[22];
    double errorplus[22];
    double errorminus[22];

    TMinuit minimizer(22);
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

    minimizer.DefineParameter(5, "A_CBO", A_CBO, 1E-5, 0, 0);
    minimizer.DefineParameter(6, "w_CBO", window_wCBO, 0.001, 0, 0);
    minimizer.DefineParameter(7, "phi_CBO", phi_CBO, 0.0001, 0, 0);

    minimizer.DefineParameter(8, "ANx2", ANx2, 0, 0, 0);
    minimizer.DefineParameter(9, "pNx2", pNx2, 0, 0, 0);

    minimizer.DefineParameter(10, "ANy1", 0, 0, 0, 0);
    minimizer.DefineParameter(11, "pNy1", 0, 0, 0, 0);
    minimizer.DefineParameter(12, "w_y" , 0, 0, 0, 0);

    minimizer.DefineParameter(13, "ANy2", guess_AVW, 0.0001, 0, 0);
    minimizer.DefineParameter(14, "pNy2", guess_pVW, 0, 0, 0);
    minimizer.DefineParameter(15, "w_vw", guess_wVW, 0, 0, 0);

    minimizer.DefineParameter(16, "LM", LM, 0.0, -0.1, 0.1); // FIX

    minimizer.DefineParameter(17, "wCBO_expCoeff", 0, 0, 0, 0);
    minimizer.DefineParameter(18, "wCBO_expOffset", 0, 0, 0, 0);
    minimizer.DefineParameter(19, "wCBO_expT", 0, 0, 0, 0);
    minimizer.DefineParameter(20, "wCBO_linCoeff", 0, 0, 0, 0);
    minimizer.DefineParameter(21, "wCBO_const", 0, 0, 0, 0);

    printf("MINUIT - FIT ONLY WIGGLE\n");
    // fix CBO parameters and only fit wiggle
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

    printf("MINUIT - FIT ONLY 1-CBO + other amplitudes \n");
    // fix wiggle and only fit CBO + As for 2-cbo and y
    minimizer.Command("RES");
    minimizer.Command("FIX 0");
    minimizer.Command("FIX 1");
    minimizer.Command("FIX 2");
    minimizer.Command("FIX 3");
    minimizer.Command("FIX 4");
    minimizer.Command("FIX 9");
    minimizer.Command("FIX 11");
    minimizer.Command("FIX 12");
    minimizer.Command("FIX 14");
    minimizer.Command("FIX 15");
    minimizer.Command("FIX 16");
    minimizer.Migrad();

    printf("MINUIT - FIT ONLY 1-CBO\n");
    // fix wiggle and only fit CBO + As for 2-cbo and y
    minimizer.Command("RES");
    minimizer.Command("FIX 0");
    minimizer.Command("FIX 1");
    minimizer.Command("FIX 2");
    minimizer.Command("FIX 3");
    minimizer.Command("FIX 4");
    minimizer.Command("FIX 8");
    minimizer.Command("FIX 9");
    minimizer.Command("FIX 11");
    minimizer.Command("FIX 12");
    minimizer.Command("FIX 13");
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

    printf("MINUIT - REFIT 5-PARAM\n");
    // fit all but 5-params and refit
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

    printf("MINUIT - 2CBO + Y\n");
    // fix all but 2-CBO and y-terms
    minimizer.Command("RES");
    minimizer.Command("FIX 0");
    minimizer.Command("FIX 1");
    minimizer.Command("FIX 2");
    minimizer.Command("FIX 3");
    minimizer.Command("FIX 4");
    minimizer.Command("FIX 5");
    minimizer.Command("FIX 6");
    minimizer.Command("FIX 7");
    minimizer.Command("FIX 16");
    minimizer.Migrad();

    printf("MINUIT - MINOS\n");
    minimizer.Command("RES");
    minimizer.Command("MIG 25000 0.01");
    minimizer.Command("MINO 25000");

    TString varname[22];
    TString chnam;
    double val, err, xlolim, xuplim;
    int iuint;
    double fmin, fedm, errdef;
    int npari, nparx, istat;
    double eplus, eminus, eparab, globc;

    for(int k=0; k<22; k++){
        minimizer.mnpout(k, chnam, val, err, xlolim, xuplim, iuint);
        par[k] = val;
        varname[k] = chnam;
        minimizer.mnerrs(k, eplus, eminus, eparab, globc);
        errorplus[k] = eplus;
        errorminus[k] = eminus;
    }

    for(int k=0; k<22; k++){
        fprintf(stderr, "%s \t %12.5e + %9.3e %9.3e\n",
                varname[k].Data(), par[k], errorplus[k], errorminus[k]);
    }
    minimizer.mnstat(fmin, fedm, errdef, npari, nparx, istat);
    int noFreeParams = minimizer.GetNumFreePars();
    fprintf(stderr, "chisq/ndf : %f / %d\n", fmin, fitrangehigh-fitrangelow - noFreeParams);
    fprintf(stderr, "reduced chisq: %f +- %f\n", fmin/(fitrangehigh-fitrangelow - noFreeParams), sqrt(2./(fitrangehigh-fitrangelow-noFreeParams)));
    minimizer.Command("SHO COR");
    double covariance[22][22];
    minimizer.mnemat(&covariance[0][0],22);

    chisq  = fmin;
    rchisq = fmin/(fitrangehigh-fitrangelow-noFreeParams);
    N0 = par[0];
    tau = par[1];
    A0 = par[2];
    phi = par[3];
    R = par[4];

    A_CBO = par[5];
    w_CBO = par[6];
    phi_CBO = par[7];

    ANx2 = par[8];
    pNx2 = par[9];

    ANy1 = par[10];
    pNy1 = par[11];
    w_y  = par[12];

    ANy2 = par[13];
    pNy2 = par[14];
    w_vw = par[15];

    LM = par[16];

    wCBO_expCoeff = par[17];
    wCBO_expOffset = par[18];
    wCBO_expT = par[19];
    wCBO_linCoeff = par[20];
    wCBO_const = par[21];

    N0err = sqrt(-errorplus[0]*errorminus[0]);
    tauerr = sqrt(-errorplus[1]*errorminus[1]);
    A0err = sqrt(-errorplus[2]*errorminus[2]);
    phierr = sqrt(-errorplus[3]*errorminus[3]);
    Rerr = sqrt(-errorplus[4]*errorminus[4]);

    A_CBOerr = sqrt(-errorplus[5]*errorminus[5]);
    w_CBOerr = sqrt(-errorplus[6]*errorminus[6]);
    phi_CBOerr = sqrt(-errorplus[7]*errorminus[7]);

    ANx2err = sqrt(-errorplus[8]*errorminus[8]);
    pNx2err = sqrt(-errorplus[9]*errorminus[9]);

    ANy1err = sqrt(-errorplus[10]*errorminus[10]);
    pNy1err = sqrt(-errorplus[11]*errorminus[11]);
    w_yerr  = sqrt(-errorplus[12]*errorminus[12]);

    ANy2err = sqrt(-errorplus[13]*errorminus[13]);
    pNy2err = sqrt(-errorplus[14]*errorminus[14]);
    w_vwerr = sqrt(-errorplus[15]*errorminus[15]);

    LMerr = sqrt(-errorplus[16]*errorminus[16]);

    wCBO_phiCBO_cov = covariance[3][4]; // covariance matrix is only floating parameters, so indices different

    fitresults->Fill();

    double sum = 0;
    for(int i=fitrangelow; i<fitrangehigh; i++){
        bestfit->Fill((i-0.5)*0.1492 + startTime, calcnu(par, i));
        double pucorr = calcnu(par, i)*calcnu(par, i)/(calcnu(par, i-0.5)*calcnu(par, i+0.5));
        sum += pucorr;
    }
    fprintf(stderr, "Mean corr %f\n", sum/nBins);
    //}; // end loop over seeds

    TFile* output = new TFile(outputFilename, "RECREATE");
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

    //TCanvas *isoResult = new TCanvas(Form("isoFit_window%i", windowNo),Form("isoFit_window%i", windowNo));
    //cboIsolateData->Draw();
    //isoFit->Draw("SAME");
    //isoResult->Write();

    cboIsolateData->Write();
    isoFit->Write();

    residua->Write();
    pulls->Write();
    fitresults->Write();  
    lambda->Write();
    bestfit->Write();
    wiggle->Write();
    output->Close();
}
