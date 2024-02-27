// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TGraphErrors.h"

// gm2 includes
#include "gm2util/blinders/Blinders.hh"
#include "CornellHistograms.hh"
#include "CornellFitAlgorithm.hh"
#include "tbb/parallel_for.h"
#include "readMinuitCommands.hh"

// cpp includes
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>

// fitting configuration
static double expcorr[nBins];		// no longer used
bool   fixtau       = false;		// fix to FR value
bool   constraintau = true;		// lifetime constraint via chisq penalty
bool   correcttau   = false; 		// no longer used
bool   asymmetry    = true;		// must be true for a-method, false for t-method
int    fitrangelow  = startBin; // TO DO: change back to startBin;
int    fitrangehigh = nBins; // TO DO: change back to nBins // 1340 for 200us
double frlifetime   = 0;

void printhelp(){
    printf(
            "Analysis Program\n"
            "Supply arguments as follows:\n"
            "-h        : Prints this help text\n"
            "-i [s]    : Input filename (required)\n"
            "-o [s]    : Output filename (required)\n"
            "-b [s]    : Boost filename (required)\n"
            "-f [s]    : Format filename (required)\n"
            "-a [i]    : Bool for asymm (1) or threshold (0)\n"
            "-c [d]    : Double for scale 0-1 to apply to c_e (default 1.0)\n"
            "-s [i]    : Int for seed (default 0)\n"
            "-n [i]    : Int for calo (default 0)\n"
          );
}

void parse_cmdline(int argc, char** argv, 
        char* &formatFilename, 
        char* &inputFilename, char* &outputFilename, char* &boostFilename,
        bool* asymmetry, double* c_e_scale, int* seed, int* calo){
    const char* const opts = "hf:o:i:b:a:c:s:n:";
    bool done = false;
    while(!done){
        const char ch = getopt(argc, argv, opts);
        switch(ch){
            case -1: done = true; break;
            case 'h': printhelp(); exit(0);
            case 'f': formatFilename = optarg; break;
            case 'i': inputFilename = optarg; break;
            case 'o': outputFilename = optarg; break;
            case 'b': boostFilename = optarg; break;
            case 'a': *asymmetry = bool(atoi(optarg)); break;
            case 'c': *c_e_scale = std::stod(optarg); break;
            case 's': *seed   = atoi(optarg); break;
            case 'n': *calo   = atoi(optarg); break;
            default: printhelp(); exit(1);
        }
    }
}

double invert(double c, double s){
    double i = atan(s/c);
    if(c < 0){
        i += M_PI;
        return i;
    }
    if(s < 0){
        i += 2*M_PI;
        return i;
    }
    return i;
}

void minuitFunction( __attribute__((unused)) int& nDim,
        __attribute__((unused)) double* gout,
        double& result, double par[],
        __attribute__((unused)) int flg){
    result = 0;
    for(int i=fitrangelow; i<=fitrangehigh; i++){
        double n = wiggle->GetBinContent(i);
        double dim[1] = {double(i)};
        double nu = calcnu(dim, par);
        //double chisq = pow(n-nu, 2)/n;
        if(asymmetry){
            result += pow(n-nu, 2)/varcorr->GetBinContent(i);
        }
        else{
            result += pow(n-nu, 2)/(wiggleraw->GetBinContent(i) + varcorr->GetBinContent(i));
        }
        /*
           double like = 0;
           if(nu > 0){
           if(n != 0)
           like = nu - n + n*log(n/nu);
           else
           like = nu;
           }
           */
        // For standard error computations
        //result += chisq;
    }
    if(constraintau and !pBinTau){
        double boostedlifetime = 2.1969811*frlifetime;
        double error = 2.1969811*0.004;
        double penalty = pow((par[1]-boostedlifetime)/error, 2.);
        result += penalty;
    }
    if(includeMopTerm and constrainMop){
        //double error = (par[27] * 2.1969811 * 0.004) / (2.1969811 * frlifetime);
        double error = (2.1969811*0.004) / (par[1]*par[1]);
        double penalty = pow(par[27]/error,2.); // diff to 0
        result += penalty; 
        //printf("par[27]: %20.15f frlifetime %20.15f error: %20.15f penalty: %20.15f\n",par[27], frlifetime, error, penalty);
    }
    //printf("result %f\n",result);
}

int main(int argc, char* argv[]){

    // get input arguments
    char* formatFilename;
    char* outputFilename;
    char* inputFilename ;
    char* boostFilename ;

    parse_cmdline(argc,argv,
            formatFilename, 
            inputFilename, outputFilename, boostFilename, 
            &asymmetry, &c_e_scale, &desiredSeed, &desiredCalo);

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
            }else if 	(line.find("pBinCBO") != std::string::npos) {
                pBinCBO = bool(atoi(line.substr(datIndex,100).c_str()));
            }else if 	(line.find("constraintau") != std::string::npos) {
                constraintau = bool(atoi(line.substr(datIndex,100).c_str()));
            }else if 	(line.find("includeMopTerm") != std::string::npos) {
                includeMopTerm = bool(atoi(line.substr(datIndex,100).c_str()));
            }else if 	(line.find("constrainMop") != std::string::npos) {
                constrainMop = bool(atoi(line.substr(datIndex,100).c_str()));
            }else if    (line.find("templatePath") != std::string::npos) {
                templatePath = line.substr(datIndex,100).c_str();
            }else if    (line.find("refFile") != std::string::npos) {
                refFile = line.substr(datIndex,100).c_str();
            }else{}
        }
    } 
    format.close();

    // open refFile (previous template iteration, or fullFits if doing templateFits0
    // save w_CBO per calo to be fixed in the template fit
    TFile *f_refFile = TFile::Open(refFile.c_str(), "READ");
    std::vector<double> ref_wCBOs(25, 0.0);
    std::vector<double> ref_Kys(25, 0.0);
    TTree *t_refFile = (TTree*)f_refFile->Get("fitresults");
    double ref_wCBO;
    double ref_Ky;
    int ref_calo;
    t_refFile->SetBranchAddress("w_CBO", &ref_wCBO);
    t_refFile->SetBranchAddress("Ky", &ref_Ky);
    t_refFile->SetBranchAddress("calo", &ref_calo);
    for (int refEntry=0; refEntry<t_refFile->GetEntries(); refEntry++){
        t_refFile->GetEntry(refEntry);
        ref_wCBOs[ref_calo] = ref_wCBO;
        ref_Kys[ref_calo] = ref_Ky;
    }

    //TO DO: stop making duplicates of caloWeights in each template file...or just make one file for both alpha/beta
    printf("loading CBO templates from templatePath = %s\n", templatePath.c_str());

    TFile *f_alpha_CBO = TFile::Open(Form("%s/templates_alpha_CBO.root", templatePath.c_str()), "READ");
    TFile *f_beta_CBO  = TFile::Open(Form("%s/templates_beta_CBO.root" , templatePath.c_str()), "READ");
    TFile *f_alpha_2CBO = TFile::Open(Form("%s/templates_alpha_2CBO.root", templatePath.c_str()), "READ");
    TFile *f_beta_2CBO  = TFile::Open(Form("%s/templates_beta_2CBO.root" , templatePath.c_str()), "READ");
    TFile *f_alpha_y = TFile::Open(Form("%s/templates_alpha_y.root", templatePath.c_str()), "READ");
    TFile *f_beta_y  = TFile::Open(Form("%s/templates_beta_y.root" , templatePath.c_str()), "READ");
    TFile *f_alpha_vw = TFile::Open(Form("%s/templates_alpha_vw.root", templatePath.c_str()), "READ");
    TFile *f_beta_vw  = TFile::Open(Form("%s/templates_beta_vw.root" , templatePath.c_str()), "READ");

    int caloNum;
    double caloWeight;

    TTree* treeOfcaloWeights = (TTree*)f_alpha_CBO->Get("caloWeights");
    treeOfcaloWeights->SetBranchAddress("calo", &caloNum);
    treeOfcaloWeights->SetBranchAddress("weight", &caloWeight);

    // if desiredCalo = 0 (i.e. calo sum) load weights as usual
    // otherwise set weight of desiredCalo==1.0 and ignore the rest
    if (desiredCalo==0){
        for (int entry=0; entry<treeOfcaloWeights->GetEntries(); entry++){
            treeOfcaloWeights->GetEntry(entry);
            caloWeights[caloNum-1] = caloWeight;
        }
    }else{
        caloWeights[desiredCalo-1] = 1.0;
    }

    for (int i=0; i<25; i++){
        TF1* temp_alpha_CBO = (TF1*)f_alpha_CBO->Get(Form("calo%i_alpha_CBO", desiredCalo));
        TF1* temp_beta_CBO  = (TF1*)f_beta_CBO ->Get(Form("calo%i_beta_CBO" , desiredCalo));
        TF1* temp_alpha_2CBO = (TF1*)f_alpha_2CBO->Get(Form("calo%i_alpha_2CBO", desiredCalo));
        TF1* temp_beta_2CBO  = (TF1*)f_beta_2CBO ->Get(Form("calo%i_beta_2CBO" , desiredCalo));
        TF1* temp_alpha_y = (TF1*)f_alpha_y->Get(Form("calo%i_alpha_y", desiredCalo));
        TF1* temp_beta_y  = (TF1*)f_beta_y ->Get(Form("calo%i_beta_y" , desiredCalo));
        TF1* temp_alpha_vw = (TF1*)f_alpha_vw->Get(Form("calo%i_alpha_vw", desiredCalo));
        TF1* temp_beta_vw  = (TF1*)f_beta_vw ->Get(Form("calo%i_beta_vw" , desiredCalo));

        alpha_CBO_TF1.push_back(*temp_alpha_CBO);
        beta_CBO_TF1 .push_back(*temp_beta_CBO); 
        alpha_2CBO_TF1.push_back(*temp_alpha_2CBO);
        beta_2CBO_TF1 .push_back(*temp_beta_2CBO); 
        alpha_y_TF1.push_back(*temp_alpha_y);
        beta_y_TF1 .push_back(*temp_beta_y); 
        alpha_vw_TF1.push_back(*temp_alpha_vw);
        beta_vw_TF1 .push_back(*temp_beta_vw); 
    }

    f_alpha_CBO->Close();
    f_beta_CBO->Close();  
    f_alpha_2CBO->Close();
    f_beta_2CBO->Close(); 
    f_alpha_y->Close();  
    f_beta_y->Close();   
    f_alpha_vw->Close();
    f_beta_vw->Close(); 

    printf("Parameters after format file: pBinTau %d pBinPhi %d pBinCBO %d constraintau %d includeMopTerm %d constrainMop %d c_e_scale %f \n", pBinTau, pBinPhi, pBinCBO, constraintau, includeMopTerm, constrainMop, c_e_scale);
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

    for(int i=fitrangelow; i<nBins; i++){
        double tt = (i-0.5)*0.1492;
        double mdl = exp(-tt/(2.1969811*frlifetime));
        double xact = 0;
        for(int j=0; j<boostgamma->GetNbinsX(); j++){
            xact += boostgamma->GetBinContent(j)/boostgamma->Integral()*exp(-tt/(2.1969811*boostgamma->GetBinCenter(j)));
        }
        expcorr[i] = xact/mdl;
    }

    //*************
    // get momentum-binned weights
    //*************

    TH1F* transform_dp_p0 = (TH1F*) boostdata->Get("transform_dp_p0");
    TH1F* transform_gamma = (TH1F*) boostdata->Get("transform_gamma");
    TH1F* transform_x     = (TH1F*) boostdata->Get("transform_x");

    int pBins = int(transform_dp_p0->GetNbinsX());
    std::pair <double,double> physicalRange (-0.0056799,0.0056085);  // in dp/p0

    // first isolate the physical range of the spectra
    // get the integral and bin numbers 
    std::vector<double> pWeights(pBins, 0.0); 

    for (int i=0; i<pBins; i++){
        double lowEdge = transform_dp_p0->GetBinLowEdge(i);
        if (physicalRange.first <= lowEdge) {
            binRange.first = i;
            break;
        }
    }
    for (int i=pBins-1; i>=0; i--){
        double upEdge = transform_dp_p0->GetBinLowEdge(i+1);
        if (physicalRange.second >= upEdge) {
            binRange.second = i;
            break;
        }
    }

    beta_0 = sqrt(1.0 - (1.0 / pow(29.3,2.0))); // this uses design gamma
    R_0    = 7112; // design R_0 in mm
    // now go through bins and store values for dp_p0, dPhi, gammas, tau, x
    for (int i=1; i<=pBins; i++){
        dp_p0 .push_back(transform_dp_p0       ->GetBinLowEdge(i));
        gammas.push_back(transform_gamma       ->GetBinLowEdge(i));
        dPhi  .push_back(-1.351 * transform_dp_p0 ->GetBinLowEdge(i)); // negative if fitting cos(wt - phi), pos if wt + phi // run6 13.51, run1 10.0
        x_e   .push_back(transform_x->GetBinLowEdge(i)); //convert from radial offset to radial position
        c_e   .push_back( c_e_scale * (2.0 * ((beta_0*beta_0)/R_0) * n 
                    * transform_x->GetBinLowEdge(i) 
                    * transform_dp_p0->GetBinLowEdge(i)) );
        printf("bin %i c_e %f\n", i, c_e.at(i-1));
    }

    printf("Physical range of dp/p0 is (%f,%f) from bins (%i,%i).\n",
            transform_dp_p0->GetBinLowEdge(binRange.first), transform_dp_p0->GetBinLowEdge(binRange.second+1),
            binRange.first, binRange.second);

    double totalContent = transform_dp_p0->Integral(binRange.first, binRange.second);
    printf("total momentum bin content %f\n", totalContent);

    // now find weights in only the physical range
    totalContent = transform_dp_p0->Integral(binRange.first, binRange.second);
    double sumWeights = 0.0;
    for (int i=1; i<=pBins; i++){
        // find magicBin 
        double lowEdge = transform_dp_p0->GetBinLowEdge(i);
        double uppEdge = transform_dp_p0->GetBinLowEdge(i+1);

        if ( lowEdge <= 0.0  and 0.0 < uppEdge){
            magicBin = i; // assumes 1-based indexing
            // update wc_0 in magic bin which we use in p-dependent wCBO; includes conversion from ns to us 
            wc_0 = (1.0 / ((transform_gamma->GetBinLowEdge(magicBin)*2.1969811) / 1000.)); 
            printf("magic bin is %i with tau %f and wc_0 %f \n", magicBin, 2.1969811*gammas[magicBin], wc_0);
        }

        // find weights according to momentum
        double thisWeight = transform_dp_p0->GetBinContent(i) / totalContent;

        // is this bin in a physical momentum bin?
        bool isPhysical = (i>=binRange.first and i<= binRange.second);

        if (pBinTau or pBinPhi){
            binWeights.push_back( isPhysical? thisWeight : 0.0 );
        }else{ // if we're not binning, set weight=1.0 in magic bin and 0.0 otherwise
            binWeights.push_back( (i==magicBin) ? 1.0 : 0.0 );
        }
        sumWeights += binWeights[i];
        weightedAvgTau += binWeights[i]*transform_gamma->GetBinLowEdge(i)*2.1969811;
        printf("bin %i content %f momentum %f isPhysical? %d weight %f\n",i,transform_dp_p0->GetBinContent(i),dp_p0[i],isPhysical,binWeights.at(binWeights.size()-1));
    }

    printf("size of dp_p0 %lu size of binWeights %lu\n", dp_p0.size(), binWeights.size());

    boostdata->Close();
    //*************
    // end get momentum-binned weights + taus
    //*************

    TH1D* bestfit = new TH1D("bestfit", "bestfit", nBins, 0.0, 650.0644);
    TTree* fitresults = new TTree("fitresults", "fitresults");

    int sidx;
    double chisq, rchisq;
    double N0, tau, A0, phi, R;
    double T_CBO, w_CBO, A_CNx1, A_SNx1, C_CBO;
    double LM, A_CNx2, A_SNx2;
    double A_CAx1, A_SAx1, A_CNy1, A_SNy1;
    double Ky, Ty, A_CNy2, A_SNy2;
    double A_Cp, A_Sp;
    double A_SNxy, A_CNxy;
    double A_ct, T_xy, w_xy;
    double N0err, tauerr, A0err, phierr, Rerr;
    double T_CBOerr, w_CBOerr, A_CNx1err, A_SNx1err, C_CBOerr;
    double LMerr, A_CNx2err, A_SNx2err;
    double A_CAx1err, A_SAx1err, A_CNy1err, A_SNy1err;
    double Kyerr, Tyerr, A_CNy2err, A_SNy2err;
    double A_Cperr, A_Sperr;
    double A_SNxyerr, A_CNxyerr;
    double A_cterr, T_xyerr, w_xyerr;
    double ANx1, ANx1err, pNx1, pNx1err;
    double ANx2, ANx2err, pNx2, pNx2err;
    double AAx1, AAx1err, pAx1, pAx1err;
    double ANy1, ANy1err, pNy1, pNy1err;
    double ANy2, ANy2err, pNy2, pNy2err;
    double Ap, Aperr, pp, pperr;
    double ANxy, ANxyerr, pNxy, pNxyerr;
    double Gamma_mop, Gamma_moperr;

    fitresults->Branch("pBinTau", &pBinTau);
    fitresults->Branch("pBinPhi", &pBinPhi);
    fitresults->Branch("pBinCBO", &pBinCBO);
    fitresults->Branch("includeMopTerm", &includeMopTerm);
    fitresults->Branch("constrainMop", &constrainMop);
    fitresults->Branch("constrainTau", &constraintau);
    fitresults->Branch("calo", &desiredCalo);
    fitresults->Branch("seed", &sidx);
    fitresults->Branch("chisq", &chisq);
    fitresults->Branch("rchisq", &rchisq);
    fitresults->Branch("N0", &N0);
    fitresults->Branch("tau", &tau);
    fitresults->Branch("A0", &A0);
    fitresults->Branch("phi", &phi);
    fitresults->Branch("R", &R);
    fitresults->Branch("C_CBO", &C_CBO);
    fitresults->Branch("T_CBO", &T_CBO);
    fitresults->Branch("w_CBO", &w_CBO);
    fitresults->Branch("A_CNx1", &A_CNx1);
    fitresults->Branch("A_SNx1", &A_SNx1);
    fitresults->Branch("LM", &LM);
    fitresults->Branch("A_CNx2", &A_CNx2);
    fitresults->Branch("A_SNx2", &A_SNx2);
    fitresults->Branch("A_CAx1", &A_CAx1);
    fitresults->Branch("A_SAx1", &A_SAx1);
    fitresults->Branch("A_CNy1", &A_CNy1);
    fitresults->Branch("A_SNy1", &A_SNy1);
    fitresults->Branch("A_SNxy", &A_SNxy);
    fitresults->Branch("A_CNxy", &A_CNxy);
    fitresults->Branch("Ky", &Ky);
    fitresults->Branch("Ty", &Ty);
    fitresults->Branch("A_CNy2", &A_CNy2);
    fitresults->Branch("A_SNy2", &A_SNy2);
    fitresults->Branch("A_Cp", &A_Cp);
    fitresults->Branch("A_Sp", &A_Sp);
    fitresults->Branch("A_ct", &A_ct);
    fitresults->Branch("T_xy", &T_xy);
    fitresults->Branch("w_xy", &w_xy);

    fitresults->Branch("N0err", &N0err);
    fitresults->Branch("tauerr", &tauerr);
    fitresults->Branch("A0err", &A0err);
    fitresults->Branch("phierr", &phierr);
    fitresults->Branch("Rerr", &Rerr);
    fitresults->Branch("C_CBOerr", &C_CBOerr);
    fitresults->Branch("T_CBOerr", &T_CBOerr);
    fitresults->Branch("w_CBOerr", &w_CBOerr);
    fitresults->Branch("A_CNx1err", &A_CNx1err);
    fitresults->Branch("A_SNx1err", &A_SNx1err);
    fitresults->Branch("LMerr", &LMerr);
    fitresults->Branch("A_CNx2err", &A_CNx2err);
    fitresults->Branch("A_SNx2err", &A_SNx2err);
    fitresults->Branch("A_CAx1err", &A_CAx1err);
    fitresults->Branch("A_SAx1err", &A_SAx1err);
    fitresults->Branch("A_CNy1err", &A_CNy1err);
    fitresults->Branch("A_SNy1err", &A_SNy1err);
    fitresults->Branch("Kyerr", &Kyerr);
    fitresults->Branch("Tyerr", &Tyerr);
    fitresults->Branch("A_CNy2err", &A_CNy2err);
    fitresults->Branch("A_SNy2err", &A_SNy2err);
    fitresults->Branch("A_Cperr", &A_Cperr);
    fitresults->Branch("A_Sperr", &A_Sperr);
    fitresults->Branch("A_SNxyerr", &A_SNxyerr);
    fitresults->Branch("A_CNxyerr", &A_CNxyerr);
    fitresults->Branch("A_cterr", &A_cterr);
    fitresults->Branch("T_xyerr", &T_xyerr);
    fitresults->Branch("w_xyerr", &w_xyerr);
    fitresults->Branch("ANx1", &ANx1);
    fitresults->Branch("Anx1err", &ANx1err);
    fitresults->Branch("pNx1", &pNx1);
    fitresults->Branch("pNx1err", &pNx1err);
    fitresults->Branch("ANx2", &ANx2);
    fitresults->Branch("ANx2err", &ANx2err);
    fitresults->Branch("pNx2", &pNx2);
    fitresults->Branch("pNx2err", &pNx2err);
    fitresults->Branch("AAx1", &AAx1);
    fitresults->Branch("AAx1err", &AAx1err);
    fitresults->Branch("pAx1", &pAx1);
    fitresults->Branch("pAx1err", &pAx1err);
    fitresults->Branch("ANy1", &ANy1);
    fitresults->Branch("ANy1err", &ANy1err);
    fitresults->Branch("pNy1", &pNy1);
    fitresults->Branch("pNy1err", &pNy1err);
    fitresults->Branch("ANy2", &ANy2);
    fitresults->Branch("ANy2err", &ANy2err);
    fitresults->Branch("pNy2", &pNy2);
    fitresults->Branch("pNy2err", &pNy2err);
    fitresults->Branch("Ap", &Ap);
    fitresults->Branch("Aperr", &Aperr);
    fitresults->Branch("pp", &pp);
    fitresults->Branch("pperr", &pperr);
    fitresults->Branch("ANxy", &ANxy);
    fitresults->Branch("ANxyerr", &ANxyerr);
    fitresults->Branch("pNxy", &pNxy);
    fitresults->Branch("pNxyerr", &pNxyerr);
    fitresults->Branch("Gamma_mop", &Gamma_mop);
    fitresults->Branch("Gamma_moperr", &Gamma_moperr);
    fitresults->Branch("c_e_scale"  , &c_e_scale);

    for(int seed=desiredSeed; seed<desiredSeed+1; seed++){
        sidx = seed;
        fprintf(stderr, "Reading data\n");
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
            if(corrtemp)
                varcorr->SetBinContent(j+1, corrtemp->GetBinContent(j+1));
        }
        double integral = 0;
        // go from bin 0 to start time bin, remember that muonloss goes from 0 - 650us
        for(int j=0; j<fitrangelow; j++){
            double ml = muonlosstemp->GetBinContent(j+1)/muonlosstemp->Integral();
            double wt = exp((j-0.5)*0.1492/64.44);
            integral += ml*wt;
        }
        for(int j=fitrangelow; j<nBins; j++){
            double ml = muonlosstemp->GetBinContent(j+1)/muonlosstemp->Integral();
            double wt = exp((j-0.5)*0.1492/64.44);
            integral += ml*wt;
            lambda->SetBinContent(j+1, integral);
            //lambdaArray[j+1] = integral;
        }
        input->Close();
        fprintf(stderr, "Done Reading data\n");

        // Now do a 28 parameter fit
        double par[29];
        double errorplus[29];
        double errorminus[29];

        TMinuit minimizer(29);
        minimizer.Command("SET PRINTOUT 3"); // change to level 1
        minimizer.Command("SET NOWARNINGS");
        minimizer.SetFCN(minuitFunction);
        //readMinuitCommands("minimizer.txt", &minimizer);
        //std::cout<<"title is: "<<minimizer.GetTitle()<<std::endl;
        minimizer.fGraphicsMode = false;
        minimizer.DefineParameter(0, "N0", 1.6*wiggle->GetBinContent(fitrangelow), 100000000, 0, 1000000000); // mistakenly put the adjustment here
        // fix tau when you have mop term unconstrained
        if(fixtau or includeMopTerm){
            minimizer.DefineParameter(1, "tau", 2.1969811*gammas[magicBin], 0, 60, 70); // constrains to tau in magic bin
            //minimizer.DefineParameter(1, "tau", weightedAvgTau, 0, 60, 70); // constrains to weighted-avg tau
        }
        else{
            minimizer.DefineParameter(1, "tau", 64.4, 0.1, 60, 70);
        }
        minimizer.DefineParameter(2, "A0", 0.37, 0.1, -1.0, 1.0);
        minimizer.DefineParameter(3, "phi", 4.208, 0.01, 2.0, 5.0);
        minimizer.DefineParameter(4, "R", -60.0, 1.0, -1000, 1000);

        minimizer.DefineParameter(5, "T_CBO", 100.0, 10.0, 1, 10000); 
        minimizer.DefineParameter(6, "w_CBO", ref_wCBOs[desiredCalo], 0.0, 2.0, 3.0);

        minimizer.DefineParameter(9, "LM", 0.001, 0.001, -0.1, 0.1);

        minimizer.DefineParameter(12,"A_CAx1", -0.2, 0.1, -1.0, 1.0);
        minimizer.DefineParameter(13,"A_SAx1", -0.3, 0.1, -1.0, 1.0);

        minimizer.DefineParameter(16,"Ky", ref_Kys[desiredCalo], 0.00, 0.0, 1.2);
        minimizer.DefineParameter(17,"Ty", 30, 10, 1.0, 10000.);

        minimizer.DefineParameter(20,"A_Cp", -0.2, 0.1, -1.0, 1.0);
        minimizer.DefineParameter(21,"A_Sp", -0.3, 0.1, -1.0, 1.0);

        minimizer.DefineParameter(22,"A_CNxy", -0.2, 0.1, -1.0, 1.0);
        minimizer.DefineParameter(23,"A_SNxy", -0.3, 0.1, -1.0, 1.0);
        minimizer.DefineParameter(25,"T_xy", 60, 10, 1, 1000);
        minimizer.DefineParameter(26,"w_xy", 11.9, 0.1, 11.5, 12.5);
        if (includeMopTerm){
            if (constrainMop){
                minimizer.DefineParameter(27,"Gamma_mop", 0.00002, 0.0000001, -0.001, 0.001);
            }else{
                minimizer.DefineParameter(27,"Gamma_mop", 0.00002, 0.0000001, -0.001, 0.001);
            }
        }else{
            minimizer.DefineParameter(27,"Gamma_mop"    , 0., 0., -   1.0,    1.0);
        }
        minimizer.DefineParameter(28, "C_CBO", 0.0, 0.0, 0, 0); // SET TO 0 TO REMOVE ENVELOPE

        printf("FITTING PHASE: 5-PARAM ONLY\n");
        minimizer.Command("SET PAR 8 0");
        minimizer.Command("SET PAR 9 0");
        minimizer.Command("SET PAR 10 0");
        minimizer.Command("SET PAR 11 0");
        minimizer.Command("SET PAR 12 0");
        minimizer.Command("SET PAR 13 0");
        minimizer.Command("SET PAR 14 0");
        minimizer.Command("SET PAR 15 0");
        minimizer.Command("SET PAR 16 0");
        minimizer.Command("SET PAR 17 0");
        minimizer.Command("SET PAR 18 0");
        minimizer.Command("SET PAR 19 0");
        minimizer.Command("SET PAR 20 0");
        minimizer.Command("SET PAR 21 0");
        minimizer.Command("SET PAR 22 0");
        minimizer.Command("SET PAR 23 0");
        minimizer.Command("SET PAR 24 0");
        minimizer.Command("SET PAR 25 0");
        minimizer.Command("SET PAR 26 0");
        minimizer.Command("SET PAR 27 0");
        minimizer.Command("SET PAR 28 0");
        minimizer.Command("SET PAR 29 0");
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
        minimizer.Command("FIX 17");
        minimizer.Command("FIX 18");
        minimizer.Command("FIX 19");
        minimizer.Command("FIX 20");
        minimizer.Command("FIX 21");
        minimizer.Command("FIX 22");
        minimizer.Command("FIX 23");
        minimizer.Command("FIX 24");
        minimizer.Command("FIX 25");
        minimizer.Command("FIX 26");
        minimizer.Command("FIX 27");
        minimizer.Command("FIX 28");
        minimizer.Command("FIX 29");
        minimizer.Migrad();

        if (pBinCBO){
            minimizer.DefineParameter(7, "A_CNx1", -0.0, 0.1, -1.0, 1.0);
            minimizer.DefineParameter(8, "A_SNx1", -0.0, 0.1, -1.0, 1.0);
            minimizer.DefineParameter(10,"A_CNx2", -0.0, 0.1, -1.0, 1.0);
            minimizer.DefineParameter(11,"A_SNx2", -0.0, 0.1, -1.0, 1.0);
            minimizer.DefineParameter(14,"A_CNy1", -0.0, 0.1, -1.0, 1.0);
            minimizer.DefineParameter(15,"A_SNy1", -0.0, 0.1, -1.0, 1.0);
            minimizer.DefineParameter(18,"A_CNy2", -0.0, 0.1, -1.0, 1.0);
            minimizer.DefineParameter(19,"A_SNy2", -0.0, 0.1, -1.0, 1.0);
            minimizer.DefineParameter(24,"A_ct", 0.0, 0.00, -0.1, 0.5);
        } else{
            minimizer.DefineParameter(7, "A_CNx1",  1.0, 0.01, -2.0, 2.0); // USE TO SCALE ALPHA/BETA
            minimizer.DefineParameter(8, "A_SNx1",  1.0, 0.0, -1.0, 1.0); // SUB WITH TF1
            minimizer.DefineParameter(10,"A_CNx2",  1.0, 0.01, -1.0, 1.0);
            minimizer.DefineParameter(11,"A_SNx2",  1.0, 0.0, -1.0, 1.0);
            minimizer.DefineParameter(14,"A_CNy1",  1.0, 0.01, -1.0, 1.0);
            minimizer.DefineParameter(15,"A_SNy1",  1.0, 0.0, -1.0, 1.0);
            minimizer.DefineParameter(18,"A_CNy2",  1.0, 0.01, -1.0, 1.0);
            minimizer.DefineParameter(19,"A_SNy2",  1.0, 0.0, -1.0, 1.0);
            minimizer.DefineParameter(24,"A_ct", 0.0, 0.0, -0.1, 0.5); // SET TO 0; SHOULD BE ABSORBED INTO ALPHA BETA
        }

        minimizer.DefineParameter(16,"Ky", ref_Kys[desiredCalo], 0.00, 0.9, 1.2);

        printf("FITTING PHASE: CBO ONLY\n");
        minimizer.Command("RES");
        minimizer.DefineParameter(7, "A_CNx1",  1.0, 0.01, -2.0, 2.0); // SUB WITH TF1
        minimizer.DefineParameter(8, "A_SNx1",  1.0, 0.0, -1.0, 1.0); // SUB WITH TF1
        minimizer.Command("FIX 1");
        minimizer.Command("FIX 2");
        minimizer.Command("FIX 3");
        minimizer.Command("FIX 4");
        minimizer.Command("FIX 5");
        minimizer.Command("FIX 10");
        minimizer.Command("FIX 13");
        minimizer.Command("FIX 14");
        minimizer.Command("FIX 15");
        minimizer.Command("FIX 16");
        minimizer.Command("FIX 17");
        minimizer.Command("FIX 18");
        minimizer.Command("FIX 19");
        minimizer.Command("FIX 20");
        minimizer.Command("FIX 21");
        minimizer.Command("FIX 22");
        minimizer.Command("FIX 23");
        minimizer.Command("FIX 24");
        minimizer.Command("FIX 25");
        minimizer.Command("FIX 26");
        minimizer.Command("FIX 27");
        minimizer.Migrad();
        minimizer.Command("RES");
        minimizer.Command("FIX 11");
        minimizer.Command("FIX 12");
        minimizer.Command("FIX 15");
        minimizer.Command("FIX 16");
        minimizer.Command("FIX 17");
        minimizer.Command("FIX 18");
        minimizer.Command("FIX 19");
        minimizer.Command("FIX 20");
        minimizer.Command("FIX 21");
        minimizer.Command("FIX 22");
        minimizer.Command("FIX 23");
        minimizer.Command("FIX 24");
        minimizer.Command("FIX 25");
        minimizer.Command("FIX 26");
        minimizer.Command("FIX 27");
        minimizer.Command("FIX 28");
        minimizer.Migrad();
        minimizer.Command("RES");
        //minimizer.Migrad();

        minimizer.Command("MIG 25000 0.01");
        //minimizer.Command("MINO 25000");

        TString varname[29];
        TString chnam;
        double val, err, xlolim, xuplim;
        int iuint;
        double fmin, fedm, errdef;
        int npari, nparx, istat;
        double eplus, eminus, eparab, globc;

        for(int k=0; k<29; k++){
            minimizer.mnpout(k, chnam, val, err, xlolim, xuplim, iuint);
            par[k] = val;
            varname[k] = chnam;
            minimizer.mnerrs(k, eplus, eminus, eparab, globc);
            errorplus[k] = eplus;
            errorminus[k] = eminus;
        }

        for(int k=0; k<29; k++){
            fprintf(stderr, "%s \t %12.5e + %9.3e %9.3e\n",
                    varname[k].Data(), par[k], errorplus[k], errorminus[k]);
        }
        minimizer.mnstat(fmin, fedm, errdef, npari, nparx, istat);
        noFreeParams = minimizer.GetNumFreePars();
        fprintf(stderr, "chisq/ndf : %f / %d\n", fmin, fitrangehigh-fitrangelow - noFreeParams);
        fprintf(stderr, "reduced chisq: %f +- %f\n", fmin/(fitrangehigh-fitrangelow - noFreeParams), sqrt(2./(fitrangehigh-fitrangelow-noFreeParams)));
        minimizer.Command("SHO COR");
        double covariance[28][28];
        minimizer.mnemat(&covariance[0][0],28);

        chisq  = fmin;
        rchisq = fmin/(fitrangehigh-fitrangelow-noFreeParams);
        N0 = par[0];
        tau = par[1];
        A0 = par[2];
        phi = par[3];
        R = par[4];
        C_CBO = par[28];
        T_CBO = par[5];
        w_CBO = par[6];
        A_CNx1 = par[7];
        A_SNx1 = par[8];
        LM = par[9];
        A_CNx2 = par[10];
        A_SNx2 = par[11];
        A_CAx1 = par[12];
        A_SAx1 = par[13];
        A_CNy1 = par[14];
        A_SNy1 = par[15];
        Ky = par[16];
        Ty = par[17];
        A_CNy2 = par[18];
        A_SNy2 = par[19];
        A_Cp = par[20];
        A_Sp = par[21];
        A_CNxy = par[22];
        A_SNxy = par[23];
        A_ct = par[24];
        T_xy = par[25];
        w_xy = par[26];
        Gamma_mop = par[27];

        N0err = sqrt(-errorplus[0]*errorminus[0]);
        tauerr = sqrt(-errorplus[1]*errorminus[1]);
        A0err = sqrt(-errorplus[2]*errorminus[2]);
        phierr = sqrt(-errorplus[3]*errorminus[3]);
        Rerr = sqrt(-errorplus[4]*errorminus[4]);
        C_CBOerr = sqrt(-errorplus[28]*errorminus[28]);
        T_CBOerr = sqrt(-errorplus[5]*errorminus[5]);
        w_CBOerr = sqrt(-errorplus[6]*errorminus[6]);
        LMerr = sqrt(-errorplus[9]*errorminus[9]);
        Kyerr = sqrt(-errorplus[16]*errorminus[16]);
        Tyerr = sqrt(-errorplus[17]*errorminus[17]);
        A_cterr = (errorplus[24]-errorminus[24])/2.;
        T_xyerr = (errorplus[25]-errorminus[25])/2.;
        w_xyerr = (errorplus[26]-errorminus[26])/2.;
        Gamma_moperr = (errorplus[27]-errorminus[27])/2.;

        A_CNx1err = sqrt(-errorplus[7]*errorminus[7]);
        A_SNx1err = sqrt(-errorplus[8]*errorminus[8]);

        ANx1 = sqrt(A_CNx1*A_CNx1 + A_SNx1*A_SNx1);
        double erracnx1 = (errorplus[7] - errorminus[7])/2.;
        double errasnx1 = (errorplus[8] - errorminus[8])/2.;
        ANx1err = sqrt(pow(erracnx1*A_CNx1, 2.) + pow(errasnx1*A_SNx1, 2.)
                + 2*A_CNx1*A_SNx1*erracnx1*errasnx1*covariance[7][8])/
            sqrt(A_CNx1*A_CNx1 + A_SNx1*A_SNx1);
        pNx1 = invert(A_CNx1, A_SNx1);
        pNx1err = sqrt(pow(A_CNx1*errasnx1, 2.) + pow(A_SNx1*erracnx1, 2.)
                + 2*A_CNx1*A_SNx1*erracnx1*errasnx1*covariance[7][8])/
            (A_CNx1*A_CNx1 + A_SNx1*A_SNx1);
        ANx2 = sqrt(A_CNx2*A_CNx2 + A_SNx2*A_SNx2);
        double erracnx2 = (errorplus[10] - errorminus[10])/2.;
        double errasnx2 = (errorplus[11] - errorminus[11])/2.;
        ANx2err = sqrt(pow(erracnx2*A_CNx2, 2.) + pow(errasnx2*A_SNx2, 2.)
                + 2*A_CNx2*A_SNx2*erracnx2*errasnx2*covariance[10][11])/
            sqrt(A_CNx2*A_CNx2 + A_SNx2*A_SNx2);
        pNx2 = invert(A_CNx2, A_SNx2);
        pNx2err = sqrt(pow(A_CNx2*errasnx2, 2.) + pow(A_SNx2*erracnx2, 2.)
                + 2*A_CNx2*A_SNx2*erracnx2*errasnx2*covariance[10][11])/
            (A_CNx2*A_CNx2 + A_SNx2*A_SNx2);

        AAx1 = sqrt(A_CAx1*A_CAx1 + A_SAx1*A_SAx1);
        double erracax1 = (errorplus[12] - errorminus[12])/2.;
        double errasax1 = (errorplus[13] - errorminus[13])/2.;
        AAx1err = sqrt(pow(erracax1*A_CAx1, 2.) + pow(errasax1*A_SAx1, 2.)
                + 2*A_CAx1*A_SAx1*erracax1*errasax1*covariance[12][13])/
            sqrt(A_CAx1*A_CAx1 + A_SAx1*A_SAx1);
        pAx1 = invert(A_CAx1, A_SAx1);
        pAx1err = sqrt(pow(A_CAx1*errasax1, 2.) + pow(A_SAx1*erracax1, 2.)
                + 2*A_CAx1*A_SAx1*erracax1*errasax1*covariance[12][13])/
            (A_CAx1*A_CAx1 + A_SAx1*A_SAx1);

        ANy1 = sqrt(A_CNy1*A_CNy1 + A_SNy1*A_SNy1);
        double erracny1 = (errorplus[14] - errorminus[14])/2.;
        double errasny1 = (errorplus[15] - errorminus[15])/2.;
        ANy1err = sqrt(pow(erracny1*A_CNy1, 2.) + pow(errasny1*A_SNy1, 2.)
                + 2*A_CNy1*A_SNy1*erracny1*errasny1*covariance[14][15])/
            sqrt(A_CNy1*A_CNy1 + A_SNy1*A_SNy1);
        pNy1 = invert(A_CNy1, A_SNy1);
        pNy1err = sqrt(pow(A_CNy1*errasny1, 2.) + pow(A_SNy1*erracny1, 2.)
                + 2*A_CNy1*A_SNy1*erracny1*errasny1*covariance[14][15])/
            (A_CNy1*A_CNy1 + A_SNy1*A_SNy1);

        ANy2 = sqrt(A_CNy2*A_CNy2 + A_SNy2*A_SNy2);
        double erracny2 = (errorplus[18] - errorminus[18])/2.;
        double errasny2 = (errorplus[19] - errorminus[19])/2.;
        ANy2err = sqrt(pow(erracny2*A_CNy2, 2.) + pow(errasny2*A_SNy2, 2.)
                + 2*A_CNy2*A_SNy2*erracny2*errasny2*covariance[18][19])/
            sqrt(A_CNy2*A_CNy2 + A_SNy2*A_SNy2);
        pNy2 = invert(A_CNy2, A_SNy2);
        pNy2err = sqrt(pow(A_CNy2*errasny2, 2.) + pow(A_SNy2*erracny2, 2.)
                + 2*A_CNy2*A_SNy2*erracny2*errasny2*covariance[18][19])/
            (A_CNy2*A_CNy2 + A_SNy2*A_SNy2);

        Ap = sqrt(A_Cp*A_Cp + A_Sp*A_Sp);
        double erracp = (errorplus[20] - errorminus[20])/2.;
        double errasp = (errorplus[21] - errorminus[21])/2.;
        Aperr = sqrt(pow(erracp*A_Cp, 2.) + pow(errasp*A_Sp, 2.)
                + 2*A_Cp*A_Sp*erracp*errasp*covariance[20][21])/
            sqrt(A_Cp*A_Cp + A_Sp*A_Sp);
        pp = invert(A_Cp, A_Sp);
        pperr = sqrt(pow(A_Cp*errasp, 2.) + pow(A_Sp*erracp, 2.)
                + 2*A_Cp*A_Sp*erracp*errasp*covariance[20][21])/
            (A_Cp*A_Cp + A_Sp*A_Sp);

        ANxy = sqrt(A_CNxy*A_CNxy + A_SNxy*A_SNxy);
        double erracxy = (errorplus[22] - errorminus[22])/2.;
        double errasxy = (errorplus[23] - errorminus[23])/2.;
        ANxyerr = sqrt(pow(erracxy*A_CNxy, 2.) + pow(errasxy*A_SNxy, 2.)
                + 2*A_CNxy*A_SNxy*erracxy*errasxy*covariance[22][23])/
            sqrt(A_CNxy*A_CNxy + A_SNxy*A_SNxy);
        pNxy = invert(A_CNxy, A_SNxy);
        pNxyerr = sqrt(pow(A_CNxy*errasxy, 2.) + pow(A_SNxy*erracxy, 2.)
                + 2*A_CNxy*A_SNxy*erracxy*errasxy*covariance[22][23])/
            (A_CNxy*A_CNxy + A_SNxy*A_SNxy);

        A_CNx2err = sqrt(-errorplus[10]*errorminus[10]);
        A_SNx2err = sqrt(-errorplus[11]*errorminus[11]);
        A_CAx1err = sqrt(-errorplus[12]*errorminus[12]);
        A_SAx1err = sqrt(-errorplus[13]*errorminus[13]);
        A_CNy1err = sqrt(-errorplus[14]*errorminus[14]);
        A_SNy1err = sqrt(-errorplus[15]*errorminus[15]);

        A_CNy2err = sqrt(-errorplus[18]*errorminus[18]);
        A_SNy2err = sqrt(-errorplus[19]*errorminus[19]);
        A_Cperr = sqrt(-errorplus[20]*errorminus[20]);
        A_Sperr = sqrt(-errorplus[21]*errorminus[21]);
        fitresults->Fill();


        double sum = 0;
        for(int i=fitrangelow; i<nBins; i++){

            double dim_i_n[1] = {double(i-0.5)};
            double dim_i_0[1] = {double(i    )};
            double dim_i_p[1] = {double(i+0.5)};

            bestfit->Fill((i-0.5)*0.1492 , calcnu(dim_i_0, par));
            double pucorr = calcnu(dim_i_0, par)*calcnu(dim_i_0, par)/(calcnu(dim_i_n, par) * calcnu(dim_i_p, par));
            sum += pucorr;
        }
        fprintf(stderr, "Mean corr %f\n", sum/nBins);

    }; // end loop over seeds

    TFile* output = new TFile(outputFilename, "RECREATE");
    output->cd();
    TH1D* residua = new TH1D("residua", "residua", nBins, 0.0, 650.0644);
    TTree* pulls = new TTree("pulls", "pulls");
    double pull;
    pulls->Branch("pull", &pull);
    fprintf(stderr, "Creating pulls tree..\n");
    for(int i=fitrangelow; i<=fitrangehigh; i++){
        double residual = wiggle->GetBinContent(i) - bestfit->GetBinContent(i);
        residua->SetBinContent(i, residual);
        pull = residual/sqrt(wiggleraw->GetBinContent(i) + varcorr->GetBinContent(i));
        pulls->Fill();
    }

    fprintf(stderr, "Writing data to file..\n");
    residua->Write();
    pulls->Write();
    fitresults->Write();  
    lambda->Write();
    bestfit->Write();
    wiggle->Write();
    fprintf(stderr, "Closing file..\n");
    output->Close();
}
