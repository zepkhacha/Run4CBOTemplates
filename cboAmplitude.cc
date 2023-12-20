#include <string>
#include <iostream>
#include <math.h>
#include <TFile.h>
#include <TStyle.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TColor.h>
#include <TAxis.h>

void printhelp(){
    printf(
           "Analysis Program\n"
           "Supply arguments as follows:\n"
           "-h        : Prints this help text\n"
           "-i [s]    : Input filename with sliding windows (required)\n"
           "-f [s]    : Input filename with full-fit result (required)\n"
           "-o [s]    : Output file label (required)\n"
           );
}

void parse_cmdline(int argc, char** argv, std::string& fullfitFilename, std::string& slidingFilename, std::string& rootOutputFilename) {
    const char* const opts = "hs:f:o:";
    bool done = false;
    while(!done){
        const char ch = getopt(argc, argv, opts);
        switch(ch){
            case -1: done = true; break;
            case 'h': printhelp(); exit(0);
            case 'f': fullfitFilename = optarg; break;
            case 's': slidingFilename = optarg; break;
            case 'o': rootOutputFilename = optarg; break;
            default: printhelp(); exit(1);
        }
    }
}

double radToDeg(double rad){
    return rad*( 360 / (2*M_PI) );
}

double radToNPi(double rad){
    return rad / M_PI;
}

double within2Pi(double init_phi){
    double phi = init_phi;
    while (phi > M_PI or phi < -M_PI){
        while (phi > M_PI){
            phi -= 2*M_PI;
        }
        while (phi < -M_PI){
            phi += 2*M_PI;
        }
    }
    return phi;
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

int main(int argc, char* argv[]){
    
    // parse command line arguments
    std::string fullFitFilename = "";
    std::string slidingFilename = "";
    std::string rootOutputFilename = "";
    parse_cmdline(argc, argv, fullFitFilename, slidingFilename, rootOutputFilename);
    
    // initializes canvases for each plot
    
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    
    TCanvas cSlidingAmplitude;
    TCanvas cSlidingFrequency;
    TCanvas cSlidingPhase;
    TCanvas cFullFitResidual;
    TCanvas cLinearResidual;
    TCanvas cSlidingVal;
    TCanvas cAmpvFreq;
    
    gStyle->SetPalette(kThermometer, 0, 0.5);
    gStyle->SetLineColor(kBlack);
    gStyle->SetMarkerColor(kBlack);
    
    bool axesInitialized = false;
    
    // add to plots for each requested run and calo
    
    TFile fullFitFile (fullFitFilename.c_str());
    TFile slidingFile (slidingFilename.c_str());
    
    TTree* fullFitResults = (TTree*) fullFitFile.Get("fitresults");
    TTree* slidingResults = (TTree*) slidingFile.Get("fitresults");
    
    // get wt-phi from fullfit
    double fullFit_R, fullFit_Rerr, fullFit_wCBO, fullFit_wCBOerr, fullFit_phi, fullFit_phierr;
    double fullFit_A_CNx1, fullFit_A_CNx1err, fullFit_A_SNx1, fullFit_A_SNx1err, fullFit_pNx1, fullFit_pNx1err;
    double fullFit_A_ct;
    
    double slidingFit_R, slidingFit_Rerr, slidingFit_wCBO, slidingFit_wCBOerr;
    double slidingFit_A_CBO, slidingFit_A_CBOerr, slidingFit_phiCBO, slidingFit_phiCBOerr;
    double slidingFit_wCBO_phiCBO_cov;
    
    double slidingFit_Ax2, slidingFit_Ay1, slidingFit_Ay2;
    double slidingFit_px2, slidingFit_py1, slidingFit_py2;
    
    double slidingFit_Ax2err, slidingFit_Ay1err, slidingFit_Ay2err;
    double slidingFit_px2err, slidingFit_py1err, slidingFit_py2err;
    
    double fitStart, fitStop;
    int windowNo;
    
    double fullchisq, slidingchisq;
    double slidingFit_w_y, slidingFit_w_vw;
    double slidingFit_w_yerr, slidingFit_w_vwerr;
    
    fullFit_pNx1 = invert(fullFit_A_CNx1, fullFit_A_SNx1);
    
    // full fit result should only have one entry in TTree
    if (fullFitResults->GetEntries()!=1){
        std::cerr<<"Full-fit file should have only 1 entry in 'fitresults' TTree.";
    }
    
    slidingResults->SetBranchAddress("chisq", &slidingchisq);
    fullFitResults->SetBranchAddress("chisq", &fullchisq);
    
    slidingResults->SetBranchAddress("windowNo", &windowNo);
    fullFitResults->SetBranchAddress("A_ct", &fullFit_A_ct);
    fullFitResults->SetBranchAddress("R", &fullFit_R);
    slidingResults->SetBranchAddress("R", &slidingFit_R);
    fullFitResults->SetBranchAddress("Rerr", &fullFit_Rerr);
    slidingResults->SetBranchAddress("Rerr", &slidingFit_Rerr);
    fullFitResults->SetBranchAddress("w_CBO", &fullFit_wCBO);
    slidingResults->SetBranchAddress("w_CBO", &slidingFit_wCBO);
    fullFitResults->SetBranchAddress("w_CBOerr", &fullFit_wCBOerr);
    slidingResults->SetBranchAddress("w_CBOerr", &slidingFit_wCBOerr);
    fullFitResults->SetBranchAddress("pNx1" , &fullFit_phi);
    slidingResults->SetBranchAddress("phi_CBO", &slidingFit_phiCBO);
    slidingResults->SetBranchAddress("phi_CBOerr", &slidingFit_phiCBOerr);
    slidingResults->SetBranchAddress("A_CBO", &slidingFit_A_CBO);
    slidingResults->SetBranchAddress("A_CBOerr", &slidingFit_A_CBOerr);
    
    slidingResults->SetBranchAddress("ANx2", &slidingFit_Ax2);
    slidingResults->SetBranchAddress("pNx2", &slidingFit_px2);
    
    slidingResults->SetBranchAddress("ANy1", &slidingFit_Ay1);
    slidingResults->SetBranchAddress("pNy1", &slidingFit_py1);
    slidingResults->SetBranchAddress("w_y", &slidingFit_w_y);
    
    slidingResults->SetBranchAddress("ANy2", &slidingFit_Ay2);
    slidingResults->SetBranchAddress("pNy2", &slidingFit_py2);
    slidingResults->SetBranchAddress("w_vw", &slidingFit_w_vw);
    
    slidingResults->SetBranchAddress("ANx2err", &slidingFit_Ax2err);
    slidingResults->SetBranchAddress("pNx2err", &slidingFit_px2err);
    
    slidingResults->SetBranchAddress("ANy1err", &slidingFit_Ay1err);
    slidingResults->SetBranchAddress("pNy1err", &slidingFit_py1err);
    slidingResults->SetBranchAddress("w_yerr", &slidingFit_w_yerr);
    
    slidingResults->SetBranchAddress("ANy2err", &slidingFit_Ay2err);
    slidingResults->SetBranchAddress("pNy2err", &slidingFit_py2err);
    slidingResults->SetBranchAddress("w_vwerr", &slidingFit_w_vwerr);
    
    fullFitResults->SetBranchAddress("pNx1err" , &fullFit_phierr);
    fullFitResults->SetBranchAddress("A_CNx1", &fullFit_A_CNx1);
    fullFitResults->SetBranchAddress("A_SNx1", &fullFit_A_SNx1);
    fullFitResults->SetBranchAddress("A_CNx1err", &fullFit_A_CNx1err);
    fullFitResults->SetBranchAddress("A_SNx1err", &fullFit_A_SNx1err);
    
    slidingResults->SetBranchAddress("wCBO_phiCBO_cov", &slidingFit_wCBO_phiCBO_cov);
    
    slidingResults->SetBranchAddress("fitStart", &fitStart);
    slidingResults->SetBranchAddress("fitStop" , &fitStop );
    
    TGraphErrors gFullFitResidual;
    TGraphErrors gLinearResidual;
    TGraphErrors gSlidingVal;
    TGraphErrors gSlidingAmplitude;
    TGraphErrors gSlidingFrequency;
    TGraphErrors gSlidingPhase;
    TGraphErrors gax2;
    TGraphErrors gpx2;
    TGraphErrors gay1;
    TGraphErrors gpy1;
    TGraphErrors gpw1;
    TGraphErrors gay2;
    TGraphErrors gpy2;
    TGraphErrors gpw2;
    
    TGraphErrors gphaseAdvance_2cbo;
    TGraphErrors gphaseAdvance_vw;

    TGraphErrors gphaseAdvanceR_2cbo;
    TGraphErrors gphaseAdvanceR_vw;
    
    TGraph gAmpvFres;
    
    gFullFitResidual.SetMarkerStyle(7);
    gLinearResidual.SetMarkerStyle(7);
    gSlidingVal.SetMarkerStyle(7);
    gSlidingAmplitude.SetMarkerStyle(7);
    gSlidingFrequency.SetMarkerStyle(7);
    gSlidingPhase.SetMarkerStyle(7);
    gAmpvFres.SetMarkerStyle(7);
    gax2.SetMarkerStyle(7);
    gpx2.SetMarkerStyle(7);
    gay1.SetMarkerStyle(7);
    gpy1.SetMarkerStyle(7);
    gpw1.SetMarkerStyle(7);
    gay2.SetMarkerStyle(7);
    gpy2.SetMarkerStyle(7);
    gpw2.SetMarkerStyle(7);
    gphaseAdvance_2cbo.SetMarkerStyle(7);
    gphaseAdvance_vw.SetMarkerStyle(7);

    gphaseAdvanceR_2cbo.SetMarkerStyle(7);
    gphaseAdvanceR_vw.SetMarkerStyle(7);
    
    gax2.SetMaximum(0.001);
    gax2.SetMinimum(-0.001);
    //gay1.SetMaximum(0.005);
    //gay2.SetMaximum(0.005);
    //gax2.SetMinimum(-0.005);
    //gay1.SetMinimum(-0.005);
    //gay2.SetMinimum(-0.005);
    gpy1.SetMaximum(M_PI);
    gpy1.SetMinimum(-M_PI);
    gpy2.SetMaximum(M_PI);
    gpy2.SetMinimum(-M_PI);
    gpx2.SetMaximum(M_PI);
    gpx2.SetMinimum(-M_PI);

    gpw2.SetMaximum(14.6);
    gpw2.SetMinimum(14.2);
    
    gSlidingAmplitude.SetMaximum(20);
    gSlidingAmplitude.SetMinimum(-1);
    
    //gLinearResidual.SetMinimum(-1.0);
    //gLinearResidual.SetMaximum(+1.5);
    
    gFullFitResidual.SetTitle("(#omega_{CBO}t - #phi_{CBO})_{window} - (#omega_{CBO}t - #phi_{CBO})_{full};Window Start [#mus];CBO Phase Residual [rad]");
    gSlidingVal.SetTitle("CBO Phase Advance from Sliding Windows;Window Start [us];(#omega_{CBO}t - #phi_{CBO}) [rad]");
    gLinearResidual.SetTitle("(#omega_{CBO}t - #phi_{CBO})_{window} - (#omega_{0}t - #phi_{0})_{late};Window Start [us];CBO Phase Residual [rad]");
    
    gSlidingAmplitude.SetTitle("A_CBO v time;Window Start [#mus]; A_{CBO} [10^{-3}]");
    gSlidingFrequency.SetTitle("w_CBO v time;Window Start [#mus]; #omega_{CBO} [MHz]");
    gSlidingPhase.SetTitle("p_CBO v time;Window Start [#mus]; #phi_{CBO} [rad]");
    gAmpvFres.SetTitle("w_CBO vs A_CBO; w_CBO [rad/s]; A_CBO [10^{-3}]");
    
    gax2.SetTitle("Ax2; Window Start [#mus]; Ax2 ");
    gpx2.SetTitle("px2; Window Start [#mus]; px2 [rad]");
    gay1.SetTitle("Ay1; Window Start [#mus]; Ay1 ");
    gpy1.SetTitle("py1; Window Start [#mus]; py1 [rad]");
    gay2.SetTitle("Ay2; Window Start [#mus]; Ay2 ");
    gpy2.SetTitle("py2; Window Start [#mus]; py2 [rad]");
    gpw1.SetTitle("w_y;  Window Start [#mus]; w_y [rad/s]");
    gpw2.SetTitle("w_vw; Window Start [#mus]; w_vw [rad/s]");
    
    gphaseAdvance_2cbo.SetTitle("2CBO Phase Advance; Window Start [#mus]; 2CBO Phase Advance [rad]");
    gphaseAdvance_vw.SetTitle("vw Phase Advance; Window Start [#mus]; vw Phase Advance [rad]");
    gphaseAdvanceR_2cbo.SetTitle("2CBO Phase Advance Residual; Window Start [#mus]; 2CBO Residual [rad]");
    gphaseAdvanceR_vw.SetTitle("vw Phase Advance Residual; Window Start [#mus]; vw Residual [rad]");
    
    // get fullFitResults
    fullFitResults->GetEntry(0);
    printf("Full fit results: \nR: %f Rerr: %f wCBO: %f phi: %f \n", fullFit_R, fullFit_Rerr, fullFit_wCBO, (fullFit_phi));
    
    double previousResidual = 0.0;
    printf("getting window results...\n");
    
    // get sliding window results
    int graphEntry = 0;
    for (int entry=0; entry<slidingResults->GetEntries(); entry++){
        printf("entry %i...\n", entry);
        slidingResults->GetEntry(entry);
        
        // skip entries with a ridiculous R -- fit didn't converge
        if (slidingFit_R > 500 or slidingFit_R < -500){
            continue;
        }
        
        // take every x-th window for now
        if (windowNo % 20 != 0){
            continue;
        }
        
        //double time = fitStart + 0.5 * (fitStop - fitStart);
        double time = fitStart ;
        
        slidingFit_phiCBO = within2Pi(slidingFit_phiCBO);
        slidingFit_px2    = within2Pi(slidingFit_px2);
        slidingFit_py1    = within2Pi(slidingFit_py1);
        slidingFit_py2    = within2Pi(slidingFit_py2);
        
        double slidingVal = (slidingFit_wCBO * time - slidingFit_phiCBO);
        //double slidingValError = std::sqrt(std::pow(slidingFit_wCBOerr * time, 2) + std::pow(slidingFit_phiCBOerr, 2) - 2 * slidingFit_wCBO_phiCBO_cov * time);
        double slidingValError = std::sqrt( std::pow(slidingFit_wCBOerr * time, 2)+std::pow(slidingFit_phiCBOerr, 2));
        
        double fullFitVal = (fullFit_wCBO * (1.0 + fullFit_A_ct*exp(-time/24.4)/time) * time - fullFit_phi);
        double residual = slidingVal - fullFitVal;
        
        if (previousResidual == 0.0) {
            previousResidual = residual;
        }
        
        double offset = std::round((residual - previousResidual) / (2*M_PI)) * (2*M_PI);
        slidingVal -= offset;
        residual -= offset;
        
        // ad-hoc shift for transition in Run-5 data
        while (residual > M_PI or residual < -M_PI) {
            while (residual > M_PI){
                slidingVal -= 2*M_PI;
                residual -= 2*M_PI;
            }
            while (residual < -M_PI){
                slidingVal += 2*M_PI;
                residual += 2*M_PI;
            }
        }
        
        previousResidual = residual;
        
        gFullFitResidual.SetPoint(graphEntry, time, residual);
        gFullFitResidual.SetPointError(graphEntry, 0, slidingValError);
        
        gSlidingVal.SetPoint(graphEntry, time, slidingVal);
        gSlidingVal.SetPointError(graphEntry, 0, slidingValError);
        
        gSlidingAmplitude.SetPoint(graphEntry, time, slidingFit_A_CBO * 1E3);
        gSlidingAmplitude.SetPointError(graphEntry, 0, slidingFit_A_CBOerr * 1E3);
        
        gSlidingFrequency.SetPoint(graphEntry, time, slidingFit_wCBO );
        gSlidingFrequency.SetPointError(graphEntry, 0, slidingFit_wCBOerr / (2*M_PI));
        
        gSlidingPhase.SetPoint(graphEntry, time, slidingFit_phiCBO);
        gSlidingPhase.SetPointError(graphEntry, 0, slidingFit_phiCBOerr);
        
        gAmpvFres.SetPoint(graphEntry, slidingFit_wCBO, slidingFit_A_CBO*1E3);
        
        gax2.SetPoint(graphEntry, time, abs(slidingFit_Ax2));
        gax2.SetPointError(graphEntry, 0, slidingFit_Ax2err);
        gpx2.SetPoint(graphEntry, time, slidingFit_px2);
        gpx2.SetPointError(graphEntry, 0, slidingFit_px2err);
        gay1.SetPoint(graphEntry, time, abs(slidingFit_Ay1));
        gay1.SetPointError(graphEntry, 0, slidingFit_Ay1err);
        gpy1.SetPoint(graphEntry, time, slidingFit_py1);
        gpy1.SetPointError(graphEntry, 0, slidingFit_py1err);
        gay2.SetPoint(graphEntry, time, abs(slidingFit_Ay2));
        gay2.SetPointError(graphEntry, 0, slidingFit_Ay2err);
        gpy2.SetPoint(graphEntry, time, slidingFit_py2);
        gpy2.SetPointError(graphEntry, 0, slidingFit_py2err);
        gpw1.SetPoint(graphEntry, time, slidingFit_w_y);
        gpw1.SetPointError(graphEntry, 0, slidingFit_w_yerr);
        gpw2.SetPoint(graphEntry, time, slidingFit_w_vw);
        gpw2.SetPointError(graphEntry, 0, slidingFit_w_vwerr);
        
        gphaseAdvance_2cbo.SetPoint(graphEntry, time, (2*slidingFit_wCBO*time - slidingFit_phiCBO));
        gphaseAdvance_vw.SetPoint(graphEntry, time, (slidingFit_w_vw*time - slidingFit_py2));
        
        graphEntry++;
    }// end for loop over window fits
    
    printf("doing linear fit...\n");
    
    // fit first time to get w_0 and phi_0
    TF1 fit ("fit", "[0]*x - [1]", gSlidingVal.GetPointX(0), gSlidingVal.GetPointX(gSlidingVal.GetN()-1));
    gSlidingVal.Fit("fit", "", "", 50, 200);
    fit.SetRange(50, 200);
    
    double w_0 = fit.GetParameter(0);
    double w_0_err = fit.GetParError(0);
    double phi_0 = fit.GetParameter(1);
    double phi_0_err = fit.GetParError(1);
    printf("FIT RESULTS: w_0 %f phi_0 %f\n", w_0, phi_0);

    // fit for vw 
    TF1 fitVW ("fitVW", "[0]*x - [1]", gphaseAdvance_vw.GetPointX(0), gphaseAdvance_vw.GetPointX(gphaseAdvance_vw.GetN()-1));
    gphaseAdvance_vw.Fit("fitVW", "", "", 50, 100);
    fitVW.SetRange(50,100);

    double lin_wVW = fitVW.GetParameter(0);
    double linwVW_err = fitVW.GetParError(0);
    double linphi_VW = fitVW.GetParameter(1);
    double linphi_VW_err = fitVW.GetParError(1);
    printf("FIT RESULTS: w_0 %f phi_0 %f\n", lin_wVW, linphi_VW);

    // fit for 2CBO
    TF1 fit2CBO ("fit2CBO", "[0]*x - [1]", gphaseAdvance_2cbo.GetPointX(0), gphaseAdvance_2cbo.GetPointX(gphaseAdvance_2cbo.GetN()-1));
    gphaseAdvance_2cbo.Fit("fit2CBO", "", "", 50, 180);
    fit2CBO.SetRange(50,180);

    double lin_w2CBO = fit2CBO.GetParameter(0);
    double linw2CBO_err = fit2CBO.GetParError(0);
    double linphi_2CBO = fit2CBO.GetParameter(1);
    double linphi_2CBO_err = fit2CBO.GetParError(1);
    printf("FIT RESULTS: w_0 %f phi_0 %f\n", lin_w2CBO, linphi_2CBO);

    // now draw residual for CBO
    for (int i=0; i<gSlidingVal.GetN(); i++) {
        gLinearResidual.SetPoint(i, gSlidingVal.GetPointX(i), gSlidingVal.GetPointY(i) - fit.Eval(gSlidingVal.GetPointX(i)));
        gLinearResidual.SetPointError(i, 0, gSlidingVal.GetErrorY(i));
    }
    
    // fit a*exp(-(x-t)/T) + b*t + c
    TF1 fitLinearCBOResidual ("fitCBORes", "[0]*exp(-(x-[1])/[2]) +[3]*x + [4]",
    gLinearResidual.GetPointX(0),
    gLinearResidual.GetPointX(gLinearResidual.GetN()-1));
    fitLinearCBOResidual.SetParameter(0,15);
    fitLinearCBOResidual.SetParameter(1,20);
    fitLinearCBOResidual.SetParameter(2,7);
    fitLinearCBOResidual.SetParLimits(3,-0.0001,0.0001);
    gLinearResidual.Fit("fitCBORes", "", "", 35, 200);
    fitLinearCBOResidual.SetRange(30,200);
    

    // now draw residual for VW
    for (int i=0; i<gphaseAdvance_vw.GetN(); i++) {
        gphaseAdvanceR_vw.SetPoint(i, gphaseAdvance_vw.GetPointX(i), 
        within2Pi(gphaseAdvance_vw.GetPointY(i) - fitVW.Eval(gphaseAdvance_vw.GetPointX(i))) );
        gphaseAdvanceR_vw.SetPointError(i, 0, gphaseAdvance_vw.GetErrorY(i));
    }

    // now draw residual for 2CBO
    for (int i=0; i<gphaseAdvance_2cbo.GetN(); i++) {
        gphaseAdvanceR_2cbo.SetPoint(i, gphaseAdvance_2cbo.GetPointX(i), 
        within2Pi(gphaseAdvance_2cbo.GetPointY(i) - fit2CBO.Eval(gphaseAdvance_2cbo.GetPointX(i))) );
        gphaseAdvanceR_2cbo.SetPointError(i, 0, gphaseAdvance_2cbo.GetErrorY(i));
    }
    
    
    // only add option 'a' for first plot, and only use palette for individual calos
    //std::string drawOption = Form("%sPE%s", axesInitialized ? "" : "A", (calo != 0 || runs.size() > 1) ? " PMC PLC" : "");
    std::string drawOption = "APE";
    axesInitialized = true;
    
    // initialize outpout files
    //TFile *f = TFile::Open(rootOutputFilename.c_str(),"RECREATE");
    
    TCanvas c;
    std::string outputFilename = Form("%s.pdf", rootOutputFilename.c_str());
    c.Print(Form("%s[", outputFilename.c_str()));
    
    gSlidingAmplitude.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gSlidingFrequency.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gSlidingPhase.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gFullFitResidual.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gSlidingVal.Draw(drawOption.c_str());
    fit.Draw("SAME");
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    //TLine line (0, 0, 240, 0);
    //line.SetLineStyle(kDotted);
    //line.Draw();
    //c.Print(outputFilename.c_str());
    
    gLinearResidual.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gAmpvFres.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gax2.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gpx2.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gay1.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gpy1.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gay2.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gpy2.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gpw1.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gpw2.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gphaseAdvance_2cbo.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    gphaseAdvance_vw.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());

    gphaseAdvanceR_vw.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());

    gphaseAdvanceR_2cbo.Draw(drawOption.c_str());
    c.Draw("Y+");
    c.Print(outputFilename.c_str());
    
    c.Print(Form("%s]", outputFilename.c_str()));
    
    //f->cd();
    //
    //gSlidingAmplitude.SetTitle("gSlidingAmplitude");
    //gSlidingAmplitude.SetName("gSlidingAmplitude");
    //gSlidingAmplitude.Write();
    //gSlidingFrequency.SetTitle("gSlidingFrequency");
    //gSlidingFrequency.SetName("gSlidingFrequency");
    //gSlidingFrequency.Write();
    //gSlidingPhase.SetTitle("gSlidingPhase");
    //gSlidingPhase.SetName("gSlidingPhase");
    //gSlidingPhase.Write();
    //
    //fullFitFile.Close();
    //slidingFile.Close();
    //f->Close();
    
} // end main
