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

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kThermometer, 0, 0.5);
    gStyle->SetLineColor(kBlack);
    gStyle->SetMarkerColor(kBlack);
    double slidingFit_R, slidingFit_Rerr, slidingFit_wCBO, slidingFit_wCBOerr;
    double slidingFit_A_CBO, slidingFit_A_CBOerr, slidingFit_phiCBO, slidingFit_phiCBOerr;
    double slidingFit_wCBO_phiCBO_cov;
    double slidingFit_Ax2, slidingFit_Ay1, slidingFit_Ay2;
    double slidingFit_px2, slidingFit_py1, slidingFit_py2;
    double slidingFit_Ax2err, slidingFit_Ay1err, slidingFit_Ay2err;
    double slidingFit_px2err, slidingFit_py1err, slidingFit_py2err;
    double fitStart, fitStop;
    int windowNo;
    int caloNum;
    double slidingchisq;
    double slidingFit_w_y, slidingFit_w_vw;
    double slidingFit_w_yerr, slidingFit_w_vwerr;

    TFile slidingFile ("/gm2data/zkhechad/cboTemplates/Run4CBOTemplates/slidingFits/run4_windowFits.root");
    // add to plots for each requested run and calo
    TTree* slidingResults = (TTree*) slidingFile.Get("fitresults");


    slidingResults->SetBranchAddress("chisq", &slidingchisq);
    slidingResults->SetBranchAddress("caloNum", &caloNum);
    slidingResults->SetBranchAddress("windowNo", &windowNo);
    slidingResults->SetBranchAddress("R", &slidingFit_R);
    slidingResults->SetBranchAddress("Rerr", &slidingFit_Rerr);
    slidingResults->SetBranchAddress("w_CBO", &slidingFit_wCBO);
    slidingResults->SetBranchAddress("w_CBOerr", &slidingFit_wCBOerr);
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
    slidingResults->SetBranchAddress("wCBO_phiCBO_cov", &slidingFit_wCBO_phiCBO_cov);
    slidingResults->SetBranchAddress("fitStart", &fitStart);
    slidingResults->SetBranchAddress("fitStop" , &fitStop );

    std::vector<TGraphErrors> gLinearResidual      ;
    std::vector<TGraphErrors> gSlidingCboPA          ;
    std::vector<TGraphErrors> gSlidingAmplitude    ;
    std::vector<TGraphErrors> gSlidingFrequency    ;
    std::vector<TGraphErrors> gSlidingPhase        ;

    std::vector<TF1> gLinearResidualFits ;
    //    std::vector<TGraphErrors> gax2                 ;
    //    std::vector<TGraphErrors> gpx2                 ;
    //    std::vector<TGraphErrors> gay1                 ;
    //    std::vector<TGraphErrors> gpy1                 ;
    //    std::vector<TGraphErrors> gpw1                 ;
    //    std::vector<TGraphErrors> gay2                 ;
    //    std::vector<TGraphErrors> gpy2                 ;
    //    std::vector<TGraphErrors> gpw2                 ;
    //    std::vector<TGraphErrors> gphaseAdvance_2cbo   ;
    //    std::vector<TGraphErrors> gphaseAdvance_vw     ;
    //    std::vector<TGraphErrors> gphaseAdvanceR_2cbo  ;
    //    std::vector<TGraphErrors> gphaseAdvanceR_vw    ;
    //    std::vector<TGraph> gAmpvFres;

    for (int i=0; i<24; i++){
        TGraphErrors temp;
        gLinearResidual.push_back(temp);
        gSlidingCboPA.push_back(temp);
        gSlidingAmplitude.push_back(temp);
        gSlidingFrequency.push_back(temp);
        gSlidingPhase.push_back(temp);

        gLinearResidual[i].SetMarkerStyle(7);
        gSlidingCboPA[i].SetMarkerStyle(7);
        gSlidingAmplitude[i].SetMarkerStyle(7);
        gSlidingAmplitude[i].SetMarkerStyle(7);
        gSlidingFrequency[i].SetMarkerStyle(7);
        gSlidingPhase[i].SetMarkerStyle(7);

        gSlidingAmplitude[i].SetMinimum(0);

    }

    //    gAmpvFres.SetMarkerStyle(7);
    //    gax2.SetMarkerStyle(7);
    //    gpx2.SetMarkerStyle(7);
    //    gay1.SetMarkerStyle(7);
    //    gpy1.SetMarkerStyle(7);
    //    gpw1.SetMarkerStyle(7);
    //    gay2.SetMarkerStyle(7);
    //    gpy2.SetMarkerStyle(7);
    //    gpw2.SetMarkerStyle(7);
    //    gphaseAdvance_2cbo.SetMarkerStyle(7);
    //    gphaseAdvance_vw.SetMarkerStyle(7);
    //
    //    gphaseAdvanceR_2cbo.SetMarkerStyle(7);
    //    gphaseAdvanceR_vw.SetMarkerStyle(7);
    //    
    //    gax2.SetMaximum(0.001);
    //    gax2.SetMinimum(-0.01);
    //    //gay1.SetMaximum(0.005);
    //    //gay2.SetMaximum(0.005);
    //    //gax2.SetMinimum(-0.005);
    //    //gay1.SetMinimum(-0.005);
    //    //gay2.SetMinimum(-0.005);
    //    gpy1.SetMaximum(M_PI);
    //    gpy1.SetMinimum(-M_PI);
    //    gpy2.SetMaximum(M_PI);
    //    gpy2.SetMinimum(-M_PI);
    //    gpx2.SetMaximum(M_PI);
    //    gpx2.SetMinimum(-M_PI);
    //
    //    gpw2.SetMaximum(14.6);
    //    gpw2.SetMinimum(14.2);
    //    
    //    gSlidingAmplitude.SetMaximum(20);
    //    gSlidingAmplitude.SetMinimum(-1);
    //    
    //    gLinearResidual.SetMinimum(-1.0);
    //    gLinearResidual.SetMaximum(+1.0);
    //    
    //    gAmpvFres.SetTitle("w_CBO vs A_CBO; w_CBO [rad/s]; A_CBO [10^{-3}]");
    //    
    //    gax2.SetTitle("Ax2; Window Start [#mus]; Ax2 ");
    //    gpx2.SetTitle("px2; Window Start [#mus]; px2 [rad]");
    //    gay1.SetTitle("Ay1; Window Start [#mus]; Ay1 ");
    //    gpy1.SetTitle("py1; Window Start [#mus]; py1 [rad]");
    //    gay2.SetTitle("Ay2; Window Start [#mus]; Ay2 ");
    //    gpy2.SetTitle("py2; Window Start [#mus]; py2 [rad]");
    //    gpw1.SetTitle("w_y;  Window Start [#mus]; w_y [rad/s]");
    //    gpw2.SetTitle("w_vw; Window Start [#mus]; w_vw [rad/s]");
    //    
    //    gphaseAdvance_2cbo.SetTitle("2CBO Phase Advance; Window Start [#mus]; 2CBO Phase Advance [rad]");
    //    gphaseAdvance_vw.SetTitle("vw Phase Advance; Window Start [#mus]; vw Phase Advance [rad]");
    //    gphaseAdvanceR_2cbo.SetTitle("2CBO Phase Advance Residual; Window Start [#mus]; 2CBO Residual [rad]");
    //    gphaseAdvanceR_vw.SetTitle("vw Phase Advance Residual; Window Start [#mus]; vw Residual [rad]");
    //    
    //    double previousResidual = 0.0;
    //    printf("getting window results...\n");

    // get sliding window results
    int graphEntry = 0;
    for (int entry=0; entry<slidingResults->GetEntries(); entry++){
        //     printf("entry %i...\n", entry);
        slidingResults->GetEntry(entry);

        // skip entries with a ridiculous R -- fit didn't converge
        if (slidingFit_R > 500 or slidingFit_R < -500){
            continue;
        }

        //     // take every x-th window for now
        //     if (windowNo % 20 != 0){
        //         continue;
        //     }

        //double time = fitStart + 0.5 * (fitStop - fitStart);
        double time = fitStart ;

        slidingFit_phiCBO = within2Pi(slidingFit_phiCBO);
        slidingFit_px2    = within2Pi(slidingFit_px2);
        slidingFit_py1    = within2Pi(slidingFit_py1);
        slidingFit_py2    = within2Pi(slidingFit_py2);

        double slidingVal = (slidingFit_wCBO * time - slidingFit_phiCBO);
        //double slidingValError = std::sqrt(std::pow(slidingFit_wCBOerr * time, 2) + std::pow(slidingFit_phiCBOerr, 2) - 2 * slidingFit_wCBO_phiCBO_cov * time);
        double slidingValError = std::sqrt( std::pow((slidingFit_wCBOerr * time)/(slidingFit_wCBO * time), 2)
                +std::pow( slidingFit_phiCBOerr/slidingFit_phiCBO, 2));

        gSlidingCboPA[caloNum-1].SetPoint(gSlidingCboPA[caloNum-1].GetN(), time, (slidingVal));
        gSlidingCboPA[caloNum-1].SetPointError(gSlidingCboPA[caloNum-1].GetN()-1, 0, slidingValError);

        gSlidingAmplitude[caloNum-1].SetPoint(gSlidingAmplitude[caloNum-1].GetN(), time, slidingFit_A_CBO * 1E3);
        gSlidingAmplitude[caloNum-1].SetPointError(gSlidingAmplitude[caloNum-1].GetN()-1, 0, slidingFit_A_CBOerr * 1E3);

        gSlidingFrequency[caloNum-1].SetPoint(gSlidingFrequency[caloNum-1].GetN(), time, slidingFit_wCBO );
        gSlidingFrequency[caloNum-1].SetPointError(gSlidingFrequency[caloNum-1].GetN()-1, 0, slidingFit_wCBOerr / (2*M_PI));

        gSlidingPhase[caloNum-1].SetPoint(gSlidingPhase[caloNum-1].GetN(), time, slidingFit_phiCBO);
        gSlidingPhase[caloNum-1].SetPointError(gSlidingPhase[caloNum-1].GetN()-1, 0, slidingFit_phiCBOerr);

        //        gAmpvFres.SetPoint(graphEntry, slidingFit_wCBO, slidingFit_A_CBO*1E3);
        //        
        //        gax2.SetPoint(graphEntry, time, abs(slidingFit_Ax2));
        //        gax2.SetPointError(graphEntry, 0, slidingFit_Ax2err);
        //        gpx2.SetPoint(graphEntry, time, slidingFit_px2);
        //        gpx2.SetPointError(graphEntry, 0, slidingFit_px2err);
        //        gay1.SetPoint(graphEntry, time, abs(slidingFit_Ay1));
        //        gay1.SetPointError(graphEntry, 0, slidingFit_Ay1err);
        //        gpy1.SetPoint(graphEntry, time, slidingFit_py1);
        //        gpy1.SetPointError(graphEntry, 0, slidingFit_py1err);
        //        gay2.SetPoint(graphEntry, time, abs(slidingFit_Ay2));
        //        gay2.SetPointError(graphEntry, 0, slidingFit_Ay2err);
        //        gpy2.SetPoint(graphEntry, time, slidingFit_py2);
        //        gpy2.SetPointError(graphEntry, 0, slidingFit_py2err);
        //        gpw1.SetPoint(graphEntry, time, slidingFit_w_y);
        //        gpw1.SetPointError(graphEntry, 0, slidingFit_w_yerr);
        //        gpw2.SetPoint(graphEntry, time, slidingFit_w_vw);
        //        gpw2.SetPointError(graphEntry, 0, slidingFit_w_vwerr);
        //        
        //        gphaseAdvance_2cbo.SetPoint(graphEntry, time, (2*slidingFit_wCBO*time - slidingFit_phiCBO));
        //        gphaseAdvance_vw.SetPoint(graphEntry, time, (slidingFit_w_vw*time - slidingFit_py2));

        graphEntry++;
    }// end for loop over window fits

    // fit first time to get w_0 and phi_0
    for (int i=0; i<24; i++){

        //TF1 fit ("fit", "[0]*x - [1]", gSlidingCboPA[i].GetPointX(0), gSlidingCboPA[i].GetPointX(gSlidingCboPA[i].GetN()-1));
        TF1 fit ("fit", "[0]*x - [1]", 50, 100);
        fit.SetParLimits(1,-M_PI,+M_PI);
        printf("FIRST CALL to FIT()\n");
        gSlidingCboPA[i].Fit("fit", "MW", "", 50, 100);
        //printf("SECOND CALL to FIT()\n");
        //gSlidingCboPA[i].Fit("fit", "M", "", 50, 100);
        gLinearResidualFits.push_back(fit);

        double w_0 = fit.GetParameter(0);
        double w_0_err = fit.GetParError(0);
        double phi_0 = fit.GetParameter(1);
        double phi_0_err = fit.GetParError(1);
        printf("CALO %i FIT RESULTS: w_0 %f phi_0 %f\n", i+1, w_0, phi_0);

        // now draw residual for CBO
        for (int j=0; j<gSlidingCboPA[i].GetN(); j++) {
            gLinearResidual[i].SetPoint(j, gSlidingCboPA[i].GetPointX(j), gSlidingCboPA[i].GetPointY(j) - fit.Eval(gSlidingCboPA[i].GetPointX(j)));
            gLinearResidual[i].SetPointError(j, 0, gSlidingCboPA[i].GetErrorY(j));
        }

    } // end loop over calos
    //
    //    // fit for vw 
    //    TF1 fitVW ("fitVW", "[0]*x - [1]", gphaseAdvance_vw.GetPointX(0), gphaseAdvance_vw.GetPointX(gphaseAdvance_vw.GetN()-1));
    //    gphaseAdvance_vw.Fit("fitVW", "", "", 50, 100);
    //    fitVW.SetRange(50,100);
    //
    //    double lin_wVW = fitVW.GetParameter(0);
    //    double linwVW_err = fitVW.GetParError(0);
    //    double linphi_VW = fitVW.GetParameter(1);
    //    double linphi_VW_err = fitVW.GetParError(1);
    //    printf("FIT RESULTS: w_0 %f phi_0 %f\n", lin_wVW, linphi_VW);
    //
    //    // fit for 2CBO
    //    TF1 fit2CBO ("fit2CBO", "[0]*x - [1]", gphaseAdvance_2cbo.GetPointX(0), gphaseAdvance_2cbo.GetPointX(gphaseAdvance_2cbo.GetN()-1));
    //    gphaseAdvance_2cbo.Fit("fit2CBO", "", "", 50, 180);
    //    fit2CBO.SetRange(50,180);
    //
    //    double lin_w2CBO = fit2CBO.GetParameter(0);
    //    double linw2CBO_err = fit2CBO.GetParError(0);
    //    double linphi_2CBO = fit2CBO.GetParameter(1);
    //    double linphi_2CBO_err = fit2CBO.GetParError(1);
    //    printf("FIT RESULTS: w_0 %f phi_0 %f\n", lin_w2CBO, linphi_2CBO);
    //    
    //    // fit a*exp(-(x-t)/T) + b*t + c
    //    TF1 fitLinearCBOResidual ("fitCBORes", "[0]*exp(-(x-[1])/[2]) +[3]*x + [4]",
    //    gLinearResidual.GetPointX(0),
    //    gLinearResidual.GetPointX(gLinearResidual.GetN()-1));
    //    fitLinearCBOResidual.SetParameter(0,15);
    //    fitLinearCBOResidual.SetParameter(1,20);
    //    fitLinearCBOResidual.SetParameter(2,7);
    //    fitLinearCBOResidual.SetParLimits(3,-0.0001,0.0001);
    //    gLinearResidual.Fit("fitCBORes", "", "", 35, 200);
    //    fitLinearCBOResidual.SetRange(30,200);
    //    
    //
    //    // now draw residual for VW
    //    for (int i=0; i<gphaseAdvance_vw.GetN(); i++) {
    //        gphaseAdvanceR_vw.SetPoint(i, gphaseAdvance_vw.GetPointX(i), 
    //        within2Pi(gphaseAdvance_vw.GetPointY(i) - fitVW.Eval(gphaseAdvance_vw.GetPointX(i))) );
    //        gphaseAdvanceR_vw.SetPointError(i, 0, gphaseAdvance_vw.GetErrorY(i));
    //    }
    //
    //    // now draw residual for 2CBO
    //    for (int i=0; i<gphaseAdvance_2cbo.GetN(); i++) {
    //        gphaseAdvanceR_2cbo.SetPoint(i, gphaseAdvance_2cbo.GetPointX(i), 
    //        within2Pi(gphaseAdvance_2cbo.GetPointY(i) - fit2CBO.Eval(gphaseAdvance_2cbo.GetPointX(i))) );
    //        gphaseAdvanceR_2cbo.SetPointError(i, 0, gphaseAdvance_2cbo.GetErrorY(i));
    //    }
    //    
    //    
    //    // only add option 'a' for first plot, and only use palette for individual calos
    //    //std::string drawOption = Form("%sPE%s", axesInitialized ? "" : "A", (calo != 0 || runs.size() > 1) ? " PMC PLC" : "");
    //    std::string drawOption = "APE";
    //    axesInitialized = true;
    //    
    //    // initialize outpout files
    //    //TFile *f = TFile::Open(rootOutputFilename.c_str(),"RECREATE");

    TCanvas c;
    std::string rootOutputFilename = "plots";
    std::string outputFilename = Form("%s.pdf", rootOutputFilename.c_str());
    c.Print(Form("%s[", outputFilename.c_str()));

    TMultiGraph mgSlidingAmplitude;
    TMultiGraph mgSlidingFrequency;
    TMultiGraph mgSlidingPhase;
    TMultiGraph mgLinearResidual;
    TMultiGraph mgSlidingCboPA;

    mgSlidingAmplitude.SetTitle("A_CBO v time;Window Start [#mus]; A_{CBO} [10^{-3}]");
    mgSlidingCboPA.SetTitle("CBO Phase Advance from Sliding Windows;Window Start [us];(#omega_{CBO}t - #phi_{CBO}) [rad]");
    mgLinearResidual.SetTitle("(#omega_{CBO}t - #phi_{CBO})_{window} - (#omega_{0}t - #phi_{0})_{late};Window Start [us];CBO Phase Residual [rad]");
    mgSlidingFrequency.SetTitle("w_CBO v time;Window Start [#mus]; #omega_{CBO} [MHz]");
    mgSlidingPhase.SetTitle("p_CBO v time;Window Start [#mus]; #phi_{CBO} [rad]");

    TLegend *lSlidingAmplitude = new TLegend(0.8,0.7,0.8,0.9);
    TLegend *lSlidingFrequency = new TLegend(0.1,0.1,0.3,0.5);
    TLegend *lSlidingPhase     = new TLegend(0.1,0.1,0.3,0.5);
    TLegend *lLinearResidual   = new TLegend(0.1,0.1,0.3,0.5);
    TLegend *lSlidingCboPA     = new TLegend(0.1,0.1,0.3,0.5);

    lSlidingAmplitude->SetTextSize(0.01);
    lSlidingFrequency->SetTextSize(0.01);
    lSlidingPhase    ->SetTextSize(0.01);
    lLinearResidual  ->SetTextSize(0.01);
    lSlidingCboPA    ->SetTextSize(0.01);

    for (unsigned int i=0; i<gSlidingAmplitude.size(); i++){
        gSlidingAmplitude[i].SetLineColor(i+1);
        gSlidingAmplitude[i].SetMarkerColor(i+1);
        mgSlidingAmplitude.Add(&gSlidingAmplitude[i]);
        lSlidingAmplitude->AddEntry(&gSlidingAmplitude[i], Form("calo%i", i+1));
    }
    c.Clear();
    mgSlidingAmplitude.SetMinimum(0);
    mgSlidingAmplitude.SetMaximum(12);
    mgSlidingAmplitude.Draw("APE");
    lSlidingAmplitude->Draw();
    c.Print(outputFilename.c_str());

    for (unsigned int i=0; i<gSlidingFrequency.size(); i++){
        gSlidingFrequency[i].SetLineColor(i+1);
        gSlidingFrequency[i].SetMarkerColor(i+1);
        mgSlidingFrequency.Add(&gSlidingFrequency[i]);
        lSlidingFrequency->AddEntry(&gSlidingFrequency[i], Form("calo%i", i+1));
    }
    c.Clear();
    mgSlidingFrequency.SetMaximum(2.46);
    mgSlidingFrequency.SetMinimum(2.2);
    mgSlidingFrequency.Draw("APE");
    //lSlidingFrequency->Draw();
    c.Print(outputFilename.c_str());

    for (unsigned int i=0; i<gSlidingPhase.size(); i++){
        gSlidingPhase[i].SetLineColor(i+1);
        gSlidingPhase[i].SetMarkerColor(i+1);
        mgSlidingPhase.Add(&gSlidingPhase[i]);
        lSlidingFrequency->AddEntry(&gSlidingPhase[i], Form("calo%i", i+1));
    }
    c.Clear();
    mgSlidingPhase.SetMaximum(2*M_PI);
    mgSlidingPhase.SetMinimum(-2*M_PI);
    mgSlidingPhase.Draw("APE");
    //lSlidingFrequency->Draw();
    c.Print(outputFilename.c_str());

    for (unsigned int i=0; i<gSlidingCboPA.size(); i++){
        gSlidingCboPA[i].SetLineColor(i+1);
        gSlidingCboPA[i].SetMarkerColor(i+1);
        gSlidingCboPA[i].SetTitle(Form("calo%i", i+1));
        mgSlidingCboPA.Add(&gSlidingCboPA[i]);
        lSlidingFrequency->AddEntry(&gSlidingCboPA[i], Form("calo%i", i+1));
        gSlidingCboPA[i].Draw("APE");
        c.Print(outputFilename.c_str());
    }
    c.Clear();
    mgSlidingCboPA.SetMaximum(1500);
    //mgSlidingCboPA.SetMinimum();
    mgSlidingCboPA.Draw("APE");
    //lSlidingFrequency->Draw();
    c.Print(outputFilename.c_str());

    for (unsigned int i=0; i<gLinearResidual.size(); i++){
        gLinearResidual[i].SetLineColor(1);
        gLinearResidual[i].SetMarkerColor(1);
        gLinearResidual[i].SetTitle(Form("Linear Residual Calo %i", i+1));
        gLinearResidual[i].SetName(Form("Linear Residual Calo %i", i+1));
        gLinearResidual[i].SetMaximum(1.5*M_PI);
        gLinearResidual[i].SetMinimum(-1.5*M_PI);
        c.Clear();
        gLinearResidual[i].Draw("APE");
        c.Print(outputFilename.c_str());

        gLinearResidual[i].SetLineColor(i+1);
        gLinearResidual[i].SetMarkerColor(i+1);
        mgLinearResidual.Add(&gLinearResidual[i]);
        lSlidingFrequency->AddEntry(&gLinearResidual[i], Form("calo%i", i+1));
    }
    c.Clear();
    mgLinearResidual.SetMaximum(1.5*M_PI);
    mgLinearResidual.SetMinimum(-1.5*M_PI);
    mgLinearResidual.Draw("AP");
    //lSlidingFrequency->Draw();
    c.Print(outputFilename.c_str());

    //    gSlidingCboPA.Draw(drawOption.c_str());
    //    fit.Draw("SAME");
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    //TLine line (0, 0, 240, 0);
    //    //line.SetLineStyle(kDotted);
    //    //line.Draw();
    //    //c.Print(outputFilename.c_str());
    //    
    //    gLinearResidual.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    gAmpvFres.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    gax2.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    gpx2.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    gay1.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    gpy1.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    gay2.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    gpy2.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    gpw1.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    gpw2.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    gphaseAdvance_2cbo.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    //    gphaseAdvance_vw.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //
    //    gphaseAdvanceR_vw.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //
    //    gphaseAdvanceR_2cbo.Draw(drawOption.c_str());
    //    c.Draw("Y+");
    //    c.Print(outputFilename.c_str());
    //    
    c.Print(Form("%s]", outputFilename.c_str()));

} // end main
