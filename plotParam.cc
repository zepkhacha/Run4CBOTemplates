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

void plotParam(
        std::string inputFilename, 
        std::string paramName, 
        std::string units, 
        double minY, 
        double maxY,
        std::string drawOption){

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kDarkRainBow, 0, 0.5);
    //gStyle->SetLineColor(kBlack);
    //gStyle->SetMarkerColor(kBlack);

    double start;
    int windowNo;
    int calo;
    double param;
    double err;

    TFile slidingFile (inputFilename.c_str());
    // add to plots for each requested run and calo
    TTree* slidingResults = (TTree*) slidingFile.Get("fitresults");

    slidingResults->SetBranchAddress("calo", &calo);
    slidingResults->SetBranchAddress("windowNo", &windowNo);
    slidingResults->SetBranchAddress(Form("%s",paramName.c_str()), &param);
    slidingResults->SetBranchAddress(Form("%s_err",paramName.c_str()), &err);
    slidingResults->SetBranchAddress("start", &start);

    std::vector<TGraphErrors> gSlidingParam    ;

    for (int i=0; i<25; i++){
        TGraphErrors temp;
        gSlidingParam.push_back(temp);
        gSlidingParam[i].SetMarkerStyle(7);
    }

    // get sliding window results
    int graphEntry = 0;
    for (int entry=0; entry<slidingResults->GetEntries(); entry++){
        slidingResults->GetEntry(entry);

        printf("calo %i entry %i param %f err %f\n", calo-1, entry ,param,err);
        //// skip if err=0
        //if (err==0.0){
        //    continue;
        //}

        //// skip entries with a ridiculous R -- fit didn't converge
        //if (slidingFit_R > 500 or slidingFit_R < -500){
        //    continue;
        //}

        double time = start ;

        gSlidingParam[calo].SetPoint(gSlidingParam[calo].GetN(), time, param );
        gSlidingParam[calo].SetPointError(gSlidingParam[calo].GetN()-1, 0, err );

        graphEntry++;
    }// end for loop over window fits

//    // fit first time to get w_0 and phi_0
//    for (int i=0; i<24; i++){
//
//        //TF1 fit ("fit", "[0]*x - [1]", gSlidingCboPA[i].GetPointX(0), gSlidingCboPA[i].GetPointX(gSlidingCboPA[i].GetN()-1));
//        TF1 fit ("fit", "[0]*x - [1]", 50, 100);
//        fit.SetParLimits(1,-M_PI,+M_PI);
//        printf("FIRST CALL to FIT()\n");
//        gSlidingCboPA[i].Fit("fit", "MW", "", 50, 100);
//        //printf("SECOND CALL to FIT()\n");
//        //gSlidingCboPA[i].Fit("fit", "M", "", 50, 100);
//        gLinearResidualFits.push_back(fit);
//
//        double w_0 = fit.GetParameter(0);
//        double w_0_err = fit.GetParError(0);
//        double phi_0 = fit.GetParameter(1);
//        double phi_0_err = fit.GetParError(1);
//        printf("CALO %i FIT RESULTS: w_0 %f phi_0 %f\n", i+1, w_0, phi_0);
//
//        // now draw residual for CBO
//        for (int j=0; j<gSlidingCboPA[i].GetN(); j++) {
//            gLinearResidual[i].SetPoint(j, gSlidingCboPA[i].GetPointX(j), gSlidingCboPA[i].GetPointY(j) - fit.Eval(gSlidingCboPA[i].GetPointX(j)));
//            gLinearResidual[i].SetPointError(j, 0, gSlidingCboPA[i].GetErrorY(j));
//        }
//
//    } // end loop over calos

    TCanvas c;
    std::string outputFilename = Form("plot_%s.pdf", paramName.c_str());
    c.Print(Form("%s[", outputFilename.c_str()));

    TMultiGraph mgSlidingParam;
    mgSlidingParam.SetTitle(
        Form("%s v time;Window Start [#mus]; %s", paramName.c_str(), units.c_str()));

    TLegend *lSlidingParam = new TLegend(0.7,0.6,0.9,0.9);

    lSlidingParam->SetTextSize(0.01);

    for (unsigned int i=0; i<gSlidingParam.size(); i++){
        //gSlidingParam[i].SetLineColor(i+1);
        //gSlidingParam[i].SetMarkerColor(i+1);
        //gSlidingParam[i].SetFillColor(i+1);
        mgSlidingParam.Add(&gSlidingParam[i]);
        lSlidingParam->AddEntry(&gSlidingParam[i], Form("calo%i", i+1));
    }

    // draw overlay
    c.Clear();
    mgSlidingParam.SetMaximum(maxY);
    mgSlidingParam.SetMinimum(minY);
    mgSlidingParam.Draw(drawOption.c_str());
    lSlidingParam->Draw();
    c.Print(outputFilename.c_str());

    // add individual calos
    for (unsigned int i=0; i<gSlidingParam.size(); i++){
        //gSlidingParam[i].SetLineColor(kBlack);
        //gSlidingParam[i].SetMarkerColor(kBlack);
        //gSlidingParam[i].SetFillColor(kGray);
        gSlidingParam[i].SetTitle(Form("calo%i",i+1));
        gSlidingParam[i].SetMaximum(maxY);
        gSlidingParam[i].SetMinimum(minY);
        gSlidingParam[i].Draw(drawOption.c_str());
        c.Print(outputFilename.c_str());
    }

    c.Print(Form("%s]", outputFilename.c_str()));

} // end main
