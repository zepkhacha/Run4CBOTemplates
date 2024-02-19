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

void makeTemplate(
        std::string inputFilename, 
        std::string weightRefFilename, // expects everything up to [calonum].root
        std::string paramName, 
        std::string units, 
        double minY, 
        double maxY,
        std::string drawOption){
    printf("weightRefFilename : %s\n", weightRefFilename.c_str());

    printf("Setting up style...\n");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kDarkRainBow, 0, 0.5);

    // declare a ttree to hold weights per calo
    TTree *caloWeights = new TTree("caloWeights", "caloWeights");
    int calo;
    double calo_N0, calo_N0_err, weight;

    caloWeights->Branch("calo", &calo);
    caloWeights->Branch("N0", &calo_N0);
    caloWeights->Branch("N0_err", &calo_N0_err);
    caloWeights->Branch("weight", &weight);

    //    // we also need N0 from calo sum result to calculate weight
    //    printf("Reading in N0 from calo0 nominal fit result...\n");
    //    TFile* calo0_weightRefFilename = TFile::Open(Form("%s%i.root", weightRefFilename.c_str(), 0), "READ");
    //    TTree* calo0_treeRef = (TTree*) calo0_weightRefFilename->Get("fitresults");
    //    double calo0_N0, calo0_N0_err;
    //
    //    calo0_treeRef->SetBranchAddress("N0", &calo0_N0);
    //    calo0_treeRef->SetBranchAddress("N0err", &calo0_N0_err);
    //    calo0_treeRef->GetEntry(0);
    //    calo0_weightRefFilename->Close();

    printf("Saving weights for each calo based on N0...\n");
    // pull the weight from the full-window fit per calo
    double sumN0= 0.0;
    for (int i=1; i<25; i++){
        TFile* calo_weightRefFilename = TFile::Open(Form("%s%i.root", weightRefFilename.c_str(), i), "READ");
        TTree* treeRef = (TTree*) calo_weightRefFilename->Get("fitresults");

        treeRef->SetBranchAddress("calo", &calo);
        treeRef->SetBranchAddress("N0", &calo_N0);
        treeRef->SetBranchAddress("N0err", &calo_N0_err);
        treeRef->GetEntry(0);
        sumN0 += calo_N0 ;
        calo_weightRefFilename->Close();
    }
    double sumWeight = 0.0;
    for (int i=1; i<25; i++){
        TFile* calo_weightRefFilename = TFile::Open(Form("%s%i.root", weightRefFilename.c_str(), i), "READ");
        TTree* treeRef = (TTree*) calo_weightRefFilename->Get("fitresults");

        treeRef->SetBranchAddress("calo", &calo);
        treeRef->SetBranchAddress("N0", &calo_N0);
        treeRef->SetBranchAddress("N0err", &calo_N0_err);
        treeRef->GetEntry(0);
        weight = calo_N0 / sumN0;
        sumWeight += weight;
        printf("i = %i calo = %i weight %f\n", i, calo, weight);
        caloWeights->Fill();
        calo_weightRefFilename->Close();
    }
    printf("sumWeight: %f\n", sumWeight);

    double fitStart;
    int windowNo;
    int caloNum;
    double param;
    double err;

    TFile slidingFile (inputFilename.c_str());
    // add to plots for each requested run and calo
    TTree* slidingResults = (TTree*) slidingFile.Get("fitresults");

    slidingResults->SetBranchAddress("calo", &caloNum);
    slidingResults->SetBranchAddress("windowNo", &windowNo);
    slidingResults->SetBranchAddress(Form("%s",paramName.c_str()), &param);
    slidingResults->SetBranchAddress(Form("%s_err",paramName.c_str()), &err);
    slidingResults->SetBranchAddress("start", &fitStart);

    std::vector<TGraphErrors> gSlidingParam    ;

    for (int i=0; i<24; i++){
        TGraphErrors temp;
        gSlidingParam.push_back(temp);
        gSlidingParam[i].SetMarkerStyle(7);
    }

    printf("Reading in sliding window values...\n");
    // get sliding window results
    int graphEntry = 0;
    for (int entry=0; entry<slidingResults->GetEntries(); entry++){
        slidingResults->GetEntry(entry);

        //printf("calo %i entry %i param %f err %f\n", caloNum-1, entry ,param,err);
        //// skip if err=0
        //if (err==0.0){
        //    continue;
        //}

        //// skip entries with a ridiculous R -- fit didn't converge
        //if (slidingFit_R > 500 or slidingFit_R < -500){
        //    continue;
        //}

        double time = fitStart ;

        gSlidingParam[caloNum-1].SetPoint(gSlidingParam[caloNum-1].GetN(), time, param );
        gSlidingParam[caloNum-1].SetPointError(gSlidingParam[caloNum-1].GetN()-1, 0, err );

        graphEntry++;
    }// end for loop over window fits

    printf("Producing root and pdf output files...\n");
    TCanvas c;
    std::string outputFilename = Form("templates_%s.pdf", paramName.c_str());
    std::string rootFilename = Form("templates_%s.root", paramName.c_str());
    c.Print(Form("%s[", outputFilename.c_str()));

    TFile *fOut = TFile::Open(rootFilename.c_str(),"RECREATE");
    caloWeights->Write();

    TMultiGraph mgSlidingParam;
    mgSlidingParam.SetTitle(
            Form("%s v time;Window Start [#mus]; %s", paramName.c_str(), units.c_str()));

    TLegend *lSlidingParam = new TLegend(0.7,0.6,0.9,0.9);

    lSlidingParam->SetTextSize(0.01);

    std::vector<TF1> fitFunction;
    for (unsigned int i=0; i<gSlidingParam.size(); i++){
        int color = gSlidingParam[i].GetLineColor();

        TF1* expModel;
        // if the parameter is CBO, then we can use double exp+c model
        if (paramName=="alpha_CBO" or paramName=="beta_CBO"){
            expModel = new TF1(Form("calo%i_%s", i+1, paramName.c_str()), 
                    "[0]*exp(-x/[1])+[2] + [3]*exp(-x/[4])+[5]", 
                    30.0, 400.0);
            expModel->SetLineColor(color);
            expModel->SetParameter(0, (gSlidingParam[i].GetPointY(0)>0 ? 0.001 : -0.001));
            expModel->SetParameter(1, 100.0);
            expModel->SetParameter(4, 10.0);
        }else{
            // create a template with exp+c model for y, vw, 2CBO
            expModel = new TF1(Form("calo%i_%s", i+1, paramName.c_str()), 
                    "[0]*exp(-x/[1])+[2]", 
                    30.0, 400.0);
            expModel->SetLineColor(color);
            expModel->SetParameter(0, (gSlidingParam[i].GetPointY(0)>0 ? 0.001 : -0.001));
            expModel->SetParameter(1, 100.0);
        }

        // now we can do a fit
        gSlidingParam[i].Fit(expModel, "ME", "", 30.0, 400.0);
        fitFunction.push_back(*expModel);
        // add to multigraph to view them overlaid
        mgSlidingParam.Add(&gSlidingParam[i]);
        lSlidingParam->AddEntry(&gSlidingParam[i], Form("calo%i", i+1));
        // write function to file
        expModel->Write();
    }

    fOut->Close();

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
        fitFunction[i].Draw("SAME");
        c.Print(outputFilename.c_str());
    }

    c.Print(Form("%s]", outputFilename.c_str()));

} // end main
