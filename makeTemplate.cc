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
    std::vector<double> weightsVector(24, 0.0);
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
        weightsVector[calo-1] = weight;
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

    for (int i=0; i<25; i++){
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

        gSlidingParam[caloNum].SetPoint(gSlidingParam[caloNum].GetN(), time, param );
        gSlidingParam[caloNum].SetPointError(gSlidingParam[caloNum].GetN()-1, 0, err );

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

    std::vector<TF1> fitFunction_short;
    std::vector<TF1> fitFunction_long;
    std::vector<TF1> fitFunction;

    for (unsigned int i=0; i<gSlidingParam.size(); i++){
        int color = gSlidingParam[i].GetLineColor();

        TF1* expShort = new TF1(Form("calo%i_%s", i, paramName.c_str()),
                "[0]*exp(-x/[1])+[2]",
                20.0,30.0);
        expShort->SetParameter(0,00.0);
        expShort->SetParameter(1,10.0);
        expShort->SetParameter(2,00.0);
        gSlidingParam[i].Fit(expShort, "MEN0", "", 20.0, 40.0);

        TF1* expLong = new TF1(Form("calo%i_%s", i, paramName.c_str()),
                "[0]*exp(-x/[1])+[2]",
                35.0,300.0);
        expLong->SetParameter(0,00.0);
        expLong->SetParameter(1,100.0);
        expLong->SetParameter(2,00.0);
        gSlidingParam[i].Fit(expLong, "MEN0", "", 30.0, 300.0);

        TF1* expModel = new TF1(Form("calo%i_%s", i, paramName.c_str()),
                "[0]*exp(-x/[1])+[2]+[3]*exp(-x/[4])",
                20.0,400.0);
        expModel->SetLineColor(kBlack);
        //expModel->SetParameter(0, (gSlidingParam[i].GetPointY(0)>0 ? 0.001 : -0.001));
        expModel->SetParameter(0, expLong ->GetParameter(0));
        expModel->SetParameter(1, expLong ->GetParameter(1));
        expModel->SetParameter(3, expShort->GetParameter(0));
        expModel->SetParameter(4, expShort->GetParameter(1));
        expModel->SetParLimits(0,-1.0,1.0);
        expModel->SetParLimits(3,-1.0,1.0);
        expModel->SetParLimits(1,5.0,1000.0);
        expModel->SetParLimits(4,0.0,10.0);
        gSlidingParam[i].Fit(expModel, "MEN0", "", 30.0, 400.0);
        double chisq_expModel = expModel->GetChisquare();
        printf("reduced chisq exp %f\n", chisq_expModel/(gSlidingParam[i].GetN()-3));

        // finally try the double exp model
        //TF1* exp2Model = new TF1(Form("calo%i_%s", i, paramName.c_str()), 
        //        "[0]*exp(-x/[1])+[2] + [3]*exp(-x/[4])+[5] + [6]*x+[7]", 
        //        30.0, 400.0);
        //exp2Model->SetLineColor(color);
        //exp2Model->SetParameter(0, (gSlidingParam[i].GetPointY(0)>0 ? 0.001 : -0.001));
        //exp2Model->SetParameter(1, 100.0);
        //exp2Model->SetParameter(4, 10.0);
        //expModel->SetParameter(6, 0.0);
        //expModel->SetParLimits(6,0.0,0.0);
        //gSlidingParam[i].Fit(exp2Model, "ME", "", 30.0, 400.0);
        //double chisq_exp2Model = exp2Model->GetChisquare();
        //printf("chisq exp2 %f\n", chisq_exp2Model);

        fitFunction.push_back(*expModel);
        fitFunction_short.push_back(*expShort);
        fitFunction_long.push_back(*expLong);
        expModel->Write();
        // add to multigraph to view them overlaid
        mgSlidingParam.Add(&gSlidingParam[i]);
        lSlidingParam->AddEntry(&gSlidingParam[i], Form("calo%i", i+1));

    }

    // now create a TF1 for the calo combination
    //TString combinationFunction = Form("%f*(%s) ",weightsVector[0],fitFunction[0].GetExpFormula().Data());
    std::string combinationFunction = Form("calo1_%s ", paramName.c_str());

    for (int i=1; i<24; i++){
        //combinationFunction += Form("+ %f*(%s) " ,weightsVector[i],fitFunction[i].GetExpFormula().Data());
        //combinationFunction += Form("+ %f*(calo%i_%s) " ,weightsVector[i],i+1,paramName.c_str());
        combinationFunction += Form("+ calo%i_%s " ,i+1,paramName.c_str());
    }

    TF1* caloCombo = new TF1(Form("caloCombo_%s", paramName.c_str()), combinationFunction.c_str());
    caloCombo->Write();

    fOut->Close();

    // draw overlay
    c.Clear();
    mgSlidingParam.SetTitle(Form("%s; time [#mus]; [arb]", paramName.c_str()));
    mgSlidingParam.SetMaximum(maxY);
    mgSlidingParam.SetMinimum(minY);
    mgSlidingParam.Draw("ALE PFC PLC");
    lSlidingParam->Draw();
    for (unsigned int i=0; i<gSlidingParam.size(); i++){
        fitFunction[i].SetLineColor(kBlack);
        fitFunction[i].SetLineStyle(8);
        fitFunction[i].Draw("SAME");
    }
    c.Print(outputFilename.c_str());

    // add individual calos
    for (unsigned int i=0; i<gSlidingParam.size(); i++){
        //gSlidingParam[i].SetLineColor(kBlack);
        //gSlidingParam[i].SetMarkerColor(kBlack);
        //gSlidingParam[i].SetFillColor(kGray);
        gSlidingParam[i].SetTitle(Form("calo%i; time [#mus]; [arb]",i));
        gSlidingParam[i].SetMaximum(maxY);
        gSlidingParam[i].SetMinimum(minY);
        gSlidingParam[i].Draw(drawOption.c_str());

        fitFunction[i].SetLineColor(kBlack);
        fitFunction_short[i].SetLineColor(kGreen);
        fitFunction_long[i].SetLineColor(kBlue);

        fitFunction[i].SetLineStyle(8);
        fitFunction_short[i].SetLineStyle(6);
        fitFunction_long[i].SetLineStyle(7);

        fitFunction[i].Draw("SAME");
        fitFunction_short[i].Draw("SAME");
        fitFunction_long[i].Draw("SAME");

        TLatex latex;
        latex.SetTextSize(0.03);
        double rchisq = fitFunction[i].GetChisquare()/(fitFunction[i].GetNDF());
        latex.DrawLatex(320, minY+0.35*(maxY-minY),Form("rchisq: %f", rchisq));
        latex.DrawLatex(320, minY+0.30*(maxY-minY),    ("p0*exp(-x/p1)+p2+p3*exp(-x/p4)"));
        latex.DrawLatex(320, minY+0.25*(maxY-minY),Form("p0: %f",fitFunction[i].GetParameter(0)));
        latex.DrawLatex(320, minY+0.20*(maxY-minY),Form("p1: %f",fitFunction[i].GetParameter(1)));
        latex.DrawLatex(320, minY+0.15*(maxY-minY),Form("p2: %f",fitFunction[i].GetParameter(2)));
        latex.DrawLatex(320, minY+0.10*(maxY-minY),Form("p3: %f",fitFunction[i].GetParameter(3)));
        latex.DrawLatex(320, minY+0.05*(maxY-minY),Form("p4: %f",fitFunction[i].GetParameter(4)));

        //printf("calo %i chisq %f rchisq %f\n", i, 
        //        fitFunction[i].GetChisquare(), 
        //        fitFunction[i].GetChisquare()/(fitFunction[i].GetNDF()));
        c.Print(outputFilename.c_str());
    }

    c.Print(Form("%s]", outputFilename.c_str()));

} // end main
