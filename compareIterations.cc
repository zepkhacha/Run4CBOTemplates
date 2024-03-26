void compareIterations(){
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kDarkRainBow, 0, 0.5);

    std::vector<std::string> filenames = {
        "fullFits/noRF_vw_pm_cbo/sBin_constraintOn_cE0.0_seed0_noRF.root",
        "templateFits0/sBin_constraintOn_cE0.0_seed0_noRF.root",
        "templateFits1/sBin_constraintOn_cE0.0_seed0_noRF.root",
        "templateFits2/sBin_constraintOn_cE0.0_seed0_noRF.root",
        "templateFits2_floatScale/sBin_constraintOn_cE0.0_seed0_noRF.root",
        "templateFits3/sBin_constraintOn_cE0.0_seed0_noRF.root",
        "templateFits3_floatBoth/sBin_constraintOn_cE0.0_seed0_noRF.root",
        "templateFits4/sBin_constraintOn_cE0.0_seed0_noRF.root",
        "templateFits5/sBin_constraintOn_cE0.0_seed0_noRF.root"
    };

    std::vector<std::vector<double>> vals_chisq ;
    std::vector<std::vector<double>> vals_rchisq;
    std::vector<std::vector<double>> vals_R     ;
    std::vector<std::vector<double>> vals_LM    ;

    std::vector<TGraphErrors> g_chisq ;
    std::vector<TGraphErrors> g_rchisq;
    std::vector<TGraphErrors> g_R     ;
    std::vector<TGraphErrors> g_LM    ;

    // param list in order { R, chisq, rcshiq, LM }
    for (int j = 0; j < 25; j++){
        std::vector<double> temp(filenames.size(), 0.0);
        vals_R     .push_back(temp);
        vals_chisq .push_back(temp);
        vals_rchisq.push_back(temp);
        vals_LM    .push_back(temp); 

        TGraphErrors g_temp;
        g_R.push_back(g_temp);
        g_chisq.push_back(g_temp);
        g_rchisq.push_back(g_temp);

        g_R[j].SetLineWidth( j==0 ? 4 : 1);
        g_chisq[j].SetLineWidth( j==0 ? 4 : 1);
        g_rchisq[j].SetLineWidth( j==0 ? 4 : 1);
    }

    for (int i = 0; i < filenames.size(); i++){

        TFile *f = TFile::Open(filenames[i].c_str(),"READ");
        TTree *t = (TTree*)f->Get("fitresults");

        int calo;
        double R, R_err;
        double chisq, rchisq;
        double LM, LM_err;

//        if (i==0){
            t->SetBranchAddress("calo", &calo);
            t->SetBranchAddress("R", &R);
            t->SetBranchAddress("Rerr", &R_err);
            t->SetBranchAddress("chisq", &chisq);
            t->SetBranchAddress("rchisq", &rchisq);
            t->SetBranchAddress("LM", &LM);
            t->SetBranchAddress("LMerr", &LM_err);
//        }else{
//            t->SetBranchAddress("calo", &calo);
//            t->SetBranchAddress("R", &R);
//            t->SetBranchAddress("R_err", &R_err);
//            t->SetBranchAddress("chisq", &chisq);
//            t->SetBranchAddress("rchisq", &rchisq);
//            t->SetBranchAddress("LM", &LM);
//            t->SetBranchAddress("LM_err", &LM_err);
//        }

        for (int entry=0; entry<t->GetEntries(); entry++){
            t->GetEntry(entry);

            vals_R     [calo][i] = R;
            vals_chisq [calo][i] = chisq;
            vals_rchisq[calo][i] = rchisq;
            vals_LM    [calo][i] = LM;

            g_R[calo].SetPoint(i, double(i), R);
            g_R[calo].SetPointError(i, 0.0, R_err);

            g_chisq[calo].SetPoint(i, double(i), chisq);
            g_rchisq[calo].SetPoint(i, double(i), rchisq);

            printf("calo %i iteration %i R %f\n", calo, i, R);
        }

        f->Close();
    }

    TCanvas *c = new TCanvas();

    TMultiGraph *mg_R      = new TMultiGraph();
    TMultiGraph *mg_chisq  = new TMultiGraph();
    TMultiGraph *mg_rchisq = new TMultiGraph();
    TMultiGraph *mg_LM     = new TMultiGraph();

    TLegend *leg_R = new TLegend(0.7,0.6,0.9,0.9);
    leg_R->SetTextSize(0.01);

    g_R[0].SetLineWidth(8);

    for (int i = 0; i < 25; i++){

        mg_R->Add(&g_R[i]);
        leg_R->AddEntry(Form("calo %i",i));

        mg_chisq->Add(&g_chisq[i]);
        mg_rchisq->Add(&g_rchisq[i]);

    }

    //mg_R->SetMaximum(-59.0);
    //mg_R->SetMaximum(-62.0);
    mg_R->SetTitle("R vs Iteration; Iteration; R[ppm]");
    mg_R->Draw("ALE PFC PLC");
    //leg_R->Draw();
    c->Print("compareIterations_R.pdf");

    mg_chisq->SetTitle("chisq vs Iteration; Iteration; chisq");
    mg_chisq->Draw("ALE PFC PLC");
    //leg_R->Draw();
    c->Print("compareIterations_chisq.pdf");

    mg_rchisq->SetTitle("Reduced chisq vs Iteration; Iteration; Reduced chisq");
    mg_rchisq->Draw("ALE PFC PLC");
    //leg_R->Draw();
    c->Print("compareIterations_rchisq.pdf");

}
