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

void checkCaloCombination(
        std::string directoryName,
        std::string frequency,
        std::string calo0_windowFitsFilename){

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kDarkRainBow, 0, 0.5);
    //gStyle->SetLineColor(kBlack);
    //gStyle->SetMarkerColor(kBlack);

    double start;
    int windowNo;
    int calo;
    std::vector<TF1> alpha;
    std::vector<TF1> beta ;
    std::vector<double> alpha_err;
    std::vector<double> beta_err;
    std::vector<double> weights(25,1.0);
    double weight;

    printf("Reading in calo0 windowFits.root...\n");
    TFile *calo0windowFits = TFile::Open(calo0_windowFitsFilename.c_str(),"READ");
    TTree *calo0fitresults = (TTree*)calo0windowFits->Get("fitresults");
    double calo0alpha, calo0beta;
    double calo0alphaerr, calo0betaerr;
    double calo0time;
    calo0fitresults->SetBranchAddress(Form("alpha_%s",frequency.c_str()),&calo0alpha);
    calo0fitresults->SetBranchAddress(Form("alpha_%s_err",frequency.c_str()),&calo0alphaerr);
    calo0fitresults->SetBranchAddress(Form("beta_%s" ,frequency.c_str()),&calo0beta );
    calo0fitresults->SetBranchAddress(Form("beta_%s_err" ,frequency.c_str()),&calo0betaerr );
    calo0fitresults->SetBranchAddress("start", &calo0time);
    TGraphErrors *g0alphaDat = new TGraphErrors();
    TGraphErrors *g0betaDat  = new TGraphErrors();

    g0alphaDat->SetLineColor(kBlue);
    g0betaDat ->SetLineColor(kBlue);

    for (int i=0; i<calo0fitresults->GetEntries(); i++){
        calo0fitresults->GetEntry(i);
        g0alphaDat->SetPoint(i, calo0time, calo0alpha);
        g0betaDat ->SetPoint(i, calo0time, calo0beta );
        g0alphaDat->SetPointError(i, 0.0, calo0alphaerr);
        g0betaDat ->SetPointError(i, 0.0, calo0betaerr );
    }

    printf("Reading in template files...\n");
    TFile *alphaFile = TFile::Open(Form("%s/templates_alpha_%s.root", directoryName.c_str(), frequency.c_str()), "READ");
    TFile *betaFile  = TFile::Open(Form("%s/templates_beta_%s.root" , directoryName.c_str(), frequency.c_str()), "READ");

    printf("...storing templates into vector<TF1> objects\n");
    // store calo templates in vector before doing combination
    for (int i=0; i<25; i++){
        TF1* alpha_temp = (TF1*)alphaFile->Get(Form("calo%i_alpha_%s", i, frequency.c_str()));
        TF1* beta_temp  = (TF1*)betaFile ->Get(Form("calo%i_beta_%s" , i, frequency.c_str()));
        alpha.push_back(*alpha_temp);
        beta .push_back(*beta_temp);
    }

    TTree *caloWeights = (TTree*)alphaFile->Get("caloWeights");
    caloWeights->SetBranchAddress("calo", &calo);
    caloWeights->SetBranchAddress("weight", &weight);

    printf("...storing calo weights into vector<double>\n");
    // fill weights into vector
    for (int i=0; i<24; i++){
        caloWeights->GetEntry(i);
        weights[calo] = weight;
    }

    TGraphErrors *g_alpha_0 = new TGraphErrors();
    TGraphErrors *g_alpha_Combo = new TGraphErrors();
    TGraphErrors *g_beta_0 = new TGraphErrors();
    TGraphErrors *g_beta_Combo = new TGraphErrors();
    TGraphErrors *g_A_0 = new TGraphErrors();
    TGraphErrors *g_A_Combo = new TGraphErrors();
    TGraphErrors *g_phi_0 = new TGraphErrors();
    TGraphErrors *g_phi_Combo = new TGraphErrors();

    printf("...iterating through time to produce graphs\n");
    for (int i=0; i<3800; i++){
        double time = (190+i)*0.1492;
        // get comparison curves from the calo-sum result
        g_alpha_0->SetPoint( i, time, alpha[0].Eval(time));
        g_beta_0 ->SetPoint( i, time, beta [0].Eval(time));
        g_A_0    ->SetPoint( i, time, sqrt( pow(alpha[0].Eval(time),2.0) + pow(beta[0].Eval(time),2.0) ));
        g_phi_0  ->SetPoint( i, time, invert( alpha[0].Eval(time),beta[0].Eval(time) ));
        // and now do the calo-combination result
        double alphaCombo = 0.0;
        double betaCombo = 0.0;
        double ACombo = 0.0;
        double phiCombo = 0.0;

        for (int j=1; j<=24; j++){
            alphaCombo += weights[j] * alpha[j].Eval(time);
            betaCombo  += weights[j] * beta [j].Eval(time);
            ACombo += weights[j] * sqrt( pow(alpha[j].Eval(time),2.0) + pow(beta[j].Eval(time),2.0) );
            phiCombo += weights[j] * invert( alpha[j].Eval(time), beta[j].Eval(time) );
        }

        g_alpha_Combo->SetPoint( g_alpha_Combo->GetN(), time, alphaCombo );
        g_beta_Combo ->SetPoint( g_beta_Combo ->GetN(), time, betaCombo );
        g_A_Combo->SetPoint( g_A_Combo->GetN(), time, ACombo );
        g_phi_Combo->SetPoint( g_phi_Combo->GetN(), time, phiCombo );
    }

    TCanvas *c = new TCanvas();
    c->Print(Form("comparison_%s.pdf[", frequency.c_str()));


    g_alpha_0->SetLineColor(kBlack);
    g_beta_0 ->SetLineColor(kBlack);
    g_A_0    ->SetLineColor(kBlack);
    g_phi_0  ->SetLineColor(kBlack);

    g_alpha_Combo->SetLineColor(kRed);
    g_beta_Combo ->SetLineColor(kRed);
    g_A_Combo    ->SetLineColor(kRed);
    g_phi_Combo  ->SetLineColor(kRed);

    printf("...creating multigraphs\n");
    TMultiGraph *mg_alpha = new TMultiGraph();
    TMultiGraph *mg_beta = new TMultiGraph();
    TMultiGraph *mg_A = new TMultiGraph();
    TMultiGraph *mg_phi = new TMultiGraph();

    printf("...setting titles\n");
    mg_alpha->SetTitle("alpha_calo0; time [#mus]; [arb]");
    mg_beta ->SetTitle("beta_calo0; time [#mus]; [arb]");
    mg_A    ->SetTitle("A_calo0; time [#mus]; [arb]");
    mg_phi  ->SetTitle("phi_calo0; time [#mus]; [rad]");
    
    mg_alpha->Add(g_alpha_0);
    mg_alpha->Add(g_alpha_Combo);
    mg_alpha->Add(g0alphaDat);
    mg_beta->Add(g_beta_Combo);
    mg_beta->Add(g_beta_0);
    mg_beta->Add(g0betaDat);
    mg_A->Add(g_A_0);
    mg_A->Add(g_A_Combo);
    mg_phi->Add(g_phi_0);
    mg_phi->Add(g_phi_Combo);

    mg_alpha->Draw("ACE");
    c->Print(Form("comparison_%s.pdf", frequency.c_str()));
    mg_beta->Draw("ACE");
    c->Print(Form("comparison_%s.pdf", frequency.c_str()));
    //mg_A->Draw("ACE");
    //c->Print(Form("comparison_%s.pdf", frequency.c_str()));
    //mg_phi->Draw("ACE");
    //c->Print(Form("comparison_%s.pdf", frequency.c_str()));

    c->Print(Form("comparison_%s.pdf]", frequency.c_str()));

}
