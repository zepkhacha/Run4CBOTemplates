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

void showDiffs(
        std::string input1Filename, 
        std::string input2Filename){

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kDarkRainBow, 0, 0.5);
    //gStyle->SetLineColor(kBlack);
    //gStyle->SetMarkerColor(kBlack);

    double start1;
    int windowNo1;
    int calo1;
    double R1, w_CBO1, A_CNx1_1, chisq1, LM1;

    double start2;
    int windowNo2;
    int calo2;
    double R2, w_CBO2, A_CNx1_2, chisq2, LM2;

    TFile fitFile1 (input1Filename.c_str());
    TFile fitFile2 (input1Filename.c_str());

    TTree* fitResults1 = (TTree*) fitFile1.Get("fitresults");
    TTree* fitResults2 = (TTree*) fitFile2.Get("fitresults");

    fitResults1->SetBranchAddress("calo", &calo1);
    fitResults1->SetBranchAddress("windowNo", &windowNo1);
    fitResults1->SetBranchAddress("start", &start1);
    fitResults1->SetBranchAddress("R", &R1);
    fitResults1->SetBranchAddress("w_CBO", &w_CBO1);
    fitResults1->SetBranchAddress("chisq", &chisq1);
    fitResults1->SetBranchAddress("LM", &LM1);
    fitResults1->SetBranchAddress("A_CNx1", &A_CNx1_1);

    fitResults2->SetBranchAddress("calo", &calo2);
    fitResults2->SetBranchAddress("windowNo", &windowNo2);
    fitResults2->SetBranchAddress("start", &start2);
    fitResults2->SetBranchAddress("R", &R2);
    fitResults2->SetBranchAddress("w_CBO", &w_CBO2);
    fitResults2->SetBranchAddress("chisq", &chisq2);
    fitResults2->SetBranchAddress("LM", &LM2);
    fitResults2->SetBranchAddress("A_CNx1", &A_CNx1_2);

    for (int entry=0; entry<fitResults1->GetEntries(); entry++){
        fitResults1->GetEntry(entry);
        fitResults2->GetEntry(entry);
        if (calo1==calo2){
            printf("calo %i diffs (file2-file1) : \n", calo1);
            printf("R : (%f - %f) = %f \n", R2, R1, R2-R1);
            printf("w_CBO : (%f - %f) = %f \n", w_CBO2, w_CBO1, w_CBO2-w_CBO1);
            printf("A_CNx1: (%f - %f) = %f \n", A_CNx1_2, A_CNx1_1, A_CNx1_2-A_CNx1_1);
            printf("chisq: (%f - %f) = %f \n", chisq2, chisq1, chisq2-chisq1);
            printf("LM: (%f - %f) = %f \n", LM2, LM1, LM2-LM1);
        }
    }


} // end main
