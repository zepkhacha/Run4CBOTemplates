void paramVparam(std::string filename, 
                std::string treeName, 
                std::string param1name, std::string param1type,
                std::string param2name, std::string param2type,
                std::string outputName){

    TFile *f = TFile::Open(filename.c_str(),"READ");
    TTree *tree = (TTree*)f->Get(treeName.c_str());
    
    int param1, param1err;
    double param2, param2err;

    tree->SetBranchAddress(Form("%s", param1name.c_str()), &param1);
    //tree->SetBranchAddress(Form("%serr", param1name.c_str()), &param1);
    tree->SetBranchAddress(Form("%s", param2name.c_str()), &param2);
    tree->SetBranchAddress(Form("%serr", param2name.c_str()), &param2err);

    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("%s vs %s; %s; %s", param1name.c_str(), param2name.c_str(), 
                                        param1name.c_str(), param2name.c_str()));
    
    for (int i=0; i<tree->GetEntries(); i++){
        tree->GetEntry(i);
        g->SetPoint(i, param1, param2);
        g->SetPointError(i, 0.0, param2err);
        printf("i %i x %f y %f\n", i, param1, param2);
    }

    TCanvas *c = new TCanvas();
    g->Draw("APE");
    c->Print(outputName.c_str());

}
