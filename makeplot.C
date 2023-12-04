void makeplot(std::string input, std::string output, std::string pdfName){

  // read in file
  TFile* f = new TFile(input.c_str(), "READ");
  // get wiggle (data) and bestfit (fit function)
  TH1D* wig = (TH1D*) f->Get("wiggle")->Clone("wig");
  TH1D* fit = (TH1D*) f->Get("bestfit")->Clone("fit");

  int startBin = 203;

  // wX holds wrapped data lines (each line is a number X)
  // fX holds the bestfit results
  TH1D* w1 = new TH1D("w1", "w1", 1000, 0, 149.2);
  TH1D* w2 = new TH1D("w2", "w2", 1000, 0, 149.2);
  TH1D* w3 = new TH1D("w3", "w3", 1000, 0, 149.2);
  TH1D* w4 = new TH1D("w4", "w4", 1000, 0, 149.2);
  TH1D* w5 = new TH1D("w5", "w5", 154, 0, 22.9768);
  TH1D* f1 = new TH1D("f1", "f1", 1000-startBin, startBin*0.1492 + 0, 149.2);
  TH1D* f2 = new TH1D("f2", "f2", 1000, 0, 149.2);
  TH1D* f3 = new TH1D("f3", "f3", 1000, 0, 149.2);
  TH1D* f4 = new TH1D("f4", "f4", 1000, 0, 149.2);
  TH1D* f5 = new TH1D("f5", "f5",  154, 0, 22.9768);

  // gX holds bestfit as a TGraph
  TGraph *g1 = new TGraph();
  TGraph *g2 = new TGraph();
  TGraph *g3 = new TGraph();
  TGraph *g4 = new TGraph();
  TGraph *g5 = new TGraph();
  // dX holds data s a TGraph
  TGraph *d1 = new TGraph();
  TGraph *d2 = new TGraph();
  TGraph *d3 = new TGraph();
  TGraph *d4 = new TGraph();
  TGraph *d5 = new TGraph();

  TMultiGraph *mg = new TMultiGraph();

  f1->SetTitle("");
  f1->GetXaxis()->SetTitle("time [us]");
  f1->GetYaxis()->SetTitle("counts");

  w1->SetMinimum(1000);
  w2->SetMinimum(1000);
  w3->SetMinimum(1000);
  w4->SetMinimum(1000);
  w5->SetMinimum(1000);
  g1->SetMinimum(1000);
  g2->SetMinimum(1000);
  g3->SetMinimum(1000);
  g4->SetMinimum(1000);
  g5->SetMinimum(1000);

  //mg->SetMinimum(1000);

  w1->SetMaximum(50000000);
  w2->SetMaximum(50000000);
  w3->SetMaximum(50000000);
  w4->SetMaximum(50000000);
  w5->SetMaximum(50000000);
  g1->SetMaximum(50000000);
  g2->SetMaximum(50000000);
  g3->SetMaximum(50000000);
  g4->SetMaximum(50000000);
  g5->SetMaximum(50000000);

  //mg->SetMaximum(50000000);

  int g1point = 0;
  // initially we had a for-loop going from i=0 to 1000
  // accessing wiggle 1 - 5000
  for(int i=0; i<1000; i++){
    w1->SetBinContent(i+1, wig->GetBinContent(i+1));
    w2->SetBinContent(i+1, wig->GetBinContent(i+1001));
    w3->SetBinContent(i+1, wig->GetBinContent(i+2001));
    w4->SetBinContent(i+1, wig->GetBinContent(i+3001));
    //if (fit->GetBinContent(i)==0){
    //   continue;
    //}
    f1->SetBinContent(i+1, fit->GetBinContent(startBin + i+1));
    f2->SetBinContent(i+1, fit->GetBinContent(i+1001));
    f3->SetBinContent(i+1, fit->GetBinContent(i+2001));
    f4->SetBinContent(i+1, fit->GetBinContent(i+3001));
    if(i<154){
      w5->SetBinContent(i+1, wig->GetBinContent(i+4001));
      f5->SetBinContent(i+1, fit->GetBinContent(i+4001));
      g5->SetPoint(i, fit->GetBinCenter(i+1), fit->GetBinContent(i+4001));
      d5->SetPoint(i, wig->GetBinCenter(i+1), wig->GetBinContent(i+4001));
    }

    if(fit->GetBinContent(i+1)>0){
       g1->SetPoint(g1point, fit->GetBinCenter(i+1), fit->GetBinContent(i+1));
       d1->SetPoint(g1point, wig->GetBinCenter(i+0), wig->GetBinContent(i+1));
       g1point++;
    }
    g2->SetPoint(i, fit->GetBinCenter(i+1), fit->GetBinContent(i+1001));
    g3->SetPoint(i, fit->GetBinCenter(i+1), fit->GetBinContent(i+2001));
    g4->SetPoint(i, fit->GetBinCenter(i+1), fit->GetBinContent(i+3001));
    d2->SetPoint(i, wig->GetBinCenter(i+1), wig->GetBinContent(i+1001));
    d3->SetPoint(i, wig->GetBinCenter(i+1), wig->GetBinContent(i+2001));
    d4->SetPoint(i, wig->GetBinCenter(i+1), wig->GetBinContent(i+3001));
  }
  f->Close();

  d1->SetLineStyle(6);
  d2->SetLineStyle(6);
  d3->SetLineStyle(6);
  d4->SetLineStyle(6);
  d5->SetLineStyle(6);
  g1->SetLineStyle(7);
  g2->SetLineStyle(7);
  g3->SetLineStyle(7);
  g4->SetLineStyle(7);
  g5->SetLineStyle(7);

  w1->SetMarkerStyle(6);
  w1->SetMarkerColor(kBlack);
  w1->SetMarkerSize(0.2);
  w2->SetMarkerStyle(6);
  w2->SetMarkerColor(kBlack);
  w2->SetMarkerSize(0.2);
  w3->SetMarkerStyle(6);
  w3->SetMarkerColor(kBlack);
  w3->SetMarkerSize(0.2);
  w4->SetMarkerStyle(6);
  w4->SetMarkerColor(kBlack);
  w4->SetMarkerSize(0.2);
  w5->SetMarkerStyle(6);
  w5->SetMarkerColor(kBlack);
  w5->SetMarkerSize(0.2);
  f1->SetLineColor(kRed);
  f1->SetLineWidth(1.0);
  f2->SetLineColor(kRed);
  f2->SetLineWidth(1.0);
  f3->SetLineColor(kRed);
  f3->SetLineWidth(1.0);
  f4->SetLineColor(kRed);
  f4->SetLineWidth(1.0);
  f5->SetLineColor(kRed);
  f5->SetLineWidth(1.0);

  g1->SetLineColor(kRed);
  g2->SetLineColor(kRed);
  g3->SetLineColor(kRed);
  g4->SetLineColor(kRed);
  g5->SetLineColor(kRed);

  d1->SetLineColor(kBlue);
  d2->SetLineColor(kBlue);
  d3->SetLineColor(kBlue);
  d4->SetLineColor(kBlue);
  d5->SetLineColor(kBlue);


  g1->SetLineWidth(2.);
  g2->SetLineWidth(2.);
  g3->SetLineWidth(2.);
  g4->SetLineWidth(2.);
  g5->SetLineWidth(2.);

  //mg->SetLineColor(kRed);
  //mg->SetLineWidth(3.);

  TCanvas* ct = new TCanvas("ct","ct", 0, 0, 600, 400);
  gPad->SetLogy();

  mg->SetTitle(";time [us]; counts");

//  printf("drawing g1-g5\n");
//  g1->Draw();
//  g2->Draw();
//  g3->Draw();
//  g4->Draw();
//  g5->Draw();

//  //g1->GetXaxis()->SetRangeUser(0,149.2);
//  //g2->GetXaxis()->SetRangeUser(0,149.2);
//  //g3->GetXaxis()->SetRangeUser(0,149.2);
//  //g4->GetXaxis()->SetRangeUser(0,149.2);
//  //g5->GetXaxis()->SetRangeUser(0,149.2);
//
//  printf("setting limits\n");
//  g1->GetXaxis()->SetLimits(0,149.2+0.1492);
//  g2->GetXaxis()->SetLimits(0,149.2+0.1492);
//  g3->GetXaxis()->SetLimits(0,149.2+0.1492);
//  g4->GetXaxis()->SetLimits(0,149.2+0.1492);
//  g5->GetXaxis()->SetLimits(0,149.2+0.1492);
//
//  printf("redrawing\n");
//  g1->Draw("AL");
//  g2->Draw("AL");
//  g3->Draw("AL");
//  g4->Draw("AL");
//  g5->Draw("AL");
//
//  ct->Update();
//

  // draw data

  // draw best-fit red curve 
  mg->Add(g1);
  mg->Add(g2);
  mg->Add(g3);
  mg->Add(g4);
  mg->Add(g5);
  mg->Add(d1);
  mg->Add(d2);
  mg->Add(d3);
  mg->Add(d4);
  mg->Add(d5);
  mg->Draw("AL");

  auto legend = new TLegend();
  legend->AddEntry(d1, "data");
  legend->AddEntry(g1, "bestfit");
  legend->Draw();

  ct->Print(pdfName.c_str());
  
//  mg->GetXaxis()->SetLimits(0,149.2+0.1492);
//  ct->Update();
//
  //exit(0); 
  //w1->Draw("SAME P");
  //w2->Draw("SAME P");
  //w3->Draw("SAME P");
  //w4->Draw("SAME P");
  //w5->Draw("SAME P");
  //
  //printf("printing PDF\n");
  //ct->Print(pdfName.c_str());
}
