#include <string>
#include <iostream>
#include <math.h>

#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TColor.h>
#include <TAxis.h>

#include </gm2data/zkhechad/cboTemplates/aMethod/fullFits/CornellFitAlgorithm.hh>

void printhelp(){
  printf(
  "Analysis Program\n"
  "Supply arguments as follows:\n"
  "-h        : Prints this help text\n"
  "-i [s]    : Input filename with sliding windows (required)\n"
  "-f [s]    : Input filename with full-fit result (required)\n"
  "-b [s]    : Boost filename (required)\n"
  "-o [s]    : Output file label (required)\n"
  );
}

void parse_cmdline(int argc, char** argv, 
  std::string& inputFilename, 
  std::string& boostFilename, 
  std::string& outputLabel){
  const char* const opts = "hi:b:o:";
  bool done = false;
  while(!done){
    const char ch = getopt(argc, argv, opts);
    switch(ch){
      case -1: done = true; break;
      case 'h': printhelp(); exit(0);
      case 'i': inputFilename = optarg; break;
      case 'b': boostFilename = optarg; break;
      case 'o': outputLabel = optarg; break;
    } 
  }
}

double radToDeg(double rad){
  return rad*( 360 / (2*M_PI) );
}

double radToNPi(double rad){
  return rad / M_PI;
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
  std::string inputFilename;
  std::string boostFilename;
  std::string outputLabel;
  parse_cmdline(argc, argv, inputFilename, boostFilename, outputLabel);

  gStyle->SetPalette(kThermometer, 0, 0.5);
  gStyle->SetLineColor(kBlack);
  gStyle->SetMarkerColor(kBlack);

  bool axesInitialized = false;

  printf("Opening %s\n", inputFilename.c_str());
  TFile *fullFitFile = TFile::Open(inputFilename.c_str(), "READ");
  TTree* fullFitResults = (TTree*) fullFitFile->Get("fitresults");

  // get best-fit histogram from full-fit to remove CBO
  TH1D* bestfit       = (TH1D*) fullFitFile->Get("bestfit");
  TH1D* wiggle        = (TH1D*) fullFitFile->Get("wiggle");

  // get all parameters from the full-fit to rerun the fit function without CBO
  int    fullFit_sidx; 
  double fullFit_chisq, fullFit_rchisq;
  double fullFit_N0, fullFit_tau, fullFit_A0, fullFit_phi, fullFit_R;
  double fullFit_T_CBO, fullFit_w_CBO, fullFit_A_CNx1, fullFit_A_SNx1, fullFit_C_CBO;
  double fullFit_LM, fullFit_A_CNx2, fullFit_A_SNx2;
  double fullFit_A_CAx1, fullFit_A_SAx1, fullFit_A_SNy1;
  double fullFit_Ky, fullFit_Ty, fullFit_A_CNy2, fullFit_A_SNy2;
  double fullFit_A_Cp, fullFit_A_Sp;
  double fullFit_A_SNxy, fullFit_A_CNxy;
  double fullFit_A_ct, fullFit_T_xy, fullFit_w_xy;
  double fullFit_A_CNy1;
  //double fullFit_ANx1, fullFit_pNx1;
  //double fullFit_ANx2, fullFit_pNx2;
  //double fullFit_AAx1, fullFit_pAx1;
  //double fullFit_ANy1, fullFit_pNy1;
  //double fullFit_ANy2, fullFit_pNy2;
  //double fullFit_Ap, fullFit_pp;
  //double fullFit_ANxy, fullFit_pNxy;
  double fullFit_Gamma_mop;
  double fullFit_c_e_scale;

  double fullFit_phiCBO; 
  double fullFit_phiX2;

  bool fullFit_pBinTau, fullFit_pBinPhi, fullFit_pBinCBO;
  bool fullFit_includeMopTerm; 
  bool fullFit_constrainMop;

  fullFit_pBinCBO = false;
  
  // full fit result should only have one entry in TTree
  if (fullFitResults->GetEntries()!=1){
    std::cerr<<"Full-fit file should have only 1 entry in 'fitresults' TTree.";
  }

  fullFitResults->SetBranchAddress("pBinTau", &fullFit_pBinTau);
  fullFitResults->SetBranchAddress("pBinPhi", &fullFit_pBinPhi);
  //fullFitResults->SetBranchAddress("pBinCBO", &fullFit_pBinCBO);
  fullFitResults->SetBranchAddress("includeMopTerm", &fullFit_includeMopTerm);
  fullFitResults->SetBranchAddress("constrainMop", &fullFit_constrainMop);
  fullFitResults->SetBranchAddress("seed", &fullFit_sidx);
  fullFitResults->SetBranchAddress("chisq", &fullFit_chisq);
  fullFitResults->SetBranchAddress("rchisq", &fullFit_rchisq);
  fullFitResults->SetBranchAddress("N0",  &fullFit_N0);
  fullFitResults->SetBranchAddress("tau", &fullFit_tau);
  fullFitResults->SetBranchAddress("A0",  &fullFit_A0);
  fullFitResults->SetBranchAddress("phi", &fullFit_phi);
  fullFitResults->SetBranchAddress("R", &fullFit_R);
  fullFitResults->SetBranchAddress("C_CBO", &fullFit_C_CBO);
  fullFitResults->SetBranchAddress("T_CBO", &fullFit_T_CBO);
  fullFitResults->SetBranchAddress("w_CBO", &fullFit_w_CBO);
  fullFitResults->SetBranchAddress("A_CNx1", &fullFit_A_CNx1);
  fullFitResults->SetBranchAddress("A_SNx1", &fullFit_A_SNx1);
  fullFitResults->SetBranchAddress("LM", &fullFit_LM);
  fullFitResults->SetBranchAddress("A_CNx2", &fullFit_A_CNx2);
  fullFitResults->SetBranchAddress("A_SNx2", &fullFit_A_SNx2);
  fullFitResults->SetBranchAddress("A_CAx1", &fullFit_A_CAx1);
  fullFitResults->SetBranchAddress("A_SAx1", &fullFit_A_SAx1);
  fullFitResults->SetBranchAddress("A_CNy1", &fullFit_A_CNy1);
  fullFitResults->SetBranchAddress("A_SNy1", &fullFit_A_SNy1);
  fullFitResults->SetBranchAddress("A_SNxy", &fullFit_A_SNxy);
  fullFitResults->SetBranchAddress("A_CNxy", &fullFit_A_CNxy);
  fullFitResults->SetBranchAddress("Ky", &fullFit_Ky);
  fullFitResults->SetBranchAddress("Ty", &fullFit_Ty);
  fullFitResults->SetBranchAddress("A_CNy2", &fullFit_A_CNy2);
  fullFitResults->SetBranchAddress("A_SNy2", &fullFit_A_SNy2);
  fullFitResults->SetBranchAddress("A_Cp", &fullFit_A_Cp);
  fullFitResults->SetBranchAddress("A_Sp", &fullFit_A_Sp);
  fullFitResults->SetBranchAddress("A_ct", &fullFit_A_ct);
  fullFitResults->SetBranchAddress("T_xy", &fullFit_T_xy);
  fullFitResults->SetBranchAddress("w_xy", &fullFit_w_xy);
  fullFitResults->SetBranchAddress("c_e_scale"  , &fullFit_c_e_scale);
  fullFitResults->SetBranchAddress("Gamma_mop", &fullFit_Gamma_mop);


  fullFitResults->SetBranchAddress("R", &fullFit_R);
  fullFitResults->SetBranchAddress("T_CBO", &fullFit_T_CBO);
  fullFitResults->SetBranchAddress("pNx1" , &fullFit_phi);
  fullFitResults->SetBranchAddress("pNx1", &fullFit_phiCBO);
  fullFitResults->SetBranchAddress("pNx2", &fullFit_phiX2);

  // we also need boost information to rerun calcnu...
  TFile* boostdata  = new TFile(boostFilename.c_str(), "READ");
  //TH1D*  boostgamma = (TH1D*) boostdata->Get("transform_gamma");
  //frlifetime = boostgamma->GetMean();
  //*************
  // get momentum-binned weights
  //*************

  TH1F* transform_dp_p0 = (TH1F*) boostdata->Get("transform_dp_p0");
  TH1F* transform_gamma = (TH1F*) boostdata->Get("transform_gamma");
  TH1F* transform_x     = (TH1F*) boostdata->Get("transform_x");

  int pBins = int(transform_dp_p0->GetNbinsX());
  std::pair <double,double> physicalRange (-0.0056799,0.0056085);  // in dp/p0
 
  // first isolate the physical range of the spectra
  // get the integral and bin numbers 
  std::vector<double> pWeights(pBins, 0.0); 

  for (int i=0; i<pBins; i++){
     double lowEdge = transform_dp_p0->GetBinLowEdge(i);
     if (physicalRange.first <= lowEdge) {
	binRange.first = i;
        break;
     }
  }
  for (int i=pBins-1; i>=0; i--){
     double upEdge = transform_dp_p0->GetBinLowEdge(i+1);
     if (physicalRange.second >= upEdge) {
	binRange.second = i;
        break;
     }
  }

  beta_0 = sqrt(1.0 - (1.0 / pow(29.3,2.0))); // this uses design gamma
  R_0    = 7112; // design R_0 in mm
  // now go through bins and store values for dp_p0, dPhi, gammas, tau, x
  for (int i=0; i<pBins; i++){
      dp_p0 .push_back(transform_dp_p0       ->GetBinLowEdge(i));
      gammas.push_back(transform_gamma       ->GetBinLowEdge(i));
      dPhi  .push_back(-1.351 * transform_dp_p0 ->GetBinLowEdge(i)); // negative if fitting cos(wt - phi), pos if wt + phi // run6 13.51, run1 10.0
      x_e   .push_back(transform_x->GetBinLowEdge(i)); //convert from radial offset to radial position
      c_e   .push_back( c_e_scale * (2.0 * ((beta_0*beta_0)/R_0) * n 
				* transform_x->GetBinLowEdge(i) 
				* transform_dp_p0->GetBinLowEdge(i)) );
  }
  double totalContent = transform_dp_p0->Integral(binRange.first, binRange.second);

  // now find weights in only the physical range
  totalContent = transform_dp_p0->Integral(binRange.first, binRange.second);
  double sumWeights = 0.0;
  for (int i=1; i<pBins; i++){
     // find magicBin 
     double lowEdge = transform_dp_p0->GetBinLowEdge(i);
     double uppEdge = transform_dp_p0->GetBinLowEdge(i+1);

     if ( lowEdge <= 0.0  and 0.0 < uppEdge){
        magicBin = i;
	// update wc_0 in magic bin which we use in p-dependent wCBO; includes conversion from ns to us 
	wc_0 = (1.0 / ((transform_gamma->GetBinLowEdge(magicBin)*2.1969811) / 1000.)); 
     }

     // find weights according to momentum
     double thisWeight = transform_dp_p0->GetBinContent(i) / totalContent;

     // is this bin in a physical momentum bin?
     bool isPhysical = (i>=binRange.first and i<= binRange.second);

     if (pBinTau or pBinPhi){
        binWeights.push_back( isPhysical? thisWeight : 0.0 );
     }else{ // if we're not binning, set weight=1.0 in magic bin and 0.0 otherwise
        binWeights.push_back( (i==magicBin) ? 1.0 : 0.0 );
     }
     sumWeights += binWeights[i];
     weightedAvgTau += binWeights[i]*transform_gamma->GetBinLowEdge(i)*2.1969811;
     //printf("bin %i content %f momentum %f isPhysical? %d weight %f\n",i,transform_dp_p0->GetBinContent(i),dp_p0[i],isPhysical,binWeights.at(binWeights.size()-1));
  }

  boostdata->Close();
  //*************
  // end get momentum-binned weights + taus
  //*************

      TH1D* gFullFitComparedWithWithoutCBODiv = new TH1D("dataDivFnc_wo_CBO", "dataDivFnc_wo_CBO", nBins, 0, 650.0644);
      TH1D* gFullFitComparedWithWithoutYDiv = new TH1D("dataDivFnc_wo_Y", "dataDivFnc_wo_Y", nBins, 0, 650.0644);

      gFullFitComparedWithWithoutCBODiv->SetMaximum(1.64);
      gFullFitComparedWithWithoutCBODiv->SetMinimum(1.56);
      gFullFitComparedWithWithoutYDiv->SetMaximum(1.7);
      gFullFitComparedWithWithoutYDiv->SetMinimum(1.5);

      // get fullFitResults
      fullFitResults->GetEntry(0);

      // now create an array with these fit results
      double par[29] = {
      fullFit_N0, fullFit_tau, fullFit_A0, fullFit_phi, fullFit_R,     
      fullFit_T_CBO, fullFit_w_CBO, fullFit_A_CNx1, fullFit_A_SNx1, 
      fullFit_LM, 
      fullFit_A_CNx2, fullFit_A_SNx2, fullFit_A_CAx1, fullFit_A_SAx1, fullFit_A_CNy1, fullFit_A_SNy1,  
      fullFit_Ky, fullFit_Ty, 
      fullFit_A_CNy2, fullFit_A_SNy2, fullFit_A_Cp, fullFit_A_Sp, fullFit_A_CNxy, fullFit_A_SNxy, 
      fullFit_A_ct, 
      fullFit_T_xy, fullFit_w_xy, 
      fullFit_Gamma_mop,
      fullFit_C_CBO};
      double noCBOpar[29] = {
      fullFit_N0, fullFit_tau, fullFit_A0, fullFit_phi, fullFit_R,     
      fullFit_T_CBO, fullFit_w_CBO, 0.0, 0.0, 
      fullFit_LM, 
      0.0, 0.0, fullFit_A_CAx1, fullFit_A_SAx1, fullFit_A_CNy1, fullFit_A_SNy1,  
      fullFit_Ky, fullFit_Ty, 
      fullFit_A_CNy2, fullFit_A_SNy2, fullFit_A_Cp, fullFit_A_Sp, fullFit_A_CNxy, fullFit_A_SNxy, 
      fullFit_A_ct, 
      fullFit_T_xy, fullFit_w_xy, 
      fullFit_Gamma_mop,
      fullFit_C_CBO};
      double noYpar[29] = {
      fullFit_N0, fullFit_tau, fullFit_A0, fullFit_phi, fullFit_R,     
      fullFit_T_CBO, fullFit_w_CBO, fullFit_A_CNx1, fullFit_A_SNx1, 
      fullFit_LM, 
      fullFit_A_CNx2, fullFit_A_SNx2, fullFit_A_CAx1, fullFit_A_SAx1, 0.0, 0.0,  
      fullFit_Ky, fullFit_Ty, 
      0.0, 0.0, fullFit_A_Cp, fullFit_A_Sp, fullFit_A_CNxy, fullFit_A_SNxy, 
      fullFit_A_ct, 
      fullFit_T_xy, fullFit_w_xy, 
      fullFit_Gamma_mop,
      fullFit_C_CBO};

      // draw wiggle-fullFit - wiggle-fullFit without CBO
      for (int point=1; point<bestfit->GetNbinsX(); point++){
          double dim[1] = {double(point)};
	  double wiggleData = wiggle->GetBinContent(point);
	  double bestfitData = bestfit->GetBinContent(point);
	  double fullFitVal_withCBO_fromFunction = calcnu(dim, par);
	  double fullFitVal_noCBO_fromFunction = calcnu(dim, noCBOpar);
	  double fullFitVal_noY_fromFunction = calcnu(dim, noYpar);
          gFullFitComparedWithWithoutCBODiv->SetBinContent(point, wiggleData / fullFitVal_noCBO_fromFunction);
	  gFullFitComparedWithWithoutYDiv->SetBinContent(point, wiggleData / fullFitVal_noY_fromFunction);
      }

      // only add option 'a' for first plot
      std::string drawOption = Form("%sPE%s", axesInitialized ? "" : "A", "PMC PLC");
      axesInitialized = true;

      // one TCanvas to hold them all
      std::string outputFilename = Form("%s", outputLabel.c_str());
      TCanvas c;
      c.cd();

      c.Clear();
      gFullFitComparedWithWithoutCBODiv->Draw();
      c.Print(Form("%s_cboIsolate.pdf", outputLabel.c_str()));

      c.Clear();
      gFullFitComparedWithWithoutYDiv->Draw();
      c.Print(Form("%s_yIsolate.pdf", outputLabel.c_str()));

      TFile *fOutput = new TFile(TString::Format("%s.root", outputLabel.c_str()),"RECREATE");
      gFullFitComparedWithWithoutCBODiv->Write();
      gFullFitComparedWithWithoutYDiv->Write();
      fOutput->Close();

} // end main
