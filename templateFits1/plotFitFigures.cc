void plotFitFigures(std::string inputFile){

// open specified file
TFile *f = TFile::Open(inputFile.c_str()); //point to fitfile
TCanvas *c = new TCanvas();
// save figures to fitFigures_[inputfilename].pdf
// first get the basename from the full path
std::string base_filename = inputFile.substr(inputFile.find_last_of("/\\") + 1);
std::string outputFilename = base_filename.substr(0, base_filename.size()-5);
c->Print(Form("fitFigures_%s.pdf[",outputFilename.c_str()));

// FIGURE 1: FFT of Residuals
// retrieve residuals histogram
TH1D  *r = (TH1D*)f->Get("residua");
// retrieve fit results for w_CBO
TTree *t = (TTree*)f->Get("fitresults");

double w_CBO, Ky;
double R;

t->SetBranchAddress("w_CBO", &w_CBO);
t->SetBranchAddress("Ky", &Ky);
t->SetBranchAddress("R", &R);
t->GetEntry(0);

double wy = Ky * w_CBO * sqrt(( (4*M_PI) / (0.1492*Ky*w_CBO) ) - 1.0);
double fy   = (wy/(2*M_PI)) * pow(10.,6.); 
double fCBO = (w_CBO/(2*M_PI)) * pow(10.,6.);
double fC   = 1.0 / (0.1492E-6); // cyclotron
double fx   = (fx + fC);
double fVO  = (fC - fy);
double fVW  = (fC - 2*fy);

printf("fy %f fVW %f\n", fy, fVW);

// initialize histogram to hold raw fft (with incorrect x-axis)
TH1* fftRaw  = 0;
fftRaw  = r->FFT(fftRaw,"MAG");

TCanvas *cRaw = new TCanvas();
fftRaw->Draw();
// draw raw histogram for viewing

//NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
// x-axis range is the number of bins (4357)
// range of function is 0 - 650.0644 us
// FFT should go up to 6.7
int    nBins    = r->GetNbinsX();
double bin_i    = 0.0;
double bin_f    = (1 / (0.1492 * 0.000001)); // same as doing x-range / function range, i.e. nBins/650.0644us
double binWidth = (bin_f-bin_i)/nBins;

// if nBins is even, take bins 1 - nBins/2
// if nBins is odd, take bins 1 - (nBins+1)/2

int adjustBins  = (nBins % 2 == 0) ? int(nBins/2) : (nBins+1)/2;
double adjust_f = (nBins % 2 == 0) ? int(bin_f/2) : (bin_f+0.1492)/2;

// initiate corrected FFT with x-axis up to 6.7MHz
TH1D * fft = new TH1D("FFT of Residuals","FFT of Residuals", adjustBins, bin_i, adjust_f);
fft->GetXaxis()->SetTitle("Frequency (Hz)");
fft->GetYaxis()->SetTitle("Power (arb. units)");

// loop through rawFFT and correct each point
for (int x=1;x<=fftRaw->GetNbinsX();x++){
    double freq_i = fftRaw->GetBinCenter(x) ; 
    double freq_f = fftRaw->GetBinCenter(x) * binWidth ; 
    int bin = (freq_f / binWidth) + 1; 
    //if(x<10){printf("freq i %f -> %f bin %i bin [%f,%f]\n", freq_i, freq_f, bin, fft->GetBinLowEdge(bin), fft->GetBinLowEdge(bin+1));}
    fft->SetBinContent(bin,fftRaw->GetBinContent(x));
}

// save output to fft_[inputfilename].pdf
c->cd();
gStyle->SetOptStat(0);
fft->Draw();

double maxY = fft->GetMaximum();

// add lines
TLine *l_fCBO = new TLine(fCBO,0.,fCBO,maxY);
TLine *l_2fCBO = new TLine(2*fCBO,0.,2*fCBO,maxY);
TLine *l_fC    = new TLine(fC  ,0.,fC  ,maxY);
TLine *l_fx   = new TLine(fx  ,0.,fx  ,maxY);
TLine *l_fy   = new TLine(fy  ,0.,fy  ,maxY); 
//TLine *l_fa   = new TLine(fa  ,0.,fa  ,maxY); 
TLine *l_fVO  = new TLine(fVO ,0.,fVO ,maxY);    
TLine *l_fVW  = new TLine(fVW ,0.,fVW ,maxY);    

l_fCBO->SetLineColor(kRed);
l_2fCBO->SetLineColor(kRed);
l_fx  ->SetLineColor(kRed);
l_fy  ->SetLineColor(kRed);
//l_fa  ->SetLineColor(kRed);
l_fCBO->SetLineColor(kRed);
l_fVO ->SetLineColor(kRed);
l_fVW ->SetLineColor(kRed);

l_fCBO->SetLineStyle(2);
l_2fCBO->SetLineStyle(2);
l_fx  ->SetLineStyle(2);
l_fy  ->SetLineStyle(2);
//l_fa  ->SetLineStyle(2);
l_fCBO->SetLineStyle(2);
l_fVO ->SetLineStyle(2);
l_fVW ->SetLineStyle(2);

l_fCBO->Draw("SAME");
l_2fCBO->Draw("SAME");
l_fx  ->Draw("SAME");
l_fy  ->Draw("SAME");
//l_fa  ->Draw("SAME");
l_fVO ->Draw("SAME");
l_fVW ->Draw("SAME");

TLatex latex;
latex.SetTextSize(0.05);
latex.SetTextAngle(90);
latex.SetTextColor(kRed);

double maxVal = 0.85*fft->GetMaximum();

latex.DrawLatex(+ fCBO, maxVal,"f_{CBO}");
latex.DrawLatex(+ 2*fCBO, maxVal,"2f_{CBO}");
latex.DrawLatex(+ fC, maxVal,"f_{C}");
latex.DrawLatex(+ fx, maxVal,"f_{x}");
latex.DrawLatex(+ fy, maxVal,"f_{y}");
//latex.DrawLatex(-50000+ fa, maxVal,"f_{a}");
latex.DrawLatex(+ fCBO, maxVal,"f_{CBO}");
latex.DrawLatex(+ fVO, maxVal,"f_{VO}");
latex.DrawLatex(+ fVW, maxVal,"f_{VW}");


c->Print(Form("fft_%s.pdf",outputFilename.c_str()));
c->Print(Form("fitFigures_%s.pdf",outputFilename.c_str()));

// FIGURE 2: Residuals
c->Clear();
r->Draw();
c->Print(Form("fitFigures_%s.pdf",outputFilename.c_str()));

// FIGURE 3: Wiggle

// close PDF
c->Print(Form("fitFigures_%s.pdf]",outputFilename.c_str()));
}
