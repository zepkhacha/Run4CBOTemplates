void plotFFT(std::string inputFile){

// open specified file
TFile *f = TFile::Open(inputFile.c_str()); //point to fitfile
// retrieve residuals histogram
TH1D  *r = (TH1D*)f->Get("residua");
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

// initiate corrected FFT with x-axis up to 6.7MHz
TH1D * fft = new TH1D("FFT of Residuals","FFT of Residuals", nBins, bin_i, bin_f);
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
std::string outputFilename = inputFile.substr(0, inputFile.size()-5);
TCanvas *c = new TCanvas();
gStyle->SetOptStat(0);
fft->Draw();
c->Print(Form("%s.pdf",outputFilename.c_str()));

}
