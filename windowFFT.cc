#include "/gm2data/zkhechad/cboTemplates/Run4CBOTemplates/slidingFits1/slidingFitFunc.hh"
#include "gm2util/blinders/Blinders.hh"
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TLatex.h>
#include <TLine.h>

void printhelp(){
    printf(
            "Supply arguments as follows:\n"
            "-h        : Prints this help text\n"
            "-i [s]    : Input filename (required)\n"
            "-o [s]    : Output filename (required)\n"
          );  
}

void parse_cmdline(int argc, char** argv, 
        char* &inputFilename, char* &fullFitFilename, char* &outputFilename){
    const char* const opts = "ho:i:f:";
    bool done = false;
    while(!done){
        const char ch = getopt(argc, argv, opts);
        switch(ch){
            case -1: done = true; break;
            case 'i': inputFilename = optarg; break;
            case 'o': outputFilename = optarg; break;
            case 'f': fullFitFilename = optarg; break;
            default: printhelp(); exit(1);
        }
    }
}

double wy(double kappa, double wcbo){
    double x = (4*M_PI)/(0.1492*kappa*wcbo) - 1.0;
    return kappa*wcbo*sqrt(x);
}

int main(int argc, char* argv[]){

    // get input arguments
    char* outputFilename;
    char* inputFilename ;
    char* fullFitFilename;
    parse_cmdline(argc,argv,inputFilename, fullFitFilename, outputFilename);

    printf("sliding input: %s\n", inputFilename);
    printf("full window input: %s\n", fullFitFilename);

    // noRF_calo1_windowFits.root
    std::string inputString = inputFilename;
    inputString = inputString.substr(inputString.find_last_of("/\\") + 1);
    // find run X
    int index = inputString.find("_");
    std::string run = inputString.substr(0,index); 
    printf("run: %s\n", run.c_str());
    inputString = inputString.substr(index,inputString.length()-index);
    index = inputString.find("_");
    // find caloXX
    inputString = inputString.substr(index+1, inputString.length() - index);
    index = inputString.find("_");
    std::string calo = inputString.substr(0, index);
    // print for checking
    printf("input %s -> %s calo %s", inputString.c_str(), run.c_str(), calo.c_str());

    const int    r_nBins  = 4357 ;
    double startTime = 0;

    printf("Opening %s\n", inputFilename);
    TFile *f = TFile::Open(inputFilename); //point to fitfile
    TTree *t = (TTree*)f->Get("fitresults");
    TFile *f_full = TFile::Open(fullFitFilename); //point to fitfile
    TTree *t_full = (TTree*)f_full->Get("fitresults");

    printf("Reading fitresults ttree\n");
    double fitStart, fitStop;
    double prevStart, prevStop;
    int windowNo;
    int windowBins;

    double w_CBO, Ky;
    double R;

    t->SetBranchAddress("start", &fitStart);
    t->SetBranchAddress("stop", &fitStop);
    t->SetBranchAddress("windowNo", &windowNo);
    t->SetBranchAddress("windowBins", &windowBins);
    //t->SetBranchAddress("w_CBO", &w_CBO);
    //t->SetBranchAddress("Ky", &Ky);
    //t->SetBranchAddress("R", &R);

    double full_wCBO, full_Ky, full_R;
    // use full fit results for the marker lines (?)
    t_full->SetBranchAddress("w_CBO", &full_wCBO);
    t_full->SetBranchAddress("Ky", &full_Ky);
    t_full->SetBranchAddress("R", &full_R);

    t_full->GetEntry(0);
    double full_wy= wy(full_Ky, full_wCBO);

    TH1D *rFFTSum = new TH1D("residua", "residua", r_nBins, startTime, 650.0644);
    TCanvas *canvas = new TCanvas();
    canvas->Print(Form("%s.pdf[", outputFilename));

    prevStop = 0.0;
    for (int entry=0; entry<t->GetEntries(); entry++){
        t->GetEntry(entry);
        //printf("doing entry %i (%f, %f)\n", entry, fitStart, fitStop);
        //int fitrangelow = int(fitStart / 0.1492);

        if (fitStart < prevStop){
            continue;
        }else{
            //printf("non overlapping entry %i (%f, %f)\n", entry, fitStart, fitStop);
            //printf("fCBO %f fVW %f fy %f \n",
            //        fCBO, fVW, fy) ;

            prevStart=fitStart;
            prevStop=fitStop;
            TFile *windowFile = TFile::Open(Form("slidingFits1/%s/%s_%s_window%04d.root",
                        run.c_str(), run.c_str(), calo.c_str(), windowNo));

            TH1D  *r = (TH1D*)windowFile->Get("residua");
            TH1D  *wiggle = (TH1D*)windowFile->Get("wiggle");
            TH1D  *bestfit = (TH1D*)windowFile->Get("bestfit");
            TH1D  *rTrimmed  = new TH1D("rTrimmed","rTrimmed",windowBins,fitStart,fitStop);
            TH1D  *windowDat = new TH1D("windowDat","windowDat",windowBins,fitStart,fitStop);
            TH1D  *windowFit = new TH1D("windowFit","windowFit",windowBins,fitStart,fitStop);

            TTree *wfitresults = (TTree*)windowFile->Get("fitresults");
            double wCBO;
            double w_y;
            double w_vw;
            double R;
            wfitresults->SetBranchAddress("w_CBO", &wCBO);
            wfitresults->SetBranchAddress("w_vw", &w_vw);
            wfitresults->SetBranchAddress("w_y", &w_y);
            wfitresults->SetBranchAddress("R", &R );
            wfitresults->GetEntry(0);

            for (int bin=1;bin<=windowBins;bin++){
                int refBin = int(fitStart/0.1492) + bin - 1;
                rTrimmed->SetBinContent(bin, r->GetBinContent(refBin));
                windowDat->SetBinContent(bin, wiggle->GetBinContent(refBin));
                windowFit->SetBinContent(bin, bestfit->GetBinContent(refBin));
            }

            TH1 * rTrimmedFFTraw = 0;
            rTrimmedFFTraw = rTrimmed->FFT(rTrimmedFFTraw,"MAG");
            rTrimmedFFTraw->SetName(Form("rTrimmedFFTraw_window%i",entry));
            //rTrimmed->Draw();
            //canvas->Print(Form("%s.pdf", outputFilename));
            //rTrimmedFFTraw->Draw();
            //canvas->Print(Form("%s.pdf", outputFilename));

            double bin_i    = 0.0;
            double bin_f    = (1 / (0.1492 * 0.000001)) ;
            //double binWidth = (bin_f-bin_i)/nBins;

            TH1D * fft = new TH1D("FFT of Residuals","FFT of Residuals", int(windowBins/2), bin_i, 0.5*bin_f);
            fft->GetXaxis()->SetTitle("Frequency (Hz)");
            fft->GetYaxis()->SetTitle("Power (arb. units)");

            for (int x=1;x<=rTrimmedFFTraw->GetNbinsX();x++){
                //double freq_i = rTrimmedFFTraw->GetBinCenter(x) ; 
                //double freq_f = rTrimmedFFTraw->GetBinCenter(x) * binWidth ; 
                //int bin = ( x ) * ( 1.0 / 0.1492 ); 
                //if(x<10){printf("freq i %f -> %f bin %i bin [%f,%f]\n", freq_i, freq_f, bin, fft->GetBinLowEdge(bin), fft->GetBinLowEdge(bin+1));}
                fft->SetBinContent(x,rTrimmedFFTraw->GetBinContent(x));
            }
            double fftMaxY = fft->GetMaximum()+2000;
            fft->SetMaximum(fftMaxY);
            fft->Draw();

            double fc = 6702413;
            double fCBO = 1000000*wCBO/ (2*M_PI);
            double fy = 1000000*w_y / (2*M_PI);
            double fx = fCBO+fc;
            double fVO = fc - fy;
            double fVW = 1000000*w_vw / (2*M_PI);
            double fa = 1000000 * (myblinders.paramToFreq(R) / (2*M_PI));

            double full_fy = 1000000*full_wy / (2*M_PI);
            double full_fVW = 6702413 - (2*full_fy);

            printf("fCBO %f fVW %f fy %f full_fy %f full_fVW %f \n", fCBO, fVW, fy, full_fy, full_fVW);
            TLine *l_fx   = new TLine(fx,0.,fx,fftMaxY); 
            TLine *l_fy   = new TLine(fy,0.,fy,fftMaxY); 
            TLine *l_fa   = new TLine(fa,0.,fa,fftMaxY); 
            TLine *l_fCBO = new TLine(fCBO,0.,fCBO,fftMaxY); 
            TLine *l_2CBO = new TLine(2*fCBO,0.,2*fCBO,fftMaxY); 
            TLine *l_fVO  = new TLine(fVO ,0.,fVO ,fftMaxY);           
            TLine *l_fVW  = new TLine(fVW ,0.,fVW ,fftMaxY);    

            TLine *l_full_fy = new TLine(full_fy, 0., full_fy, fftMaxY);
            TLine *l_full_fVW = new TLine(full_fVW, 0., full_fVW, fftMaxY);

            l_fx  ->SetLineColor(kRed);
            l_fy  ->SetLineColor(kRed);
            l_fa  ->SetLineColor(kRed);
            l_fCBO->SetLineColor(kRed);
            l_2CBO->SetLineColor(kRed);
            l_fVO ->SetLineColor(kRed);
            l_fVW ->SetLineColor(kRed);

            l_full_fy->SetLineColor(kBlue);
            l_full_fVW->SetLineColor(kBlue);

            l_fx  ->SetLineStyle(2);
            l_fy  ->SetLineStyle(2);
            l_fa  ->SetLineStyle(2);
            l_fCBO->SetLineStyle(2);
            l_2CBO->SetLineStyle(2);
            l_fVO ->SetLineStyle(2);
            l_fVW ->SetLineStyle(2);

            l_full_fy->SetLineStyle(3);
            l_full_fVW->SetLineStyle(3);

            //l_fx  ->Draw("SAME");
            l_fy  ->Draw("SAME");
            l_fa  ->Draw("SAME");
            l_fCBO->Draw("SAME");
            l_2CBO->Draw("SAME");
            //l_fVO ->Draw("SAME");
            l_fVW ->Draw("SAME");

            l_full_fy->Draw("SAME");
            l_full_fVW->Draw("SAME");

            TLatex latex;
            latex.SetTextSize(0.03);
            latex.SetTextAngle(90);
            latex.SetTextColor(kRed);
            //latex.DrawLatex(-50000+ fx, 0.7*fftMaxY,"f_{x}");
            latex.DrawLatex(-50000+ fCBO, 0.7*fftMaxY,"f_{CBO}");
            latex.DrawLatex(-50000+ 2*fCBO, 0.7*fftMaxY,"f_{2CBO}");
            latex.DrawLatex(-50000+ fa  , 0.7*fftMaxY,"f_{a}");
            latex.DrawLatex(-20000+ fy, 0.5*fftMaxY,"f_{y}");
            latex.DrawLatex(-50000+ fVW, 0.7*fftMaxY,"f_{VW}");
            //latex.DrawLatex(-50000+ fVO, 0.7*fftMaxY,"f_{VO}");

            latex.DrawLatex(-50000+ full_fy, 0.7*fftMaxY,"(full) f_{y}");
            latex.DrawLatex(-50000+ full_fVW, 0.7*fftMaxY,"(full) f_{VW}");

            canvas->Print(Form("%s.pdf", outputFilename));

            windowDat->Draw();
            windowFit->SetLineColor(kRed);
            windowFit->Draw("SAME");
            canvas->Print(Form("%s.pdf", outputFilename));
        }
    }

    canvas->Print(Form("%s.pdf]", outputFilename));

}
