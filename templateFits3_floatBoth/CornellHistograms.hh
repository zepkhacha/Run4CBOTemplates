#include <cmath>
#include "HistogramBase.hh"
#include <map>

class CornellHistograms: public HistogramBase {

  public:

    static constexpr int nEnergyBins = 80;
    static constexpr int energyBinWidth = 50; // MeV

    static constexpr int endSpectrum = 10000; // MeV
    static constexpr int spectrumBinWidth = 10; // MeV
    static const     int nSpectrumBins = std::round(endSpectrum / spectrumBinWidth);

    static constexpr int nCalos = 24;
    static constexpr int nBunches = 8; // cycle of 16 is really 2 cycles of the same 8 bunch configurations

    static constexpr int nSeeds = 100;
    static constexpr int nTimeOfFlightBins = 520;

    static constexpr double cyclotronPeriod = 0.1492; // microseconds
    static constexpr double endTime = 650.0644; // microseconds
    static const     int nTimeBins = std::round(endTime / cyclotronPeriod);

    static constexpr double usPerClockTick = 0.00125;

    static constexpr double shadowWindowDelay = 9.0; // in ct
    static constexpr double artificialDeadTime = 6.0; // for muonloss, in ct

    // variables used in lost-muon histogram
    double deadtime;
    int prevCalo;
    int prevFill;
    int prevSubrun;
    int prevRun;

    // histograms structured as vector[energy][calo]
    // note that internally calo index 0 stands for calo sum
    std::vector<TH2I*> fitHist;
    std::vector<TH2F*> pileup2;
    std::vector<TH2F*> pileup3;
    // std::vector<TH2I*> doubleN;
    // std::vector<TH2I*> doubleDeltaN;
    std::vector<TH2I*> pileupN;
    std::vector<TH2I*> pileupNDiffSq;

    // vector[energy][bunch]
    std::vector<TH2I*> bunchFitHist;
    std::vector<TH2F*> bunchPileup2;
    std::vector<TH2F*> bunchPileup3;
    // std::vector<TH2I*> bunchDoubleN;
    // std::vector<TH2I*> bunchDoubleDeltaN;
    std::vector<TH2I*> bunchPileupN;
    std::vector<TH2I*> bunchPileupNDiffSq;

    // vector[calo][bunch]
    std::vector<std::vector<TH1F*>> singleEnergySpectrum;
    std::vector<std::vector<TH1F*>> doubleEnergySpectrum;
    std::vector<std::vector<TH1F*>> tripleEnergySpectrum;

    // muon loss histogram
    TH1D* muonLoss;

    // forward declarations of constructor and inherited methods
    CornellHistograms();
    void bookHistograms(int seedIndex, int skimIndex) override;
    void fillSinglesHistograms(PositronData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) override;
    void fillDoublesHistograms(PileupData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) override;
    void fillTriplesHistograms(PileupData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) override;
    void fillLostMuonHistograms(LostMuonData& entry, LostMuonInput& lmInput, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) override;
    void writeHistograms(TFile* outputFile, int seedIndex) override;
    // void loadHistograms(TDirectory* inputFile, int seedIndex) override;
    // void addHistograms(HistogramBase* other, int seedIndex) override;

  private:

    // forward declarations of private helper functions
    template<class T> T* makeTimeEnergyHistogram(std::string name);
    TH1D* makeTimeHistogram(std::string name);
    TH1F* makeEnergyHistogram(std::string name);
    void renameWriteReset(TH1* histogram, std::string name);

};
