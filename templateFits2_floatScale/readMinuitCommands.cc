#include <iostream>
#include <fstream>
#include <TMinuit.h>

void myFCN(int& npar, double* deriv, double& f, double par[], int flag) {
    // Example of a simple quadratic function
    f = (par[0] - 2) * (par[0] - 2) + (par[1] - 3) * (par[1] - 3);
}

void readMinuitCommands(const char* filename, TMinuit* minuit) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        minuit->Command(line.c_str());
    }

    file.close();
}

int main() {
    TMinuit* minuit = new TMinuit;
    minuit->SetPrintLevel(0); // Suppress TMinuit printouts
    minuit->SetFCN(myFCN); // Set your fitting function
    std::cout<<"title before commands: "<<minuit->GetTitle()<<std::endl;

    // Load commands from external file
    readMinuitCommands("commands.txt", minuit);
    std::cout<<"title after commands: "<<minuit->GetTitle()<<std::endl;

    // Perform fitting
    minuit->Migrad();

    delete minuit;
    return 0;
}
