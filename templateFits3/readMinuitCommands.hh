#include <iostream>
#include <fstream>
#include <TMinuit.h>

void readMinuitCommands(const char* filename, TMinuit *minuit) {
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
