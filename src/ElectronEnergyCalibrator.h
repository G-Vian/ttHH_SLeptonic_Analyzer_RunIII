#ifndef TTHHANALYZER_TRIGGER_H
#define TTHHANALYZER_TRIGGER_H

#include "ElectronEnergyCalibrator.h"
#include <vector>
#include <string>

class ttHHanalyzer {
public:
    ttHHanalyzer(const std::string& inputFile,
                 eventBuffer* buffer,
                 float lumi,
                 bool isMC,
                 std::string year,
                 std::string DataOrMC,
                 std::string sampleName);

    void analyze(event* thisEvent);

private:
    ElectronEnergyCalibrator calibrator;

    // Exemplo de vetores (se quiser usar diretos do evento):
    std::vector<float> Electron_pt;
    std::vector<float> Electron_eta;
    std::vector<float> Electron_r9;
    std::vector<int>   Electron_seedGain;

    bool _isMC;
    std::string _year;
    std::string _DataOrMC;
};

#endif
