#ifndef ELECTRON_ENERGY_CALIBRATOR_H
#define ELECTRON_ENERGY_CALIBRATOR_H

#include <vector>
#include <string>
#include <random>
#include <memory>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    // Constructor: recebe caminho para JSON, tipo de dado (Data ou MC) e ano
    ElectronEnergyCalibrator(const std::string& jsonPath, const std::string& dataOrMC, int year);

    // Aplica calibração nos elétrons
    void applyElectronCalibration(
        std::vector<float>& Electron_pt,
        const std::vector<float>& Electron_eta,
        const std::vector<float>& Electron_r9,
        const std::vector<int>& Electron_seedGain,
        unsigned int runNumber,
        long long eventNumber,
        int isMC
    );

private:
    std::unique_ptr<correction::CorrectionSet> cset;
    std::string _dataOrMC;
    int _year;
    mutable std::mt19937 rng; // gerador aleatório para smearing MC
};

#endif
