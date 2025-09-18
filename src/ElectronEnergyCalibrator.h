#pragma once

#include <vector>
#include <string>
#include <random>
#include <memory>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    // Construtor: ano e tipo de dado ("DATA" ou "MC")
    ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC);

    // Aplica calibração de pT
    void applyElectronCalibration(
        std::vector<float>& Electron_pt,
        const std::vector<float>& Electron_eta,
        const std::vector<float>& Electron_r9,
        const std::vector<int>& Electron_seedGain,
        unsigned int runNumber,
        long long eventNumber,
        bool isMC
    );

private:
    std::unique_ptr<correction::CorrectionSet> cset; // Corrections carregadas
    std::string _year;
    std::string _dataOrMC;
    mutable std::mt19937 rng; // Para smearing MC

    std::string getElectronJSONPath() const;
    std::string getCompoundName(bool isMC) const;
};
