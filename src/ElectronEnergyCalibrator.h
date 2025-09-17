#pragma once

#include <vector>
#include <string>
#include <random>
#include <memory>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    // Construtor: define ano e Data/MC
    ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC);

    // Aplica a calibração em elétrons
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
    std::unique_ptr<correction::CorrectionSet> cset;
    std::string _dataOrMC;
    std::string _year;
    mutable std::mt19937 rng;

    // Métodos auxiliares
    std::string getElectronJSONPath() const;        // retorna o caminho ET-dependent correto
    std::string getCompoundName(bool isMC) const;   // retorna o nome do compound dentro do JSON
};
