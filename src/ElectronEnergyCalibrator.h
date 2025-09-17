#pragma once

#include <vector>
#include <string>
#include <random>
#include <memory>
#include <correction.h>
#include <stdexcept>

class ElectronEnergyCalibrator {
public:
    // Construtor: dataOrMC = "DATA" ou "MC", year = "2022", "2023", "2024", etc.
    ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC);

    // Aplica calibração nos vetores de elétrons
    void applyElectronCalibration(
        std::vector<float>& Electron_pt,
        const std::vector<float>& Electron_eta,
        const std::vector<float>& Electron_r9,
        const std::vector<int>& Electron_seedGain,
        unsigned int runNumber,
        long long eventNumber,
        int isMC // 1 para MC (smearing), 0 para DATA (scale)
    );

private:
    std::unique_ptr<correction::CorrectionSet> cset;
    std::string _dataOrMC;
    std::string _year;          // agora é string
    mutable std::mt19937 rng;

    // Retorna o caminho absoluto do JSON dependendo do ano e data/MC
    std::string getElectronJSONPath() const;
};


