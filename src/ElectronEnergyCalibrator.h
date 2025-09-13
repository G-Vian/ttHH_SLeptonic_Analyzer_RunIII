#ifndef ELECTRON_ENERGY_CALIBRATOR_H
#define ELECTRON_ENERGY_CALIBRATOR_H

#include <vector>
#include <string>
#include <random>
#include <memory>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    // Construtor: dataOrMC = "DATA" ou "MC", year = 2022, 2023, 2024
    ElectronEnergyCalibrator(const std::string& dataOrMC, int year);

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
    int _year;
    mutable std::mt19937 rng;

    std::string getElectronJSONPath() const;
};

#endif
