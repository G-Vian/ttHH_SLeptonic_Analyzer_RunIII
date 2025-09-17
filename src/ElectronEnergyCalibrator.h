#ifndef ELECTRON_ENERGY_CALIBRATOR_H
#define ELECTRON_ENERGY_CALIBRATOR_H

#include <vector>
#include <string>
#include <random>
#include <memory>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC);

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
    std::string _year;
    std::string _dataOrMC;
    std::unique_ptr<correction::CorrectionSet> cset;
    mutable std::mt19937 rng;

    std::string getElectronJSONPath() const;
    std::string getCompoundName() const;
};

#endif
