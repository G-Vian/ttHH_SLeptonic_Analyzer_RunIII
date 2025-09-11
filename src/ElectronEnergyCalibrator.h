#ifndef ELECTRON_ENERGY_CALIBRATOR_H
#define ELECTRON_ENERGY_CALIBRATOR_H

#include <vector>
#include <string>
#include <random>
#include <fstream>
#include <memory>
#include <correction.h>
#include "ttHHanalyzer_trigger.h"  // para LOG_INTERVAL, _DataOrMC, _year etc.

class ElectronEnergyCalibrator {
public:
    ElectronEnergyCalibrator(const std::string& jsonPath, const std::string& dataOrMC, int year);
    void applyElectronCalibration(
        std::vector<float>& Electron_pt,
        const std::vector<float>& Electron_eta,
        const std::vector<float>& Electron_r9,
        const std::vector<int>& Electron_seedGain,
        unsigned int runNumber,
        long long eventNumber,
        std::ofstream& logFile
    );

private:
    std::unique_ptr<correction::CorrectionSet> cset;
    std::string _dataOrMC;
    int _year;
    mutable std::mt19937 rng;
};

#endif
