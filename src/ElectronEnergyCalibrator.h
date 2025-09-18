#pragma once

#include <vector>
#include <string>
#include <memory>
#include <random>
#include <variant>
#include <iostream>

#include "correction.h"

class ElectronEnergyCalibrator {
public:
    ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC);

    // Função principal para aplicar a calibração
    void calibrateElectrons(
        std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains,
        int runNumber
    );

    float getMin(const std::string& var) const;
    float getMax(const std::string& var) const;

private:
 //   std::string _year;
 //   std::string _dataOrMC;
    std::unique_ptr<correction::CorrectionSet> cset;
    std::mt19937 rng;

    std::string getElectronJSONPath() const;
};
