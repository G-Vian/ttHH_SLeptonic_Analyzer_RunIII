#pragma once

#include <string>
#include <vector>
#include <memory>
#include <variant>
#include <random>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC);

    // Aplica a calibração nos vetores de elétrons
    void calibrateElectrons(
        std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains,
        int runNumber
    );

    // Recupera os limites válidos do JSON
    float getMin(const std::string& variable) const;
    float getMax(const std::string& variable) const;

private:
    std::string _year;
    std::string _dataOrMC;

    std::unique_ptr<correction::CorrectionSet> cset;
    std::mt19937 rng;

    std::string getElectronJSONPath() const;
};
