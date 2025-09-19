#pragma once

#include <vector>
#include <string>
#include <memory>
#include <random>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    ElectronEnergyCalibrator(const std::string& DataOrMC, const std::string& year);

    void calibrateElectrons(std::vector<float>& pts,
                            const std::vector<float>& etas,
                            const std::vector<float>& r9s,
                            const std::vector<int>& gains,
                            int runNumber);

private:
    std::string _year;
    std::string _DataOrMC;
    std::unique_ptr<correction::CorrectionSet> cset;

    // Gerador para smearing gaussiano
    std::mt19937 rng;

    // Funções auxiliares para pegar limites de variáveis
    float getMin(const std::string& var) const;
    float getMax(const std::string& var) const;
};
