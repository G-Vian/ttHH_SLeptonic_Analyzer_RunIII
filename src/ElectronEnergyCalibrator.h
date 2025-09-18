#pragma once

#include <vector>
#include <string>
#include <random>
#include <memory>
#include <iostream>
#include <cmath>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    // Construtor: recebe ano e se é MC ou DATA
    ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC);

    // Função principal que aplica correção dependendo se é DATA ou MC
    void calibrateElectrons(
        std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains,
        int runNumber
    );

private:
    std::string _year;
    std::string _dataOrMC;
    std::mt19937 rng;
    correction::CorrectionSet cset;

    std::string getElectronJSONPath() const;
};
