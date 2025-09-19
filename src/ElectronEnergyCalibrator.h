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
    // Construtor com parâmetros
    ElectronEnergyCalibrator(const std::string& year, const std::string& DataOrMC);

    // Função principal para aplicar a calibração
    void calibrateElectrons(
        std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains,
        int runNumber
    );

    // Limites de segurança
    float getMin(const correction::Variable& var) const;
    float getMax(const correction::Variable& var) const;

private:
    std::string _year;
    std::string _DataOrMC;
    std::unique_ptr<correction::CorrectionSet> cset;
    std::mt19937 rng;

    // Caminho do JSON
    std::string getElectronJSONPath() const;
};
