#pragma once

#include <vector>
#include <string>
#include <random>
#include <memory>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    // Construtor: dataOrMC = "DATA" ou "MC"
    ElectronEnergyCalibrator(const std::string& dataOrMC, const std::string& year);

    // Calibra energia de um elétron (pT, eta, isData)
    float calibrateElectron(float pt, float eta, bool isData);

private:
    std::string _dataOrMC;
    std::string _year;

    // CorrectionSet carregado a partir do JSON
    std::unique_ptr<correction::CorrectionSet> cset;

    // RNG para smearing
    std::mt19937 rng;

    // Caminho do JSON (privado)
    std::string getElectronJSONPath() const;

    // Funções auxiliares para limites de variáveis
    float getMin(const correction::Variable& var) const;
    float getMax(const correction::Variable& var) const;
};
