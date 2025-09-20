#pragma once

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <cmath>
#include "ttHHanalyzer_trigger.h"

// Precisamos declarar o objectLep antes de usar
// Ele está definido em ttHHanalyzer_trigger.h
#include "ttHHanalyzer_trigger.h"

class ElectronEnergyCalibrator {
public:
    // Construtor
    ElectronEnergyCalibrator(const std::string& dataOrMC, const std::string& year);

    // Método principal para calibrar elétrons
    // Modifica os pTs diretamente dentro do vetor de objetos
    void calibrateElectrons(
        std::vector<objectLep*>& electrons,
        unsigned int runNumber,
        const std::string& syst // "central", "up", "down"
    );

    // Métodos auxiliares
    float getMin(const std::string& var) const;
    float getMax(const std::string& var) const;

private:
    std::string _dataOrMC; // "DATA" ou "MC"
    std::string _year;     // ex: "2022", "2023", "2024"

    // Aqui você pode armazenar os evaluators ou calibradores reais
    // Exemplo: std::unique_ptr<correction::CorrectionSet> correction;
};
