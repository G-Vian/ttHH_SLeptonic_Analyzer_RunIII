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

    // Função para aplicar incertezas sistemáticas (opcional - para uso futuro)
    void calibrateElectronsWithSystematics(
        std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains,
        int runNumber,
        const std::string& systematic = "nominal"  // "nominal", "scale_up", "scale_down", "smear_up", "smear_down"
    );

    // Limites de segurança
    float getMin(const std::string& varName) const;
    float getMax(const std::string& varName) const;

    // Função auxiliar para verificar se o calibrador está configurado para DATA ou MC
    bool isData() const { return _DataOrMC == "DATA"; }
    bool isMC() const { return _DataOrMC == "MC"; }

private:
    std::string _year;
    std::string _DataOrMC;
    std::unique_ptr<correction::CorrectionSet> cset;
    std::mt19937 rng;

    // Caminho do JSON
    std::string getElectronJSONPath() const;
    
    // Helpers para obter nomes dos corrections
    std::string getScaleCompoundName() const;
    std::string getSmearEvaluatorName() const;
};
