#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include "correction.h"  // sua implementação do CorrectionSet

class ElectronEnergyCalibrator {
public:
    // Construtor: recebe o caminho para o JSON de correções
    ElectronEnergyCalibrator(const std::string& path_json);

    // ========================
    // Aplicar scale correction (DATA)
    // ========================
    // pts: vetor de pT dos elétrons (será modificado)
    // etas: vetor de eta dos elétrons
    // r9s: vetor de r9 dos elétrons
    // gains: vetor de seedGain dos elétrons
    // runNumber/eventNumber: usados para o JSON
    void applyElectronScaleCorrection(
        std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains,
        int runNumber,
        int eventNumber
    );

    // ========================
    // Obter smearing (MC)
    // ========================
    // pts, etas, r9s, gains: vetores de elétrons
    // Retorna vetor de fatores de smearing nominais
    std::vector<float> getElectronSmearing(
        const std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains
    );

private:
    CorrectionSet correctionSet;  // objeto que carrega o JSON de correções
};
