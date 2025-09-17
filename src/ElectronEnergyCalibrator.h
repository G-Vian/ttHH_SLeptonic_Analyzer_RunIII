#pragma once

#include <vector>
#include <string>
#include <random>
#include <memory>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    // Construtor: agora recebe apenas o ano e se é Data ou MC
    ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC);

    // Aplica a calibração de pT dos elétrons
    // - Electron_pt, Electron_eta, Electron_r9, Electron_seedGain: vetores de elétrons do evento
    // - runNumber, eventNumber: informações do evento
    // - isMC: true se o evento é MC, false se é dado
    void applyElectronCalibration(
        std::vector<float>& Electron_pt,
        const std::vector<float>& Electron_eta,
        const std::vector<float>& Electron_r9,
        const std::vector<int>& Electron_seedGain,
        unsigned int runNumber,
        long long eventNumber,
        bool isMC
    );

private:
    std::unique_ptr<correction::CorrectionSet> cset; // Corrections carregadas
    std::string _dataOrMC;   // "DATA" ou "MC"
    std::string _year;       // Ano do dataset
    mutable std::mt19937 rng; // Para geração de números aleatórios se necessário

    // Retorna o path correto do JSON ET-dependent, baseado no ano e se é Data ou MC
    std::string getElectronJSONPath() const;

    // Retorna o nome do compound dentro do JSON dependendo se é MC ou DATA
    std::string getCompoundName(bool isMC) const;
};
