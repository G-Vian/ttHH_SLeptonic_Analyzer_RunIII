#ifndef ELECTRON_ENERGY_CALIBRATOR_H
#define ELECTRON_ENERGY_CALIBRATOR_H

#include <vector>
#include <string>
#include <random>
#include <memory>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    // Construtor: agora recebe ano como string e DataOrMC
    ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC);

    // Aplica calibração em vetores de elétrons
    void applyElectronCalibration(
        std::vector<float>& Electron_pt,
        const std::vector<float>& Electron_eta,
        const std::vector<float>& Electron_r9,
        const std::vector<int>& Electron_seedGain,
        unsigned int runNumber,
        long long eventNumber,
        int isMC // flag para diferenciar MC (smearing) de DATA (scale)
    );

private:
    std::unique_ptr<correction::CorrectionSet> cset;
    std::string _dataOrMC;
    std::string _year;
    mutable std::mt19937 rng;

    // Retorna o caminho do JSON ET-dependent apropriado
    std::string getElectronJSONPath() const;

    // Retorna o nome do compound usado internamente (para logs ou uso futuro)
    std::string getCompoundName() const;
};

#endif
