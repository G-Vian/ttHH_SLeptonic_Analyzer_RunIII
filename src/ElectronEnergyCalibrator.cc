#include "ElectronEnergyCalibrator.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC)
    : _year(year), _DataOrMC(dataOrMC), rng(std::random_device{}())
{
    std::string jsonPath = getElectronJSONPath();
    try {
        cset = std::make_unique<correction::CorrectionSet>(correction::CorrectionSet::from_file(jsonPath));
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Falha ao carregar JSON do calibrador: " << e.what() << std::endl;
        throw;
    }
}

std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    // Ajuste o caminho para o JSON conforme seu setup
    if (_DataOrMC == "MC") {
        return "ElectronCalibration_MC_" + _year + ".json";
    } else {
        return "ElectronCalibration_DATA_" + _year + ".json";
    }
}

float ElectronEnergyCalibrator::getMin(const std::string& var) const {
    if (var == "pt") return 5.f;
    if (var == "ScEta") return -2.5f;
    if (var == "r9") return 0.f;
    if (var == "seedGain") return 0;
    throw std::runtime_error("getMin: variável desconhecida " + var);
}

float ElectronEnergyCalibrator::getMax(const std::string& var) const {
    if (var == "pt") return 500.f;
    if (var == "ScEta") return 2.5f;
    if (var == "r9") return 1.0f;
    if (var == "seedGain") return 12;
    throw std::runtime_error("getMax: variável desconhecida " + var);
}

void ElectronEnergyCalibrator::calibrateElectrons(
    std::vector<float>& pts,
    const std::vector<float>& etas,
    const std::vector<float>& r9s,
    const std::vector<int>& gains,
    int runNumber
) {
    if (!cset) throw std::runtime_error("calibrateElectrons: calibrador não inicializado!");

    std::string syst = "nominal"; // ⚠️ string obrigatória para correctionlib

    for (size_t i = 0; i < pts.size(); ++i) {
        try {
            auto corr = (*cset)["EleEnergyCorrection"]; // nome do correction no JSON
            pts[i] = corr->evaluate({syst, _year, pts[i], etas[i], r9s[i], gains[i]});
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] Falha ao calibrar elétron " << i << ": " << e.what() << std::endl;
        }
    }
}
