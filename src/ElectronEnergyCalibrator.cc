#include "ElectronEnergyCalibrator.h"
#include <cmath>
#include <iostream>
#include <variant>
#include <map>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year,
                                                   const std::string& dataOrMC)
    : _year(year), _DataOrMC(dataOrMC), rng(std::random_device{}()) 
{
    std::string jsonPath = getElectronJSONPath();
    try {
        cset = std::make_unique<correction::CorrectionSet>(correction::CorrectionSet::from_file(jsonPath));
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Falha ao carregar JSON: " << jsonPath << " | " << e.what() << std::endl;
        cset = nullptr;
    }
}

std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    std::string basePath = "data/ElectronSF/";

    if (_DataOrMC.empty() || _year.empty()) {
        std::cerr << "[WARNING] _DataOrMC ou _year vazio. Usando arquivo default." << std::endl;
        return basePath + "ElectronSF_default.json";
    }

    // Arquivo esperado: ElectronSF_<year>_<DataOrMC>.json
    std::string fileName = "ElectronSF_" + _year + "_" + _DataOrMC + ".json";
    return basePath + fileName;
}

float ElectronEnergyCalibrator::getMin(const std::string& var) const {
    // Ajuste conforme os limites desejados ou lidos do JSON
    if (var == "pt") return 5.0f;
    if (var == "ScEta") return -2.5f;
    if (var == "r9") return 0.0f;
    if (var == "seedGain") return 0;
    return 0.0f;
}

float ElectronEnergyCalibrator::getMax(const std::string& var) const {
    if (var == "pt") return 200.0f;
    if (var == "ScEta") return 2.5f;
    if (var == "r9") return 1.0f;
    if (var == "seedGain") return 12;
    return 0.0f;
}

void ElectronEnergyCalibrator::calibrateElectrons(
    std::vector<float>& pts,
    const std::vector<float>& etas,
    const std::vector<float>& r9s,
    const std::vector<int>& gains,
    int runNumber
) {
    if (!cset) {
        std::cerr << "[ERROR] Calibrator não inicializado, pulando calibração." << std::endl;
        return;
    }

    auto eleCorr = cset->at("ElectronCalibration");
    if (!eleCorr) {
        std::cerr << "[ERROR] Não encontrou correção 'ElectronCalibration' no JSON" << std::endl;
        return;
    }

    for (size_t i = 0; i < pts.size(); i++) {
        std::map<std::string, std::variant<std::string,double,int>> inputs;
        inputs["pt"] = static_cast<double>(pts[i]);
        inputs["eta"] = static_cast<double>(etas[i]);
        inputs["r9"] = static_cast<double>(r9s[i]);
        inputs["gain"] = gains[i];
        inputs["syst"] = std::string("nominal"); // ⚠️ obrigatório ser string

        try {
            double corr = eleCorr->evaluate(inputs);
            pts[i] *= static_cast<float>(corr);
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] Falha ao avaliar correção do elétron " << i
                      << ": " << e.what() << std::endl;
        }
    }
}
