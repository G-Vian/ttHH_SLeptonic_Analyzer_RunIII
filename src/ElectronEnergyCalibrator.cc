#include "ElectronEnergyCalibrator.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC)
    : _year(year), _DataOrMC(dataOrMC), rng(std::random_device{}())
{
    // Caminho do JSON dependendo do ano e DATA/MC
    std::string jsonPath = getElectronJSONPath();

    try {
        // Inicializa o CorrectionSet corretamente
        cset = std::unique_ptr<correction::CorrectionSet>(
            new correction::CorrectionSet(correction::CorrectionSet::from_file(jsonPath))
        );
        std::cout << "[INFO] Carregado ElectronEnergyCalibrator: " 
                  << (_DataOrMC.empty() ? "UNKNOWN" : _DataOrMC)
                  << ", ano: " << (_year.empty() ? "UNKNOWN" : _year)
                  << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Falha ao carregar JSON do ElectronEnergyCalibrator: "
                  << e.what() << std::endl;
        cset = nullptr;
    }
}

std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    // Ajuste este caminho conforme sua organização de arquivos
    std::string baseDir = "/eos/user/g/gvian/ElectronCalibrations/";
    std::string type = (_DataOrMC == "DATA") ? "data" : "mc";

    if (_year == "2024") return baseDir + type + "_2024.json";
    if (_year == "2023") return baseDir + type + "_2023.json";
    if (_year == "2022") return baseDir + type + "_2022.json";

    throw std::runtime_error("Ano desconhecido para ElectronEnergyCalibrator: " + _year);
}

float ElectronEnergyCalibrator::getMin(const std::string& var) const {
    if (var == "pt") return 5.0;
    if (var == "ScEta") return -2.5;
    if (var == "r9") return 0.0;
    if (var == "seedGain") return 0;
    throw std::runtime_error("getMin: variável desconhecida: " + var);
}

float ElectronEnergyCalibrator::getMax(const std::string& var) const {
    if (var == "pt") return 2000.0;
    if (var == "ScEta") return 2.5;
    if (var == "r9") return 1.0;
    if (var == "seedGain") return 12;
    throw std::runtime_error("getMax: variável desconhecida: " + var);
}

void ElectronEnergyCalibrator::calibrateElectrons(
    std::vector<float>& pts,
    const std::vector<float>& etas,
    const std::vector<float>& r9s,
    const std::vector<int>& gains,
    int runNumber
) {
    if (!cset) {
        std::cerr << "[ERROR] calibrateElectrons chamado sem JSON carregado!" << std::endl;
        return;
    }

    if (pts.size() != etas.size() || pts.size() != r9s.size() || pts.size() != gains.size()) {
        throw std::runtime_error("Vetores de elétrons com tamanhos inconsistentes!");
    }

    for (size_t i = 0; i < pts.size(); i++) {
        // Preparar inputs para correctionlib
        std::map<std::string, std::variant<std::string, double, int>> inputs;

        // Correção de tipo: syst sempre string
        inputs["syst"] = std::string("Nominal");
        inputs["pt"] = static_cast<double>(pts[i]);
        inputs["eta"] = static_cast<double>(etas[i]);
        inputs["r9"] = static_cast<double>(r9s[i]);
        inputs["gain"] = static_cast<int>(gains[i]);
        inputs["run"] = runNumber;

        try {
            auto eleCorr = cset->at("EleEnergyCorr");
            double corr = eleCorr->evaluate(inputs);
            pts[i] *= corr;
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] Calibração do elétron " << i << " falhou: " 
                      << e.what() << std::endl;
        }
    }
}
