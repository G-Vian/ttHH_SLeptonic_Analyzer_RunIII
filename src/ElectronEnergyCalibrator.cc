#include "ElectronEnergyCalibrator.h"
#include <iostream>
#include <cmath>
#include <random>
#include <map>
#include <variant>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <memory>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC)
    : _year(year), _DataOrMC(dataOrMC), rng(std::random_device{}())
{
    std::string jsonPath = getElectronJSONPath();

    try {
        // Carregar JSON diretamente como CorrectionSet
        cset = correction::CorrectionSet::from_file(jsonPath);
        std::cout << "[INFO] ElectronEnergyCalibrator carregado: "
                  << (_DataOrMC.empty() ? "UNKNOWN" : _DataOrMC)
                  << ", ano: " << (_year.empty() ? "UNKNOWN" : _year) << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Falha ao carregar JSON: " << e.what() << std::endl;
    }
}

// Função que retorna o caminho correto do JSON dependendo do ano e MC/DATA
std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    std::string basePath = "/eos/user/g/gvian/corrections/electron/";
    std::string filename;

    if (_DataOrMC == "MC") {
        if (_year == "2022") filename = "ElectronEnergyCorrections_MC_2022.json";
        else if (_year == "2023") filename = "ElectronEnergyCorrections_MC_2023.json";
        else if (_year == "2024") filename = "ElectronEnergyCorrections_MC_2024.json";
    } else if (_DataOrMC == "DATA") {
        if (_year == "2022") filename = "ElectronEnergyCorrections_DATA_2022.json";
        else if (_year == "2023") filename = "ElectronEnergyCorrections_DATA_2023.json";
        else if (_year == "2024") filename = "ElectronEnergyCorrections_DATA_2024.json";
    }

    if (filename.empty()) {
        throw std::runtime_error("JSON para ElectronEnergyCalibrator não encontrado para "
                                 + _DataOrMC + " " + _year);
    }

    return basePath + filename;
}

// Função de calibração de elétrons
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
        std::map<std::string, std::variant<std::string, double, int>> inputs;

        // Correção de tipo: syst sempre string
        inputs["syst"] = std::string("Nominal");
        inputs["pt"] = static_cast<double>(pts[i]);
        inputs["eta"] = static_cast<double>(etas[i]);
        inputs["r9"] = static_cast<double>(r9s[i]);
        inputs["gain"] = static_cast<int>(gains[i]);
        inputs["run"] = runNumber;

        try {
            const auto& eleCorr = cset->at("EleEnergyCorr");
            double corr = eleCorr->evaluate(inputs);
            pts[i] *= corr;
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] Calibração do elétron " << i << " falhou: "
                      << e.what() << std::endl;
        }
    }
}
