#include "ElectronEnergyCalibrator.h"
#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <algorithm>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& DataOrMC)
    : _year(year), _DataOrMC(DataOrMC)
{
    // Inicializa o gerador de números aleatórios
    std::random_device rd;
    rng = std::mt19937(rd());

    // Recupera o caminho do JSON correto
    std::string jsonPath = getElectronJSONPath();
    try {
        // Inicialização correta do CorrectionSet
        cset = correction::CorrectionSet::from_file(jsonPath);

        if (!cset) {
            throw std::runtime_error("Falha ao carregar CorrectionSet do arquivo: " + jsonPath);
        }

        std::cout << "[INFO] ElectronEnergyCalibrator inicializado com sucesso." << std::endl;
        std::cout << "       Year: " << _year << " | Type: " << _DataOrMC << std::endl;
        std::cout << "       JSON: " << jsonPath << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Não foi possível inicializar CorrectionSet: " << e.what() << std::endl;
        throw;
    }
}

// Aplica a calibração aos elétrons
void ElectronEnergyCalibrator::calibrateElectrons(
    std::vector<float>& pts,
    const std::vector<float>& etas,
    const std::vector<float>& r9s,
    const std::vector<int>& gains,
    int runNumber
)
{
    if (!cset) {
        throw std::runtime_error("CorrectionSet não inicializado!");
    }

    if (pts.size() != etas.size() || pts.size() != r9s.size() || pts.size() != gains.size()) {
        throw std::runtime_error("Vetores de entrada com tamanhos diferentes.");
    }

    try {
        // Itera sobre todos os elétrons
        for (size_t i = 0; i < pts.size(); i++) {
            // Cria o objeto de correção do elétron
            auto eleCorr = (*cset)["electron_pt"];
            
            // Monta o mapa de variáveis esperado pelo CorrectionSet
            std::map<std::string, std::variant<std::string, double, int>> inputs;
            inputs["pt"]      = pts[i];
            inputs["eta"]     = etas[i];
            inputs["r9"]      = r9s[i];
            inputs["gain"]    = gains[i];
            inputs["year"]    = _year;        // ⚠️ string
            inputs["syst"]    = "nominal";    // ⚠️ string
            inputs["isMC"]    = (_DataOrMC == "MC" ? 1 : 0);

            // Aplica a correção
            double corr = eleCorr->evaluate(inputs);

            // Atualiza pt do elétron
            pts[i] *= corr;
        }

    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Falha na calibração dos elétrons: " << e.what() << std::endl;
        throw;
    }
}

// Retorna valor mínimo de uma variável
float ElectronEnergyCalibrator::getMin(const std::string& var) const {
    if (var == "pt") return 5.0;
    if (var == "ScEta") return -2.5;
    if (var == "r9") return 0.0;
    if (var == "seedGain") return 0;
    throw std::runtime_error("Variável desconhecida no getMin: " + var);
}

// Retorna valor máximo de uma variável
float ElectronEnergyCalibrator::getMax(const std::string& var) const {
    if (var == "pt") return 1000.0;
    if (var == "ScEta") return 2.5;
    if (var == "r9") return 1.0;
    if (var == "seedGain") return 12;
    throw std::runtime_error("Variável desconhecida no getMax: " + var);
}

// Retorna caminho do JSON com base no ano e tipo
std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    std::string baseDir = "/caminho/para/JSONs/"; // ajuste conforme seu setup
    std::string fileName;

    if (_DataOrMC == "MC") {
        if (_year == "2022") fileName = "EleCalib_2022_MC.json";
        else if (_year == "2023") fileName = "EleCalib_2023_MC.json";
        else fileName = "EleCalib_default_MC.json";
    } else {
        if (_year == "2022") fileName = "EleCalib_2022_DATA.json";
        else if (_year == "2023") fileName = "EleCalib_2023_DATA.json";
        else fileName = "EleCalib_default_DATA.json";
    }

    return baseDir + fileName;
}
