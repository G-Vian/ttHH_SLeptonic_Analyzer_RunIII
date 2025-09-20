#include "ElectronEnergyCalibrator.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept> // Para std::runtime_error

// ------------------------
// Construtor
// ------------------------
ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& DataOrMC)
    : _year(year), _DataOrMC(DataOrMC), rng(std::random_device{}())
{
    std::cerr << "[DEBUG] Inicializando ElectronEnergyCalibrator: "
              << _DataOrMC << " / " << _year << std::endl;

    std::string jsonPath = getElectronJSONPath();
    std::cerr << "[DEBUG] JSON path: " << jsonPath << std::endl;
    cset = correction::CorrectionSet::from_file(jsonPath);
}

// ------------------------
// Caminho do JSON
// ------------------------
std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    // Usando os caminhos corretos e mais detalhados para cada era
    if (_year == "2022") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoBCD/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2022EE") { 
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoE+PromptFG/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2023") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23C/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2023B") { 
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2024") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/SS/electronSS_EtDependent_v1.json.gz";
    } else {
        throw std::runtime_error("Year not supported for electron corrections!");
    }
}

// ------------------------
// Calibração principal (Nominal)
// ------------------------
void ElectronEnergyCalibrator::calibrateElectrons(
        std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains,
        int runNumber
) {
    if (pts.empty()) {
        return;
    }

    try {
        for (size_t i = 0; i < pts.size(); ++i) {
            float pt_original = pts[i];
            float eta = etas[i];
            float r9 = r9s[i];
            int gain = gains[i];
            float absEta = std::abs(eta);

            std::vector<std::variant<int, double, std::string>> args;

            if (_DataOrMC == "DATA") {
                // ========================
                // DATA: Scale correction
                // ========================
                args.emplace_back(std::string("scale"));
                args.emplace_back(static_cast<double>(runNumber));
                args.emplace_back(static_cast<double>(eta));
                args.emplace_back(static_cast<double>(r9));
                
                // Para 2024, a dependência redundante de abs(eta) foi removida
                if (_year != "2024") {
                    args.emplace_back(static_cast<double>(absEta));
                }
                
                args.emplace_back(static_cast<double>(pt_original));
                args.emplace_back(static_cast<double>(gain));

                auto scale_corr = cset->compound().at(getScaleCompoundName());
                double scale = scale_corr->evaluate(args);
                pts[i] *= static_cast<float>(scale);

            } else { // MC
                // ========================
                // MC: Smearing correction
                // ========================
                args.emplace_back(std::string("smear"));
                args.emplace_back(static_cast<double>(pt_original));
                args.emplace_back(static_cast<double>(r9));
                args.emplace_back(static_cast<double>(absEta));

                double smear = cset->at(getSmearEvaluatorName())->evaluate(args);
                
                // Aplica smearing estocástico
                std::normal_distribution<float> gauss(0.0, 1.0);
                float random_number = gauss(rng);
                pts[i] *= (1.0f + static_cast<float>(smear) * random_number);
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Exceção na calibração: " << e.what() << std::endl;
        // Opcional: reverter as alterações em caso de erro para não propagar pts parcialmente corrigidos.
    }
}

// ------------------------
// Limites de segurança
// ------------------------
float ElectronEnergyCalibrator::getMin(const std::string& varName) const {
    // Limite de pT ajustado para 15 GeV conforme a recomendação da documentação
    if (varName == "pt") return 15.0f;
    if (varName == "ScEta") return -2.5f;
    if (varName == "r9") return 0.0f;
    if (varName == "seedGain") return 0;
    return 0.0f;
}

float ElectronEnergyCalibrator::getMax(const std::string& varName) const {
    if (varName == "pt") return 1000.0f;
    if (varName == "ScEta") return 2.5f;
    if (varName == "r9") return 1.5f;
    if (varName == "seedGain") return 12;
    return 1.0f;
}

// ------------------------
// Helpers para nomes dos corrections
// ------------------------
std::string ElectronEnergyCalibrator::getScaleCompoundName() const {
    if (_year == "2022") return "EGMScale_Compound_Ele_2022preEE";
    if (_year == "2022EE") return "EGMScale_Compound_Ele_2022postEE";
    if (_year == "2023") return "EGMScale_Compound_Ele_2023preBPIX";
    if (_year == "2023B") return "EGMScale_Compound_Ele_2023postBPIX";
    if (_year == "2024") return "EGMScale_Compound_Ele_2024";
    throw std::runtime_error("Invalid year for scale compound name: " + _year);
}

std::string ElectronEnergyCalibrator::getSmearEvaluatorName() const {
    if (_year == "2022") return "EGMSmearAndSyst_ElePTsplit_2022preEE";
    if (_year == "2022EE") return "EGMSmearAndSyst_ElePTsplit_2022postEE";
    if (_year == "2023") return "EGMSmearAndSyst_ElePTsplit_2023preBPIX";
    if (_year == "2023B") return "EGMSmearAndSyst_ElePTsplit_2023postBPIX";
    if (_year == "2024") return "EGMSmearAndSyst_ElePTsplit_2024";
    throw std::runtime_error("Invalid year for smear evaluator name: " + _year);
}

// ------------------------
// Calibração com Sistemáticas (opcional)
// ------------------------
void ElectronEnergyCalibrator::calibrateElectronsWithSystematics(
    std::vector<float>& pts,
    const std::vector<float>& etas,
    const std::vector<float>& r9s,
    const std::vector<int>& gains,
    int runNumber,
    const std::string& systematic
) {
    if (pts.empty() || systematic == "nominal") {
        // Se for nominal, chama a função principal
        calibrateElectrons(pts, etas, r9s, gains, runNumber);
        return;
    }

    // Sistemáticas são aplicadas apenas no MC
    if (isData()) return;

    try {
        for (size_t i = 0; i < pts.size(); ++i) {
            float pt_original = pts[i];
            float eta = etas[i];
            float r9 = r9s[i];
            float absEta = std::abs(eta);

            std::vector<std::variant<int, double, std::string>> args;
            std::string eval_type;

            if (systematic == "scale_up" || systematic == "scale_down") {
                eval_type = "escale"; // Incerteza da escala
            } else if (systematic == "smear_up" || systematic == "smear_down") {
                eval_type = "esmear"; // Incerteza do smearing
            } else {
                continue; // Sistemática desconhecida
            }
            
            args.emplace_back(eval_type);
            args.emplace_back(static_cast<double>(pt_original));
            args.emplace_back(static_cast<double>(r9));
            args.emplace_back(static_cast<double>(absEta));

            double uncertainty = cset->at(getSmearEvaluatorName())->evaluate(args);

            if (eval_type == "esmear") {
                // Para smear_up/down, precisamos do smear nominal primeiro
                std::vector<std::variant<int, double, std::string>> nominal_args;
                nominal_args.emplace_back(std::string("smear"));
                nominal_args.emplace_back(static_cast<double>(pt_original));
                nominal_args.emplace_back(static_cast<double>(r9));
                nominal_args.emplace_back(static_cast<double>(absEta));
                double nominal_smear = cset->at(getSmearEvaluatorName())->evaluate(nominal_args);

                std::normal_distribution<float> gauss(0.0, 1.0);
                float random_number = gauss(rng);
                float total_smear = (systematic == "smear_up") ? (nominal_smear + uncertainty) : (nominal_smear - uncertainty);
                pts[i] *= (1.0f + total_smear * random_number);

            } else if (eval_type == "escale") {
                // Primeiro aplicamos o smearing nominal
                std::vector<std::variant<int, double, std::string>> nominal_args;
                nominal_args.emplace_back(std::string("smear"));
                nominal_args.emplace_back(static_cast<double>(pt_original));
                nominal_args.emplace_back(static_cast<double>(r9));
                nominal_args.emplace_back(static_cast<double>(absEta));
                double nominal_smear = cset->at(getSmearEvaluatorName())->evaluate(nominal_args);
                std::normal_distribution<float> gauss(0.0, 1.0);
                pts[i] *= (1.0f + nominal_smear * gauss(rng));

                // E então aplicamos a incerteza de escala sobre o pT já "smired"
                float scale_factor = (systematic == "scale_up") ? (1.0f + uncertainty) : (1.0f - uncertainty);
                pts[i] *= scale_factor;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Exceção na calibração com sistemáticas: " << e.what() << std::endl;
    }
}
