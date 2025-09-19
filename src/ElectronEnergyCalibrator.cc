#include "ElectronEnergyCalibrator.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

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
    if (_year == "2022") {
        return (_DataOrMC == "MC")
            ? "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/MC/electronSS_EtDependent.json.gz"
            : "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2022EE") { // pós EE
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoE+PromptFG/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2023") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2023B") { // pós BPIX
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023B/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2024") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/SS/electronSS_EtDependent_v1.json.gz";
    } else {
        throw std::runtime_error("Year not supported for electron corrections!");
    }
}

// ------------------------
// Calibração principal
// ------------------------
void ElectronEnergyCalibrator::calibrateElectrons(
        std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains,
        int runNumber
) {
    if (pts.empty()) {
        std::cerr << "[DEBUG] Nenhum elétron para calibrar." << std::endl;
        return;
    }

    auto checkArgTypes = [](const std::vector<std::variant<int,double,std::string>>& args) -> bool {
        // Checa se todos os elementos estão no tipo esperado: int/double/string
        for (auto& a : args) {
            if (!std::holds_alternative<int>(a) &&
                !std::holds_alternative<double>(a) &&
                !std::holds_alternative<std::string>(a)) {
                return false;
            }
        }
        return true;
    };

    try {
        for (size_t i = 0; i < pts.size(); ++i) {

            float pt = pts[i];
            float eta = etas[i];
            float r9 = r9s[i];
            int gain = gains[i];
            float absEta = std::abs(eta);

            std::cerr << "[DEBUG] Eletrón " << i << " pt/eta/r9/gain: "
                      << pt << "/" << eta << "/" << r9 << "/" << gain << std::endl;

            std::vector<std::variant<int,double,std::string>> args;

            if (_DataOrMC == "DATA") {
                std::string syst = "central";
                args = {syst, runNumber, eta, r9, absEta, pt, gain};
            } else {
                std::string syst = "central";
                args = {syst, pt, r9, absEta};
            }

            // ---------- Validação de tipos ----------
            if (!checkArgTypes(args)) {
                std::cerr << "[ERROR] Argumentos para evaluate() têm tipos inválidos!" << std::endl;
                continue;
            }

            // ---------- Debug ----------
            std::cerr << "[DEBUG] Args para evaluate: ";
            for (auto& a : args) std::visit([](auto&& val){ std::cerr << val << " "; }, a);
            std::cerr << std::endl;

            // ---------- Calibração ----------
            if (_DataOrMC == "DATA") {
                auto scale_corr = cset->compound().at(
                    (_year=="2022") ? "EGMScale_Compound_Ele_2022preEE" :
                    (_year=="2022EE") ? "EGMScale_Compound_Ele_2022postEE" :
                    (_year=="2023") ? "EGMScale_Compound_Ele_2023preBPIX" :
                    (_year=="2023B") ? "EGMScale_Compound_Ele_2023postBPIX" :
                    "EGMScale_Compound_Ele_2024"
                );
                double scale = scale_corr->evaluate(args);
                pts[i] *= static_cast<float>(scale);
            } else {
                std::normal_distribution<float> gauss(0.0, 1.0);
                double smear = cset->at(
                    (_year=="2022") ? "EGMSmearAndSyst_ElePTsplit_2022preEE" :
                    (_year=="2022EE") ? "EGMSmearAndSyst_ElePTsplit_2022postEE" :
                    (_year=="2023") ? "EGMSmearAndSyst_ElePTsplit_2023preBPIX" :
                    (_year=="2023B") ? "EGMSmearAndSyst_ElePTsplit_2023postBPIX" :
                    "EGMSmearAndSyst_ElePTsplit_2024"
                )->evaluate(args);
                pts[i] *= 1.0f + smear * gauss(rng);
            }

            std::cerr << "[DEBUG] Electron[" << i << "] Pt antes/depois: "
                      << pt << " / " << pts[i] << std::endl;
        }
    } catch (const std::out_of_range& e) {
        std::cerr << "[ERROR] Correção falhou: " << e.what() << std::endl;
    }
}

// ------------------------
// Limites de segurança
// ------------------------
float ElectronEnergyCalibrator::getMin(const std::string& varName) const {
    if (varName == "pt") return 5.0f;
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
