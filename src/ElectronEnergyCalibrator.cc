#include "ElectronEnergyCalibrator.h"
#include <algorithm>
#include <cmath>
#include <iostream>

// ------------------------
// Construtor
// ------------------------
ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC)
    : _year(year), _dataOrMC(dataOrMC), rng(std::random_device{}())
{
    std::string jsonPath = getElectronJSONPath();
    cset = correction::CorrectionSet::from_file(jsonPath);
}

// ------------------------
// Caminho do JSON
// ------------------------
std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    if (_year == "2022") {
        return (_dataOrMC == "MC")
            ? "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/MC/electronSS_EtDependent.json.gz"
            : "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2022EE") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoE+PromptFG/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2023") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2023B") {
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
    if (pts.empty()) return;

    try {
        for (size_t i = 0; i < pts.size(); ++i) {

            float pt = pts[i];
            float eta = etas[i];
            float r9 = r9s[i];
            int gain = gains[i];
            float absEta = std::abs(eta);

            std::vector<std::variant<int,double,std::string>> args;

            if (_dataOrMC == "DATA") {
                // ========================
                // DATA: Compound correction
                // ========================
                std::string syst = "nominal"; // sistemática padrão
                args = { syst, runNumber, eta, r9, absEta, pt, gain };

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
                // ========================
                // MC: Smearing + sistemática
                // ========================
                std::string syst = "nominal"; // sempre string

                auto smear_corr = cset->at(
                    (_year=="2022") ? "EGMSmearAndSyst_ElePTsplit_2022preEE" :
                    (_year=="2022EE") ? "EGMSmearAndSyst_ElePTsplit_2022postEE" :
                    (_year=="2023") ? "EGMSmearAndSyst_ElePTsplit_2023preBPIX" :
                    (_year=="2023B") ? "EGMSmearAndSyst_ElePTsplit_2023postBPIX" :
                    "EGMSmearAndSyst_ElePTsplit_2024"
                );

                // ⚠️ Ajuste: string syst deve vir primeiro
                args = { syst, pt, r9, absEta };

                double smear = smear_corr->evaluate(args);

                std::normal_distribution<float> gauss(0.0, 1.0);
                pts[i] *= 1.0f + smear * gauss(rng);
            }
        }
    } catch (const std::out_of_range& e) {
        std::cerr << "[ERROR] Correção falhou: " << e.what() << std::endl;
    }
}

// ------------------------
// Limites de segurança
// ------------------------
float ElectronEnergyCalibrator::getMin(const std::string& var) const {
    if (var=="pt") return 5.0f;
    if (var=="ScEta") return -2.5f;
    if (var=="r9") return 0.0f;
    if (var=="seedGain") return 0;
    return 0.0f;
}

float ElectronEnergyCalibrator::getMax(const std::string& var) const {
    if (var=="pt") return 1000.0f;
    if (var=="ScEta") return 2.5f;
    if (var=="r9") return 1.5f;
    if (var=="seedGain") return 12;
    return 1.0f;
}
