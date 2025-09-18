#include "ElectronEnergyCalibrator.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC)
    : _year(year), _dataOrMC(dataOrMC), rng(std::random_device{}())
{
    std::string jsonPath = getElectronJSONPath();
    cset = std::make_unique<correction::CorrectionSet>(correction::CorrectionSet::from_file(jsonPath));
}

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
        return (_dataOrMC == "MC")
            ? "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/MC_B/electronSS_EtDependent.json.gz"
            : "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/SS_B/electronSS_EtDependent.json.gz";
    } else if (_year == "2024") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/SS/electronSS_EtDependent_v1.json.gz";
    } else {
        throw std::runtime_error("Year not supported for electron corrections!");
    }
}

float ElectronEnergyCalibrator::getMin(const std::string& var) const {
    // Exemplos de limites (ajuste conforme JSON)
    if (var == "pt") return 5.0f;
    if (var == "ScEta") return -2.5f;
    if (var == "r9") return 0.0f;
    if (var == "seedGain") return 0;
    return 0.0f;
}

float ElectronEnergyCalibrator::getMax(const std::string& var) const {
    if (var == "pt") return 1000.0f;
    if (var == "ScEta") return 2.5f;
    if (var == "r9") return 1.0f;
    if (var == "seedGain") return 12;
    return 0.0f;
}

void ElectronEnergyCalibrator::calibrateElectrons(std::vector<float>& pts,
                                                  const std::vector<float>& etas,
                                                  const std::vector<float>& r9s,
                                                  const std::vector<int>& gains,
                                                  int runNumber)
{
    if (pts.empty()) return;

    for (size_t i = 0; i < pts.size(); i++) {
        float pt_clamped   = std::clamp(pts[i], getMin("pt"), getMax("pt"));
        float eta_clamped  = std::clamp(etas[i], getMin("ScEta"), getMax("ScEta"));
        float r9_clamped   = std::clamp(r9s[i], getMin("r9"), getMax("r9"));
        int gain_clamped   = std::clamp(gains[i], (int)getMin("seedGain"), (int)getMax("seedGain"));

        try {
            if (_dataOrMC == "DATA") {
                // Scale correction para DATA
                auto scale_corr = cset->compound.at(
                    (_year == "2022") ? "EGMScale_Compound_Ele_2022preEE" :
                    (_year == "2023") ? "EGMScale_Compound_Ele_2023preBPIX" :
                    (_year == "2023B") ? "EGMScale_Compound_Ele_2023B" :
                    "EGMScale_Compound_Ele_2024"
                );
                std::vector<std::variant<int, double, std::string>> inputs;
                if (_year == "2024") {
                    inputs = {runNumber, eta_clamped, r9_clamped, pt_clamped, gain_clamped};
                } else {
                    inputs = {runNumber, eta_clamped, r9_clamped, std::abs(eta_clamped), pt_clamped, gain_clamped};
                }
                double scale = scale_corr->evaluate(inputs);
                pts[i] *= scale;
            } else {
                // Smearing para MC
                auto smear_corr = cset->at(
                    (_year == "2022") ? "EGMSmearAndSyst_ElePTsplit_2022preEE" :
                    (_year == "2023") ? "EGMSmearAndSyst_ElePTsplit_2023preBPIX" :
                    (_year == "2023B") ? "EGMSmearAndSyst_ElePTsplit_2023B" :
                    "EGMSmearAndSyst_ElePTsplit_2024"
                );
                std::vector<std::variant<int, double, std::string>> inputs = {pt_clamped, r9_clamped, std::abs(eta_clamped)};
                double smear = smear_corr->evaluate(inputs);

                std::normal_distribution<double> gauss(0.0, 1.0);
                double rnd = gauss(rng);
                pts[i] *= (1.0 + smear * rnd);
            }
        } catch (const std::out_of_range& e) {
            std::cerr << "[ERROR] Correção falhou: " << e.what() << std::endl;
        }
    }
}
