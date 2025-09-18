#include "ElectronEnergyCalibrator.h"
#include <algorithm>
#include <iostream>
#include <cmath>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC)
    : _year(year), _DataOrMC(dataOrMC), rng(std::random_device{}())
{
    std::string jsonPath = getElectronJSONPath();
    cset = std::make_unique<correction::CorrectionSet>(correction::CorrectionSet::from_file(jsonPath));
}

std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    if (_year == "2022") {
        return (_DataOrMC == "MC")
            ? "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/MC/electronSS_EtDependent.json.gz"
            : "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2022EE") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoE+PromptFG/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2023") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2024") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/SS/electronSS_EtDependent_v1.json.gz";
    } else {
        throw std::runtime_error("Year not supported for electron corrections!");
    }
}

void ElectronEnergyCalibrator::calibrateElectrons(
    std::vector<float>& pts,
    const std::vector<float>& etas,
    const std::vector<float>& r9s,
    const std::vector<int>& gains,
    int runNumber,
    int eventNumber
)
{
    if (pts.empty()) return;

    for (size_t i = 0; i < pts.size(); ++i) {
        try {
            float pt_clamped   = pts[i];   // opcional: aplicar limites se desejar
            float eta_clamped  = etas[i];
            float r9_clamped   = r9s[i];
            int gain_clamped   = gains[i];

            float absEta = std::abs(eta_clamped);

            if (_DataOrMC == "DATA") {
                auto scale_corr = cset->compound().at("EGMScale_Compound_Ele_2022preEE");
                float scale = scale_corr->evaluate("scale", runNumber, eta_clamped, r9_clamped, absEta, pt_clamped, gain_clamped);
                pts[i] *= scale;
            } else if (_DataOrMC == "MC") {
                auto smear_corr = cset->at("EGMSmearAndSyst_ElePTsplit_2022preEE");
                float smear = smear_corr->evaluate("smear", pt_clamped, r9_clamped, absEta);

                std::normal_distribution<float> dist(0.0, 1.0);
                float rnd = dist(rng);

                pts[i] *= (1.0f + smear * rnd);
            }

        } catch (const std::out_of_range& e) {
            std::cerr << "[ERROR] Calibração falhou: " << e.what() << std::endl;
        }
    }
}
