#include "ElectronEnergyCalibrator.h"
#include <random>
#include <stdexcept>

// ========================
// Construtor
// ========================
ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC)
    : _year(year), _dataOrMC(dataOrMC), rng(std::random_device{}())
{
    std::string jsonPath = getElectronJSONPath();
    cset = correction::CorrectionSet::from_file(jsonPath);
}

// ========================
// Obter caminho do JSON
// ========================
std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    if (_year == "2022") {
        return (_dataOrMC == "MC")
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

// ========================
// Calibração de elétrons (DATA ou MC)
// ========================
void ElectronEnergyCalibrator::calibrateElectrons(
    std::vector<float>& pts,
    const std::vector<float>& etas,
    const std::vector<float>& r9s,
    const std::vector<int>& gains,
    int runNumber
)
{
    for (size_t i = 0; i < pts.size(); i++) {
        try {
            float absEta = std::abs(etas[i]);

            if (_dataOrMC != "MC") {
                // ========================
                // DATA: aplica scale correction
                // ========================
                float scale = cset.compound.at("EGMScale_Compound_Ele_2022preEE")
                                   .evaluate("scale", "Data", runNumber, etas[i], r9s[i], absEta, pts[i], gains[i]);
                pts[i] *= scale;
            } else {
                // ========================
                // MC: aplica smearing
                // ========================
                float smear = cset.at("EGMSmearAndSyst_ElePTsplit_2022preEE")
                                   .evaluate("smear", pts[i], r9s[i], absEta);
                std::normal_distribution<float> gauss(0.0, 1.0);
                float rnd = gauss(rng);
                pts[i] *= (1.0f + smear * rnd);
            }

        } catch (const std::out_of_range& e) {
            std::cerr << "[ERROR] Calibração falhou: map::at inválido (verifique pt, eta, r9 ou seedGain)" << std::endl;
        }
    }
}
