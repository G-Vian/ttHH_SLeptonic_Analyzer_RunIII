#include "ElectronEnergyCalibrator.h"
#include <stdexcept>
#include <cmath>
#include <random>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC)
    : _year(year), _dataOrMC(dataOrMC), rng(std::random_device{}())
{
    std::string jsonPath = getElectronJSONPath();
    cset = correction::CorrectionSet::from_file(jsonPath);
}

std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    if (_year == "2022") {
        return (_dataOrMC == "MC")
            ? "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/MC/electronSS_EtDependent.json.gz"
            : "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/SS/electronSS_EtDependent.json.gz";
    } 
    else if (_year == "2022EE") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoE+PromptFG/SS/electronSS_EtDependent.json.gz";
    } 
    else if (_year == "2023") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/SS/electronSS_EtDependent.json.gz";
    } 
    else if (_year == "2024") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/SS/electronSS_EtDependent_v1.json.gz";
    }
    else {
        throw std::runtime_error("Year not supported for electron corrections!");
    }
}

std::string ElectronEnergyCalibrator::getCompoundName(bool isMC) const {
    if (_year == "2022" || _year == "2022EE") {
        return isMC ? "EGMSmearAndSyst_ElePT_2022" : "EGMScale_Compound_Ele_2022preEE";
    } 
    else if (_year == "2023") {
        return isMC ? "EGMSmearAndSyst_ElePT_2023" : "EGMScale_Compound_Ele_2023";
    } 
    else if (_year == "2024") {
        return isMC ? "EGMSmearAndSyst_ElePT_2024" : "EGMScale_Compound_Ele_2024";
    } 
    else {
        throw std::runtime_error("Year not supported for compound name!");
    }
}
void ElectronEnergyCalibrator::applyElectronCalibration(
    std::vector<float>& Electron_pt,
    const std::vector<float>& Electron_eta,
    const std::vector<float>& Electron_r9,
    const std::vector<int>& Electron_seedGain,
    unsigned int runNumber,
    long long eventNumber,
    bool isMC
) {
    if (Electron_pt.size() != Electron_eta.size() ||
        Electron_pt.size() != Electron_r9.size() ||
        Electron_pt.size() != Electron_seedGain.size()) 
    {
        throw std::runtime_error("Electron vectors have different sizes!");
    }

    for (size_t i = 0; i < Electron_pt.size(); ++i) {
        double pt  = Electron_pt[i];
        double eta = Electron_eta[i];
        double r9  = Electron_r9[i];
        int gain   = Electron_seedGain[i];

        double newPt = pt;

        // Pega o nome do compound correto
        std::string compoundName = getCompoundName(isMC);
        auto compound = (*cset).compound().at(compoundName);

        if (isMC) {
            // MC: aplica smearing
            std::vector<correction::Variable::Type> args = {pt, r9, std::abs(eta)};
            newPt = pt * compound->evaluate(args);
        } else {
            // Data: aplica scale
            std::vector<correction::Variable::Type> args = {
                static_cast<int>(runNumber),
                eta,
                r9,
                pt,
                static_cast<double>(gain)
            };
            newPt = compound->evaluate(args);
        }

        Electron_pt[i] = static_cast<float>(newPt);
    }
}
