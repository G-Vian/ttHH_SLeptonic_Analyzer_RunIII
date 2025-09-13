#include "ElectronEnergyCalibrator.h"
#include <stdexcept>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& dataOrMC)
    : _year(year), _dataOrMC(dataOrMC), rng(std::random_device{}())
{
    std::string jsonPath = getElectronJSONPath();
    cset = correction::CorrectionSet::from_file(jsonPath);
}

std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    if (_year == "2022") {
        return (_dataOrMC == "MC") 
            ? "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/MC/electronSS_EtDependent.json"
            : "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/SS/electronSS_EtDependent.json";
    } 
    else if (_year == "2022EE") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoE+PromptFG/SS/electronSS_EtDependent.json";
    } 
    else if (_year == "2023") {
        return (_dataOrMC == "MC")
            ? "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/MC/electronSS_EtDependent.json"
            : "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/SS/electronSS_EtDependent.json";
    } 
    else if (_year == "2023B" || _year == "2024") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/SS/electronSS_EtDependent.json";
    } 
    else {
        throw std::runtime_error("Year not supported for electron corrections!");
    }
}

void ElectronEnergyCalibrator::applyElectronCalibration(
    std::vector<float>& Electron_pt,
    const std::vector<float>& Electron_eta,
    const std::vector<float>& Electron_r9,
    const std::vector<int>& Electron_seedGain,
    unsigned int runNumber,
    long long eventNumber,
    int isMC
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

        if (isMC) {
            auto smearCorr = (*cset)["ElectronSmear"];
            std::vector<correction::Variable::Type> args = {pt, r9, std::abs(eta)};
            newPt = pt * smearCorr->evaluate(args);
        } else {
            auto scaleCorr = (*cset)["ElectronScale"];
            std::vector<correction::Variable::Type> args = {
                static_cast<int>(runNumber),
                eta,
                r9,
                pt,
                static_cast<double>(gain)
            };
            newPt = scaleCorr->evaluate(args);
        }

        Electron_pt[i] = static_cast<float>(newPt);
    }
}
