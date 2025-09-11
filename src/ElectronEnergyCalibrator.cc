#include "ElectronEnergyCalibrator.h"
#include <cmath>
#include <stdexcept>
#include <variant>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(
    const std::string& jsonPath,
    const std::string& dataOrMC,
    int year
) : _dataOrMC(dataOrMC), _year(year), rng(std::random_device{}())
{
    cset = correction::CorrectionSet::from_file(jsonPath);
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
        double pt   = Electron_pt[i];
        double eta  = Electron_eta[i];
        double r9   = Electron_r9[i];
        int gain    = Electron_seedGain[i];

        double newPt = pt;

        if (isMC) {
            auto smearCorr = cset->at("ElectronSmear");
            std::vector<std::variant<int,double,std::string>> args;
            args.emplace_back(pt);
            args.emplace_back(r9);
            args.emplace_back(std::abs(eta));

            newPt = pt * correction::any_cast<double>(smearCorr->evaluate(args));

        } else {
            auto scaleCorr = cset->at("ElectronScale");
            std::vector<std::variant<int,double,std::string>> args;
            args.emplace_back(runNumber);
            args.emplace_back(eta);
            args.emplace_back(r9);
            args.emplace_back(pt);
            args.emplace_back(static_cast<double>(gain));

            newPt = correction::any_cast<double>(scaleCorr->evaluate(args));
        }

        Electron_pt[i] = static_cast<float>(newPt);
    }
}
