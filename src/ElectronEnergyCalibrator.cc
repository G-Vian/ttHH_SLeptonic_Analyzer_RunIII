#include "ElectronEnergyCalibrator.h"
#include <cmath>
#include <stdexcept>

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
        std::vector<correction::Variable::Type> args;
        args.push_back(pt);
        args.push_back(r9);
        args.push_back(std::abs(eta));
        newPt = pt * smearCorr->evaluate(args);
    
    } else {
        auto scaleCorr = cset->at("ElectronScale");
        std::vector<correction::Variable::Type> args;
        args.push_back(static_cast<int>(runNumber)); // explicit int
        args.push_back(eta);
        args.push_back(r9);
        args.push_back(pt);
        args.push_back(static_cast<double>(gain));   // explicit double
        newPt = scaleCorr->evaluate(args);
    }

        Electron_pt[i] = static_cast<float>(newPt);
    }
}
