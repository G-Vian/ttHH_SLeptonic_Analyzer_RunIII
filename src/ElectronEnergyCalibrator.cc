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
            // Correctionlib aceita std::vector<correction::Variable>
            std::vector<correction::Variable> args = {pt, r9, std::abs(eta)};
            double scale = smearCorr->evaluate(args); // já retorna double
            newPt = pt * scale;

        } else {
            auto scaleCorr = cset->at("ElectronScale");
            std::vector<correction::Variable> args = {runNumber, eta, r9, pt, static_cast<double>(gain)};
            newPt = scaleCorr->evaluate(args); // já retorna double
        }

        Electron_pt[i] = static_cast<float>(newPt);
    }
}
