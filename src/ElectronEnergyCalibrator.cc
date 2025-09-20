#include "ElectronEnergyCalibrator.h"
#include <cmath>
#include <fstream>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year,
                                                   const std::string& DataOrMC)
    : _year(year), _DataOrMC(DataOrMC), rng(std::random_device{}())
{
    std::string jsonPath = getElectronJSONPath();
    cset = correction::CorrectionSet::from_file(jsonPath);

    // Scale para dados, Smear para MC
    _scaleEvaluator = cset->at("electron_scale");
    _smearEvaluator = cset->at("electron_smear");
}

std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    return "./ElectronScale_Smear_" + _year + ".json";  // ajuste conforme seu arquivo
}

float ElectronEnergyCalibrator::clamp(float val, float minVal, float maxVal) const {
    return std::max(minVal, std::min(maxVal, val));
}

void ElectronEnergyCalibrator::calibrateElectrons(std::vector<objectLep*>& electrons,
                                                  unsigned int runNumber,
                                                  const std::string& syst)
{
    for (auto* ele : electrons) {
        float pt  = clamp(ele->pt(), 5.f, 1000.f);
        float eta = clamp(ele->eta(), -2.5f, 2.5f);
        float r9  = clamp(ele->r9(), 0.f, 1.5f);
        float gain = clamp(ele->gain(), 0.f, 12.f);

        float newPt = pt;

        if (_DataOrMC == "DATA") {
            newPt = _scaleEvaluator->evaluate({pt, eta, r9, gain, syst});
        } else {
            newPt = _smearEvaluator->evaluate({"smear", pt, eta, r9, gain, syst});
        }

        ele->setCalibratedPt(newPt);
    }
}
