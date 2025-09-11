#include "ElectronEnergyCalibrator.h"

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& jsonPath, const std::string& dataOrMC, int year)
    : _dataOrMC(dataOrMC), _year(year) 
{
    cset = correction::CorrectionSet::from_file(jsonPath);
}

void ElectronEnergyCalibrator::applyElectronCalibration(
    std::vector<float>& Electron_pt,
    const std::vector<float>& Electron_eta,
    const std::vector<float>& Electron_r9,
    const std::vector<int>& Electron_seedGain,
    unsigned int runNumber,
    long long eventNumber
) {
    if (Electron_pt.empty()) return;

    for (size_t i = 0; i < Electron_pt.size(); ++i) {
        float pt  = Electron_pt[i];
        float eta = Electron_eta[i];
        float r9  = Electron_r9[i];
        int gain  = Electron_seedGain[i];

        if (pt < 15.0) continue;  // Correções não são válidas abaixo de ~15 GeV

        if (_dataOrMC == "MC") {
            auto smearCorr = cset->at("EGMSmearAndSyst_ElePTsplit_2023preBPIX");
            float rho = smearCorr.evaluate({"rho", pt, r9, std::abs(eta)});

            std::mt19937 rng_event(eventNumber); // usa event como seed
            std::normal_distribution<float> gaus(1.0, rho);
            float smear = gaus(rng_event);

            Electron_pt[i] *= smear;
        } 
        else if (_dataOrMC == "Data") {
            auto scaleCorr = cset->at("EGMScale_Compound_Ele_2023preBPIX");
            float scale = scaleCorr.evaluate({"scale", runNumber, eta, r9, pt, (double)gain});

            Electron_pt[i] *= scale;
        }
    }
}
