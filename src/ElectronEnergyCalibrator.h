#pragma once
#include "ttHHanalyzer_trigger.h"  // Para objectLep

#include <vector>
#include <string>
#include <random>
#include <iostream>
#include "correction.h"             // correctionlib

class ElectronEnergyCalibrator {
public:
    ElectronEnergyCalibrator(const std::string& year, const std::string& DataOrMC);

    // Função principal de calibração
    void calibrateElectrons(std::vector<objectLep*>& electrons,
                            unsigned int runNumber,
                            const std::string& syst = "central");

private:
    std::string _year;
    std::string _DataOrMC;

    // correctionlib objects
    std::unique_ptr<correction::CorrectionSet> cset;
    std::shared_ptr<correction::Correction> _scaleEvaluator;
    std::shared_ptr<correction::Correction> _smearEvaluator;

    std::mt19937 rng;

    // Helpers
    std::string getElectronJSONPath() const;
    float clamp(float val, float minVal, float maxVal) const;
};
