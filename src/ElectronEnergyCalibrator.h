#pragma once

#include <vector>
#include <string>
#include <memory>
#include <random>

// Aqui você garante que objectLep está definido
#include "ttHHanalyzer_trigger.h"

class ElectronEnergyCalibrator {
public:
    // Construtor: define se é DATA ou MC e o ano
    ElectronEnergyCalibrator(const std::string& DataOrMC, const std::string& year);

    // Método para calibrar os elétrons
    void calibrateElectrons(std::vector<objectLep*>& electrons,
                            unsigned int runNumber,
                            const std::string& syst);

private:
    std::string _DataOrMC;
    std::string _year;
};

