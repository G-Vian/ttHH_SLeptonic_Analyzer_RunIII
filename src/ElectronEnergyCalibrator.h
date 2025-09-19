#pragma once

#include <vector>
#include <string>
#include <memory>
#include <correction.h>

class ElectronEnergyCalibrator {
public:
    // Construtor
    ElectronEnergyCalibrator(const std::string& year, const std::string& DataOrMC);

    // Função de calibração de elétrons
    void calibrateElectrons(
        std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains,
        const std::vector<int>& runs
    );

private:
    std::string _year;
    std::string _DataOrMC;

    // CorrectionSet do correctionlib
    std::unique_ptr<correction::CorrectionSet> cset;

    // Retorna o caminho correto do JSON de acordo com ano e tipo
    std::string getJSONPath() const;
};
