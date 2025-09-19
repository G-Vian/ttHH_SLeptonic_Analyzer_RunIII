#pragma once

#include <vector>
#include <string>
#include <memory>
#include <correction.h>

// Classe para calibrar energia de elétrons usando correctionlib
class ElectronEnergyCalibrator {
public:
    // Construtor recebe ano e tipo de dado ("DATA" ou "MC")
    ElectronEnergyCalibrator(const std::string& year, const std::string& DataOrMC);

    // Função principal para calibrar elétrons
    void calibrateElectrons(std::vector<float>& pts,
                            const std::vector<float>& etas,
                            const std::vector<float>& r9s,
                            const std::vector<int>& gains,
                            int nEle);

private:
    std::string _year;
    std::string _DataOrMC;

    // CorrectionSet carregado a partir do JSON
    // 🔑 ESTE É O ÚNICO MEMBRO NECESSÁRIO
    std::unique_ptr<correction::CorrectionSet> cset;

    // 🔴 NÃO declare eleCorr aqui! Ele é buscado dentro do .cc com:
    // auto eleCorr = cset->at("ElectronEnergyCorrection");

    // Retorna o caminho do JSON correto para ano e tipo (DATA/MC)
    std::string getJSONPath() const;
};
