#include "ElectronEnergyCalibrator.h"
#include <iostream>
#include <cmath>
#include <random>

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& DataOrMC, const std::string& year)
    : _year(year), _DataOrMC(DataOrMC), rng(std::random_device{}()) 
{
    try {
        std::cerr << "[DEBUG] Inicializando ElectronEnergyCalibrator para " 
                  << _DataOrMC << " / " << _year << std::endl;

        cset = std::make_unique<correction::CorrectionSet>(
            correction::CorrectionSet::from_file("electron_corrections.json"));

        std::cerr << "[DEBUG] CorrectionSet carregado com sucesso." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Falha ao carregar CorrectionSet: " << e.what() << std::endl;
    }
}

void ElectronEnergyCalibrator::calibrateElectrons(
    std::vector<float>& pts,
    const std::vector<float>& etas,
    const std::vector<float>& r9s,
    const std::vector<int>& gains,
    int runNumber)
{
    if (!cset) {
        std::cerr << "[ERROR] CorrectionSet não carregado!" << std::endl;
        return;
    }

    correction::Correction::Ref corrObject = cset->at("electron_calibration");

    for (size_t i = 0; i < pts.size(); ++i) {
        float pt = pts[i];
        float eta = etas[i];
        float r9 = r9s[i];
        int gain = gains[i];

        std::cerr << "[DEBUG] Eletrón " << i 
                  << " pt/eta/r9/gain: " << pt << "/" << eta << "/" << r9 << "/" << gain 
                  << std::endl;

        float corr = 1.0;

        try {
            if (_DataOrMC == "DATA") {
                std::vector<correction::Variable::Type> args = {
                    correction::Variable::Type(runNumber),
                    correction::Variable::Type(eta),
                    correction::Variable::Type(r9),
                    correction::Variable::Type(gain),
                    correction::Variable::Type("central")
                };

                std::cerr << "[DEBUG][DATA] Chamando evaluate com args: "
                          << "run=" << runNumber << " eta=" << eta 
                          << " r9=" << r9 << " gain=" << gain 
                          << " syst=central" << std::endl;

                corr = corrObject->evaluate(args);

            } else if (_DataOrMC == "MC") {
                std::vector<correction::Variable::Type> args = {
                    correction::Variable::Type(pt),
                    correction::Variable::Type(r9),
                    correction::Variable::Type(fabs(eta))
                };

                std::cerr << "[DEBUG][MC] Chamando evaluate com args: "
                          << "pt=" << pt << " r9=" << r9 
                          << " absEta=" << fabs(eta) << std::endl;

                corr = corrObject->evaluate(args);

                // Aplica smearing gaussiano no MC
                std::normal_distribution<float> gauss(0.0, 1.0);
                float smear = 0.01; // desvio relativo fixo (ajustar conforme necessidade)
                corr = 1.0 + smear * gauss(rng);

                std::cerr << "[DEBUG][MC] Correção com smearing aplicado = " << corr << std::endl;
            }

            std::cerr << "[DEBUG] Correção final aplicada = " << corr << std::endl;
            pts[i] *= corr;

        } catch (const std::exception& e) {
            std::cerr << "[ERROR] Calibração falhou: " << e.what() << std::endl;
            continue;
        }

        std::cerr << "[DEBUG] Electron[" << i << "] Pt antes/depois: "
                  << pt << " / " << pts[i] << std::endl;
    }

    std::cerr << "[DEBUG] Calibrator type: " << _DataOrMC 
              << ", Year: " << _year << std::endl;
}

// Funções auxiliares de range
float ElectronEnergyCalibrator::getMin(const std::string& var) const {
    if (var == "eta") return -2.5;
    if (var == "pt") return 5.0;
    if (var == "r9") return 0.0;
    return 0.0;
}

float ElectronEnergyCalibrator::getMax(const std::string& var) const {
    if (var == "eta") return 2.5;
    if (var == "pt") return 500.0;
    if (var == "r9") return 1.0;
    return 0.0;
}
