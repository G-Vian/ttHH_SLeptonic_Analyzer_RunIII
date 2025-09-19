#include "ElectronEnergyCalibrator.h"
#include <iostream>
#include <iomanip> // para std::setprecision

void ElectronEnergyCalibrator::calibrateElectrons(
    std::vector<float>& pts,
    const std::vector<float>& etas,
    const std::vector<float>& r9s,
    const std::vector<int>& gains,
    int nEle)
{
    if (!eleCorr) {
        std::cerr << "[ERROR] Correção não inicializada!" << std::endl;
        return;
    }

    std::cerr << "[DEBUG] Iniciando calibrateElectrons()" << std::endl;
    std::cerr << "  Calibrator type: " << _DataOrMC << std::endl;
    std::cerr << "  Year: " << _year << std::endl;
    std::cerr << "  nEle = " << nEle << std::endl;

    for (int i = 0; i < nEle; i++) {
        double pt = pts[i];
        double eta = etas[i];
        double r9 = r9s[i];
        int gain = gains[i];

        std::cerr << "\n[DEBUG] --- Eletrón " << i << " ---" << std::endl;
        std::cerr << "  pt   = " << std::setprecision(6) << pt << std::endl;
        std::cerr << "  eta  = " << std::setprecision(6) << eta << std::endl;
        std::cerr << "  r9   = " << std::setprecision(6) << r9 << std::endl;
        std::cerr << "  gain = " << gain << std::endl;

        try {
            double corr = 1.0;

            if (_DataOrMC == "DATA") {
                // Para DATA: string "central" + pt + r9 + abs(eta)
                std::vector<correction::Variable::Type> args;
                args.emplace_back(std::string("central"));
                args.emplace_back(pt);
                args.emplace_back(r9);
                args.emplace_back(std::abs(eta));

                std::cerr << "  [ARGS-DATA] ";
                for (auto& a : args) {
                    if (std::holds_alternative<int>(a)) std::cerr << std::get<int>(a) << " ";
                    else if (std::holds_alternative<double>(a)) std::cerr << std::get<double>(a) << " ";
                    else if (std::holds_alternative<std::string>(a)) std::cerr << "\"" << std::get<std::string>(a) << "\" ";
                }
                std::cerr << std::endl;

                corr = eleCorr->evaluate(args);

            } else if (_DataOrMC == "MC") {
                // Para MC: somente floats pt, r9, absEta
                std::vector<correction::Variable::Type> args;
                args.emplace_back(pt);
                args.emplace_back(r9);
                args.emplace_back(std::abs(eta));

                std::cerr << "  [ARGS-MC] ";
                for (auto& a : args) {
                    if (std::holds_alternative<int>(a)) std::cerr << std::get<int>(a) << " ";
                    else if (std::holds_alternative<double>(a)) std::cerr << std::get<double>(a) << " ";
                    else if (std::holds_alternative<std::string>(a)) std::cerr << "\"" << std::get<std::string>(a) << "\" ";
                }
                std::cerr << std::endl;

                corr = eleCorr->evaluate(args);
            }

            double pt_before = pts[i];
            double pt_after = pt_before * corr;

            std::cerr << "  Correction factor = " << corr << std::endl;
            std::cerr << "  Pt antes/depois   = " << pt_before << " / " << pt_after << std::endl;

            pts[i] = pt_after;

        } catch (const std::exception& e) {
            std::cerr << "[ERROR] Calibração falhou para elétron " << i
                      << ": " << e.what() << std::endl;
            std::cerr << "  --> Inputs: pt=" << pt
                      << " r9=" << r9
                      << " eta=" << eta
                      << " gain=" << gain
                      << " DataOrMC=" << _DataOrMC
                      << " Year=" << _year << std::endl;
        }
    }
}

