#include "ElectronEnergyCalibrator.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include "ttHHanalyzer_trigger.h"   // <-- aqui está objectLep

// ------------------------
// Construtor
// ------------------------
ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& DataOrMC)
    : _year(year), _DataOrMC(DataOrMC), rng(std::random_device{}())
{
    std::cerr << "[DEBUG] Inicializando ElectronEnergyCalibrator: "
              << _DataOrMC << " / " << _year << std::endl;

    std::string jsonPath = getElectronJSONPath();
    std::cerr << "[DEBUG] JSON path: " << jsonPath << std::endl;
    cset = correction::CorrectionSet::from_file(jsonPath);
}

// ------------------------
// Caminho do JSON
// ------------------------
std::string ElectronEnergyCalibrator::getElectronJSONPath() const {
    if (_year == "2022") {
        return (_DataOrMC == "MC")
            ? "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/MC/electronSS_EtDependent.json.gz"
            : "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2022EE") { // pós EE
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoE+PromptFG/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2023") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2023B") { // pós BPIX
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023B/SS/electronSS_EtDependent.json.gz";
    } else if (_year == "2024") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/SS/electronSS_EtDependent_v1.json.gz";
    } else {
        throw std::runtime_error("Year not supported for electron corrections!");
    }
}

// ------------------------
// Calibração principal
// ------------------------
void ElectronEnergyCalibrator::calibrateElectrons(
    std::vector<objectLep*>& electrons, 
    unsigned int runNumber, 
    const std::string& syst // pode ser "central", "up", "down"
) {
    try {
        std::cout << "[INFO] Iniciando calibração de elétrons. Número de elétrons: " 
                  << electrons.size() << std::endl;
        std::cout << "      _DataOrMC = " << _DataOrMC 
                  << " | Year = " << _year 
                  << " | RunNumber = " << runNumber << std::endl;

        for (size_t i = 0; i < electrons.size(); ++i) {
            auto* ele = electrons[i];

            float pt  = clamp(ele->pt(), 5.f, 1000.f);
            float eta = clamp(ele->eta(), -2.5f, 2.5f);
            float r9  = clamp(ele->r9(), 0.f, 1.5f);
            float gain = clamp(ele->gain(), 0.f, 12.f);

            std::cout << "[DEBUG] Eletrón " << i 
                      << " original: pt=" << ele->pt() 
                      << " eta=" << ele->eta() 
                      << " r9=" << ele->r9() 
                      << " gain=" << ele->gain() 
                      << ", clamped: pt=" << pt 
                      << " eta=" << eta 
                      << " r9=" << r9 
                      << " gain=" << gain << std::endl;

            float newPt = pt;

            // ===== DATA =====
            if (_DataOrMC == "DATA") {
                if (_year == "2024") {
                    // 5 args: "scale", run, eta, r9, pt, gain
                    std::cout << "[DEBUG] Args DATA 2024: "
                              << runNumber << " " << eta << " " << r9 
                              << " " << pt << " " << gain << std::endl;

                    newPt = _scaleEvaluator->evaluate(
                        "scale", runNumber, eta, r9, pt, gain
                    );

                } else {
                    // 6 args: "scale", run, eta, r9, abs(eta), pt, gain
                    std::cout << "[DEBUG] Args DATA <=2023: "
                              << runNumber << " " << eta << " " << r9 
                              << " " << std::abs(eta) << " " << pt << " " << gain << std::endl;

                    newPt = _scaleEvaluator->evaluate(
                        "scale", runNumber, eta, r9, std::abs(eta), pt, gain
                    );
                }
            }

            // ===== MC =====
            else if (_DataOrMC == "MC") {
                // Sempre 4 args: "smear", pt, r9, abs(eta)
                std::cout << "[DEBUG] Args MC: " 
                          << pt << " " << r9 << " " << std::abs(eta) << std::endl;

                newPt = _smearEvaluator->evaluate("smear", pt, r9, std::abs(eta));

                // Se quiser aplicar incertezas:
                if (syst == "up") {
                    float unc = _smearEvaluator->evaluate("esmear", pt, r9, std::abs(eta));
                    newPt += unc;
                    std::cout << "[DEBUG] Aplicando sistemático UP: +" << unc << std::endl;
                } else if (syst == "down") {
                    float unc = _smearEvaluator->evaluate("esmear", pt, r9, std::abs(eta));
                    newPt -= unc;
                    std::cout << "[DEBUG] Aplicando sistemático DOWN: -" << unc << std::endl;
                }
            }

            else {
                std::cerr << "[ERROR] Tipo de calibração desconhecido: " << _DataOrMC << std::endl;
            }

            std::cout << "[DEBUG] Electron[" << i << "] Pt antes/depois: " 
                      << ele->pt() << " / " << newPt << std::endl;

            // Atualiza o valor calibrado no objeto
            ele->setCalibratedPt(newPt);
        }
    } 
    catch (const std::exception& e) {
        std::cerr << "[ERROR] Exceção inesperada na calibração: " 
                  << e.what() << std::endl;
    }
}


// ------------------------
// Limites de segurança
// ------------------------
float ElectronEnergyCalibrator::getMin(const std::string& varName) const {
    if (varName == "pt") return 5.0f;
    if (varName == "ScEta") return -2.5f;
    if (varName == "r9") return 0.0f;
    if (varName == "seedGain") return 0;
    return 0.0f;
}

float ElectronEnergyCalibrator::getMax(const std::string& varName) const {
    if (varName == "pt") return 1000.0f;
    if (varName == "ScEta") return 2.5f;
    if (varName == "r9") return 1.5f;
    if (varName == "seedGain") return 12;
    return 1.0f;
}
