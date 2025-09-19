#include "ElectronEnergyCalibrator.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>

// ------------------------
// Construtor
// ------------------------
ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year,
                                                   const std::string& DataOrMC)
    : _year(year), _DataOrMC(DataOrMC), rng(std::random_device{}()) {
    std::string jsonPath = getJSONPath();
    std::cerr << "[DEBUG] Inicializando Calibrator: Year=" << _year
              << ", Type=" << _DataOrMC << std::endl;
    std::cerr << "[DEBUG] Carregando JSON: " << jsonPath << std::endl;

    cset = correction::CorrectionSet::from_file(jsonPath);
    if (!cset) {
        throw std::runtime_error("Falha ao carregar CorrectionSet de " + jsonPath);
    }
}

// ------------------------
// Caminho do JSON
// ------------------------
std::string ElectronEnergyCalibrator::getJSONPath() const {
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
        throw std::runtime_error("Ano não suportado: " + _year);
    }
}

// ------------------------
// Calibração principal
// ------------------------
void ElectronEnergyCalibrator::calibrateElectrons(
        std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains,
        int runNumber
) {
    if (pts.empty()) {
        std::cerr << "[DEBUG] Nenhum elétron no evento" << std::endl;
        return;
    }

    try {
        for (size_t i = 0; i < pts.size(); ++i) {
            float pt = pts[i];
            float eta = etas[i];
            float r9 = r9s[i];
            int gain = gains[i];
            float absEta = std::abs(eta);

            std::cerr << "[DEBUG] Eletrón " << i
                      << " pt/eta/r9/gain: " << pt << "/" << eta
                      << "/" << r9 << "/" << gain << std::endl;

            double corr = 1.0;

            if (_DataOrMC == "DATA") {
                // ========================
                // DATA: scale corrections
                // ========================
                std::string syst = "central";
                std::vector<std::variant<int,double,std::string>> args;
                args.emplace_back(syst);
                args.emplace_back(runNumber);
                args.emplace_back(eta);
                args.emplace_back(r9);
                args.emplace_back(absEta);
                args.emplace_back(pt);
                args.emplace_back(gain);

                std::cerr << "  [ARGS DATA] syst=" << syst
                          << " run=" << runNumber
                          << " eta=" << eta
                          << " r9=" << r9
                          << " absEta=" << absEta
                          << " pt=" << pt
                          << " gain=" << gain << std::endl;

                auto scale_corr = cset->compound().at(
                    (_year=="2022") ? "EGMScale_Compound_Ele_2022preEE" :
                    (_year=="2022EE") ? "EGMScale_Compound_Ele_2022postEE" :
                    (_year=="2023") ? "EGMScale_Compound_Ele_2023preBPIX" :
                    (_year=="2023B") ? "EGMScale_Compound_Ele_2023postBPIX" :
                    "EGMScale_Compound_Ele_2024"
                );

                corr = scale_corr->evaluate(args);

            } else {
                // ========================
                // MC: Smearing
                // ========================
                std::string syst = "central";
                std::vector<std::variant<int,double,std::string>> args;
                args.emplace_back(pt);
                args.emplace_back(r9);
                args.emplace_back(absEta);
                args.emplace_back(syst);

                std::cerr << "  [ARGS MC] pt=" << pt
                          << " r9=" << r9
                          << " absEta=" << absEta
                          << " syst=" << syst << std::endl;

                auto smear_corr = cset->at(
                    (_year=="2022") ? "EGMSmearAndSyst_ElePTsplit_2022preEE" :
                    (_year=="2022EE") ? "EGMSmearAndSyst_ElePTsplit_2022postEE" :
                    (_year=="2023") ? "EGMSmearAndSyst_ElePTsplit_2023preBPIX" :
                    (_year=="2023B") ? "EGMSmearAndSyst_ElePTsplit_2023postBPIX" :
                    "EGMSmearAndSyst_ElePTsplit_2024"
                );

                std::normal_distribution<float> gauss(0.0, 1.0);
                double smear = smear_corr->evaluate(args);
                corr = 1.0 + smear * gauss(rng);
            }

            std::cerr << "  -> Correção aplicada: " << corr << std::endl;
            float oldPt = pts[i];
            pts[i] *= static_cast<float>(corr);
            std::cerr << "  -> Pt antes/depois: " << oldPt << " / " << pts[i] << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Calibração falhou: " << e.what() << std::endl;
    }
}

// ------------------------
// Limites de segurança
// ------------------------
float ElectronEnergyCalibrator::getMin(const std::string& var) const {
    if (var=="pt") return 5.0f;
    if (var=="ScEta") return -2.5f;
    if (var=="r9") return 0.0f;
    if (var=="seedGain") return 0;
    return 0.0f;
}

float ElectronEnergyCalibrator::getMax(const std::string& var) const {
    if (var=="pt") return 1000.0f;
    if (var=="ScEta") return 2.5f;
    if (var=="r9") return 1.5f;
    if (var=="seedGain") return 12;
    return 1.0f;
}
