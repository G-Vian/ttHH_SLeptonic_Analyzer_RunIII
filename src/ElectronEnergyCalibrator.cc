#include "ElectronEnergyCalibrator.h"
#include <algorithm>
#include <cmath>
#include <iostream>

// ------------------------
// Construtor
// ------------------------
ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& DataOrMC, const std::string& year)
    : _year(year), _DataOrMC(DataOrMC), rng(std::random_device{}())
{
    try {
        std::cerr << "[DEBUG] Inicializando ElectronEnergyCalibrator para "
                  << _DataOrMC << " / " << _year << std::endl;

        // ✅ usar diretamente o retorno de from_file
        cset = correction::CorrectionSet::from_file(getElectronJSONPath());

        std::cerr << "[DEBUG] CorrectionSet carregado com sucesso." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Falha ao carregar CorrectionSet: " << e.what() << std::endl;
    }
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
        std::vector<float>& pts,
        const std::vector<float>& etas,
        const std::vector<float>& r9s,
        const std::vector<int>& gains,
        int runNumber
) {
    if (pts.empty()) return;

    try {
        for (size_t i = 0; i < pts.size(); ++i) {
            float pt   = pts[i];
            float eta  = etas[i];
            float r9   = r9s[i];
            int gain   = gains[i];
            float absEta = std::abs(eta);

            std::cerr << "[DEBUG] Eletrón " << i
                      << " pt/eta/r9/gain: "
                      << pt << " / " << eta << " / "
                      << r9 << " / " << gain << std::endl;

            std::vector<std::variant<int,double,std::string>> args;
            double corr = 1.0;

            if (_DataOrMC == "DATA") {
                // ========================
                // DATA
                // ========================
                std::string syst = "central";

                args.emplace_back(syst);
                args.emplace_back(runNumber);
                args.emplace_back(eta);
                args.emplace_back(r9);
                args.emplace_back(absEta);
                args.emplace_back(pt);
                args.emplace_back(gain);

                std::cerr << "[DEBUG][DATA] Args:";
                for (auto& a : args) {
                    if (std::holds_alternative<int>(a))
                        std::cerr << " int=" << std::get<int>(a);
                    else if (std::holds_alternative<double>(a))
                        std::cerr << " double=" << std::get<double>(a);
                    else
                        std::cerr << " string=" << std::get<std::string>(a);
                }
                std::cerr << std::endl;

                auto scale_corr = cset->compound().at(
                    (_year=="2022")  ? "EGMScale_Compound_Ele_2022preEE" :
                    (_year=="2022EE")? "EGMScale_Compound_Ele_2022postEE" :
                    (_year=="2023")  ? "EGMScale_Compound_Ele_2023preBPIX" :
                    (_year=="2023B") ? "EGMScale_Compound_Ele_2023postBPIX" :
                                       "EGMScale_Compound_Ele_2024"
                );

                corr = scale_corr->evaluate(args);
                std::cerr << "[DEBUG][DATA] scale factor = " << corr << std::endl;

            } else {
                // ========================
                // MC
                // ========================
                args.emplace_back(pt);
                args.emplace_back(r9);
                args.emplace_back(absEta);

                std::cerr << "[DEBUG][MC] Args:";
                for (auto& a : args) {
                    if (std::holds_alternative<int>(a))
                        std::cerr << " int=" << std::get<int>(a);
                    else if (std::holds_alternative<double>(a))
                        std::cerr << " double=" << std::get<double>(a);
                    else
                        std::cerr << " string=" << std::get<std::string>(a);
                }
                std::cerr << std::endl;

                auto smear_corr = cset->at(
                    (_year=="2022")  ? "EGMSmearAndSyst_ElePTsplit_2022preEE" :
                    (_year=="2022EE")? "EGMSmearAndSyst_ElePTsplit_2022postEE" :
                    (_year=="2023")  ? "EGMSmearAndSyst_ElePTsplit_2023preBPIX" :
                    (_year=="2023B") ? "EGMSmearAndSyst_ElePTsplit_2023postBPIX" :
                                       "EGMSmearAndSyst_ElePTsplit_2024"
                );

                double smear = smear_corr->evaluate(args);
                std::normal_distribution<float> gauss(0.0, 1.0);
                corr = 1.0 + smear * gauss(rng);

                std::cerr << "[DEBUG][MC] smear = " << smear
                          << " → applied corr = " << corr << std::endl;
            }

            float pt_before = pts[i];
            pts[i] *= static_cast<float>(corr);
            std::cerr << "[DEBUG] Electron[" << i << "] Pt antes/depois: "
                      << pt_before << " / " << pts[i] << std::endl;
        }
    } catch (const std::out_of_range& e) {
        std::cerr << "[ERROR] Correção falhou: " << e.what() << std::endl;
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
