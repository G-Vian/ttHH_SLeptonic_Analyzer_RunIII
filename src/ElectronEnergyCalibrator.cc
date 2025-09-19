#include "ElectronEnergyCalibrator.h"
#include <stdexcept>
#include <iostream>

// Construtor
ElectronEnergyCalibrator::ElectronEnergyCalibrator(const std::string& year, const std::string& DataOrMC)
    : _year(year), _DataOrMC(DataOrMC)
{
    std::string jsonPath = getJSONPath();
    try {
        cset = std::make_unique<correction::CorrectionSet>(correction::CorrectionSet::from_file(jsonPath));
        std::cout << "[ElectronEnergyCalibrator] Loaded JSON: " << jsonPath << std::endl;
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("[ElectronEnergyCalibrator] Failed to load JSON: ") + e.what());
    }
}

// Função para retornar o caminho correto do JSON
std::string ElectronEnergyCalibrator::getJSONPath() const {
    if (_year == "2022") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoBCD/SS/electronSS_EtDependent.json.gz";
    }
    if (_year == "2022EE") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoE+PromptFG/SS/electronSS_EtDependent.json.gz";
    }
    if (_year == "2023") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23C/SS/electronSS_EtDependent.json.gz";
    }
    if (_year == "2023B") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/SS/electronSS_EtDependent.json.gz";
    }
    if (_year == "2024") {
        return "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/SS/electronSS_EtDependent_v1.json.gz";
    }
    throw std::runtime_error("[ElectronEnergyCalibrator] Ano não reconhecido: " + _year);
}

// Função principal de calibração
void ElectronEnergyCalibrator::calibrateElectrons(
    std::vector<float>& pts,
    const std::vector<float>& etas,
    const std::vector<float>& r9s,
    const std::vector<int>& gains,
    const std::vector<int>& runs
) {
    if (!cset) throw std::runtime_error("[ElectronEnergyCalibrator] CorrectionSet not initialized!");

    for (size_t i = 0; i < pts.size(); ++i) {
        float pt_before = pts[i];
        double Et = pt_before; // usamos pT como proxy de energia

        const auto& eleCorr = cset->at("EGMScale_Compound_Ele_" + _year); // ajustável por JSON
        std::string key;
        if (_DataOrMC == "DATA") {
            key = "scale";  // apenas escala para DATA
        } else {
            key = "escale"; // smearing / scale MC
        }

        // Monta vetor de inputs na ordem correta exigida pelo JSON
        std::vector<correction::Variable::Type> inputs;
        inputs.push_back(key);                   // "scale" ou "escale"
        inputs.push_back(runs[i]);               // run
        inputs.push_back(etas[i]);               // eta ou supercluster eta
        inputs.push_back(r9s[i]);                // r9
        inputs.push_back(Et);                    // Et
        inputs.push_back(static_cast<double>(gains[i])); // seed crystal gain

        try {
            double corr = eleCorr.evaluate(inputs);
            pts[i] *= corr;
            std::cout << "[DEBUG] Electron[" << i << "] Pt antes/depois: "
                      << pt_before << " / " << pts[i] << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] Calibração falhou para o elétron " << i
                      << ": " << e.what() << std::endl;
        }
    }
}
