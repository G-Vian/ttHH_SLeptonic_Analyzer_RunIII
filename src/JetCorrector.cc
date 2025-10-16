#include "src/JetCorrector.h"
#include <iostream>
#include <stdexcept>
#include <string>

JetCorrector::JetCorrector(const std::string& year, const std::string& dataOrMC) {
    std::string base_path = "/afs/cern.ch/user/g/gvian/JME_corrections/";
    std::string json_path;

    if (year == "2022") {
        json_path = base_path + "Run3-22CDSep23-Summer22-NanoAODv12/jet_jerc.json.gz";
    } else if (year == "2022EE") {
        json_path = base_path + "Run3-22EFGSep23-Summer22EE-NanoAODv12/jet_jerc.json.gz";
    } else if (year == "2023") {
        json_path = base_path + "Run3-23CSep23-Summer23-NanoAODv12/jet_jerc.json.gz";
    } else if (year == "2023B") {
        json_path = base_path + "Run3-23DSep23-Summer23BPix-NanoAODv12/jet_jerc.json.gz";
    } else if (year == "2024") {
        json_path = base_path + "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.json.gz";
    } else {
        throw std::runtime_error("[JetCorrector] Ano não suportado: " + year);
    }

    std::cout << "[JetCorrector] Carregando correções de: " << json_path << std::endl;
    cset_ = correction::CorrectionSet::from_file(json_path);

    // Define a pilha de correções JEC com base em Data/MC
    isMC_ = (dataOrMC == "MC");
    if (isMC_) {
        // Os nomes exatos podem variar um pouco entre os arquivos, ajuste se necessário
        if (year == "2024") {
             jec_stack_names_ = {"Summer24_V1_MC_L1FastJet_AK4PFPuppi", "Summer24_V1_MC_L2Relative_AK4PFPuppi", "Summer24_V1_MC_L3Absolute_AK4PFPuppi"};
        } else {
             jec_stack_names_ = {"Summer22_22Sep2023_V3_MC_L1FastJet_AK4PFPuppi", "Summer22_22Sep2023_V3_MC_L2Relative_AK4PFPuppi", "Summer22_22Sep2023_V3_MC_L3Absolute_AK4PFPuppi"};
        }
    } else { // DATA
         if (year == "2024") {
             jec_stack_names_ = {"Summer24_V1_DATA_L1FastJet_AK4PFPuppi", "Summer24_V1_DATA_L2Relative_AK4PFPuppi", "Summer24_V1_DATA_L3Absolute_AK4PFPuppi", "Summer24_V1_DATA_L2L3Residual_AK4PFPuppi"};
         } else if (year == "2022") { // RunCD
            jec_stack_names_ = {"Summer22_22Sep2023_RunCD_V3_DATA_L1FastJet_AK4PFPuppi", "Summer22_22Sep2023_RunCD_V3_DATA_L2Relative_AK4PFPuppi", "Summer22_22Sep2023_RunCD_V3_DATA_L3Absolute_AK4PFPuppi", "Summer22_22Sep2023_RunCD_V3_DATA_L2L3Residual_AK4PFPuppi"};
         } // Adicione else if para 2022EE, 2023, etc., se os nomes forem diferentes.
    }
    
    // Carrega os nomes das fontes de incerteza (exemplos, pode haver mais)
    // No seu código principal, você pode usar getUncertaintySources() para pegar a lista completa.
    if (isMC_) {
        // Exemplo de como pegar algumas fontes de incerteza
        // Os nomes podem mudar entre os anos, verifique seu JSON.
         uncertainty_sources_.push_back("FlavorQCD");
         uncertainty_sources_.push_back("RelativeBal");
         uncertainty_sources_.push_back("HF");
         // ... adicione outras fontes conforme necessário
    }
}

TLorentzVector JetCorrector::getCorrectedP4(const pat::Jet& rawJet, double rho, const std::string& systematic) const {
    double pt_raw = rawJet.pt;
    double eta = rawJet.eta;
    double area = rawJet.area;

    // --- Passo 1: Aplicar a correção JEC nominal ---
    double pt_nominal = pt_raw;
    for (const auto& level_name : jec_stack_names_) {
        auto evaluator = cset_->at(level_name);
        double factor = (level_name.find("L1FastJet") != std::string::npos) ?
                        evaluator->evaluate({area, eta, pt_nominal, rho}) :
                        evaluator->evaluate({eta, pt_nominal});
        pt_nominal *= factor;
    }

    double pt_final = pt_nominal;

    // --- Passo 2: Aplicar a variação sistemática, se solicitada ---
    if (systematic != "nominal") {
        size_t pos_up = systematic.find("_up");
        size_t pos_down = systematic.find("_down");

        if (pos_up != std::string::npos || pos_down != std::string::npos) {
            size_t pos_end = (pos_up != std::string::npos) ? pos_up : pos_down;
            std::string source_name = systematic.substr(0, pos_end);
            
            // O nome no JSON pode ser mais longo, ex: "Summer22..._Regrouped_FlavorQCD_..."
            // Esta é uma busca simplificada, pode precisar de ajuste.
            std::string unc_evaluator_name;
            for(auto const& corr : cset_->corrections()){
                if(corr.name().find(source_name) != std::string::npos && corr.name().find("_AK4PFPuppi") != std::string::npos){
                    unc_evaluator_name = corr.name();
                    break;
                }
            }

            if (!unc_evaluator_name.empty()) {
                double unc_value = cset_->at(unc_evaluator_name)->evaluate({eta, pt_nominal});
                double shift = (pos_up != std::string::npos) ? (1.0 + unc_value) : (1.0 - unc_value);
                pt_final = pt_nominal * shift;
            } else {
                 std::cerr << "[JetCorrector] AVISO: Fonte de incerteza '" << source_name << "' não encontrada." << std::endl;
            }
        }
    }
    
    TLorentzVector corrected_p4;
    corrected_p4.SetPtEtaPhiM(pt_final, rawJet.eta, rawJet.phi, rawJet.mass * (pt_final / pt_raw) );
    return corrected_p4;
}

const std::vector<std::string>& JetCorrector::getUncertaintySources() const {
    return uncertainty_sources_;
}
