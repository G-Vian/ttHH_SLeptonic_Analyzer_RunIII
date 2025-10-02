#include "ttHHanalyzer_trigger.h" // Inclua seu header principal

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>

// --- DEFINIÇÃO DAS VARIÁVEIS GLOBAIS DE SF ---
json muonTrigSFJson;
json muonHighPtTrigSFJson;
json muonLowPtIDSFJson;
json muonMediumPtIDSFJson;
json muonHighPtIDSFJson;
json muonIsoSFJson;

/**
 * @brief Função auxiliar para carregar, tratar e fazer o parse de um arquivo JSON.
 */
json loadSFJson(const TString& filePath) {
    if (filePath.IsNull()) {
        std::cerr << "[loadSFJson] ERRO: O caminho do arquivo (filePath) está vazio!" << std::endl;
        return json();
    }
    std::ifstream input(filePath.Data());
    if (!input.is_open()) {
        std::cerr << "[loadSFJson] ERRO: Falha ao abrir o arquivo: " << filePath.Data() << std::endl;
        return json();
    }
    std::string json_str((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    input.close();
    size_t pos = 0;
    while ((pos = json_str.find("Infinity", pos)) != std::string::npos) {
        json_str.replace(pos, 8, "1e10");
    }
    try {
        return json::parse(json_str);
    } catch (const std::exception& e) {
        std::cerr << "[loadSFJson] ERRO: Falha no parse do JSON do arquivo " << filePath.Data() << ": " << e.what() << std::endl;
        return json();
    }
}

// =================================================================================
//          FUNÇÕES PARA MUON TRIGGER SCALE FACTORS
// =================================================================================

void ttHHanalyzer::initMuonHLTriggerSF() {
    TString sfMediumPtFilePath;
    TString sfHighPtFilePath;
    TString localDir = "/afs/cern.ch/user/g/gvian/muon_SF/muonefficiencies/Run3/";

    if (_year == "2022") {
        sfMediumPtFilePath = localDir + "2022/2022_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2022_eta_pt_schemaV2.json";
        sfHighPtFilePath = localDir + "2022/2022_HighPt/ScaleFactors_Muon_highPt_RECO_2022_schemaV2.json";
    } else if (_year == "2022EE") {
        sfMediumPtFilePath = localDir + "2022_EE/2022_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2022_EE_eta_pt_schemaV2.json";
        sfHighPtFilePath = localDir + "2022_EE/2022_HighPt/ScaleFactors_Muon_highPt_HLT_2022_EE_schemaV2.json";
    } else if (_year == "2023") {
        sfMediumPtFilePath = localDir + "2023/2023_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2023_eta_pt_schemaV2.json";
        sfHighPtFilePath = localDir + "2023/2023_HighPt/ScaleFactors_Muon_highPt_HLT_2023_schemaV2.json";
    } else if (_year == "2023B") {
        sfMediumPtFilePath = localDir + "2023_BPix/2023_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2023_BPix_eta_pt_schemaV2.json";
        sfHighPtFilePath = localDir + "2023_BPix/2023_HighPt/ScaleFactors_Muon_highPt_HLT_2023_BPix_schemaV2.json";
    } else if (_year == "2024") {
        sfMediumPtFilePath = localDir + "2023_BPix/2023_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2023_BPix_eta_pt_schemaV2.json";
        sfHighPtFilePath = localDir + "2024/2024_HighPt/ScaleFactors_Muon_highPt_HLT_2024_schemaV2.json";
    } else {
        std::cerr << "[initMuonHLTriggerSF] ERRO: Ano não suportado para SF de múons: " << _year << std::endl;
        return;
    }
    
    muonTrigSFJson = loadSFJson(sfMediumPtFilePath);
    muonHighPtTrigSFJson = loadSFJson(sfHighPtFilePath);

    // Lógica dos histogramas, se aplicável
    if (h_sf_muon_vs_pt) { delete h_sf_muon_vs_pt; h_sf_muon_vs_pt = nullptr; }
    // ... (resto da sua limpeza de histogramas) ...

    h_sf_muon_vs_pt = new TH1F("h_sf_muon_vs_pt", "Muon SF vs pT;Muon pT [GeV];SF", 100, 0, 700);
    // ... (resto da sua criação de histogramas) ...

    for (auto& h : {h_sf_muon_vs_pt, h_sf_muon_vs_eta /*, etc */}) {
        if (h) h->SetDirectory(0);
    }
}

float ttHHanalyzer::getMuonTrigSF(float eta, float pt) {
    const json* sfJson = nullptr;
    std::string correctionName;
    float eta_for_lookup;

    if (pt > 200.0) {
        sfJson = &muonHighPtTrigSFJson;
        correctionName = "NUM_HLT_DEN_TrkHighPtTightRelIsoProbes";
        eta_for_lookup = fabs(eta);
    } else {
        sfJson = &muonTrigSFJson;
        correctionName = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight";
        eta_for_lookup = eta;
    }

    if (!sfJson || sfJson->empty()) return 1.0;

    const json* correction = nullptr;
    for (const auto& corr : (*sfJson)["corrections"]) {
        if (corr.contains("name") && corr["name"] == correctionName) {
            correction = &corr;
            break;
        }
    }
    if (!correction) return 1.0;

    const auto& data_eta = (*correction)["data"];
    if (!data_eta.contains("edges") || !data_eta.contains("content")) return 1.0;

    const std::vector<float> eta_edges = data_eta["edges"].get<std::vector<float>>();
    int eta_bin = -1;
    for (size_t i = 0; i < eta_edges.size() - 1; ++i) {
        if (eta_for_lookup >= eta_edges[i] && eta_for_lookup < eta_edges[i+1]) {
            eta_bin = i;
            break;
        }
    }
    if (eta_bin == -1) return 1.0;

    const auto& data_pt = data_eta["content"][eta_bin];
    if (!data_pt.contains("edges") || !data_pt.contains("content")) return 1.0;

    const std::vector<float> pt_edges = data_pt["edges"].get<std::vector<float>>();
    int pt_bin = -1;
    for (size_t i = 0; i < pt_edges.size() - 1; ++i) {
        if (pt >= pt_edges[i] && pt < pt_edges[i+1]) {
            pt_bin = i;
            break;
        }
    }
    if (pt_bin == -1) return 1.0;

    const auto& categories = data_pt["content"][pt_bin]["content"];
    for (const auto& entry : categories) {
        if (entry.contains("key") && entry["key"] == "nominal" && entry.contains("value")) {
            return entry["value"].get<float>();
        }
    }
    
    return 1.0;
}

// =================================================================================
//          FUNÇÕES PARA MUON ID SCALE FACTORS
// =================================================================================

void ttHHanalyzer::initMuonIDSF() {
    TString sfLowPtFilePath, sfMediumPtFilePath, sfHighPtFilePath;
    TString localDir = "/afs/cern.ch/user/g/gvian/muon_SF/muonefficiencies/Run3/";

    if (_year == "2022") {
        sfLowPtFilePath    = localDir + "2022/2022_Jpsi/ScaleFactors_Muon_Jpsi_ID_2022_schemaV2.json";
        sfMediumPtFilePath = localDir + "2022/2022_Z/ScaleFactors_Muon_Z_ID_ISO_2022_schemaV2.json";
        sfHighPtFilePath   = localDir + "2022/2022_HighPt/ScaleFactors_Muon_highPt_IDISO_2022_schemaV2.json";
    } else if (_year == "2022EE") {
        sfLowPtFilePath    = localDir + "2022_EE/2022_Jpsi/ScaleFactors_Muon_Jpsi_ID_2022_EE_schemaV2.json";
        sfMediumPtFilePath = localDir + "2022_EE/2022_Z/ScaleFactors_Muon_Z_ID_ISO_2022_EE_schemaV2.json";
        sfHighPtFilePath   = localDir + "2022_EE/2022_HighPt/ScaleFactors_Muon_highPt_IDISO_2022_EE_schemaV2.json";
    } else if (_year == "2023") {
        sfLowPtFilePath    = localDir + "2023/2023_JPsi/ScaleFactors_Muon_Jpsi_ID_2023_schemaV2.json";
        sfMediumPtFilePath = localDir + "2023/2023_Z/ScaleFactors_Muon_Z_ID_ISO_2023_schemaV2.json";
        sfHighPtFilePath   = localDir + "2023/2023_HighPt/ScaleFactors_Muon_highPt_IDISO_2023_schemaV2.json";
    } else if (_year == "2023B") {
        sfLowPtFilePath    = localDir + "2023_BPix/2023_JPsi/ScaleFactors_Muon_Jpsi_ID_2023_schemaV2.json";
        sfMediumPtFilePath = localDir + "2023_BPix/2023_Z/ScaleFactors_Muon_Z_ID_ISO_2023_BPix_schemaV2.json";
        sfHighPtFilePath   = localDir + "2023_BPix/2023_HighPt/ScaleFactors_Muon_highPt_IDISO_2023_BPix_schemaV2.json";
    } else if (_year == "2024") {
        sfLowPtFilePath    = localDir + "2024/2024_JPsi/ScaleFactors_Muon_Jpsi_ID_2024_schemaV2.json";
        sfMediumPtFilePath = localDir + "2024/2024_Z/ScaleFactors_Muon_ID_ISO_2024_schemaV2.json";
        sfHighPtFilePath   = localDir + "2024/2024_HighPt/ScaleFactors_Muon_highPt_IDISO_2024_schemaV2.json";
    }
    
    muonLowPtIDSFJson = loadSFJson(sfLowPtFilePath);
    muonMediumPtIDSFJson = loadSFJson(sfMediumPtFilePath);
    muonHighPtIDSFJson = loadSFJson(sfHighPtFilePath);
}

float ttHHanalyzer::getMuonIDSF(float eta, float pt) {
    const json* sfJson = nullptr;
    std::string correctionName;
    float eta_for_lookup;

    if (pt < 30.0) {
        sfJson = &muonLowPtIDSFJson;
        correctionName = "NUM_TightID_DEN_TrackerMuons";
        eta_for_lookup = fabs(eta);
    } else if (pt < 200.0) {
        sfJson = &muonMediumPtIDSFJson;
        correctionName = "NUM_TightID_DEN_TrackerMuons";
        if (_year == "2023" || _year == "2023B" || _year == "2024") {
            eta_for_lookup = eta;
        } else {
            eta_for_lookup = fabs(eta);
        }
    } else {
        sfJson = &muonHighPtIDSFJson;
        correctionName = "NUM_TightID_DEN_GlobalMuonProbes";
        eta_for_lookup = fabs(eta);
    }
    
    if (!sfJson || sfJson->empty()) return 1.0;

    const json* correction = nullptr;
    for (const auto& corr : (*sfJson)["corrections"]) {
        if (corr.contains("name") && corr["name"] == correctionName) {
            correction = &corr;
            break;
        }
    }
    if (!correction) return 1.0;

    const auto& data_eta = (*correction)["data"];
    if (!data_eta.contains("edges") || !data_eta.contains("content")) return 1.0;
    
    const std::vector<float> eta_edges = data_eta["edges"].get<std::vector<float>>();
    int eta_bin = -1;
    for (size_t i = 0; i < eta_edges.size() - 1; ++i) {
        if (eta_for_lookup >= eta_edges[i] && eta_for_lookup < eta_edges[i+1]) {
            eta_bin = i;
            break;
        }
    }
    if (eta_bin == -1) return 1.0;
    
    const auto& data_pt = data_eta["content"][eta_bin];
    if (!data_pt.contains("edges") || !data_pt.contains("content")) return 1.0;
    
    const std::vector<float> pt_edges = data_pt["edges"].get<std::vector<float>>();
    int pt_bin = -1;
    for (size_t i = 0; i < pt_edges.size() - 1; ++i) {
        if (pt >= pt_edges[i] && pt < pt_edges[i+1]) {
            pt_bin = i;
            break;
        }
    }
    if (pt_bin == -1) return 1.0;
    
    const auto& categories = data_pt["content"][pt_bin]["content"];
    for (const auto& entry : categories) {
        if (entry.contains("key") && entry["key"] == "nominal" && entry.contains("value")) {
            return entry["value"].get<float>();
        }
    }    
    return 1.0;
}

// =================================================================================
//          FUNÇÕES PARA MUON ISO SCALE FACTORS
// =================================================================================

void ttHHanalyzer::initMuonIsoSF() {
    TString sfFilePath;
    TString baseDir = "/afs/cern.ch/user/g/gvian/muon_SF/muonefficiencies/Run3/";

    if (_year == "2022") {
        sfFilePath = baseDir + "2022/2022_Z/ScaleFactors_Muon_Z_ID_ISO_2022_schemaV2.json";
    } else if (_year == "2022EE") {
        sfFilePath = baseDir + "2022_EE/2022_Z/ScaleFactors_Muon_Z_ID_ISO_2022_EE_schemaV2.json";
    } else if (_year == "2023") {
        sfFilePath = baseDir + "2023/2023_Z/ScaleFactors_Muon_Z_ID_ISO_2023_schemaV2.json";
    } else if (_year == "2023B") {
        sfFilePath = baseDir + "2023_BPix/2023_Z/ScaleFactors_Muon_Z_ID_ISO_2023_BPix_schemaV2.json";
    } else if (_year == "2024") {
        sfFilePath = baseDir + "2024/2024_Z/ScaleFactors_Muon_ID_ISO_2024_schemaV2.json";
    } else {
        return;
    }
    
    muonIsoSFJson = loadSFJson(sfFilePath);
}

float ttHHanalyzer::getMuonIsoSF(float eta, float pt) {
    const json* sfJson = &muonIsoSFJson;
    std::string correctionName = "NUM_TightPFIso_DEN_TightID";
    float eta_for_lookup = fabs(eta);

    if (!sfJson || sfJson->empty()) return 1.0;

    const json* correction = nullptr;
    for (const auto& corr : (*sfJson)["corrections"]) {
        if (corr.contains("name") && corr["name"] == correctionName) {
            correction = &corr;
            break;
        }
    }
    if (!correction) return 1.0;

    const auto& data_eta = (*correction)["data"];
    if (!data_eta.contains("edges") || !data_eta.contains("content")) return 1.0;

    const std::vector<float> eta_edges = data_eta["edges"].get<std::vector<float>>();
    int eta_bin = -1;
    for (size_t i = 0; i < eta_edges.size() - 1; ++i) {
        if (eta_for_lookup >= eta_edges[i] && eta_for_lookup < eta_edges[i+1]) {
            eta_bin = i;
            break;
        }
    }
    if (eta_bin == -1) return 1.0;

    const auto& data_pt = data_eta["content"][eta_bin];
    if (!data_pt.contains("edges") || !data_pt.contains("content")) return 1.0;

    const std::vector<float> pt_edges = data_pt["edges"].get<std::vector<float>>();
    int pt_bin = -1;
    for (size_t i = 0; i < pt_edges.size() - 1; ++i) {
        if (pt >= pt_edges[i] && pt < pt_edges[i+1]) {
            pt_bin = i;
            break;
        }
    }
    if (pt_bin == -1) return 1.0;

    const auto& categories = data_pt["content"][pt_bin]["content"];
    for (const auto& entry : categories) {
        if (entry.contains("key") && entry["key"] == "nominal" && entry.contains("value")) {
            return entry["value"].get<float>();
        }
    }
    
    return 1.0;
}
