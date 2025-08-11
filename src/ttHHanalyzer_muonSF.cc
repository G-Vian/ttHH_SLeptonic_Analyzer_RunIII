// ttHHanalyzer_muonSF.cc
#include "ttHHanalyzer_trigger.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;
#include <cstdlib>/// this is for MUON trigger SF  (TSFmu)
#include "json.hpp"// this is for MUON trigger SF (TSFmu)
using json = nlohmann::json;  /// this is for MUON trigger SF  (TSFmu)
json muonTrigSFJson; /// this is for MUON trigger SF  (TSFmu)



void ttHHanalyzer::initMuonHLTriggerSF() {
    TString repoPath = "muonefficiencies";
    TString sfFilePath;


    if (_year == "2022") {
        sfFilePath = repoPath + "/Run3/2022/2022_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2022_eta_pt_schemaV2.json";
    } else if (_year == "2022EE") {
        sfFilePath = repoPath + "/Run3/2022_EE/2022_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2022_EE_eta_pt_schemaV2.json";
    } else if (_year == "2023") {
        sfFilePath = repoPath + "/Run3/2023/2023_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2023_eta_pt_schemaV2.json";
    } else if (_year == "2023B") {
        sfFilePath = repoPath + "/Run3/2023_BPix/2023_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2023_BPix_eta_pt_schemaV2.json";
    } else if (_year == "2024") {
        sfFilePath = repoPath + "/Run3/2023_BPix/2023_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2023_BPix_eta_pt_schemaV2.json";
    } else {
        std::cerr << "Ano não suportado para SF de muons: " << _year << std::endl;
        return;
    }

    std::ifstream input(sfFilePath.Data());
    if (!input.is_open()) {
        std::cerr << "Erro ao abrir arquivo de SF: " << sfFilePath << std::endl;
        return;
    }

    std::string json_str((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    input.close();

    size_t pos = 0;
    while ((pos = json_str.find("Infinity", pos)) != std::string::npos) {
        json_str.replace(pos, 8, "1e10");
        pos += 4;
    }

    try {
        muonTrigSFJson = json::parse(json_str);
    } catch (const std::exception& e) {
        std::cerr << "Erro ao ler JSON: " << e.what() << std::endl;
        return;
    }

    if (h_sf_muon_vs_pt) { delete h_sf_muon_vs_pt; h_sf_muon_vs_pt = nullptr; }
    if (h_sf_muon_vs_eta) { delete h_sf_muon_vs_eta; h_sf_muon_vs_eta = nullptr; }
    if (h_sf_muon_vs_pt_sum) { delete h_sf_muon_vs_pt_sum; h_sf_muon_vs_pt_sum = nullptr; }
    if (h_sf_muon_vs_pt_count) { delete h_sf_muon_vs_pt_count; h_sf_muon_vs_pt_count = nullptr; }
    if (h_sf_muon_vs_eta_sum) { delete h_sf_muon_vs_eta_sum; h_sf_muon_vs_eta_sum = nullptr; }
    if (h_sf_muon_vs_eta_count) { delete h_sf_muon_vs_eta_count; h_sf_muon_vs_eta_count = nullptr; }
    if (h_sf_muon_vs_pt_avg) { delete h_sf_muon_vs_pt_avg; h_sf_muon_vs_pt_avg = nullptr; }
    if (h_sf_muon_vs_eta_avg) { delete h_sf_muon_vs_eta_avg; h_sf_muon_vs_eta_avg = nullptr; }

    h_sf_muon_vs_pt        = new TH1F("h_sf_muon_vs_pt", "Muon SF vs pT;Muon pT [GeV];SF", 100, 0, 700);
    h_sf_muon_vs_eta       = new TH1F("h_sf_muon_vs_eta", "Muon SF vs Eta;Muon #eta;SF", 100, -5, 5);
    h_sf_muon_vs_pt_sum    = new TH1F("h_sf_muon_vs_pt_sum", "Sum SF vs pT", 100, 0, 700);
    h_sf_muon_vs_pt_count  = new TH1F("h_sf_muon_vs_pt_count", "Count SF vs pT", 100, 0, 700);
    h_sf_muon_vs_eta_sum   = new TH1F("h_sf_muon_vs_eta_sum", "Sum SF vs eta", 100, -5, 5);
    h_sf_muon_vs_eta_count = new TH1F("h_sf_muon_vs_eta_count", "Count SF vs eta", 100, -5, 5);
    h_sf_muon_vs_pt_avg    = new TH1F("h_sf_muon_vs_pt_avg", "Avg SF vs pT", 100, 0, 700);
    h_sf_muon_vs_eta_avg   = new TH1F("h_sf_muon_vs_eta_avg", "Avg SF vs eta", 100, -5, 5);

    std::vector<TH1*> hists = {
        h_sf_muon_vs_pt, h_sf_muon_vs_eta,
        h_sf_muon_vs_pt_sum, h_sf_muon_vs_pt_count,
        h_sf_muon_vs_eta_sum, h_sf_muon_vs_eta_count,
        h_sf_muon_vs_pt_avg, h_sf_muon_vs_eta_avg
    };
    for (auto& h : hists) {
        if (h) h->SetDirectory(0);
    }

    std::cout << "SF de muons carregado com sucesso!" << std::endl;
}

float ttHHanalyzer::getMuonTrigSF(float eta, float pt) {
    if (muonTrigSFJson.empty()) return 1.0;

    const json* correction = nullptr;
    for (const auto& corr : muonTrigSFJson["corrections"]) {
        if (corr.contains("name") && corr["name"] == "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight") {
            correction = &corr;
            break;
        }
    }
    if (!correction) {
        std::cerr << "Correção NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight não encontrada no JSON!" << std::endl;
        return 1.0;
    }

    const auto& data_eta = (*correction)["data"];
    if (!data_eta.contains("edges") || !data_eta.contains("content")) {
        std::cerr << "Formato inválido na parte data do JSON." << std::endl;
        return 1.0;
    }

    const std::vector<float> eta_edges = data_eta["edges"].get<std::vector<float>>();
    int eta_bin = -1;
    for (size_t i = 0; i < eta_edges.size() - 1; ++i) {
        if (eta >= eta_edges[i] && eta < eta_edges[i+1]) {
            eta_bin = i;
            break;
        }
    }
    if (eta_bin == -1) {
        if (eta == eta_edges.back()) eta_bin = int(eta_edges.size()) - 2;
        else return 1.0;
    }

    const auto& data_pt = data_eta["content"][eta_bin];
    if (!data_pt.contains("edges") || !data_pt.contains("content")) {
        std::cerr << "Formato inválido na parte pt do JSON." << std::endl;
        return 1.0;
    }

    const std::vector<float> pt_edges = data_pt["edges"].get<std::vector<float>>();
    int pt_bin = -1;
    for (size_t i = 0; i < pt_edges.size() - 1; ++i) {
        if (pt >= pt_edges[i] && pt < pt_edges[i+1]) {
            pt_bin = i;
            break;
        }
    }
    if (pt_bin == -1) {
        if (pt == pt_edges.back()) pt_bin = int(pt_edges.size()) - 2;
        else return 1.0;
    }

    const auto& categories = data_pt["content"][pt_bin]["content"];
    for (const auto& entry : categories) {
        if (entry.contains("key") && entry["key"] == "nominal" && entry.contains("value")) {
            return entry["value"].get<float>();
        }
    }

    return 1.0;
}
