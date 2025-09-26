#include "ttHHanalyzer_trigger.h" // Inclua seu header principal

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>

// --- DEFINIÇÃO DAS VARIÁVEIS GLOBAIS DE SF ---
// A definição real (alocação de memória) das variáveis acontece aqui.
json muonTrigSFJson;
json muonHighPtTrigSFJson;
json muonLowPtIDSFJson;
json muonMediumPtIDSFJson;
json muonHighPtIDSFJson;
/**
 * @brief Função auxiliar para carregar, tratar e fazer o parse de um arquivo JSON.
 */
json loadSFJson(const TString& filePath) {
    std::ifstream input(filePath.Data());
    if (!input.is_open()) {
        std::cerr << "[loadSFJson] ERRO: Não foi possível abrir o arquivo de SF: " << filePath << std::endl;
        return json();
    }
    std::string json_str((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    input.close();
    size_t pos = 0;
    while ((pos = json_str.find("Infinity", pos)) != std::string::npos) {
        json_str.replace(pos, 8, "1e10");
        pos += 4;
    }
    try {
        return json::parse(json_str);
    } catch (const std::exception& e) {
        std::cerr << "[loadSFJson] ERRO ao fazer o parse do JSON de " << filePath << ": " << e.what() << std::endl;
        return json();
    }
}

/**
 * @brief Inicializa os objetos JSON com os fatores de escala para múons e os histogramas.
 */
void ttHHanalyzer::initMuonHLTriggerSF() {
    TString sfMediumPtFilePath;
    TString sfHighPtFilePath;
    TString localDir = "/afs/cern.ch/user/g/gvian/muon_SF/muonefficiencies/Run3/";

    // Define os caminhos dos arquivos JSON
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

    // Carrega os arquivos JSON
    std::cout << "[INFO] Carregando SF de Múon para pT Médio: " << sfMediumPtFilePath << std::endl;
    muonTrigSFJson = loadSFJson(sfMediumPtFilePath);
    if (!muonTrigSFJson.empty()) {
        std::cout << "[INFO] SF de Múon para pT Médio carregado com sucesso." << std::endl;
    }

    std::cout << "[INFO] Carregando SF de Múon para pT Alto: " << sfHighPtFilePath << std::endl;
    muonHighPtTrigSFJson = loadSFJson(sfHighPtFilePath);
    if (!muonHighPtTrigSFJson.empty()) {
        std::cout << "[INFO] SF de Múon para pT Alto carregado com sucesso." << std::endl;
    }

    // ==========================================
    // Limpa histos antigos
    // ==========================================
    if (h_sf_muon_vs_pt)       { delete h_sf_muon_vs_pt;       h_sf_muon_vs_pt = nullptr; }
    if (h_sf_muon_vs_eta)      { delete h_sf_muon_vs_eta;      h_sf_muon_vs_eta = nullptr; }
    if (h_sf_muon_vs_pt_sum)   { delete h_sf_muon_vs_pt_sum;   h_sf_muon_vs_pt_sum = nullptr; }
    if (h_sf_muon_vs_pt_count) { delete h_sf_muon_vs_pt_count; h_sf_muon_vs_pt_count = nullptr; }
    if (h_sf_muon_vs_eta_sum)  { delete h_sf_muon_vs_eta_sum;  h_sf_muon_vs_eta_sum = nullptr; }
    if (h_sf_muon_vs_eta_count){ delete h_sf_muon_vs_eta_count;h_sf_muon_vs_eta_count = nullptr; }
    if (h_sf_muon_vs_pt_avg)   { delete h_sf_muon_vs_pt_avg;   h_sf_muon_vs_pt_avg = nullptr; }
    if (h_sf_muon_vs_eta_avg)  { delete h_sf_muon_vs_eta_avg;  h_sf_muon_vs_eta_avg = nullptr; }

    // ==========================================
    // Cria histos novos
    // ==========================================
    h_sf_muon_vs_pt       = new TH1F("h_sf_muon_vs_pt", "Muon SF vs pT;Muon pT [GeV];SF", 100, 0, 700);
    h_sf_muon_vs_eta      = new TH1F("h_sf_muon_vs_eta", "Muon SF vs Eta;Muon #eta;SF", 100, -5, 5);
    h_sf_muon_vs_pt_sum   = new TH1F("h_sf_muon_vs_pt_sum", "Sum SF vs pT", 100, 0, 700);
    h_sf_muon_vs_pt_count = new TH1F("h_sf_muon_vs_pt_count", "Count SF vs pT", 100, 0, 700);
    h_sf_muon_vs_eta_sum  = new TH1F("h_sf_muon_vs_eta_sum", "Sum SF vs eta", 100, -5, 5);
    h_sf_muon_vs_eta_count= new TH1F("h_sf_muon_vs_eta_count", "Count SF vs eta", 100, -5, 5);
    h_sf_muon_vs_pt_avg   = new TH1F("h_sf_muon_vs_pt_avg", "Avg SF vs pT", 100, 0, 700);
    h_sf_muon_vs_eta_avg  = new TH1F("h_sf_muon_vs_eta_avg", "Avg SF vs eta", 100, -5, 5);

    // ==========================================
    // Sumw2 para incertezas
    // ==========================================
    if (h_sf_muon_vs_pt_sum)   h_sf_muon_vs_pt_sum->Sumw2();
    if (h_sf_muon_vs_pt_count) h_sf_muon_vs_pt_count->Sumw2();
    if (h_sf_muon_vs_eta_sum)  h_sf_muon_vs_eta_sum->Sumw2();
    if (h_sf_muon_vs_eta_count)h_sf_muon_vs_eta_count->Sumw2();

    // ==========================================
    // Protege histos contra "auto-delete" do ROOT
    // ==========================================
    std::vector<TH1*> hists = {
        h_sf_muon_vs_pt, h_sf_muon_vs_eta,
        h_sf_muon_vs_pt_sum, h_sf_muon_vs_pt_count,
        h_sf_muon_vs_eta_sum, h_sf_muon_vs_eta_count,
        h_sf_muon_vs_pt_avg, h_sf_muon_vs_eta_avg
    };
    for (auto& h : hists) {
        if (h) h->SetDirectory(0);
    }
}

/**
 * @brief Obtém o fator de escala (SF) de trigger para um dado múon.
 * A função seleciona a correção de médio/alto pT e trata o eta (normal/absoluto) internamente.
 * @param eta O eta do múon (pode ser positivo ou negativo).
 * @param pt O momento transverso (pT) do múon.
 * @return O valor do fator de escala. Retorna 1.0 se não encontrado.
 */
float ttHHanalyzer::getMuonTrigSF(float eta, float pt) {
    const json* sfJson = nullptr;
    std::string correctionName;
    float eta_for_lookup; // Variável para guardar o eta a ser usado na busca

    // --- Lógica de Seleção ---
    if (pt > 200.0) {
        sfJson = &muonHighPtTrigSFJson;
        correctionName = "NUM_HLT_DEN_TrkHighPtTightRelIsoProbes";
        // Para o SF de alto pT, o JSON é binnado em 'abseta', então usamos o valor absoluto.
        eta_for_lookup = fabs(eta);
    } else {
        sfJson = &muonTrigSFJson;
        correctionName = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight";
        // Para o SF de pT intermediário, usamos o valor original de eta.
        eta_for_lookup = eta;
    }

    if (!sfJson || sfJson->empty()) {
        return 1.0;
    }

    // Lógica para extrair o SF do JSON
    const json* correction = nullptr;
    for (const auto& corr : (*sfJson)["corrections"]) {
        if (corr.contains("name") && corr["name"] == correctionName) {
            correction = &corr;
            break;
        }
    }
    if (!correction) {
        std::cerr << "AVISO: Correção '" << correctionName << "' não encontrada no JSON para pt=" << pt << " e eta=" << eta << "!" << std::endl;
        return 1.0;
    }

    const auto& data_eta = (*correction)["data"];
    if (!data_eta.contains("edges") || !data_eta.contains("content")) {
        std::cerr << "ERRO: Formato inválido na parte 'data' do JSON." << std::endl;
        return 1.0;
    }

    // Encontra o bin de eta usando a variável 'eta_for_lookup'
    const std::vector<float> eta_edges = data_eta["edges"].get<std::vector<float>>();
    int eta_bin = -1;
    for (size_t i = 0; i < eta_edges.size() - 1; ++i) {
        if (eta_for_lookup >= eta_edges[i] && eta_for_lookup < eta_edges[i+1]) {
            eta_bin = i;
            break;
        }
    }
    if (eta_bin == -1) { // Trata o caso de estar exatamente na borda superior
        if (eta_for_lookup == eta_edges.back()) eta_bin = int(eta_edges.size()) - 2;
        else return 1.0;
    }

    const auto& data_pt = data_eta["content"][eta_bin];
    if (!data_pt.contains("edges") || !data_pt.contains("content")) {
        std::cerr << "ERRO: Formato inválido na parte 'pt' do JSON." << std::endl;
        return 1.0;
    }

    // Encontra o bin de pT
    const std::vector<float> pt_edges = data_pt["edges"].get<std::vector<float>>();
    int pt_bin = -1;
    for (size_t i = 0; i < pt_edges.size() - 1; ++i) {
        if (pt >= pt_edges[i] && pt < pt_edges[i+1]) {
            pt_bin = i;
            break;
        }
    }
    if (pt_bin == -1) { // Trata o caso de estar na borda superior ou acima
        if (pt >= pt_edges.back()) pt_bin = int(pt_edges.size()) - 2;
        else return 1.0;
    }

    // Extrai o valor nominal do SF
    const auto& categories = data_pt["content"][pt_bin]["content"];
    for (const auto& entry : categories) {
        if (entry.contains("key") && entry["key"] == "nominal" && entry.contains("value")) {
            return entry["value"].get<float>();
        }
    }

    return 1.0; // Retorno padrão caso algo falhe
}

////////////////// ID SF for muons /////////////////////////////////////
/**
 * @brief Inicializa os SFs de ID para múons, carregando os arquivos JSON para
 * as três faixas de pT (Low, Medium, High).
 */
void ttHHanalyzer::initMuonIDSF() {
    TString sfLowPtFilePath, sfMediumPtFilePath, sfHighPtFilePath;
    TString localDir = "/afs/cern.ch/user/g/gvian/muonefficiencies/Run3/";

    // Define os caminhos dos arquivos para cada ano
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
    } else {
        std::cerr << "[initMuonIDSF] ERRO: Ano não suportado para SF de ID de múons: " << _year << std::endl;
        return;
    }

    // Carrega os três arquivos JSON
    std::cout << "[INFO] Carregando SF de ID de Múon (Low pT): " << sfLowPtFilePath << std::endl;
    muonLowPtIDSFJson = loadSFJson(sfLowPtFilePath);

    std::cout << "[INFO] Carregando SF de ID de Múon (Medium pT): " << sfMediumPtFilePath << std::endl;
    muonMediumPtIDSFJson = loadSFJson(sfMediumPtFilePath);

    std::cout << "[INFO] Carregando SF de ID de Múon (High pT): " << sfHighPtFilePath << std::endl;
    muonHighPtIDSFJson = loadSFJson(sfHighPtFilePath);
}

/**
 * @brief Obtém o SF de ID para um múon, selecionando o arquivo e o nome da 
 * correção corretos para as faixas de pT baixo, médio e alto.
 * A função também trata a dependência do ano para o uso de 'eta' vs 'abseta'.
 */
float ttHHanalyzer::getMuonIDSF(float eta, float pt) {
    const json* sfJson = nullptr;
    std::string correctionName;
    float eta_for_lookup;

    // ===================================================================
    // Seleciona o arquivo JSON, o nome da correção e o tipo de eta (eta/abseta)
    // de acordo com a faixa de pT do múon.
    // ===================================================================
    if (pt < 30.0) { // --- CASO DE PT BAIXO ---
        sfJson = &muonLowPtIDSFJson;
        correctionName = "NUM_TightID_DEN_TrackerMuons";
        eta_for_lookup = fabs(eta); // PT baixo sempre usa abseta
    
    } else if (pt < 200.0) { // --- CASO DE PT MÉDIO ---
        sfJson = &muonMediumPtIDSFJson;
 //       correctionName = "NUM_TightPFIso_DEN_TightID";
         correctionName = "NUM_TightID_DEN_TrackerMuons";
        // Lógica dependente do ano para pT médio, conforme sua verificação
        if (_year == "2023" || _year == "2023B" || _year == "2024") {
            eta_for_lookup = eta; // Anos recentes usam eta
        } else { // 2022 e 2022EE
            eta_for_lookup = fabs(eta); // Anos mais antigos usam abseta
        }

    } else { // --- CASO DE PT ALTO (pt >= 200.0) ---
        sfJson = &muonHighPtIDSFJson;
        correctionName = "NUM_TightPFIso_DEN_TightID";
        eta_for_lookup = fabs(eta); // PT alto sempre usa abseta
    }

    // ===================================================================
    // Lógica padrão para buscar o valor no arquivo JSON selecionado
    // ===================================================================
    if (!sfJson || sfJson->empty()) return 1.0;

    const json* correction = nullptr;
    for (const auto& corr : (*sfJson)["corrections"]) {
        if (corr.contains("name") && corr["name"] == correctionName) {
            correction = &corr;
            break;
        }
    }
    if (!correction) {
        std::cerr << "AVISO: Correção de ID '" << correctionName << "' não encontrada no JSON para pt=" << pt << "!" << std::endl;
        return 1.0;
    }

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
    if (eta_bin == -1) {
        if (eta_for_lookup == eta_edges.back()) eta_bin = int(eta_edges.size()) - 2;
        else return 1.0;
    }

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
    if (pt_bin == -1) {
        if (pt >= pt_edges.back()) pt_bin = int(pt_edges.size()) - 2;
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
