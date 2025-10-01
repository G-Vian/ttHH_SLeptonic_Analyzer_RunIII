#include "ttHHanalyzer_trigger.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include <iostream>
#include <vector>

using namespace std;

// =================================================================================
//          FUNÇÕES PARA ELECTRON SCALE FACTORS
// =================================================================================

// ==========================
// Electron Trigger SF (MODIFICADO PARA USAR O MÉTODO PADRÃO DO ROOT)
// ==========================
void ttHHanalyzer::initTriggerSF() {
    TString sfFilePath;

    if (_year == "2022") {
        sfFilePath = "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoBCD/tnpEleHLT/HLT_SF_Ele30_MVAiso90ID/egammaEffi.txt_EGM2D.root";
    } else if (_year == "2022EE") {
        sfFilePath = "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoE+PromptFG/tnpEleHLT/HLT_SF_Ele30_MVAiso90ID/egammaEffi.txt_EGM2D.root";
    } else if (_year == "2023") {
        sfFilePath = "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23C/tnpEleHLT/HLT_SF_Ele30_MVAiso90ID/egammaEffi.txt_EGM2D.root";
    } else if (_year == "2023B" || _year == "2024") {
        sfFilePath = "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/tnpEleHLT/HLT_SF_Ele30_MVAiso90ID/egammaEffi.txt_EGM2D.root";
    } else {
        std::cerr << "[initTriggerSF] Ano não suportado para electron SF: " << _year << std::endl;
        return;
    }
    
    // Limpa ponteiros antigos
    if (h2_eleTrigSF) { delete h2_eleTrigSF; h2_eleTrigSF = nullptr; }
    if (h2_eleTrigSF_unc) { delete h2_eleTrigSF_unc; h2_eleTrigSF_unc = nullptr; }

    TFile* tempFile = TFile::Open(sfFilePath, "READ");
    if (!tempFile || tempFile->IsZombie()) {
        std::cerr << "[initTriggerSF] Erro ao abrir arquivo de electron SF: " << sfFilePath << std::endl;
        if (tempFile) delete tempFile;
        return;
    }

    // Carrega o SF principal usando Clone() e SetDirectory(0) para evitar crashes
    if (TH2F* tempSF = dynamic_cast<TH2F*>(tempFile->Get("EGamma_SF2D"))) {
        h2_eleTrigSF = (TH2F*)tempSF->Clone("h2_eleTrigSF");
        h2_eleTrigSF->SetDirectory(0);
    } else {
        std::cerr << "[initTriggerSF] AVISO: EGamma_SF2D não encontrado em " << sfFilePath << "!" << std::endl;
    }
    
    // Carrega a incerteza (opcional, mas boa prática)
    // O nome do histograma de incerteza pode variar, verifique o seu arquivo
    if (TH2F* tempUnc = dynamic_cast<TH2F*>(tempFile->Get("sys"))) { // O nome pode ser "unc" ou "err"
        h2_eleTrigSF_unc = (TH2F*)tempUnc->Clone("h2_eleTrigSF_unc");
        h2_eleTrigSF_unc->SetDirectory(0);
    }

    tempFile->Close();
    delete tempFile;

    std::cout << "[initTriggerSF] SF de elétrons (Trigger) carregado com sucesso (método ROOT)!" << std::endl;
}

float ttHHanalyzer::getEleTrigSF(float eta, float pt, float& sf_unc) {
    if (!h2_eleTrigSF) {
        sf_unc = 0.0f;
        return 1.0f;
    }

    int binX = h2_eleTrigSF->GetXaxis()->FindBin(eta);
    int binY = h2_eleTrigSF->GetYaxis()->FindBin(pt);

    // Proteção contra valores fora da faixa
    if (binX == 0 || binX > h2_eleTrigSF->GetNbinsX() || binY == 0 || binY > h2_eleTrigSF->GetNbinsY()) {
        sf_unc = 0.0f;
        return 1.0f;
    }
    
    float sf = h2_eleTrigSF->GetBinContent(binX, binY);
    sf_unc   = h2_eleTrigSF_unc ? h2_eleTrigSF_unc->GetBinContent(binX, binY) : 0.0f;

    return (sf > 0) ? sf : 1.0f;
}


// ==========================
// Electron Reco SF
// ==========================
void ttHHanalyzer::initRecoSF() {
    auto safeDelete = [](auto*& ptr) { if (ptr) { delete ptr; ptr = nullptr; } };
    safeDelete(h2_eleRecoSF_low);
    safeDelete(h2_eleRecoSF_mid);
    safeDelete(h2_eleRecoSF_high);

    auto loadReco = [&](const TString& path, const TString& tag) -> TH2F* {
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "[initRecoSF] Erro ao abrir arquivo: " << path << std::endl;
            if (f) delete f;
            return nullptr;
        }
        TH2F* h = nullptr;
        if (auto* temp = dynamic_cast<TH2F*>(f->Get("EGamma_SF2D"))) {
            h = (TH2F*)temp->Clone("h2_" + tag);
            h->SetDirectory(0);
        } else {
            std::cerr << "[initRecoSF] EGamma_SF2D não encontrado em " << path << std::endl;
        }
        f->Close();
        delete f;
        return h;
    };

    if (_year == "2024") {
        h2_eleRecoSF_low  = nullptr;
        h2_eleRecoSF_mid  = loadReco("/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/EleReco/midPt/egammaEffi.txt_EGM2D.root", "eleRecoSF_mid_2024");
        h2_eleRecoSF_high = loadReco("/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/EleReco/highPt/egammaEffi.txt_EGM2D.root", "eleRecoSF_high_2024");
    } else if (_year == "2023") {
        h2_eleRecoSF_low  = loadReco("/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/SF_prompt_2023_19012024/lowpT/Run3_2023C_New_lowpT_mergeEta_Added_symmetrizationsystEta_29052024/passingRECO/egammaEffi.txt_EGM2D.root", "eleRecoSF_low_2023");
        h2_eleRecoSF_mid  = loadReco("/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/SF_prompt_2023_19012024/midpT/Run3_2023C_New_midpT2_eta/passingRECO/egammaEffi.txt_EGM2D.root", "eleRecoSF_mid_2023");
        h2_eleRecoSF_high = loadReco("/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/SF_prompt_2023_19012024/highpT/Run3_2023C_New_highpt1_eta/passingRECO/egammaEffi.txt_EGM2D.root", "eleRecoSF_high_2023");
    } // ... (o resto do seu código para outros anos)
    
    std::cout << "[initRecoSF] SF de Reco de elétrons carregado para ano " << _year << std::endl;
}

float ttHHanalyzer::getEleRecoSF(float eta, float pt, float& sf_unc) {
    TH2F* h = nullptr;
    if (_year == "2024") {
        if (pt >= 20.f && pt < 75.f)   h = h2_eleRecoSF_mid;
        else if (pt >= 75.f)           h = h2_eleRecoSF_high;
    } else {
        if (pt > 10.f && pt < 20.f)    h = h2_eleRecoSF_low;
        else if (pt >= 20.f && pt < 75.f)   h = h2_eleRecoSF_mid;
        else if (pt >= 75.f)           h = h2_eleRecoSF_high;
    }
    if (!h) {
        sf_unc = 0.f;
        return 1.f;
    }
    int binX = h->GetXaxis()->FindBin(eta);
    int binY = h->GetYaxis()->FindBin(pt);
    float sf = h->GetBinContent(binX, binY);
    sf_unc   = h->GetBinError(binX, binY);
    if (sf <= 0.f) {
        sf = 1.f;
        sf_unc = 0.f;
    }
    return sf;
}

// ==========================
// Electron ID SF (MVA ISO WP90)
// ==========================
void ttHHanalyzer::initEleIDSF() {
    auto safeDelete = [](auto*& ptr){ if(ptr){ delete ptr; ptr=nullptr; } };
    safeDelete(h2_eleIDSF_2023B_Hole);
    safeDelete(h2_eleIDSF_2023B_NoHole);
    safeDelete(h2_eleIDSF_others);

    auto loadSF = [&](const TString& path, const TString& tag) -> TH2F* {
        TFile* f = TFile::Open(path, "READ");
        if(!f || f->IsZombie()){
            std::cerr << "[initEleIDSF] Erro ao abrir arquivo: " << path << std::endl;
            if(f) delete f;
            return nullptr;
        }
        TH2F* h = nullptr;
        if(auto* temp = dynamic_cast<TH2F*>(f->Get("EGamma_SF2D"))){
            h = (TH2F*)temp->Clone("h2_" + tag);
            h->SetDirectory(0);
        } else {
            std::cerr << "[initEleIDSF] EGamma_SF2D não encontrado em " << path << std::endl;
        }
        f->Close();
        delete f;
        return h;
    };

    if(_year == "2024"){
        h2_eleIDSF_others = loadSF("/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/EleID/passingMVA122Xwp90isoV1/merged_EGamma_SF2D_wp90iso.root", "eleIDSF_2024");
    } // ... (o resto do seu código para outros anos)

    std::cout << "[initEleIDSF] SF de Electron ID (MVA ISO WP90) carregado para ano " << _year << std::endl;
}

float ttHHanalyzer::getEleIDSF(float eta, float phi, float pt, float& sf_unc){
    TH2F* h = nullptr;
    if(_year == "2023B"){
        if(eta > -1.5 && eta < 0 && phi > -1.2 && phi < -0.8){
            h = h2_eleIDSF_2023B_Hole;
        } else {
            h = h2_eleIDSF_2023B_NoHole;
        }
    } else {
        h = h2_eleIDSF_others;
    }
    if(!h){
        sf_unc = 0.f;
        return 1.f;
    }
    int binX = h->GetXaxis()->FindBin(eta);
    int binY = h->GetYaxis()->FindBin(pt);
    float sf = h->GetBinContent(binX, binY);
    sf_unc = h->GetBinError(binX, binY);
    if(sf <= 0.f){
        sf = 1.f;
        sf_unc = 0.f;
    }
    return sf;
}
