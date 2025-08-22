#include "ttHHanalyzer_trigger.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include "TH2.h" // Trigger SF for electron (TSFel)
#include <vector>

using namespace std;

// ==========================
// Electron Trigger SF
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
        // usar mesmo arquivo até sair oficial de 2024
        sfFilePath = "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/tnpEleHLT/HLT_SF_Ele30_MVAiso90ID/egammaEffi.txt_EGM2D.root";
    } else {
        std::cerr << "[initTriggerSF] Ano não suportado para electron SF: " << _year << std::endl;
        return;
    }

    TFile* tempFile = TFile::Open(sfFilePath, "READ");
    if (!tempFile || tempFile->IsZombie()) {
        std::cerr << "[initTriggerSF] Erro ao abrir arquivo de electron SF: " << sfFilePath << std::endl;
        if (tempFile) delete tempFile;
        return;
    }

    // Safe delete
    auto safeDelete = [](auto*& ptr) { if (ptr) { delete ptr; ptr = nullptr; } };
    safeDelete(h_sf_vs_pt);    safeDelete(h_sf_vs_eta);
    safeDelete(h_effMC_vs_pt); safeDelete(h_effMC_vs_eta);
    safeDelete(h2_effMC);      safeDelete(h2_eleTrigSF); safeDelete(h2_eleTrigSF_unc);

    safeDelete(h_sf_vs_pt_sum);    safeDelete(h_sf_vs_pt_count);
    safeDelete(h_sf_vs_eta_sum);   safeDelete(h_sf_vs_eta_count);
    safeDelete(h_effMC_vs_pt_sum); safeDelete(h_effMC_vs_pt_count);
    safeDelete(h_effMC_vs_eta_sum); safeDelete(h_effMC_vs_eta_count);

    safeDelete(h_sf_vs_pt_avg);    safeDelete(h_sf_vs_eta_avg);
    safeDelete(h_effMC_vs_pt_avg); safeDelete(h_effMC_vs_eta_avg);

    // Cria histogramas auxiliares (apenas para monitorar Trigger; Reco não cria)
    h_sf_vs_pt        = new TH1F("h_sf_vs_pt", "Electron SF vs pT;Electron pT [GeV];SF", 100, 0, 700);
    h_sf_vs_eta       = new TH1F("h_sf_vs_eta", "Electron SF vs Eta;Electron #eta;SF", 100, -5, 5);
    h_effMC_vs_pt     = new TH1F("h_effMC_vs_pt", "Electron EffMC vs pT;Electron pT [GeV];Eff.", 100, 0, 700);
    h_effMC_vs_eta    = new TH1F("h_effMC_vs_eta", "Electron EffMC vs Eta;Electron #eta;Eff.", 100, -5, 5);

    h_sf_vs_pt_sum    = new TH1F("h_sf_vs_pt_sum", "Sum Electron SF vs pT", 100, 0, 700);
    h_sf_vs_pt_count  = new TH1F("h_sf_vs_pt_count", "Count Electron SF vs pT", 100, 0, 700);
    h_sf_vs_eta_sum   = new TH1F("h_sf_vs_eta_sum", "Sum Electron SF vs eta", 100, -5, 5);
    h_sf_vs_eta_count = new TH1F("h_sf_vs_eta_count", "Count Electron SF vs eta", 100, -5, 5);

    h_effMC_vs_pt_sum    = new TH1F("h_effMC_vs_pt_sum", "Sum Electron Eff vs pT", 100, 0, 700);
    h_effMC_vs_pt_count  = new TH1F("h_effMC_vs_pt_count", "Count Electron Eff vs pT", 100, 0, 700);
    h_effMC_vs_eta_sum   = new TH1F("h_effMC_vs_eta_sum", "Sum Electron Eff vs eta", 100, -5, 5);
    h_effMC_vs_eta_count = new TH1F("h_effMC_vs_eta_count", "Count Electron Eff vs eta", 100, -5, 5);

    for (auto* h : {h_sf_vs_pt_sum, h_sf_vs_pt_count, h_sf_vs_eta_sum, h_sf_vs_eta_count,
                    h_effMC_vs_pt_sum, h_effMC_vs_pt_count, h_effMC_vs_eta_sum, h_effMC_vs_eta_count}) {
        if (h) h->Sumw2();
    }

    std::vector<TH1*> hists = {
        h_sf_vs_pt, h_sf_vs_eta, h_effMC_vs_pt, h_effMC_vs_eta,
        h_sf_vs_pt_sum, h_sf_vs_pt_count, h_sf_vs_eta_sum, h_sf_vs_eta_count,
        h_effMC_vs_pt_sum, h_effMC_vs_pt_count, h_effMC_vs_eta_sum, h_effMC_vs_eta_count
    };
    for (auto& h : hists) if (h) h->SetDirectory(0);

    if (TH2F* tempSF = dynamic_cast<TH2F*>(tempFile->Get("EGamma_SF2D"))) {
        h2_eleTrigSF = (TH2F*)tempSF->Clone("h2_eleTrigSF");
        h2_eleTrigSF->SetDirectory(0);
    } else {
        std::cerr << "[initTriggerSF] WARNING: EGamma_SF2D não encontrado!" << std::endl;
    }

    TString uncHistName = (_DataOrMC == "Data") ? "statData" : (_DataOrMC == "MC") ? "statMC" : "";
    if (uncHistName != "") {
        if (TH2F* tempUnc = dynamic_cast<TH2F*>(tempFile->Get(uncHistName))) {
            h2_eleTrigSF_unc = (TH2F*)tempUnc->Clone("h2_eleTrigSF_unc");
            h2_eleTrigSF_unc->SetDirectory(0);
        } else {
            std::cerr << "[initTriggerSF] WARNING: histograma " << uncHistName << " não encontrado!" << std::endl;
        }
    }

    if (TH2F* tempEffMC = dynamic_cast<TH2F*>(tempFile->Get("EGamma_EffMC2D"))) {
        h2_effMC = (TH2F*)tempEffMC->Clone("h2_effMC");
        h2_effMC->SetDirectory(0);
    }

    tempFile->Close();
    delete tempFile;
    gROOT->cd();

    std::cout << "[initTriggerSF] SF de elétrons (Trigger) carregado com sucesso!" << std::endl;
}

float ttHHanalyzer::getEleTrigSF(float eta, float pt, float& sf_unc) {
    if (!h2_eleTrigSF) {
        std::cerr << "[getEleTrigSF] WARNING: h2_eleTrigSF não inicializado!" << std::endl;
        sf_unc = 0.;
        return 1.;
    }

    int binX = h2_eleTrigSF->GetXaxis()->FindBin(eta);
    int binY = h2_eleTrigSF->GetYaxis()->FindBin(pt);

    float sf = h2_eleTrigSF->GetBinContent(binX, binY);
    sf_unc   = h2_eleTrigSF_unc ? h2_eleTrigSF_unc->GetBinContent(binX, binY) : 0.;

    return sf;
}

// ==========================
// Electron Reco SF (sem histos auxiliares)
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
        gROOT->cd();
        return h;
    };

    // Carrega conforme ano
    if (_year == "2024") {
        h2_eleRecoSF_mid = loadReco(
            "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/EleReco/midPt/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_mid_2024"
        );
    }
    else if (_year == "2023") {
        h2_eleRecoSF_low  = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/SF_prompt_2023_19012024/lowpT/Run3_2023C_New_lowpT_mergeEta_Added_symmetrizationsystEta_29052024/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_low_2023"
        );
        h2_eleRecoSF_mid  = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/SF_prompt_2023_19012024/midpT/Run3_2023C_New_midpT2_eta/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_mid_2023"
        );
        h2_eleRecoSF_high = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/SF_prompt_2023_19012024/highpT/Run3_2023C_New_highpt1_eta/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_high_2023"
        );
    }
    else if (_year == "2023B") {
        h2_eleRecoSF_low  = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/SF_prompt_2023_19012024/lowpT/Run3_2023D_New_lowpT_mergeEta_Added_symmetrizationsystEta_29052024/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_low_2023B"
        );
        h2_eleRecoSF_mid  = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/SF_prompt_2023_19012024/midpT/Run3_2023D_New_midpT_eta2/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_mid_2023B"
        );
        h2_eleRecoSF_high = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/SF_prompt_2023_19012024/highpT/Run3_2023D_New_highpt_eta2/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_high_2023B"
        );
    }
    else if (_year == "2022") {
        h2_eleRecoSF_low  = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/New_SF_19122023/lowpT/Run3_2022BCD_PreEEMC/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_low_2022"
        );
        h2_eleRecoSF_mid  = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/New_SF_19122023/midpT/Run3_2022BCD_New_midpT7/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_mid_2022"
        );
        h2_eleRecoSF_high = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/New_SF_19122023/highpT/Run3_2022BCD_New_highpt6/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_high_2022"
        );
    }
    else if (_year == "2022EE") {
        h2_eleRecoSF_low  = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/New_SF_19122023/lowpT/Run3_2022EFG_PostEEMC/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_low_2022EE"
        );
        h2_eleRecoSF_mid  = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/New_SF_19122023/midpT/Run3_2022EFG_New_midpT5/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_mid_2022EE"
        );
        h2_eleRecoSF_high = loadReco(
            "/eos/cms/store/group/phys_egamma/validation/web/Run3_egm_reco_SF/New_SF_19122023/highpT/Run3_2022EFG_New_highpt5/passingRECO/egammaEffi.txt_EGM2D.root",
            "eleRecoSF_high_2022EE"
        );
    }
    else {
        std::cerr << "[initRecoSF] Ano não suportado: " << _year << std::endl;
    }

    std::cout << "[initRecoSF] SF de Reco de elétrons carregado para ano " << _year << std::endl;
}

float ttHHanalyzer::getEleRecoSF(float eta, float pt, float& sf_unc) {
    TH2F* h = nullptr;

    // Seleção de hist por pT conforme ano
    if (_year == "2024") {
        // Disponível apenas midPt 20 < pT < 75 GeV
        if (pt > 20.f && pt < 75.f) {
            h = h2_eleRecoSF_mid;
        }
    } else {
        if (pt > 10.f && pt < 20.f)         h = h2_eleRecoSF_low;
        else if (pt >= 20.f && pt < 75.f)   h = h2_eleRecoSF_mid;
        else if (pt >= 75.f)                h = h2_eleRecoSF_high;
    }

    if (!h) {
        std::cerr << "[getEleRecoSF] Nenhum histograma correspondente para pT=" << pt
                  << " GeV e ano=" << _year << ". Retornando SF=1." << std::endl;
        sf_unc = 0.f;
        return 1.f;
    }

    // Faz lookup no TH2
    int binX = h->GetXaxis()->FindBin(eta);
    int binY = h->GetYaxis()->FindBin(pt);

    float sf = h->GetBinContent(binX, binY);
    // Usa o erro do próprio hist como incerteza (usualmente contém estatística; syst pode não estar embutida)
    sf_unc   = h->GetBinError(binX, binY);

    // Proteção contra bins vazios/zeros (opcional)
    if (sf <= 0.f) {
        std::cerr << "[getEleRecoSF] Bin vazio ou SF<=0 (eta=" << eta << ", pt=" << pt
                  << "). Retornando SF=1." << std::endl;
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
        gROOT->cd();
        return h;
    };

    if(_year == "2024"){
        h2_eleIDSF_others = loadSF(
            "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/EleID/passingMVA122Xwp90isoV1/merged_EGamma_SF2D_wp90iso.root",
            "eleIDSF_2024"
        );
    }
    else if(_year == "2023"){
        h2_eleIDSF_others = loadSF(
            "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23C/tnpEleID/passingMVA122Xwp90isoV1/egammaEffi.txt_EGM2D.root",
            "eleIDSF_2023"
        );
    }
    else if(_year == "2023B"){
        h2_eleIDSF_2023B_NoHole = loadSF(
            "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/tnpEleID/passingMVA122Xwp90isoV1/NoHole_egammaEffi.txt_EGM2D.root",
            "eleIDSF_2023B_NoHole"
        );
        h2_eleIDSF_2023B_Hole = loadSF(
            "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/tnpEleID/passingMVA122Xwp90isoV1/Hole_egammaEffi.txt_EGM2D.root",
            "eleIDSF_2023B_Hole"
        );
    }
    else if(_year == "2022"){
        h2_eleIDSF_others = loadSF(
            "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoBCD/tnpEleID/passingMVA122Xwp90isoV1/egammaEffi.txt_EGM2D.root",
            "eleIDSF_2022"
        );
    }
    else if(_year == "2022EE"){
        h2_eleIDSF_others = loadSF(
            "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoE+PromptFG/tnpEleID/passingMVA122Xwp90isoV1/egammaEffi.txt_EGM2D.root",
            "eleIDSF_2022EE"
        );
    }
    else{
        std::cerr << "[initEleIDSF] Ano não suportado: " << _year << std::endl;
    }

    std::cout << "[initEleIDSF] SF de Electron ID (MVA ISO WP90) carregado para ano " << _year << std::endl;
}

float ttHHanalyzer::getEleIDSF(float eta, float phi, float pt, float& sf_unc){
    TH2F* h = nullptr;

    if(_year == "2023B"){
        // Caso especial: -1.5 < eta < 0 e -1.2 < phi < -0.8 usa Hole
        if(eta > -1.5 && eta < 0 && phi > -1.2 && phi < -0.8){
            h = h2_eleIDSF_2023B_Hole;
        } else {
            h = h2_eleIDSF_2023B_NoHole;
        }
    } else {
        h = h2_eleIDSF_others;
    }

    if(!h){
        std::cerr << "[getEleIDSF] Nenhum histograma disponível para SF, retornando 1." << std::endl;
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


