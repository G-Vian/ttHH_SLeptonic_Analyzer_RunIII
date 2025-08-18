#include "ttHHanalyzer_trigger.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include "TH2.h" // Trigger SF for electron (TSFel)
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

    // ============================
    // Safe delete de anteriores
    auto safeDelete = [](auto*& ptr) { if (ptr) { delete ptr; ptr = nullptr; } };

    safeDelete(h_sf_vs_pt);   safeDelete(h_sf_vs_eta);
    safeDelete(h_effMC_vs_pt); safeDelete(h_effMC_vs_eta);
    safeDelete(h2_effMC);      safeDelete(h2_eleTrigSF); safeDelete(h2_eleTrigSF_unc);

    safeDelete(h_sf_vs_pt_sum);   safeDelete(h_sf_vs_pt_count);
    safeDelete(h_sf_vs_eta_sum);  safeDelete(h_sf_vs_eta_count);

    safeDelete(h_effMC_vs_pt_sum); safeDelete(h_effMC_vs_pt_count);
    safeDelete(h_effMC_vs_eta_sum); safeDelete(h_effMC_vs_eta_count);

    safeDelete(h_sf_vs_pt_avg);   safeDelete(h_sf_vs_eta_avg);
    safeDelete(h_effMC_vs_pt_avg); safeDelete(h_effMC_vs_eta_avg);

    // ============================
    // Cria histogramas novos
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

    // ============================
    // Ativa propagação de erros
    for (auto* h : {h_sf_vs_pt_sum, h_sf_vs_pt_count, h_sf_vs_eta_sum, h_sf_vs_eta_count,
                    h_effMC_vs_pt_sum, h_effMC_vs_pt_count, h_effMC_vs_eta_sum, h_effMC_vs_eta_count}) {
        if (h) h->Sumw2();
    }

    // SetDirectory(0) para todos
    std::vector<TH1*> hists = {
        h_sf_vs_pt, h_sf_vs_eta, h_effMC_vs_pt, h_effMC_vs_eta,
        h_sf_vs_pt_sum, h_sf_vs_pt_count, h_sf_vs_eta_sum, h_sf_vs_eta_count,
        h_effMC_vs_pt_sum, h_effMC_vs_pt_count, h_effMC_vs_eta_sum, h_effMC_vs_eta_count
    };
    for (auto& h : hists) if (h) h->SetDirectory(0);

    // ============================
    // SF central
    if (TH2F* tempSF = dynamic_cast<TH2F*>(tempFile->Get("EGamma_SF2D"))) {
        h2_eleTrigSF = (TH2F*)tempSF->Clone("h2_eleTrigSF");
        h2_eleTrigSF->SetDirectory(0);
    } else {
        std::cerr << "[initTriggerSF] WARNING: histograma EGamma_SF2D não encontrado no arquivo!" << std::endl;
    }

    // Incerteza
    TString uncHistName = (_DataOrMC == "Data") ? "statData" : (_DataOrMC == "MC") ? "statMC" : "";
    if (uncHistName != "") {
        if (TH2F* tempUnc = dynamic_cast<TH2F*>(tempFile->Get(uncHistName))) {
            h2_eleTrigSF_unc = (TH2F*)tempUnc->Clone("h2_eleTrigSF_unc");
            h2_eleTrigSF_unc->SetDirectory(0);
        } else {
            std::cerr << "[initTriggerSF] WARNING: histograma " << uncHistName << " não encontrado!" << std::endl;
        }
    }

    // Eficiência MC
    if (TH2F* tempEffMC = dynamic_cast<TH2F*>(tempFile->Get("EGamma_EffMC2D"))) {
        h2_effMC = (TH2F*)tempEffMC->Clone("h2_effMC");
        h2_effMC->SetDirectory(0);
    }

    tempFile->Close();
    delete tempFile;
    gROOT->cd();

    std::cout << "[initTriggerSF] SF de elétrons carregado com sucesso!" << std::endl;
}

// ==========================
// Getter com incerteza
// ==========================
float ttHHanalyzer::getEleTrigSF(float eta, float pt, float& sf_unc) {
    if (!h2_eleTrigSF) {
        std::cerr << "[getEleTrigSF] WARNING: h2_eleTrigSF não inicializado!" << std::endl;
        sf_unc = 0.;
        return 1.;
    }

    int binX = h2_eleTrigSF->GetXaxis()->FindBin(eta);
    int binY = h2_eleTrigSF->GetYaxis()->FindBin(pt);

    float sf = h2_eleTrigSF->GetBinContent(binX, binY);

    if (h2_eleTrigSF_unc) {
        sf_unc = h2_eleTrigSF_unc->GetBinContent(binX, binY);
    } else {
        sf_unc = 0.;
    }

    return sf;
}





