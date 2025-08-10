#include "ttHHanalyzer_trigger.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include "TH2.h"//Trigger SF for electron  (TSFel)
using namespace std;



// ==========================
// Electron Trigger SF
// ==========================
void ttHHanalyzer::initTriggerSF() {
    TString sfFilePath;
    sf_summary_log_file.open("SF_summary_log.txt");
    if (!sf_summary_log_file.is_open()) {
        std::cerr << "Erro ao abrir arquivo SF_summary_log.txt para escrita!" << std::endl;
    }
    if (_year == "2022") {
        sfFilePath = "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoBCD/tnpEleHLT/HLT_SF_Ele30_MVAiso90ID/egammaEffi.txt_EGM2D.root";
    } else if (_year == "2022EE") {
        sfFilePath = "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2022/ForRe-recoE+PromptFG/tnpEleHLT/HLT_SF_Ele30_MVAiso90ID/egammaEffi.txt_EGM2D.root";
    } else if (_year == "2023") {
        sfFilePath = "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23C/tnpEleHLT/HLT_SF_Ele30_MVAiso90ID/egammaEffi.txt_EGM2D.root";
    } else if (_year == "2023B") {
        sfFilePath = "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/tnpEleHLT/HLT_SF_Ele30_MVAiso90ID/egammaEffi.txt_EGM2D.root";
    } else if (_year == "2024") {
        sfFilePath = "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2023/ForPrompt23D/tnpEleHLT/HLT_SF_Ele30_MVAiso90ID/egammaEffi.txt_EGM2D.root";
    } else {
        std::cerr << "Unknown year for electron trigger SF! Year: " << _year << std::endl;
        return;
    }

    TFile* tempFile = TFile::Open(sfFilePath, "READ");
    if (!tempFile || tempFile->IsZombie()) {
        std::cerr << "Failed to open electron trigger SF file: " << sfFilePath << std::endl;
        if (tempFile) delete tempFile;
        return;
    }

    // Deleta anteriores, se existirem
    if (h_sf_vs_pt) delete h_sf_vs_pt;
    if (h_sf_vs_eta) delete h_sf_vs_eta;
    if (h_effMC_vs_pt) delete h_effMC_vs_pt;
    if (h_effMC_vs_eta) delete h_effMC_vs_eta;
    if (h2_effMC) delete h2_effMC;
    if (h2_eleTrigSF) delete h2_eleTrigSF;
    if (h2_eleTrigSF_unc) delete h2_eleTrigSF_unc;
    if (h_sf_vs_pt_sum) delete h_sf_vs_pt_sum;
    if (h_sf_vs_pt_count) delete h_sf_vs_pt_count;
    if (h_sf_vs_eta_sum) delete h_sf_vs_eta_sum;
    if (h_sf_vs_eta_count) delete h_sf_vs_eta_count;
    if (h_effMC_vs_pt_sum) delete h_effMC_vs_pt_sum;
    if (h_effMC_vs_pt_count) delete h_effMC_vs_pt_count;
    if (h_effMC_vs_eta_sum) delete h_effMC_vs_eta_sum;
    if (h_effMC_vs_eta_count) delete h_effMC_vs_eta_count;
    if (h_sf_vs_pt_avg) delete h_sf_vs_pt_avg;
    if (h_sf_vs_eta_avg) delete h_sf_vs_eta_avg;
    if (h_effMC_vs_pt_avg) delete h_effMC_vs_pt_avg;
    if (h_effMC_vs_eta_avg) delete h_effMC_vs_eta_avg;

    h_sf_vs_pt = new TH1F("h_sf_vs_pt", "SF vs pT;Electron pT [GeV];SF", 100000, 0, 700);
    h_sf_vs_eta = new TH1F("h_sf_vs_eta", "SF vs Eta;Electron #eta;SF", 100000, -5, 5);
    h_effMC_vs_pt = new TH1F("h_effMC_vs_pt", "EffMC vs pT;Electron pT [GeV];Eff.", 100000, 0, 700);
    h_effMC_vs_eta = new TH1F("h_effMC_vs_eta", "EffMC vs Eta;Electron #eta;Eff.", 100000, -5, 5);
    h_sf_vs_pt_sum = new TH1F("h_sf_vs_pt_sum", "Sum SF vs pT", 100000, 0, 700);
    h_sf_vs_pt_count = new TH1F("h_sf_vs_pt_count", "Count SF vs pT", 100000, 0, 700);
    h_sf_vs_eta_sum = new TH1F("h_sf_vs_eta_sum", "Sum SF vs eta", 100000, -5, 5);
    h_sf_vs_eta_count = new TH1F("h_sf_vs_eta_count", "Count SF vs eta", 100000, -5, 5);
    h_effMC_vs_pt_sum = new TH1F("h_effMC_vs_pt_sum", "Sum Eff vs pT", 100000, 0, 700);
    h_effMC_vs_pt_count = new TH1F("h_effMC_vs_pt_count", "Count Eff vs pT", 100000, 0, 700);
    h_effMC_vs_eta_sum = new TH1F("h_effMC_vs_eta_sum", "Sum Eff vs eta", 100000, -5, 5);
    h_effMC_vs_eta_count = new TH1F("h_effMC_vs_eta_count", "Count Eff vs eta", 100000, -5, 5);

    // SetDirectory(0)
    std::vector<TH1*> hists = {
        h_sf_vs_pt, h_sf_vs_eta, h_effMC_vs_pt, h_effMC_vs_eta,
        h_sf_vs_pt_sum, h_sf_vs_pt_count, h_sf_vs_eta_sum, h_sf_vs_eta_count,
        h_effMC_vs_pt_sum, h_effMC_vs_pt_count, h_effMC_vs_eta_sum, h_effMC_vs_eta_count
    };
    for (auto& h : hists) {
        if (h) h->SetDirectory(0);
    }

    // SF central
    TH2F* tempSF = dynamic_cast<TH2F*>(tempFile->Get("EGamma_SF2D"));
    if (tempSF) {
        h2_eleTrigSF = (TH2F*)tempSF->Clone("h2_eleTrigSF");
        h2_eleTrigSF->SetDirectory(0);
    }

    // Incerteza
    TString uncHistName = (_DataOrMC == "Data") ? "statData" : (_DataOrMC == "MC") ? "statMC" : "";
    if (uncHistName != "") {
        TH2F* tempUnc = dynamic_cast<TH2F*>(tempFile->Get(uncHistName));
        if (tempUnc) {
            h2_eleTrigSF_unc = (TH2F*)tempUnc->Clone("h2_eleTrigSF_unc");
            h2_eleTrigSF_unc->SetDirectory(0);
        }
    }

    // Eficiência MC
    TH2F* tempEffMC = dynamic_cast<TH2F*>(tempFile->Get("EGamma_EffMC2D"));
    if (tempEffMC) {
        h2_effMC = (TH2F*)tempEffMC->Clone("h2_effMC");
        h2_effMC->SetDirectory(0);
    }

    tempFile->Close();
    delete tempFile;
    gROOT->cd();
}

float ttHHanalyzer::getEleTrigSF(float eta, float pt, float& sf_unc) {
    if (!h2_eleTrigSF || !h2_eleTrigSF_unc) {
        sf_unc = 0.;
        return 1.;
    }

    int binX = h2_eleTrigSF->GetXaxis()->FindBin(eta);
    int binY = h2_eleTrigSF->GetYaxis()->FindBin(pt);

    float sf = h2_eleTrigSF->GetBinContent(binX, binY);
    sf_unc = h2_eleTrigSF_unc->GetBinContent(binX, binY);

    return sf;
}





