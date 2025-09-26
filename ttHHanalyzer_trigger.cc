#include "tnm.h"
#include <cmath> 
#include <algorithm>
#include <vector>
#include <map>
//#include <filesystem> // Para criar diretórios
#include <TSystem.h> // ROOT utilitário para mkdir
#include "TVector3.h"
#include "ttHHanalyzer_trigger.h"
#include <iostream>
#include <fstream>
#include <iomanip> //Para usar setw e manipuladores de fluxo --> pro txt do cutflow ficar bonito
#include "TH2.h"//Trigger SF for electron  (TSFel)
using namespace std;
#include <cstdlib>/// this is for MUON trigger SF  (TSFmu)
//#include "json.hpp"// this is for MUON trigger SF (TSFmu)
//using json = nlohmann::json;  /// this is for MUON trigger SF  (TSFmu)
//json muonTrigSFJson; /// this is for MUON trigger SF  (TSFmu)


////Log of selection and SF///////////
static std::unique_ptr<std::ofstream> event_log_file;
static std::unique_ptr<std::ofstream> sf_log_file;
static std::unique_ptr<std::ofstream> sf_summary_log_file;

//static int event_counter = 0;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////cuts//////////////////

map<std::string, float> cut { 
    {"nJets", 5} // nJets higher than 
    , {"nLeptons", 1} // nLepton equals to (Single-Leptonic)
    //    , {"nVetoLeptons", 0} // nVetoLepton equals to
    , {"HT", 0}
    , {"MET", 20} // MET higher than
    , {"nbJets", 4} //nBjets higher than
    , {"jetPt", 30} // jet pT higher than
    , {"leadElePt", 30} // leadElectron pT higher than
    , {"leadMuonPt", 29} // leadMuon pT higher than
    , {"subLeadElePt", 15} // subLeadElectron pT higher than
    , {"subLeadMuonPt", 15} // subLeadMuon pT higher than
    // , {"vetoLepPt", 15} // lepton pT higher than
    , {"boostedJetPt", 10} // boostedJet pT higher than
    , {"hadHiggsPt", 20} // hadronic Higgs pT higher than
    , {"jetEta", 2.6} // jet eta higher than
    , {"eleEta", 2.5} // electron eta higher than
    , {"muonEta", 2.4} // muon eta higher than
    , {"boostedJetEta", 2.5} // boostedJet eta higher than
    , {"muonIso", 0.15} // muon isolation less than
    , {"eleIso", 0.15}  // ele isolation less than (we set = to the muon Isolation)
    , {"jetID", 6}   // pass tight and tightLepVeto ID
    , {"jetPUid", 4}   // pass loose cut fail tight and medium
    , {"bTagDisc", 0.80} // this is not the same as BTagDisc, with capital B. bTagDisc in this case is just a cut for the boosted jets.
    , {"trigger", 1} // trigger
    , {"filter", 1} // MET filter
    , {"pv", 1}}; // primary vertex  

/////////////////////////////////
//////////_DO NOT ERASE THIS_////
void ttHHanalyzer::performAnalysis(){
    loop(noSys, false);}
////////////////////////////////

void ttHHanalyzer::loop(sysName sysType, bool up) {

    // Cria pasta com o nome do processo, se não existir
    if (_sampleName == "nothing") {
        std::cerr << "[ERROR] SampleName is not defined before log initialization!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    TString logDir = _sampleName;
    if (gSystem->AccessPathName(logDir)) {   // Se não existe
        if (gSystem->mkdir(logDir, kTRUE) != 0) { // kTRUE permite criar recursivamente
            std::cerr << "[ERROR] Failed to create directory: " << logDir << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    // Inicializa os arquivos de log, se ainda não foram criados
    if (!event_log_file || !sf_log_file) {
      //  std::string log1_name = logDir + "/event_selection_log_" + _sampleName + ".txt";
      //  std::string log2_name = logDir + "/log_trigger_sf_" + _sampleName + ".txt";
		TString log1_name = logDir + "/event_selection_log_" + _sampleName + ".txt";
		TString log2_name = logDir + "/log_applied_sf_" + _sampleName + ".txt";

 //       event_log_file = std::make_unique<std::ofstream>(log1_name);
  //      sf_log_file    = std::make_unique<std::ofstream>(log2_name);
		event_log_file = std::make_unique<std::ofstream>(log1_name.Data());
		sf_log_file    = std::make_unique<std::ofstream>(log2_name.Data());



		
        if (!event_log_file->is_open() || !sf_log_file->is_open()) {
            std::cerr << "[ERROR] Failed to open log files for writing!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    if (!sf_summary_log_file) {
  //      std::string sf_summary_log_name = logDir + "/sf_summary_log_" + _sampleName + ".txt";
  //      sf_summary_log_file = std::make_unique<std::ofstream>(sf_summary_log_name);
			TString sf_summary_log_name = logDir + "/sf_summary_log_" + _sampleName + ".txt";
			sf_summary_log_file = std::make_unique<std::ofstream>(sf_summary_log_name.Data());
        if (!sf_summary_log_file->is_open()) {
            std::cerr << "[ERROR] Failed to open sf summary log file for writing!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }


    int nevents = _ev->size();

    // helper para imprimir o cutflow (periodicamente e no final)
    auto printCutflowSummary = [&](const char* tag){
        if (!event_log_file || !event_log_file->is_open()) return;

        auto get = [&](const char* k){ return static_cast<double>(cutflow[k]); };
        auto safe_eff = [&](double num, double den){
            return (den > 0.0) ? (100.0 * num / den) : 0.0;
        };

        std::ofstream& log = *event_log_file;
        log << "==== Cutflow summary with efficiencies"
            << (tag ? std::string(" (") + tag + ")" : std::string())
            << " ====" << std::endl;

        log << std::left << std::setw(20) << "No cut"
            << ": " << std::right << std::setw(10) << static_cast<long long>(get("noCut")) << std::endl;

        log << std::left << std::setw(20) << "Trigger"
            << ": " << std::right << std::setw(10) << static_cast<long long>(get("nHLTrigger"))
            << " | Abs eff: " << std::setw(6) << std::fixed << std::setprecision(2)
            << safe_eff(get("nHLTrigger"), get("noCut")) << "%"
            << " | Seq eff: N/A" << std::endl;

        log << std::left << std::setw(20) << "Elec_Trigger"
            << ": " << std::right << std::setw(10) << static_cast<long long>(get("Elec_Trigger"))
            << " | Abs eff: " << std::setw(6) << safe_eff(get("Elec_Trigger"), get("noCut")) << "%"
            << " | Seq eff: N/A" << std::endl;

        log << std::left << std::setw(20) << "Muon_Trigger"
            << ": " << std::right << std::setw(10) << static_cast<long long>(get("Muon_Trigger"))
            << " | Abs eff: " << std::setw(6) << safe_eff(get("Muon_Trigger"), get("noCut")) << "%"
            << " | Seq eff: N/A" << std::endl;

        log << std::left << std::setw(20) << "Filters"
            << ": " << std::right << std::setw(10) << static_cast<long long>(get("nFilter"))
            << " | Abs eff: " << std::setw(6) << safe_eff(get("nFilter"), get("noCut")) << "%"
            << " | Seq eff: " << std::setw(6) << safe_eff(get("nFilter"), get("nHLTrigger")) << "%" << std::endl;

        log << std::left << std::setw(20) << "Good PV"
            << ": " << std::right << std::setw(10) << static_cast<long long>(get("nPV"))
            << " | Abs eff: " << std::setw(6) << safe_eff(get("nPV"), get("noCut")) << "%"
            << " | Seq eff: " << std::setw(6) << safe_eff(get("nPV"), get("nFilter")) << "%" << std::endl;

        log << std::left << std::setw(20) << "nJets > 5"
            << ": " << std::right << std::setw(10) << static_cast<long long>(get("njets>5"))
            << " | Abs eff: " << std::setw(6) << safe_eff(get("njets>5"), get("noCut")) << "%"
            << " | Seq eff: " << std::setw(6) << safe_eff(get("njets>5"), get("nPV")) << "%" << std::endl;

        log << std::left << std::setw(20) << "nbJets > 4"
            << ": " << std::right << std::setw(10) << static_cast<long long>(get("nbjets>4"))
            << " | Abs eff: " << std::setw(6) << safe_eff(get("nbjets>4"), get("noCut")) << "%"
            << " | Seq eff: " << std::setw(6) << safe_eff(get("nbjets>4"), get("njets>5")) << "%" << std::endl;

        log << std::left << std::setw(20) << "nLeptons==1"
            << ": " << std::right << std::setw(10) << static_cast<long long>(get("nlepton==1"))
            << " | Abs eff: " << std::setw(6) << safe_eff(get("nlepton==1"), get("noCut")) << "%"
            << " | Seq eff: " << std::setw(6) << safe_eff(get("nlepton==1"), get("nbjets>4")) << "%" << std::endl;

        log << std::left << std::setw(20) << "MET > 20"
            << ": " << std::right << std::setw(10) << static_cast<long long>(get("MET>20"))
            << " | Abs eff: " << std::setw(6) << safe_eff(get("MET>20"), get("noCut")) << "%"
            << " | Seq eff: " << std::setw(6) << safe_eff(get("MET>20"), get("nlepton==1")) << "%" << std::endl;

        log << std::endl;
        log.flush();
    };

    std::cout << std::endl;
    std::cout << "[WARNING] This analyzer commented out [ \"WTF\" log ] in the header, Please check if you want!!!" << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
    std::cout << "Before start, Let's check the analysis information" << std::endl;
    std::cout << "Run Year    ----> [  " << _year << "  ]" << std::endl;
    std::cout << "Data or MC  ----> [  " << _DataOrMC << "  ]" << std::endl;
    std::cout << "Sample Name ----> [  " << _sampleName << "  ]" << std::endl;

    std::string checklist = "[ tnm.cc ] & [ analyzer header ] & [ main ] & [ analyzer constructor ]";
    bool exitFlag = false;

    if (_year == "nothing") {
        std::cout << "[ERROR] year is not defined, Please check the " << checklist << std::endl;
        exitFlag = true;
    }
    if (_DataOrMC == "nothing") {
        std::cout << "[ERROR] Whether Data or MC is not defined, Please check the " << checklist << std::endl;
        exitFlag = true;
    }
    if (_sampleName == "nothing") {
        std::cout << "[ERROR] SampleName is not defined, Please check the " << checklist << std::endl;
        exitFlag = true;
    }

    std::cout << "--------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    if (exitFlag) std::exit(EXIT_FAILURE);

    std::string analysisInfo = _year + ", " + _DataOrMC + ", " + _sampleName;

    for (int entry = 0; entry < nevents; entry++) {
        _entryInLoop = entry;

        event *currentEvent = new event;
        _ev->read(entry);

        process(currentEvent, sysType, up);

    // Log periódico do evento, independente de objetos
    if (sf_log_file && sf_log_file->is_open() && (_entryInLoop % LOG_INTERVAL == 0)) {
        (*sf_log_file) << "=== Entry " << _entryInLoop << " ===\n";
        (*sf_log_file) << "Number of selected electrons: " << currentEvent->getSelElectrons()->size() << "\n";
        (*sf_log_file) << "Number of selected muons: " << currentEvent->getSelMuons()->size() << "\n";
    }

		
        if (entry % LOG_INTERVAL == 0) {
            std::cout << "[INFO] Processed events of " << analysisInfo << ": " << entry << std::endl;
            currentEvent->summarize();
        }
auto logSFSummary = [&](std::ofstream* logFile) {
    if (!logFile || !logFile->is_open()) return;

    *logFile << "\n========= SF SUMMARY =========\n";
    *logFile << "Total Electrons Processed: " << this->total_electrons_processed << "\n";
    *logFile << "Total Muons Processed: " << this->total_muons_processed << "\n";

    // Lambda para imprimir médias "on-the-fly"
    auto logAverage = [&](TH1F* sumH, TH1F* countH, const char* type, const char* var) {
        if (!sumH || !countH) return;
        *logFile << "\n[" << type << " Avg per " << var << " Bin]\n";
        for (int i = 1; i <= sumH->GetNbinsX(); ++i) {
            double sum   = sumH->GetBinContent(i);
            double sum2  = std::pow(sumH->GetBinError(i), 2); // erro acumulado
            double count = countH->GetBinContent(i);
            if (count > 0) {
                double avg = sum / count;              // média correta
                double err = std::sqrt(sum2) / count; // aproximação do erro da média
                *logFile << "Bin " << i << " (~" << sumH->GetBinCenter(i) << "): "
                         << avg << " ± " << err << " [entries: " << count << "]\n";
            } else {
                *logFile << "Bin " << i << " (~" << sumH->GetBinCenter(i) << "): "
                         << "0 [entries: 0]\n";
            }
        }
    };

    // Agora chama para todos os histogramas de SF/eff
    logAverage(h_sf_vs_pt_sum,       h_sf_vs_pt_count,       "Electron SF", "pt");
    logAverage(h_sf_vs_eta_sum,      h_sf_vs_eta_count,      "Electron SF", "eta");
    logAverage(h_effMC_vs_pt_sum,    h_effMC_vs_pt_count,    "Electron EffMC", "pT");
    logAverage(h_effMC_vs_eta_sum,   h_effMC_vs_eta_count,   "Electron EffMC", "η");
    logAverage(h_sf_muon_vs_pt_sum,  h_sf_muon_vs_pt_count,  "Muon SF", "pt");
    logAverage(h_sf_muon_vs_eta_sum, h_sf_muon_vs_eta_count, "Muon SF", "eta");
};

// Chamar periodicamente dentro do loop
if (_entryInLoop % LOG_INTERVAL == 0) {
    logSFSummary(this->sf_summary_log_file.get());
}



		

        events.push_back(currentEvent);
    }

    // garante um resumo final SEMPRE
    printCutflowSummary("final");

    writeHistos();
    writeTree();

    for (const auto &x : cutflow) {
        std::cout << x.first << ": " << x.second << std::endl;
    }

    hCutFlow->Write();
    hCutFlow_w->Write();
}/// fim do método loop

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ttHHanalyzer::createObjects(event * thisEvent, sysName sysType, bool up){

    ////Log of selection///
//    _entryInLoop++;
    bool doLog = (_entryInLoop % LOG_INTERVAL == 0);
    if (doLog) (*event_log_file) << "==== Evento " << _entryInLoop << " ====" << std::endl;

    cutflow["noCut"]+=1;
    hCutFlow->Fill("noCut",1);
    hCutFlow_w->Fill("noCut",_weight);

    _ev->fillObjects();

    // Trigger paths for SL channel
    if(_year == "2022" || _year == "2022EE" || _year == "2023" || _year == "2023B" || _year == "2024"){
        thisEvent->setFilter(_ev->Flag_goodVertices && _ev->Flag_globalSuperTightHalo2016Filter &&
                             _ev->Flag_EcalDeadCellTriggerPrimitiveFilter && _ev->Flag_BadPFMuonFilter &&
                             _ev->Flag_BadPFMuonDzFilter && _ev->Flag_hfNoisyHitsFilter &&
                             _ev->Flag_eeBadScFilter && _ev->Flag_ecalBadCalibFilter);
        thisEvent->setTrigger(_ev->HLT_Ele30_WPTight_Gsf || _ev->HLT_IsoMu24 || _ev->HLT_Mu50);
    }

    // Primary vertex
    thisEvent->setPV(_ev->PV_npvsGood > 0);

    // References to objects
    auto &genPart = _ev->GenPart;
    auto &jet = _ev->Jet;
    auto &muonT = _ev->Muon;
    auto &ele = _ev->Electron;
    auto &boostedJet = _ev->FatJet;

    objectGenPart * currentGenPart;
    objectBoostedJet * currentBoostedJet;
    objectJet * currentJet;
    objectLep * currentMuon;
    objectLep * currentEle;

    int nVetoMuons = 0, nVetoEle = 0;

    objectMET * MET = new objectMET(_ev->PuppiMET_pt, 0, _ev->PuppiMET_phi, 0);
    thisEvent->setMET(MET);

    ///Test with logs ///
    // === Boosted Jets ===
    int nBoostedJets = 0;
    int nHadronicHiggs = 0;
    for (int i = 0; i < boostedJet.size(); i++) {
        currentBoostedJet = new objectBoostedJet(boostedJet[i].pt, boostedJet[i].eta, boostedJet[i].phi, boostedJet[i].mass);
        currentBoostedJet->softDropMass = boostedJet[i].msoftdrop;

        if (currentBoostedJet->getp4()->Pt() > cut["boostedJetPt"] &&
            fabs(currentBoostedJet->getp4()->Eta()) < fabs(cut["boostedJetEta"])) {

            thisEvent->selectBoostedJet(currentBoostedJet);
            nBoostedJets++;

            if (currentBoostedJet->getp4()->Pt() > cut["hadHiggsPt"] && boostedJet[i].particleNet_HbbvsQCD > cut["bTagDisc"]) {
                thisEvent->selectHadronicHiggs(currentBoostedJet);
                nHadronicHiggs++;
            }
        }
    }
    if(doLog) (*event_log_file) << "Boosted jets selecionados: " << nBoostedJets << ", Hadronic Higgs: " << nHadronicHiggs << std::endl;
	
/////////////////////////////////////////////////////////////////////////////
// Bloco 1: CALIBRAÇÃO DE ELÉTRONS E ATUALIZAÇÃO DO MET
/////////////////////////////////////////////////////////////////////////////

std::vector<float> Electron_pt_before;
std::vector<float> Electron_pt_after;
std::vector<size_t> valid_indices;
std::vector<float> pts, etas_for_calib, r9s;
std::vector<int> gains;

// Limpa os vetores membro no início de cada evento
_final_electron_pts_before_calib.clear();
_final_electron_pts_after_calib.clear();

if (!ele.empty()) {

    // Recuperar limites de segurança do calibrador
    float pt_min   = calibrator.getMin("pt");
    float pt_max   = calibrator.getMax("pt");
    float eta_min  = calibrator.getMin("ScEta");
    float eta_max  = calibrator.getMax("ScEta");

    // Garante que o run number é 1 para amostras MC
    int run_for_calibrator = _runNumber;
    if (calibrator.isMC()) {
        run_for_calibrator = 1;
    }

    // Prepara os vetores com as variáveis dos elétrons para a calibração
    for (size_t i = 0; i < ele.size(); ++i) {
        float scEta = ele[i].deltaEtaSC + ele[i].eta;
        float pt_clamped   = std::min(std::max(ele[i].pt, pt_min), pt_max);
        float scEta_clamped = std::min(std::max(scEta, eta_min), eta_max);

        if (std::isfinite(pt_clamped) && std::isfinite(scEta_clamped)) {
            pts.push_back(pt_clamped);
            etas_for_calib.push_back(scEta_clamped);
            r9s.push_back(ele[i].r9);
            gains.push_back(ele[i].seedGain);
            Electron_pt_before.push_back(ele[i].pt);
            valid_indices.push_back(i);
        }
    }

    // Se existirem elétrons válidos, aplica a calibração
    if (!pts.empty()) {
        calibrator.calibrateElectrons(pts, etas_for_calib, r9s, gains, run_for_calibrator);
        Electron_pt_after = pts;

        // Atualiza a coleção original 'ele' com os pTs calibrados
        for (size_t j = 0; j < valid_indices.size(); ++j) {
            size_t original_idx = valid_indices[j];
            ele[original_idx].pt = pts[j];
        }

        // Recalcula o MET com base na diferença do pt dos elétrons 
		//futuramenbte adicionar a contribuição do PT CALIBRADO dos MUONS 
		//e também adicionar a conbtribuição da calibração dos muons e eletrons para :
	    //hleptonHT (Soma escalar do pT dos léptons): Se houver elétrons, esta soma mudará.
		//hST (Soma escalar do pT de todos os objetos + MET): Esta variável é afetada por duas vias
		
        if (MET) {
            float met_px = MET->getp4()->Px();
            float met_py = MET->getp4()->Py();
            for (size_t j = 0; j < valid_indices.size(); ++j) {
                size_t idx = valid_indices[j];
                float old_px = Electron_pt_before[j] * cos(ele[idx].phi);
                float old_py = Electron_pt_before[j] * sin(ele[idx].phi);
                float new_px = Electron_pt_after[j]  * cos(ele[idx].phi);
                float new_py = Electron_pt_after[j]  * sin(ele[idx].phi);
                met_px += (old_px - new_px);
                met_py += (old_py - new_py);
            }
            MET->getp4()->SetPxPyPzE(met_px, met_py, MET->getp4()->Pz(), sqrt(met_px*met_px + met_py*met_py));
        }

        // ===================================================================
        //  SALVA OS VETORES DE PT NAS VARIÁVEIS MEMBRO DA CLASSE
        // ===================================================================
        _final_electron_pts_before_calib = Electron_pt_before;
        _final_electron_pts_after_calib  = Electron_pt_after;
    }
}
////////////////////////////////////////////////////////////////
// Bloco 2: SELEÇÃO DE LÉPTONS E VERIFICAÇÃO FINAL
////////////////////////////////////////////////////////////////

bool thereIsALeadLepton = false;
int nLeadingMuons = 0;
int nLeadingElectrons = 0;

for (size_t i = 0; i < muonT.size(); i++) {
    if (cut.count("muonEta") && cut.count("muonIso") && cut.count("leadMuonPt")) {
        if (fabs(muonT[i].eta) < cut["muonEta"] && muonT[i].tightId && muonT[i].pfRelIso04_all < cut["muonIso"] && muonT[i].pt > cut["leadMuonPt"]) {
            thereIsALeadLepton = true; nLeadingMuons++; break;
        }
    }
}
if (!thereIsALeadLepton) {
    for (size_t i = 0; i < ele.size(); i++) {
        if (cut.count("eleEta") && cut.count("leadElePt")) {
            if ((fabs(ele[i].deltaEtaSC + ele[i].eta) < 1.4442 || fabs(ele[i].deltaEtaSC + ele[i].eta) > 1.5660) && fabs(ele[i].eta) < cut["eleEta"] && ele[i].mvaIso_WP90 && ele[i].pt > cut["leadElePt"]) {
                thereIsALeadLepton = true; nLeadingElectrons++; break;
            }
        }
    }
}

if (doLog) (*event_log_file) << "Lepton líder encontrado? " << (thereIsALeadLepton ? "SIM" : "NÃO") << " | Muons líderes: " << nLeadingMuons << " | Elétrons líderes: " << nLeadingElectrons << std::endl;

int nSubMuons = 0;
int nSubEles = 0;
if (thereIsALeadLepton) {
    for (size_t i = 0; i < muonT.size(); i++) {
        if (cut.count("muonEta") && cut.count("muonIso") && cut.count("subLeadMuonPt")) {
            if (fabs(muonT[i].eta) < cut["muonEta"] && muonT[i].tightId && muonT[i].pfRelIso04_all < cut["muonIso"] && muonT[i].pt > cut["subLeadMuonPt"]) {
                std::cout << ">> DIAGNOSTICO: Criando Muon objeto com pt=" << muonT[i].pt << std::endl;
                currentMuon = new objectLep(muonT[i].pt, muonT[i].eta, muonT[i].phi, 0.);
                std::cout << ">> DIAGNOSTICO: Objeto Muon RECEM-CRIADO tem pt=" << currentMuon->getp4()->Pt() << std::endl;
                thisEvent->selectMuon(currentMuon);
                nSubMuons++;
            }
        }
    }

    for (size_t i = 0; i < ele.size(); i++) {
        if (cut.count("eleEta") && cut.count("subLeadElePt")) {
            if ((fabs(ele[i].deltaEtaSC + ele[i].eta) < 1.4442 || fabs(ele[i].deltaEtaSC + ele[i].eta) > 1.5660) && fabs(ele[i].eta) < cut["eleEta"] && ele[i].mvaIso_WP90 && ele[i].pt > cut["subLeadElePt"]) {
                std::cout << ">> DIAGNOSTICO: Criando Electron objeto com pt=" << ele[i].pt << " (calibrado)" << std::endl;
                currentEle = new objectLep(ele[i].pt, ele[i].eta, ele[i].phi, 0.);
                // ===================================================================
                //  NOVO LOG DE DIAGNÓSTICO FINAL
                // ===================================================================
                std::cout << ">> DIAGNOSTICO: Objeto Electron RECEM-CRIADO tem pt=" << currentEle->getp4()->Pt() << std::endl;
                thisEvent->selectEle(currentEle);
                nSubEles++;
            }
        }
    }
}

if (doLog) (*event_log_file) << "Sub-leading leptons: Muons: " << nSubMuons << ", Electrons: " << nSubEles << std::endl;

std::cout << "--- VERIFICACAO DO LEPTON LIDER NO FINAL DO EVENTO ---" << std::endl;
objectLep* finalLepton = nullptr;
std::string lepType = "Nenhum";
std::vector<objectLep*>* selMuons = thisEvent->getSelMuons();
std::vector<objectLep*>* selEles = thisEvent->getSelElectrons();

if (nLeadingMuons > 0 && selMuons && !selMuons->empty()) {
    finalLepton = selMuons->at(0); lepType = "Muon";
} else if (nLeadingElectrons > 0 && selEles && !selEles->empty()) {
    finalLepton = selEles->at(0); lepType = "Electron";
}

if (finalLepton) {
    std::cout << "Lepton Lider Selecionado (" << lepType << "): pt = " << finalLepton->getp4()->Pt() << ", eta = " << finalLepton->getp4()->Eta() << std::endl;
} else {
    std::cout << "Nenhum lepton valido foi selecionado e armazenado no thisEvent." << std::endl;
}


	
	
    // === Jets ===
    int nJets = 0, nLightJets = 0, nLooseBJets = 0, nMediumBJets = 0;
    for (int i = 0; i < jet.size(); i++) {
        currentJet = new objectJet(jet[i].pt, jet[i].eta, jet[i].phi, jet[i].mass);
        currentJet->bTagCSV = jet[i].btagUParTAK4B;

        float pt = jet[i].pt;
        float eta = jet[i].eta;
        float nhf = jet[i].neHEF;
        float nemf = jet[i].neEmEF;
        float chf = jet[i].chHEF;
        float chemf = jet[i].chEmEF;
        float muf = jet[i].muEF;
        int chm = jet[i].chMultiplicity;
        int numConst = jet[i].nConstituents;
        int numNeutral = jet[i].neMultiplicity;
        float absEta = fabs(eta);
        bool passJetID = false;

        if (absEta <= 2.6) passJetID = (chemf < 0.8 && chm > 0 && chf > 0.01 && numConst > 1 && nemf < 0.9 && muf < 0.8 && nhf < 0.99);
        else if (absEta > 2.6 && absEta <= 2.7) passJetID = (chemf < 0.8 && nemf < 0.99 && muf < 0.8 && nhf < 0.9);
        else if (absEta > 2.7 && absEta <= 3.0) passJetID = (nhf < 0.99);
        else if (absEta > 3.0) passJetID = (nemf < 0.4 && numNeutral >= 2);

        bool passPUJetID = true;

        if(currentJet->getp4()->Pt() > cut["jetPt"] && fabs(currentJet->getp4()->Eta()) < abs(cut["jetEta"]) && passJetID && passPUJetID){
            thisEvent->selectJet(currentJet);
            nJets++;
            if(jet[i].btagUParTAK4B <= currentJet->getValbTagLoose(_year)) { thisEvent->selectLightJet(currentJet); nLightJets++; }
            else if(jet[i].btagUParTAK4B > currentJet->getValbTagMedium(_year)) { thisEvent->selectbJet(currentJet); nMediumBJets++; }
            if(jet[i].btagUParTAK4B > currentJet->getValbTagLoose(_year)) { thisEvent->selectLoosebJet(currentJet); nLooseBJets++; }
        }
    }
    if(doLog) (*event_log_file) << "Jets selecionados: " << nJets << " | Light: " << nLightJets << " | Loose b-jets: " << nLooseBJets << " | Medium b-jets: " << nMediumBJets << std::endl;

    // === GenPart (b-quarks com Higgs ou Top mãe) ===
    int nGenBJets = 0, nFromHiggs = 0, nFromTop = 0;
    for (int i = 0; i < genPart.size(); i++){
        if ((abs(genPart[i].pdgId) == 5) && (genPart[i].statusFlags & 256)){
            currentGenPart = new objectGenPart(genPart[i].pt, genPart[i].eta, genPart[i].phi, genPart[i].mass);
            currentGenPart->hasHiggsMother = false;
            currentGenPart->hasTopMother = false;
            int motherInd = genPart[i].genPartIdxMother;
            while (motherInd > 1 && motherInd < genPart.size() && motherInd >= 0){
                int pdg = abs(genPart[motherInd].pdgId);
                if ((pdg == 25 || pdg == 6) && (genPart[motherInd].statusFlags & 256)){
                    if(pdg==25) currentGenPart->hasHiggsMother = true;
                    if(pdg==6) currentGenPart->hasTopMother = true;
                    break;
                }
                motherInd = genPart[motherInd].genPartIdxMother;
            }
            thisEvent->selectGenPart(currentGenPart);
            nGenBJets++;
            if(currentGenPart->hasHiggsMother) nFromHiggs++;
            if(currentGenPart->hasTopMother) nFromTop++;
        }
    }
    if(doLog) {
        (*event_log_file) << "Gen b-quarks selecionados: " << nGenBJets
                          << " | Com mãe Higgs: " << nFromHiggs
                          << " | Com mãe Top: " << nFromTop << std::endl << std::endl;
    }
}//fim do método createObjects
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ttHHanalyzer::selectObjects(event *thisEvent) {
    // Trigger cut
bool doLog = (_entryInLoop % LOG_INTERVAL == 0);

if (cut["trigger"] > 0) {
    if (!thisEvent->getTriggerAccept()) {
        if (doLog) (*event_log_file) << "Evento rejeitado pelo trigger.\n";
        return false;
    } else {
        cutflow["nHLTrigger"] += 1;
        hCutFlow->Fill("nHLTrigger", 1);
        hCutFlow_w->Fill("nHLTrigger", _weight);

        if (_ev->HLT_IsoMu24) {
            if (doLog) (*event_log_file) << "Trigger múon OK.\n";
            cutflow["Muon_Trigger"] += 1;
            hCutFlow->Fill("Muon_Trigger", 1);
            hCutFlow_w->Fill("Muon_Trigger", _weight);
        } 
        else if (_ev->HLT_Ele30_WPTight_Gsf) {
            if (doLog) (*event_log_file) << "Trigger elétron OK.\n";
            cutflow["Elec_Trigger"] += 1;
            hCutFlow->Fill("Elec_Trigger", 1);
            hCutFlow_w->Fill("Elec_Trigger", _weight);
        }
    }
}

// Filter cut
if (cut["filter"] > 0) {
    if (!thisEvent->getMETFilter()) {
        if (doLog) (*event_log_file) << "Evento rejeitado pelo Filters.\n";
        return false;
    } else {
        if (doLog) (*event_log_file) << "Filters OK.\n";
        cutflow["nFilter"] += 1;
        hCutFlow->Fill("nFilter", 1);
        hCutFlow_w->Fill("nFilter", _weight);
    }
}

// Primary vertex cut
if (doLog) (*event_log_file) << "PV_npvsGood=" << _ev->PV_npvsGood << "\n";
if (cut["pv"] > 0) {
    if (!thisEvent->getPVvalue()) {
        if (doLog) (*event_log_file) << "Evento rejeitado pelo PV.\n";
        return false;
    } else {
        if (doLog) (*event_log_file) << "PV OK.\n";
        cutflow["nPV"] += 1;
        hCutFlow->Fill("nPV", 1);
        hCutFlow_w->Fill("nPV", _weight);
    }
}

// Jet multiplicity
if (thisEvent->getnSelJet() < cut["nJets"]) {
    if (doLog) (*event_log_file) << "Evento rejeitado: nJets = " << thisEvent->getnSelJet() 
                                 << " (mínimo requerido = " << cut["nJets"] << ")\n";
    return false;
}
cutflow["njets>5"] += 1;
hCutFlow->Fill("njets>5", 1);
hCutFlow_w->Fill("njets>5", _weight);
if (doLog) (*event_log_file) << "Jets selecionados = " << thisEvent->getnSelJet() << "\n";

// b-jet multiplicity
if (thisEvent->getnbJet() < cut["nbJets"]) {
    if (doLog) (*event_log_file) << "Evento rejeitado: nbJets = " << thisEvent->getnbJet() 
                                 << " (mínimo requerido = " << cut["nbJets"] << ")\n";
    return false;
}
cutflow["nbjets>4"] += 1;
hCutFlow->Fill("nbjets>4", 1);
hCutFlow_w->Fill("nbjets>4", _weight);
if (doLog) (*event_log_file) << "b-Jets selecionados = " << thisEvent->getnbJet() << "\n";

// Lepton multiplicity
if (thisEvent->getnSelLepton() != cut["nLeptons"]) {
    if (doLog) (*event_log_file) << "Evento rejeitado: nLeptons = " << thisEvent->getnSelLepton()
                                 << " (esperado = " << cut["nLeptons"] << ")\n";
    return false;
}
cutflow["nlepton==1"] += 1;
hCutFlow->Fill("nlepton==1", 1);
hCutFlow_w->Fill("nlepton==1", _weight);
if (doLog) (*event_log_file) << "Leptons selecionados = " << thisEvent->getnSelLepton() << "\n";

thisEvent->getStatsComb(thisEvent->getSelJets(), thisEvent->getSelLeptons(), ljetStat);
thisEvent->getStatsComb(thisEvent->getSelbJets(), thisEvent->getSelLeptons(), lbjetStat);

// Leading lepton cuts
if (thisEvent->getSelElectrons()->size() == 1 && thisEvent->getSelMuons()->size() == 0) {
    auto ele = thisEvent->getSelElectrons()->at(0);
    if (!(ele->getp4()->Pt() >= cut["leadElePt"] && fabs(ele->getp4()->Eta()) <= cut["eleEta"])) {
        if (doLog) (*event_log_file) << "Evento rejeitado: elétron líder pT=" << ele->getp4()->Pt() 
                                     << " eta=" << ele->getp4()->Eta() << "\n";
        return false;
    }
    if (doLog) (*event_log_file) << "Elétron líder OK: pT=" << ele->getp4()->Pt()
                                 << " eta=" << ele->getp4()->Eta() << "\n";
} 
else if (thisEvent->getSelElectrons()->size() == 0 && thisEvent->getSelMuons()->size() == 1) {
    auto mu = thisEvent->getSelMuons()->at(0);
    if (!(mu->getp4()->Pt() >= cut["leadMuonPt"] && fabs(mu->getp4()->Eta()) <= cut["muonEta"])) {
        if (doLog) (*event_log_file) << "Evento rejeitado: múon líder pT=" << mu->getp4()->Pt() 
                                     << " eta=" << mu->getp4()->Eta() << "\n";
        return false;
    }
    if (doLog) (*event_log_file) << "Múon líder OK: pT=" << mu->getp4()->Pt()
                                 << " eta=" << mu->getp4()->Eta() << "\n";
}

// MET cut
if (!(thisEvent->getMET()->getp4()->Pt() > cut["MET"])) {
    if (doLog) (*event_log_file) << "Evento rejeitado: MET=" << thisEvent->getMET()->getp4()->Pt()
                                 << " (mínimo requerido = " << cut["MET"] << ")\n";
    return false;
}
cutflow["MET>20"] += 1;
hCutFlow->Fill("MET>20", 1);
hCutFlow_w->Fill("MET>20", _weight);
if (doLog) (*event_log_file) << "MET OK: " << thisEvent->getMET()->getp4()->Pt() << "\n";

return true;
}
/// end of selection
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void ttHHanalyzer::motherReco(const TLorentzVector & dPar1p4,const TLorentzVector & dPar2p4, const float mother1mass, float & _minChi2,float & _bbMassMin1){
    float bbMass1, chi2;
    bbMass1 = (dPar1p4+dPar2p4).M();
    chi2 = pow((bbMass1 - mother1mass),2)/pow((dPar1p4.Pt()+dPar2p4.Pt())/2.*0.03,0.5);
    if(_minChi2 > chi2){
	_minChi2      = chi2;
	_bbMassMin1   = bbMass1;
	_bpTHiggs1    = (dPar1p4+dPar2p4).Pt();
    }
} 


void ttHHanalyzer::diMotherReco(const TLorentzVector & dPar1p4,const TLorentzVector & dPar2p4,const TLorentzVector & dPar3p4,const TLorentzVector & dPar4p4, const float mother1mass, const float  mother2mass, float & _minChi2,float & _bbMassMin1, float & _bbMassMin2){
    float bbMass1, bbMass2, chi2;
    bbMass1 = (dPar1p4+dPar2p4).M();
    bbMass2 = (dPar3p4+dPar4p4).M();
    chi2 = pow((bbMass1 - mother1mass),2)/pow((dPar1p4.Pt()+dPar2p4.Pt())/2.*0.03,0.5) + pow((bbMass2 - mother2mass),2)/pow((dPar3p4.Pt()+dPar4p4.Pt())/2.*0.03,0.5);
    if(_minChi2 > chi2){
	_minChi2      = chi2;
	_bbMassMin1   = bbMass1;
	_bbMassMin2   = bbMass2;
	_bpTHiggs1    = (dPar1p4+dPar2p4).Pt();
	_bpTHiggs2    = (dPar3p4+dPar4p4).Pt();
    }
} 
void ttHHanalyzer::analyze(event *thisEvent) {
    // ========================
    // Pega vetores de leptons do evento
    // ========================
    std::vector<objectLep*>* selElePtr = thisEvent->getSelElectrons();
    std::vector<objectLep*>* selMuPtr  = thisEvent->getSelMuons();

    std::vector<objectLep*> selectedElectrons;
    std::vector<objectLep*> selectedMuons;

    if (selElePtr) selectedElectrons = *selElePtr;
    if (selMuPtr)  selectedMuons     = *selMuPtr;

    // As variáveis de pT "before" e "after" agora são lidas das variáveis membro
    // que foram salvas no passo anterior. Não precisamos recriá-las aqui.

    bool isMC = (_DataOrMC == "MC");
    float weight_before_SFs = _weight;
    float triggerSF = 1.0, recoSF = 1.0, idSF = 1.0, totalSFUnc = 0.0;
    float trigSF_unc = 0.0, recoSF_unc = 0.0, idSF_unc = 0.0;

    if (!selectedElectrons.empty()) {
        objectLep* ele = selectedElectrons[0];
        triggerSF = getEleTrigSF(ele->getp4()->Eta(), ele->getp4()->Pt(), trigSF_unc);
        recoSF    = getEleRecoSF(ele->getp4()->Eta(), ele->getp4()->Pt(), recoSF_unc);
        idSF      = getEleIDSF(ele->getp4()->Eta(), ele->getp4()->Phi(), ele->getp4()->Pt(), idSF_unc);
        totalSFUnc = std::sqrt(trigSF_unc*trigSF_unc + recoSF_unc*recoSF_unc + idSF_unc*idSF_unc);
        _weight *= (triggerSF * recoSF * idSF);
    } else if (!selectedMuons.empty()) {
        objectLep* mu = selectedMuons[0];
        
        // Pega os SFs de Trigger e de ID para o múon
        triggerSF = getMuonTrigSF(mu->getp4()->Eta(), mu->getp4()->Pt());
        idSF      = getMuonIDSF(mu->getp4()->Eta(), mu->getp4()->Pt()); // <-- NOVA LINHA

        totalSFUnc = 0.0; // Incertezas para múons ainda não implementadas
        
        // Aplica ambos os SFs ao peso do evento
        _weight *= (triggerSF * idSF); // <-- LINHA MODIFICADA
    }

    // ========================
    // Log detalhado (versão corrigida)
    // ========================
    bool hasLeadLepton = (!selectedElectrons.empty() || !selectedMuons.empty());
    if (sf_log_file && sf_log_file->is_open() && (_entryInLoop % LOG_INTERVAL == 0 && hasLeadLepton)) {
        (*sf_log_file) << "=== Entry " << _entryInLoop << " ===\n";
        size_t nEle = selectedElectrons.size();
        size_t nMu  = selectedMuons.size();

        (*sf_log_file) << "Number of selected electrons: " << nEle << "\n";
        (*sf_log_file) << "Number of selected muons: "     << nMu << "\n";

        if (nEle > 0) {
            // Garante que não vamos ler além dos limites dos vetores
            size_t nLog = std::min({nEle, _final_electron_pts_before_calib.size(), _final_electron_pts_after_calib.size()});
            for (size_t i = 0; i < nLog; ++i) {
                objectLep* ele = selectedElectrons[i];
                if (!ele) continue;

                // USA AS VARIÁVEIS MEMBRO CORRETAS PARA O LOG
                float pt_before = _final_electron_pts_before_calib[i];
                float pt_after  = _final_electron_pts_after_calib[i];

                (*sf_log_file) << "Electron " << i
                               << " | η = " << ele->getp4()->Eta()
                               << ", φ = " << ele->getp4()->Phi()
                               << ", pT before calib = " << pt_before
                               << ", pT after calib  = " << pt_after
                               << " | Calibration applied = " << (isMC ? "Smearing (MC)" : "Scale Correction (Data)")
                               << " | Trigger SF = " << triggerSF
                               << " | Reco SF = " << recoSF
                               << " | ID SF = " << idSF
                               << " | Weight before = " << weight_before_SFs
                               << " | Weight after = " << _weight
                               << " | Total SF Uncertainty = " << totalSFUnc
                               << "\n";
            }
        } else if (nMu > 0) {
            objectLep* mu = selectedMuons[0];
            if (mu) {
                (*sf_log_file) << "Muon | η = " << mu->getp4()->Eta()
                               << ", pT = " << mu->getp4()->Pt()
                               << " | Trigger SF = " << triggerSF      // <-- LOG MODIFICADO
                               << " | ID SF = " << idSF                 // <-- NOVA LINHA NO LOG
                               << " | Weight before = " << weight_before_SFs
                               << " | Weight after = " << _weight
                               << "\n";
            }
        }
    }

    _entryInLoop++;

///////////////////////////////////////

	
    std::vector<objectJet*>* bJetsInv = thisEvent->getSelbJets(); 
    std::vector<objectJet*>* lbJetsInv = thisEvent->getLoosebJets(); 
    std::vector<objectJet*>* jetsInv = thisEvent->getSelJets(); 

    std::vector<TVector3> vectorsJet, vectorsBjet;
    // Event Shape Calculation & genbjet matching for mother particle
    for(int k = 0; k < jetsInv->size(); k++){
	vectorsJet.push_back(jetsInv->at(k)->getp4()->Vect());
    }
    for(int m = 0; m < bJetsInv->size(); m++){
	vectorsBjet.push_back(bJetsInv->at(m)->getp4()->Vect());
	if(thisEvent->getnGenPart() < 1) continue;
	bJetsInv->at(m)->matchedtoHiggs = false;	    	
	for(auto genParticle: *thisEvent->getGenParts()){
	    float dR = bJetsInv->at(m)->getp4()->DeltaR( *genParticle->getp4());
	    if(genParticle->hasHiggsMother == true && dR < 0.8){
	      if(dR > genParticle->dRmatched) {
		  //std::cout << "this was matched to a closer particle before" << std::endl;
	      }else if( genParticle->matched){
		  //std::cout << "this was matched before" << std::endl;
	      }
	      bJetsInv->at(m)->matchedtoHiggs   = true;	
	      bJetsInv->at(m)->matchedtoHiggsdR = dR;	
	      genParticle->matched = true;
	      genParticle->dRmatched = dR;
	      break;
	    }
	} 
	//	std::cout << bJetsInv->at(m)->matchedtoHiggsb << std::endl;  
    }

    _minChi2Higgs  = cLargeValue;
    _minChi2Z      = cLargeValue;
    _minChi2HiggsZ = cLargeValue;
    _minChi2SHiggsNotMatched = cLargeValue;
    _minChi2SHiggsMatched = cLargeValue;
    _minChi2HHNotMatched = cLargeValue;
    _minChi2HHMatched = cLargeValue;
    _bbMassMinSHiggsMatched = -1;
    _bbMassMinSHiggsNotMatched = -1;
    _bbMassMinHH1Matched = -1;
    _bbMassMinHH1NotMatched = -1;
    _bbMassMinHH2Matched = -1;
    _bbMassMinHH2NotMatched = -1;

    float tempminChi2 = cLargeValue, tmpMassMin1HiggsZ = 0., tmpMassMin2HiggsZ = 0.;
    float tempMinChi2SHiggs = cLargeValue, tempMinChi2SHiggs_r = cLargeValue, tmpMassMinSHiggs = 0.;
    float tempMinChi2SHiggsMatched = cLargeValue, tempMinChi2SHiggsMatched_r = cLargeValue, tmpMassMinSHiggsMatched = 0.;
    float tempMinChi2SHiggsNotMatched = cLargeValue, tempMinChi2SHiggsNotMatched_r = cLargeValue, tmpMassMinSHiggsNotMatched = 0.;
    //extract H
    for( int ibjet1 = 0; ibjet1 < bJetsInv->size(); ibjet1++){
	tempMinChi2SHiggs = cLargeValue;
	bJetsInv->at(ibjet1)->minChiHiggsIndex = -1;
	for( int ibjet2 = 1; ibjet2 < bJetsInv->size(); ibjet2++){
	    if( ibjet1 == ibjet2) continue;	   
	    tempMinChi2SHiggs_r = tempMinChi2SHiggs;
	    motherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4(),cHiggsMass, tempMinChi2SHiggs, tmpMassMinSHiggs);
	    if(tempMinChi2SHiggs_r > tempMinChi2SHiggs){
		bJetsInv->at(ibjet1)->minChiHiggsIndex = ibjet2;
	    }
	    if(bJetsInv->at(ibjet1)->matchedtoHiggs == true && bJetsInv->at(ibjet2)->matchedtoHiggs == true){
		motherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4(),cHiggsMass, _minChi2SHiggsMatched, _bbMassMinSHiggsMatched);		
	    } else {
		motherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4(),cHiggsMass, _minChi2SHiggsNotMatched, _bbMassMinSHiggsNotMatched);		
	    }
	}
	bJetsInv->at(ibjet1)->minChiHiggs = tempMinChi2SHiggs;
    }


    // HH & ZZ reco : 4 medium b jet case

    if(thisEvent->getnbJet() >  3){
    	for( int ibjet1 = 0; ibjet1 < bJetsInv->size(); ibjet1++){
   	    for( int ibjet2 = ibjet1+1; ibjet2 < bJetsInv->size(); ibjet2++){
    		if( ibjet1 == ibjet2) continue;
    		for( int ibjet3 = 1; ibjet3 < bJetsInv->size(); ibjet3++){
    		    if(ibjet1 == ibjet3 || ibjet2 == ibjet3) continue;
    		    for( int ibjet4 = ibjet3+1; ibjet4 < bJetsInv->size(); ibjet4++){
    			if(ibjet1 == ibjet4 || ibjet2 == ibjet4 || ibjet3 == ibjet4 ) continue;		
			if(bJetsInv->at(ibjet1)->matchedtoHiggs == true && bJetsInv->at(ibjet2)->matchedtoHiggs == true && bJetsInv->at(ibjet3)->matchedtoHiggs == true && bJetsInv->at(ibjet4)->matchedtoHiggs == true){
			    diMotherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
					 , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4() 
					 , cHiggsMass, cHiggsMass, _minChi2HHMatched, _bbMassMinHH1Matched, _bbMassMinHH2Matched);
			} else {
			    diMotherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
					 , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4() 
					 , cHiggsMass, cHiggsMass, _minChi2HHNotMatched, _bbMassMinHH1NotMatched, _bbMassMinHH2NotMatched);
			}
			diMotherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4() 
				     , cHiggsMass, cHiggsMass, _minChi2Higgs, _bbMassMin1Higgs, _bbMassMin2Higgs);
			diMotherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cZMass, cZMass, _minChi2Z, _bbMassMin1Z, _bbMassMin2Z);  
			diMotherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4() //ZH 
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cHiggsMass, cZMass, _minChi2HiggsZ, _bbMassMin1HiggsZ, _bbMassMin2HiggsZ); 
		    }
		}
	    }
	}
	// HH & ZZ reco : 3 medium + 1 loose b jet case
    } else if(thisEvent->getnbJet() == 3 && thisEvent->getnbLooseJet() > 3){
	for( int ibjet1 = 0; ibjet1 < lbJetsInv->size(); ibjet1++){
	    for( int ibjet2 = 0; ibjet2 < bJetsInv->size(); ibjet2++){
		if( lbJetsInv->at(ibjet1) == bJetsInv->at(ibjet2)) continue;
		for( int ibjet3 = 0; ibjet3 < bJetsInv->size(); ibjet3++){
		    if(lbJetsInv->at(ibjet1) == bJetsInv->at(ibjet3) || ibjet2 == ibjet3) continue;
		    for( int ibjet4 = ibjet3+1; ibjet4 < bJetsInv->size(); ibjet4++){
			if(lbJetsInv->at(ibjet1) == bJetsInv->at(ibjet4) || ibjet2 == ibjet4 || ibjet3 == ibjet4 ) continue;		
			if( bJetsInv->at(ibjet2)->matchedtoHiggs == true && bJetsInv->at(ibjet3)->matchedtoHiggs == true && bJetsInv->at(ibjet4)->matchedtoHiggs == true){
			    diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
					 , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4() 
					 , cHiggsMass, cHiggsMass, _minChi2HHMatched, _bbMassMinHH1Matched, _bbMassMinHH2Matched);
			} else {
			    diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
					 , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4() 
					 , cHiggsMass, cHiggsMass, _minChi2HHNotMatched, _bbMassMinHH1NotMatched, _bbMassMinHH2NotMatched);
			}

			diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cHiggsMass, cHiggsMass, _minChi2Higgs, _bbMassMin1Higgs, _bbMassMin2Higgs); 
			diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cZMass, cZMass, _minChi2Z, _bbMassMin1Z, _bbMassMin2Z); 
			// ZH combinatorics 
			diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4() // ZH H(lbb)Z(bb)
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cHiggsMass, cZMass, _minChi2HiggsZ, _bbMassMin1HiggsZ, _bbMassMin2HiggsZ);  
			tempminChi2 = _minChi2HiggsZ; tmpMassMin1HiggsZ = _bbMassMin1HiggsZ; tmpMassMin2HiggsZ = _bbMassMin2HiggsZ;
			diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4() // ZH Z(lbb)H(bb)
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cZMass, cHiggsMass, _minChi2HiggsZ, _bbMassMin2HiggsZ, _bbMassMin1HiggsZ); 
			// pick the lowest minChi2 for the two combinatorics
			if(tempminChi2 < _minChi2HiggsZ){
			    _minChi2HiggsZ = tempminChi2; _bbMassMin1HiggsZ = tmpMassMin1HiggsZ; _bbMassMin2HiggsZ = tmpMassMin2HiggsZ;
			}
		    }
		}
	    }
	}
    } else if(thisEvent->getnbJet() == 3){   	// HH & ZZ reco : 3 medium b jet + 1 jet case 
	for( int ijet1 = 0; ijet1 < jetsInv->size(); ijet1++){
	    for( int ibjet2 = 0; ibjet2 < bJetsInv->size(); ibjet2++){
		if( jetsInv->at(ijet1) == bJetsInv->at(ibjet2)) continue;
		for( int ibjet3 = 0; ibjet3 < bJetsInv->size(); ibjet3++){
		    if(jetsInv->at(ijet1) == bJetsInv->at(ibjet3) || ibjet2 == ibjet3) continue;
		    for( int ibjet4 = ibjet3+1; ibjet4 < bJetsInv->size(); ibjet4++){
			if(jetsInv->at(ijet1) == bJetsInv->at(ibjet4) || ibjet2 == ibjet4 || ibjet3 == ibjet4 ) continue;		
			diMotherReco(*jetsInv->at(ijet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cHiggsMass, cHiggsMass, _minChi2Higgs, _bbMassMin1Higgs, _bbMassMin2Higgs);
			diMotherReco(*jetsInv->at(ijet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cZMass, cZMass, _minChi2Z, _bbMassMin1Z, _bbMassMin2Z);
			// ZH combinatorics 
			diMotherReco(*jetsInv->at(ijet1)->getp4(), *bJetsInv->at(ibjet2)->getp4() // ZH H(jb)Z(bb)
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cHiggsMass, cZMass, _minChi2HiggsZ, _bbMassMin1HiggsZ, _bbMassMin2HiggsZ);
			tempminChi2 = _minChi2HiggsZ; tmpMassMin1HiggsZ = _bbMassMin1HiggsZ; tmpMassMin2HiggsZ = _bbMassMin2HiggsZ;
			diMotherReco(*jetsInv->at(ijet1)->getp4(), *bJetsInv->at(ibjet2)->getp4() // ZH Z(jb)H(bb)
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cZMass, cHiggsMass, _minChi2HiggsZ, _bbMassMin2HiggsZ, _bbMassMin1HiggsZ);
			// pick the lowest minChi2 for the two combinatorics
			if(tempminChi2 < _minChi2HiggsZ){
			    _minChi2HiggsZ = tempminChi2; _bbMassMin1HiggsZ = tmpMassMin1HiggsZ; _bbMassMin2HiggsZ = tmpMassMin2HiggsZ;
			}
		    }
		}
	    }
	}
    } 


    //    std::map<std::string, float> testVars = HypoComb->GetBestPermutation(getLepP4(thisEvent),getJetP4(thisEvent),getJetCSV(thisEvent),*(thisEvent->getMET()->getp4()));
    //    HypoComb.GetBestPermutation(getLepP4(thisEvent),getJetP4(thisEvent),getJetCSV(thisEvent),*(thisEvent->getMET()->getp4()));
    //    std::cout<< "BLR: " << testVars["Evt_blr"] << std::endl;

    thisEvent->eventShapeJet  = new EventShape(vectorsJet);
    thisEvent->eventShapeBjet  = new EventShape(vectorsBjet);
}

///(TSFel) (IDSFel) (TSFel)
void ttHHanalyzer::process(event* thisEvent, sysName sysType, bool up) {
    // Reset do peso antes de aplicar SFs
    _weight = _initialWeight;

    // Criação de objetos para o evento
    createObjects(thisEvent, sysType, up);
    if (!selectObjects(thisEvent)) return;

    // Inicializar SFs de elétrons (Trigger, Reco, ID)
    static bool eleSFInitialized = false;
    if (!eleSFInitialized) {
        initTriggerSF();  // Trigger SF
        initRecoSF();     // Reco SF
        initEleIDSF();      // ID SF (MVA ISO WP90)
        eleSFInitialized = true;
    }

    // Inicializar SFs de múons (Trigger e ID/ISO)
    static bool muonSFInitialized = false;
    if (!muonSFInitialized) {
        initMuonHLTriggerSF(); // Carrega SFs de Trigger do Múon
        initMuonIDSF();        // <-- NOVA LINHA: Carrega SFs de ID/ISO do Múon
        muonSFInitialized = true;
    }
    // Análise do evento: agora já aplica todos os SFs (Trigger * Reco * ID) e loga os valores
    analyze(thisEvent);

    // Preenchimento de histogramas e árvore
    fillHistos(thisEvent);
    fillTree(thisEvent);
}


void ttHHanalyzer::fillHistos(event * thisEvent){
static long long total_electrons_processed = 0;
static long long total_muons_processed = 0;
      
	int bin = 1;
	for (auto &x : cutflow) {
	    hCutFlow->GetXaxis()->SetBinLabel(bin, x.first.c_str());
	    hCutFlow_w->GetXaxis()->SetBinLabel(bin, x.first.c_str());
	    bin++;
	}
// /////////////////////////// Electron Trigger SF ///////////////////////////  (TSFel)
// ======================== Configurações de log ========================

auto logEvent = [&](const char* type, float eta, float pt, int binX, int binY, float val, bool out_of_range = false) {
    if (this->sf_summary_log_file && this->sf_summary_log_file->is_open()) {
        *(this->sf_summary_log_file) << "[" << type << "] eta: " << eta << ", pt: " << pt << "\n";
        *(this->sf_summary_log_file) << " → BinX: " << binX << ", BinY: " << binY
                                     << ", Value: " << val
                                     << (out_of_range ? " *** OUT OF RANGE ***" : "") << "\n";
    }
};

// ======================== Elétrons ========================
auto electrons = thisEvent->getSelElectrons();
if (electrons && !electrons->empty()) {
    std::vector<std::tuple<float,float,int,int,float,bool>> ele_log; // guardar info de log
    for (objectLep* ele : *electrons) {
        float pt = ele->getp4()->Pt();
        float eta = ele->getp4()->Eta();
        this->total_electrons_processed++;

        // Encontrar bin do SF
        int binX = h2_eleTrigSF->GetXaxis()->FindBin(eta);
        int binY = h2_eleTrigSF->GetYaxis()->FindBin(pt);
        float sf_val = 0;
        bool out_of_range = false;

        if (binX < 1 || binX > h2_eleTrigSF->GetNbinsX() ||
            binY < 1 || binY > h2_eleTrigSF->GetNbinsY()) {
            out_of_range = true;
            sf_val = 0; // ou 1.0
        } else {
            sf_val = h2_eleTrigSF->GetBinContent(binX, binY);
        }

        // Guardar info para log
        ele_log.push_back({eta, pt, binX, binY, sf_val, out_of_range});

        // Preencher histogramas
        if (!out_of_range) {
            h_sf_vs_pt_sum->Fill(pt, sf_val);
            h_sf_vs_pt_count->Fill(pt, 1);
            h_sf_vs_eta_sum->Fill(eta, sf_val);
            h_sf_vs_eta_count->Fill(eta, 1);
        }

        // Eficiência MC
        int binX_eff = h2_effMC->GetXaxis()->FindBin(eta);
        int binY_eff = h2_effMC->GetYaxis()->FindBin(pt);
        if (binX_eff >= 1 && binX_eff <= h2_effMC->GetNbinsX() &&
            binY_eff >= 1 && binY_eff <= h2_effMC->GetNbinsY()) {
            float eff_val = h2_effMC->GetBinContent(binX_eff, binY_eff);
            h_effMC_vs_pt_sum->Fill(pt, eff_val);
            h_effMC_vs_pt_count->Fill(pt, 1);
            h_effMC_vs_eta_sum->Fill(eta, eff_val);
            h_effMC_vs_eta_count->Fill(eta, 1);
        }
    }

    // Logar apenas eventos 10000, 20000, 30000...
    if (_entryInLoop % LOG_INTERVAL == 0) {
        for (auto& t : ele_log) {
            float eta, pt, val;
            int binX, binY;
            bool out_of_range;
            std::tie(eta, pt, binX, binY, val, out_of_range) = t;
            logEvent("Electron", eta, pt, binX, binY, val, out_of_range);
        }
    }
}

// ======================== Múons ========================
auto muons = thisEvent->getSelMuons();
if (muons && !muons->empty()) {
    std::vector<std::tuple<float,float,int,int,float,bool>> mu_log; // guardar info de log
    for (objectLep* mu : *muons) {
        float pt = mu->getp4()->Pt();
        float eta = mu->getp4()->Eta();
        this->total_muons_processed++;

        float sf_val = getMuonTrigSF(eta, pt);
        int eta_bin = -1, pt_bin = -1;

        // Encontrar bins do JSON (ATENÇÃO muonTrigSFJson SÓ ESTÁ VERIFICANDO OS SFs dos muons de PT médio!)
        if (!muonTrigSFJson.empty()) {
            const json* correction = nullptr;
            for (const auto& corr : muonTrigSFJson["corrections"]) {
                if (corr.contains("name") && corr["name"] == "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight") {
                    correction = &corr; break;
                }
            }
            if (correction) {
                const auto& data_eta = (*correction)["data"];
                const std::vector<float> eta_edges = data_eta["edges"].get<std::vector<float>>();
                for (size_t i = 0; i < eta_edges.size() - 1; ++i) {
                    if (eta >= eta_edges[i] && eta < eta_edges[i+1]) { eta_bin = i; break; }
                }
                if (eta_bin == -1 && eta == eta_edges.back()) eta_bin = int(eta_edges.size()) - 2;

                if (eta_bin != -1) {
                    const auto& data_pt = data_eta["content"][eta_bin];
                    const std::vector<float> pt_edges = data_pt["edges"].get<std::vector<float>>();
                    for (size_t i = 0; i < pt_edges.size() - 1; ++i) {
                        if (pt >= pt_edges[i] && pt < pt_edges[i+1]) { pt_bin = i; break; }
                    }
                    if (pt_bin == -1 && pt == pt_edges.back()) pt_bin = int(pt_edges.size()) - 2;
                }
            }
        }

        bool out_of_range = (eta_bin < 0 || pt_bin < 0);
        mu_log.push_back({eta, pt, eta_bin, pt_bin, sf_val, out_of_range});

        if (h_sf_muon_vs_pt_sum)    h_sf_muon_vs_pt_sum->Fill(pt, sf_val);
        if (h_sf_muon_vs_pt_count)  h_sf_muon_vs_pt_count->Fill(pt, 1);
        if (h_sf_muon_vs_eta_sum)   h_sf_muon_vs_eta_sum->Fill(eta, sf_val);
        if (h_sf_muon_vs_eta_count) h_sf_muon_vs_eta_count->Fill(eta, 1);
    }

    // Logar apenas eventos 10000, 20000, 30000...
    if (_entryInLoop % LOG_INTERVAL == 0) {
        for (auto& t : mu_log) {
            float eta, pt, val;
            int binX, binY;
            bool out_of_range;
            std::tie(eta, pt, binX, binY, val, out_of_range) = t;

            if (!out_of_range) {
                *(this->sf_summary_log_file) << "[Muon] eta: " << eta << ", pt: " << pt
                                             << " → SF JSON eta_bin: " << binX
                                             << ", pt_bin: " << binY
                                             << ", SF value: " << val << "\n";
            } else {
                *(this->sf_summary_log_file) << "[Muon] eta: " << eta << ", pt: " << pt
                                             << " → WARNING: Muon out of JSON bin range, SF fallback 1.0\n";
            }
        }
    }
}


	

    thisEvent->getCentrality(thisEvent->getSelJets(), thisEvent->getSelbJets(), jbjetCent);
    thisEvent->getCentrality(thisEvent->getSelJets(), thisEvent->getSelLeptons(), jlepCent);
    thisEvent->getStats(thisEvent->getSelJets(), jetStat);
    thisEvent->getStats(thisEvent->getSelbJets(), bjetStat);
    thisEvent->getStatsComb(thisEvent->getSelJets(), thisEvent->getSelbJets(), bjStat);
    thisEvent->getMaxPTComb(thisEvent->getSelJets(), thisEvent->getSelbJets(), jbbMaxs);
    thisEvent->getMaxPTSame(thisEvent->getSelJets(), jjjMaxs);

    thisEvent->getStatsComb(thisEvent->getSelJets(), thisEvent->getSelLeptons(), ljetStat);
    thisEvent->getStatsComb(thisEvent->getSelbJets(), thisEvent->getSelLeptons(), lbjetStat);
    
    thisEvent->getFoxWolfram(thisEvent->getSelJets(), jetFoxWolfMom);
    thisEvent->getFoxWolfram(thisEvent->getSelbJets(), bjetFoxWolfMom);
 

    //    std::cout << "Number of Hadronic Higgs: " << thisEvent->getnHadronicHiggs() << std::endl;

    hjetNumber->Fill(thisEvent->getnSelJet(),_weight*thisEvent->getbTagSys());
    hHadronicHiggsNumber->Fill(thisEvent->getnHadronicHiggs(),_weight*thisEvent->getbTagSys());
    hBjetNumber->Fill(thisEvent->getnbJet(),_weight*thisEvent->getbTagSys());
    hLightJetNumber->Fill(thisEvent->getnLightJet(),_weight*thisEvent->getbTagSys());

    hjetAverageMass->Fill(thisEvent->getSumSelJetMass()/(float)thisEvent->getnSelJet(),_weight*thisEvent->getbTagSys());
    hHadronicHiggsAverageMass->Fill(thisEvent->getSumSelHadronicHiggsMass()/(float)thisEvent->getnHadronicHiggs(),_weight*thisEvent->getbTagSys());
    hBjetAverageMass->Fill(thisEvent->getSumSelbJetMass()/(float)thisEvent->getnbJet(),_weight*thisEvent->getbTagSys());
    hLightJetAverageMass->Fill(thisEvent->getSumSelLightJetMass()/(float)thisEvent->getnLightJet(),_weight*thisEvent->getbTagSys());
    hBjetAverageMassSqr->Fill((thisEvent->getSumSelbJetMass()*thisEvent->getSumSelbJetMass())/(float)thisEvent->getnbJet(), _weight*thisEvent->getbTagSys());

    if(thisEvent->getnHadronicHiggs() > 0){ 
	hHadronicHiggsSoftDropMass1->Fill(thisEvent->getSelHadronicHiggses()->at(0)->softDropMass,_weight*thisEvent->getbTagSys());
    }

    if(thisEvent->getnHadronicHiggs() > 1){ 
	hHadronicHiggsSoftDropMass2->Fill(thisEvent->getSelHadronicHiggses()->at(1)->softDropMass,_weight*thisEvent->getbTagSys());
    }


    hjetHT->Fill(thisEvent->getSumSelJetScalarpT(),_weight*thisEvent->getbTagSys());
    hBjetHT->Fill(thisEvent->getSumSelbJetScalarpT(),_weight*thisEvent->getbTagSys());
    hHadronicHiggsHT->Fill(thisEvent->getSumSelHadronicHiggsScalarpT(),_weight*thisEvent->getbTagSys());
    hLightJetHT->Fill(thisEvent->getSumSelLightJetScalarpT(),_weight*thisEvent->getbTagSys());
    hAvgDeltaRjj->Fill(jetStat.meandR,_weight*thisEvent->getbTagSys());
    hAvgDeltaEtajj->Fill(jetStat.meandEta,_weight*thisEvent->getbTagSys());
    hminDeltaRjj->Fill(jetStat.mindR,_weight*thisEvent->getbTagSys());
    hminDeltaRMassjj->Fill(jetStat.mindRMass,_weight*thisEvent->getbTagSys());
    hminDeltaRpTjj->Fill(jetStat.mindRpT,_weight*thisEvent->getbTagSys());
    hAvgDeltaRbb->Fill(bjetStat.meandR,_weight*thisEvent->getbTagSys());
    hAvgDeltaRbj->Fill(bjStat.meandR,_weight*thisEvent->getbTagSys());
    hAvgDeltaEtabj->Fill(bjStat.meandEta,_weight*thisEvent->getbTagSys());
    hminDeltaRbj->Fill(bjStat.mindR,_weight*thisEvent->getbTagSys());
    hminDeltaRMassbj->Fill(bjStat.mindRMass,_weight*thisEvent->getbTagSys());
    hminDeltaRpTbj->Fill(bjStat.mindRpT,_weight*thisEvent->getbTagSys());    
    hAvgDeltaEtabb->Fill(bjetStat.meandEta,_weight*thisEvent->getbTagSys());
    hminDeltaRbb->Fill(bjetStat.mindR,_weight*thisEvent->getbTagSys());
    hminDeltaRMassbb->Fill(bjetStat.mindRMass,_weight*thisEvent->getbTagSys());
    hminDeltaRpTbb->Fill(bjetStat.mindRpT,_weight*thisEvent->getbTagSys());
    hmaxDeltaEtabb->Fill(bjetStat.maxdEta,_weight*thisEvent->getbTagSys());
    hmaxDeltaEtajj->Fill(jetStat.maxdEta,_weight*thisEvent->getbTagSys());
    hAvgDeltaEtabj->Fill(bjStat.meandEta,_weight*thisEvent->getbTagSys());
    hmaxDeltaEtabj->Fill(bjStat.maxdEta,_weight*thisEvent->getbTagSys());
    hmaxPTmassjbb->Fill(jbbMaxs.maxPTmass, _weight*thisEvent->getbTagSys());
    hmaxPTmassjjj->Fill(jjjMaxs.maxPTmass, _weight*thisEvent->getbTagSys());

    hInvMassHSingleMatched->Fill(_bbMassMinSHiggsMatched,_weight*thisEvent->getbTagSys());
    hInvMassHSingleNotMatched->Fill(_bbMassMinSHiggsNotMatched,_weight*thisEvent->getbTagSys());
    hChi2HiggsSingleMatched->Fill(_minChi2SHiggsMatched,_weight*thisEvent->getbTagSys());
    hChi2HiggsSingleNotMatched->Fill(_minChi2SHiggsNotMatched,_weight*thisEvent->getbTagSys());

    hInvMassHH1Matched->Fill(_bbMassMinHH1Matched,_weight*thisEvent->getbTagSys());
    hInvMassHH1NotMatched->Fill(_bbMassMinHH1NotMatched,_weight*thisEvent->getbTagSys());
    hInvMassHH2Matched->Fill(_bbMassMinHH2Matched,_weight*thisEvent->getbTagSys());
    hInvMassHH2NotMatched->Fill(_bbMassMinHH2NotMatched,_weight*thisEvent->getbTagSys());
    hChi2HHMatched->Fill(_minChi2HHMatched,_weight*thisEvent->getbTagSys());
    hChi2HHNotMatched->Fill(_minChi2HHNotMatched,_weight*thisEvent->getbTagSys());


    hInvMassH1->Fill(_bbMassMin1Higgs,_weight*thisEvent->getbTagSys());
    hInvMassH2->Fill(_bbMassMin2Higgs,_weight*thisEvent->getbTagSys());
    hPTH1->Fill(_bpTHiggs1,_weight*thisEvent->getbTagSys());
    hPTH2->Fill(_bpTHiggs2,_weight*thisEvent->getbTagSys());
    hInvMassZ1->Fill(_bbMassMin1Z,_weight*thisEvent->getbTagSys());
    hInvMassZ2->Fill(_bbMassMin2Z,_weight*thisEvent->getbTagSys());
    hInvMassHZ1->Fill(_bbMassMin1HiggsZ,_weight*thisEvent->getbTagSys());
    hInvMassHZ2->Fill(_bbMassMin2HiggsZ,_weight*thisEvent->getbTagSys());
    if(fabs(_bbMassMin1Higgs-cHiggsMass) < fabs(_bbMassMin2Higgs-cHiggsMass)){
	hInvMassH1mChi->Fill(_bbMassMin1Higgs,_weight*thisEvent->getbTagSys());
	hInvMassH2mChi->Fill(_bbMassMin2Higgs,_weight*thisEvent->getbTagSys());
    } else  {
	hInvMassH2mChi->Fill(_bbMassMin1Higgs,_weight*thisEvent->getbTagSys());
	hInvMassH1mChi->Fill(_bbMassMin2Higgs,_weight*thisEvent->getbTagSys());
    }
    hChi2Higgs->Fill(_minChi2Higgs,_weight*thisEvent->getbTagSys());
    hChi2Z->Fill(_minChi2Z,_weight*thisEvent->getbTagSys());
    hChi2HiggsZ->Fill(_minChi2HiggsZ,_weight*thisEvent->getbTagSys());

    
    hmet->Fill(thisEvent->getMET()->getp4()->Pt(),_weight*thisEvent->getbTagSys());
    //    hmetPhi->Fill(thisEvent->getMET()->getp4()->Phi(),_weight*thisEvent->getbTagSys());
    // hmetEta->Fill(thisEvent->getMET()->getp4()->Eta(),_weight*thisEvent->getbTagSys());
    
    hAplanarity->Fill(thisEvent->eventShapeJet->getAplanarity(), _weight*thisEvent->getbTagSys());
    hSphericity->Fill(thisEvent->eventShapeJet->getSphericity(), _weight*thisEvent->getbTagSys());
    hTransSphericity->Fill(thisEvent->eventShapeJet->getTransSphericity(), _weight*thisEvent->getbTagSys());
    hCvalue->Fill(thisEvent->eventShapeJet->getC(), _weight*thisEvent->getbTagSys());
    hDvalue->Fill(thisEvent->eventShapeJet->getD(), _weight*thisEvent->getbTagSys());
    hCentralityjb->Fill(jbjetCent.centrality, _weight*thisEvent->getbTagSys());    
    hCentralityjl->Fill(jlepCent.centrality, _weight*thisEvent->getbTagSys());    

    hH0->Fill(jetFoxWolfMom.h0, _weight*thisEvent->getbTagSys());
    hH1->Fill(jetFoxWolfMom.h1, _weight*thisEvent->getbTagSys());
    hH2->Fill(jetFoxWolfMom.h2, _weight*thisEvent->getbTagSys());
    hH3->Fill(jetFoxWolfMom.h3, _weight*thisEvent->getbTagSys());
    hH4->Fill(jetFoxWolfMom.h4, _weight*thisEvent->getbTagSys());
    hR1->Fill(jetFoxWolfMom.r1, _weight*thisEvent->getbTagSys());
    hR2->Fill(jetFoxWolfMom.r2, _weight*thisEvent->getbTagSys());
    hR3->Fill(jetFoxWolfMom.r3, _weight*thisEvent->getbTagSys());
    hR4->Fill(jetFoxWolfMom.r4, _weight*thisEvent->getbTagSys()); 


    hBjetH0->Fill(bjetFoxWolfMom.h0, _weight*thisEvent->getbTagSys());
    hBjetH1->Fill(bjetFoxWolfMom.h1, _weight*thisEvent->getbTagSys());
    hBjetH2->Fill(bjetFoxWolfMom.h2, _weight*thisEvent->getbTagSys());
    hBjetH3->Fill(bjetFoxWolfMom.h3, _weight*thisEvent->getbTagSys());
    hBjetH4->Fill(bjetFoxWolfMom.h4, _weight*thisEvent->getbTagSys());
    hBjetR1->Fill(bjetFoxWolfMom.r1, _weight*thisEvent->getbTagSys());
    hBjetR2->Fill(bjetFoxWolfMom.r2, _weight*thisEvent->getbTagSys());
    hBjetR3->Fill(bjetFoxWolfMom.r3, _weight*thisEvent->getbTagSys());
    hBjetR4->Fill(bjetFoxWolfMom.r4, _weight*thisEvent->getbTagSys()); 


    hBjetAplanarity->Fill(thisEvent->eventShapeBjet->getAplanarity(), _weight*thisEvent->getbTagSys());
    hBjetSphericity->Fill(thisEvent->eventShapeBjet->getSphericity(), _weight*thisEvent->getbTagSys());
    hBjetTransSphericity->Fill(thisEvent->eventShapeBjet->getTransSphericity(), _weight*thisEvent->getbTagSys());
    hBjetCvalue->Fill(thisEvent->eventShapeBjet->getC(), _weight*thisEvent->getbTagSys());
    hBjetDvalue->Fill(thisEvent->eventShapeBjet->getD(), _weight*thisEvent->getbTagSys());


    for(int ih=0; ih < thisEvent->getnSelJet() && ih < nHistsJets; ih++){
	hjetsPTs.at(ih)->Fill(thisEvent->getSelJets()->at(ih)->getp4()->Pt(),_weight*thisEvent->getbTagSys());
	hjetsEtas.at(ih)->Fill(thisEvent->getSelJets()->at(ih)->getp4()->Eta(),_weight*thisEvent->getbTagSys());
	hjetsBTagDisc.at(ih)->Fill(getJetCSV(thisEvent).at(ih),_weight*thisEvent->getbTagSys());
    }

    for(int ih=0; ih < thisEvent->getnbJet() && ih < nHistsbJets; ih++){
	hbjetsPTs.at(ih)->Fill(thisEvent->getSelbJets()->at(ih)->getp4()->Pt(),_weight*thisEvent->getbTagSys());
	hbjetsEtas.at(ih)->Fill(thisEvent->getSelbJets()->at(ih)->getp4()->Eta(),_weight*thisEvent->getbTagSys());
	hbjetsBTagDisc.at(ih)->Fill(getbJetCSV(thisEvent).at(ih),_weight*thisEvent->getbTagSys());
    }


    for(int ih=0; ih < thisEvent->getnLightJet() && ih < nHistsLightJets; ih++){
	hLightJetsPTs.at(ih)->Fill(thisEvent->getSelLightJets()->at(ih)->getp4()->Pt(),_weight*thisEvent->getbTagSys());
	hLightJetsEtas.at(ih)->Fill(thisEvent->getSelLightJets()->at(ih)->getp4()->Eta(),_weight*thisEvent->getbTagSys());
	hLightJetsBTagDisc.at(ih)->Fill(getlightJetCSV(thisEvent).at(ih),_weight*thisEvent->getbTagSys());
    }

    hleptonNumber->Fill(thisEvent->getnSelLepton(),_weight*thisEvent->getbTagSys());
    hElecNumber->Fill(thisEvent->getnSelElectron(),_weight*thisEvent->getbTagSys());
    hMuonNumber->Fill(thisEvent->getnSelMuon(),_weight*thisEvent->getbTagSys());

//    if(thisEvent->getnSelMuon() == 2){
//	hDiMuonMass->Fill(thisEvent->getSelMuonsMass(),_weight*thisEvent->getbTagSys());
//	hDiMuonPT->Fill(thisEvent->getSelMuonsPT(),_weight*thisEvent->getbTagSys());
//	hDiMuonEta->Fill(thisEvent->getSelMuonsEta(),_weight*thisEvent->getbTagSys());
//    }

//    if(thisEvent->getnSelElectron() == 2){
//	hDiElectronMass->Fill(thisEvent->getSelElectronsMass(),_weight*thisEvent->getbTagSys());
//	hDiElectronPT->Fill(thisEvent->getSelElectronsPT(),_weight*thisEvent->getbTagSys());
//	hDiElectronEta->Fill(thisEvent->getSelElectronsEta(),_weight*thisEvent->getbTagSys());
//    }

    hleptonHT->Fill(thisEvent->getSelLeptonHT(),_weight*thisEvent->getbTagSys());
    hST->Fill(thisEvent->getSelLeptonST(),_weight*thisEvent->getbTagSys());
    hLeptonPT1->Fill(thisEvent->getSelLeptons()->at(0)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
    hLeptonEta1->Fill(thisEvent->getSelLeptons()->at(0)->getp4()->Eta(), _weight*thisEvent->getbTagSys());
   // hLeptonPT2->Fill(thisEvent->getSelLeptons()->at(1)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
  //  hLeptonEta2->Fill(thisEvent->getSelLeptons()->at(1)->getp4()->Eta(), _weight*thisEvent->getbTagSys());


    if(thisEvent->getnSelMuon() > 0){
	hMuonPT1->Fill(thisEvent->getSelMuons()->at(0)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
	hMuonEta1->Fill(thisEvent->getSelMuons()->at(0)->getp4()->Eta(), _weight*thisEvent->getbTagSys());
    }
    
    if(thisEvent->getnSelElectron() > 0){
	hElePT1->Fill(thisEvent->getSelElectrons()->at(0)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
	hEleEta1->Fill(thisEvent->getSelElectrons()->at(0)->getp4()->Eta(), _weight*thisEvent->getbTagSys());
    }
  /*  
    if(thisEvent->getnSelMuon() > 1){
	hMuonPT2->Fill(thisEvent->getSelMuons()->at(1)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
	hMuonEta2->Fill(thisEvent->getSelMuons()->at(1)->getp4()->Eta(), _weight*thisEvent->getbTagSys());
    }
    
    if(thisEvent->getnSelElectron() > 1){
	hElePT2->Fill(thisEvent->getSelElectrons()->at(1)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
	hEleEta2->Fill(thisEvent->getSelElectrons()->at(1)->getp4()->Eta(), _weight*thisEvent->getbTagSys());
    } 
*/
    hLepCharge1->Fill(thisEvent->getSelLeptons()->at(0)->charge, _weight*thisEvent->getbTagSys());
//    hLepCharge2->Fill(thisEvent->getSelLeptons()->at(1)->charge, _weight*thisEvent->getbTagSys());
    
}



void ttHHanalyzer::writeHistos() {


// ========= Cálculo das médias com erro correto =========
auto calcAverage = [](TH1F* sumH, TH1F* countH, TH1F*& avgH, const char* name) {
    if (!sumH || !countH) return;
    if (avgH) delete avgH;
    avgH = (TH1F*) sumH->Clone(name);
    avgH->Reset();
    for (int i = 1; i <= sumH->GetNbinsX(); ++i) {
        double sum = sumH->GetBinContent(i);
        double sum2 = sumH->GetBinError(i) * sumH->GetBinError(i);
        double count = countH->GetBinContent(i);
        if (count > 0) {
            avgH->SetBinContent(i, sum / count);
            avgH->SetBinError(i, std::sqrt(sum2) / count);
        } else {
            avgH->SetBinContent(i, 0);
            avgH->SetBinError(i, 0);
        }
    }
    avgH->SetDirectory(0);
};

// Calcular médias
calcAverage(h_sf_vs_pt_sum, h_sf_vs_pt_count, h_sf_vs_pt_avg, "h_sf_vs_pt_avg");
calcAverage(h_sf_vs_eta_sum, h_sf_vs_eta_count, h_sf_vs_eta_avg, "h_sf_vs_eta_avg");
calcAverage(h_effMC_vs_pt_sum, h_effMC_vs_pt_count, h_effMC_vs_pt_avg, "h_effMC_vs_pt_avg");
calcAverage(h_effMC_vs_eta_sum, h_effMC_vs_eta_count, h_effMC_vs_eta_avg, "h_effMC_vs_eta_avg");
calcAverage(h_sf_muon_vs_pt_sum, h_sf_muon_vs_pt_count, h_sf_muon_vs_pt_avg, "h_sf_muon_vs_pt_avg");
calcAverage(h_sf_muon_vs_eta_sum, h_sf_muon_vs_eta_count, h_sf_muon_vs_eta_avg, "h_sf_muon_vs_eta_avg");

// ========= Função de log de resumo =========
	/*
auto logSFSummary = [&](std::ofstream* logFile) {
    if (!logFile || !logFile->is_open()) return;

    *logFile << "\n========= SF SUMMARY =========\n";
    *logFile << "Total Electrons Processed: " << this->total_electrons_processed << "\n";
    *logFile << "Total Muons Processed: " << this->total_muons_processed << "\n";

    auto logAverage = [&](TH1F* avgH, TH1F* countH, const char* type, const char* var) {
        if (!avgH || !countH) return;
        *logFile << "\n[" << type << " Avg per " << var << " Bin]\n";
        for (int i = 1; i <= avgH->GetNbinsX(); ++i) {
            double count = countH->GetBinContent(i);
            if (count > 0) {
                double avg = avgH->GetBinContent(i);
                *logFile << "Bin " << i << " (~" << avgH->GetBinCenter(i) << "): "
                         << avg << " [entries: " << count << "]\n";
            }
        }
    };


    // Exemplo de uso do helper logAverage
    logAverage(h_sf_vs_pt_sum, h_sf_vs_pt_count, "Electron SF", "pt");
    logAverage(h_sf_vs_eta_sum, h_sf_vs_eta_count, "Electron SF", "eta");    
	logAverage(h_effMC_vs_pt_avg, h_effMC_vs_pt_count, "Electron EffMC", "pT");
    logAverage(h_effMC_vs_eta_avg, h_effMC_vs_eta_count, "Electron EffMC", "η");
    logAverage(h_sf_muon_vs_pt_sum, h_sf_muon_vs_pt_count, "Muon SF", "pt");
    logAverage(h_sf_muon_vs_eta_sum, h_sf_muon_vs_eta_count, "Muon SF", "eta");
};


	*/

	
	
    _of->file->cd();
    _histoDirs.at(0)->cd();
    for(int ih=0; ih<nHistsJets; ih++){
	hjetsPTs.at(ih)->Write();
	hjetsEtas.at(ih)->Write();
	hjetsBTagDisc.at(ih)->Write();
    }
    for(int ih=0; ih<nHistsbJets; ih++){
	hbjetsPTs.at(ih)->Write();
	hbjetsEtas.at(ih)->Write();
	hbjetsBTagDisc.at(ih)->Write();
    }

    for(int ih=0; ih<nHistsLightJets; ih++){
	hLightJetsPTs.at(ih)->Write();
	hLightJetsEtas.at(ih)->Write();
	hLightJetsBTagDisc.at(ih)->Write();
    }
	
    hInvMassHSingleMatched->Write();
    hInvMassHSingleNotMatched->Write();
    hChi2HiggsSingleNotMatched->Write();
    hChi2HiggsSingleMatched->Write();
    hInvMassHH1Matched->Write();
    hInvMassHH1NotMatched->Write();
    hInvMassHH2Matched->Write();
    hInvMassHH2NotMatched->Write();
    hChi2HHNotMatched->Write();
    hChi2HHMatched->Write();

    hjetNumber->Write();
    hBjetNumber->Write();
    hHadronicHiggsNumber->Write();
    hLightJetNumber->Write();
    hjetAverageMass->Write();
    hBjetAverageMass->Write();
    hHadronicHiggsAverageMass->Write();
    hLightJetAverageMass->Write();
    hBjetAverageMassSqr->Write();
    hHadronicHiggsSoftDropMass1->Write();
    hHadronicHiggsSoftDropMass2->Write();
    hjetHT->Write();
    hBjetHT->Write();
    hHadronicHiggsHT->Write();
    hLightJetHT->Write();
    hAvgDeltaRjj->Write();
    hminDeltaRjj->Write();
    hminDeltaRMassjj->Write();
    hminDeltaRpTjj->Write();
    hAvgDeltaRbb->Write();
    hAvgDeltaEtajj->Write();
    hAvgDeltaEtabb->Write();
    hminDeltaRbb->Write();
    hminDeltaRMassbb->Write();
    hminDeltaRpTbb->Write();
    hmaxDeltaEtabb->Write();
    hmaxDeltaEtajj->Write();
    hAvgDeltaRbj->Write();
    hAvgDeltaEtabj->Write();
    hminDeltaRbj->Write();
    hminDeltaRMassbj->Write();
    hminDeltaRpTbj->Write();
    hmaxDeltaEtabj->Write();
    hmaxPTmassjbb->Write();
    hmaxPTmassjjj->Write();

    hPTH1->Write();
    hPTH2->Write();
    hInvMassH1->Write();
    hInvMassH2->Write();
    hInvMassH1mChi->Write();
    hInvMassH2mChi->Write();
    hInvMassHZ1->Write();
    hInvMassHZ2->Write();
    hInvMassZ1->Write();
    hInvMassZ2->Write();
    hChi2Higgs->Write();
    hChi2HiggsZ->Write();
    hChi2Z->Write();


    hmet->Write();
    // hmetPhi->Write();
    // hmetEta->Write();

    hAplanarity->Write();
    hSphericity->Write();
    hTransSphericity->Write();
    hCvalue->Write();
    hDvalue->Write();
    hCentralityjb->Write();
    hCentralityjl->Write();

    hH0->Write();
    hH1->Write();
    hH2->Write();
    hH3->Write();
    hH4->Write();
    hR1->Write();
    hR2->Write();
    hR3->Write();
    hR4->Write(); 

    hBjetH0->Write();
    hBjetH1->Write();
    hBjetH2->Write();
    hBjetH3->Write();
    hBjetH4->Write();
    hBjetR1->Write();
    hBjetR2->Write();
    hBjetR3->Write();
    hBjetR4->Write(); 

    hBjetAplanarity->Write();
    hBjetSphericity->Write();
    hBjetTransSphericity->Write();
    hBjetCvalue->Write();
    hBjetDvalue->Write();

    _histoDirs.at(1)->cd();

// ===SF Trigger Electron===  (TSFel)
if (h_sf_vs_pt)          h_sf_vs_pt->Write();
if (h_sf_vs_eta)         h_sf_vs_eta->Write();
if (h_effMC_vs_pt)       h_effMC_vs_pt->Write();
if (h_effMC_vs_eta)      h_effMC_vs_eta->Write();

// Escreve os histogramas de média
if (h_sf_vs_pt_avg)     h_sf_vs_pt_avg->Write();
if (h_sf_vs_eta_avg)    h_sf_vs_eta_avg->Write();
if (h_effMC_vs_pt_avg)  h_effMC_vs_pt_avg->Write();
if (h_effMC_vs_eta_avg) h_effMC_vs_eta_avg->Write();


//Histograms for muons
// === SF Trigger Muon ===  (TSFmu)
if (h_sf_muon_vs_pt)      h_sf_muon_vs_pt->Write();
if (h_sf_muon_vs_eta)     h_sf_muon_vs_eta->Write();
if (h_sf_muon_vs_pt_avg)  h_sf_muon_vs_pt_avg->Write();
if (h_sf_muon_vs_eta_avg) h_sf_muon_vs_eta_avg->Write();

	
    hLepCharge1->Write();
  //  hLepCharge2->Write();

    hleptonNumber->Write();
    hElecNumber->Write();
    hMuonNumber->Write();
   // hDiMuonMass->Write();
   // hDiElectronMass->Write();
    hDiMuonPT->Write();
    hDiElectronPT->Write();
    hDiMuonEta->Write();
    hDiElectronEta->Write();
    hleptonHT->Write();
    hST->Write();
    hLeptonEta1->Write();
    hLeptonPT1->Write();
 //   hLeptonEta2->Write();
 //   hLeptonPT2->Write();
    
    hMuonEta1->Write();
    hMuonPT1->Write();
    hEleEta1->Write();
    hElePT1->Write();
//    hMuonEta2->Write();
//    hMuonPT2->Write();
//    hEleEta2->Write();
//    hElePT2->Write(); 
if (sf_summary_log_file && sf_summary_log_file->is_open()) {
    sf_summary_log_file->close();
}
	
}
////end of writeHistos

void ttHHanalyzer::fillTree(event * thisEvent){

   
    jetPT1 = thisEvent->getSelJets()->at(0)->getp4()->Pt();
    jetPT2 = thisEvent->getSelJets()->at(1)->getp4()->Pt();
    jetPT3 = thisEvent->getSelJets()->at(2)->getp4()->Pt();
    jetPT4 = thisEvent->getSelJets()->at(3)->getp4()->Pt();
    jetEta1 = thisEvent->getSelJets()->at(0)->getp4()->Eta();
    jetEta2 = thisEvent->getSelJets()->at(1)->getp4()->Eta();
    jetEta3 = thisEvent->getSelJets()->at(2)->getp4()->Eta();
    jetEta4 = thisEvent->getSelJets()->at(3)->getp4()->Eta();
    jetBTagDisc1 = getJetCSV(thisEvent).at(0);
    jetBTagDisc2 = getJetCSV(thisEvent).at(1);
    jetBTagDisc3 = getJetCSV(thisEvent).at(2);
    jetBTagDisc4 = getJetCSV(thisEvent).at(3);


    bjetPT1 = thisEvent->getSelbJets()->at(0)->getp4()->Pt();
    bjetPT2 = thisEvent->getSelbJets()->at(1)->getp4()->Pt();
    bjetPT3 = thisEvent->getSelbJets()->at(2)->getp4()->Pt();
    bjetEta1 = thisEvent->getSelbJets()->at(0)->getp4()->Eta();
    bjetEta2 = thisEvent->getSelbJets()->at(1)->getp4()->Eta();
    bjetEta3 = thisEvent->getSelbJets()->at(2)->getp4()->Eta();
    bjetBTagDisc1 = getbJetCSV(thisEvent).at(0);
    bjetBTagDisc2 = getbJetCSV(thisEvent).at(1);
    bjetBTagDisc3 = getbJetCSV(thisEvent).at(2);
    bbjetHiggsMatched1 = thisEvent->getSelbJets()->at(0)->matchedtoHiggs;
    bbjetHiggsMatched2 = thisEvent->getSelbJets()->at(1)->matchedtoHiggs;
    bbjetHiggsMatched3 = thisEvent->getSelbJets()->at(2)->matchedtoHiggs;
    bbjetHiggsMatcheddR1 = thisEvent->getSelbJets()->at(0)->matchedtoHiggsdR;
    bbjetHiggsMatcheddR2 = thisEvent->getSelbJets()->at(1)->matchedtoHiggsdR;
    bbjetHiggsMatcheddR3 = thisEvent->getSelbJets()->at(2)->matchedtoHiggsdR;
    bbjetMinChiHiggsIndex1 = thisEvent->getSelbJets()->at(0)->minChiHiggsIndex;
    bbjetMinChiHiggsIndex2 = thisEvent->getSelbJets()->at(1)->minChiHiggsIndex;
    bbjetMinChiHiggsIndex3 = thisEvent->getSelbJets()->at(2)->minChiHiggsIndex;


    if(thisEvent->getnSelJet() > 4){
    	jetPT5 = thisEvent->getSelJets()->at(4)->getp4()->Pt();
	jetEta5 = thisEvent->getSelJets()->at(4)->getp4()->Eta();
	jetBTagDisc5 = getJetCSV(thisEvent).at(4);
    } else{
	jetPT5 = -6;
	jetEta5 = -6;
	jetBTagDisc5 = -6;
    }
    
    if(thisEvent->getnSelJet() > 5){
	jetPT6 = thisEvent->getSelJets()->at(5)->getp4()->Pt();
	jetEta6 = thisEvent->getSelJets()->at(5)->getp4()->Eta();
	jetBTagDisc6 = getJetCSV(thisEvent).at(5);
    } else{
	jetPT6 = -6;
	jetEta6 = -6;
	jetBTagDisc6 = -6;
    }


    if(thisEvent->getnSelJet() > 6){
	jetPT7 = thisEvent->getSelJets()->at(6)->getp4()->Pt();
	jetEta7 = thisEvent->getSelJets()->at(6)->getp4()->Eta();
	jetBTagDisc7 = getJetCSV(thisEvent).at(6);
    } else{
	jetPT7 = -6;
	jetEta7 = -6;
	jetBTagDisc7 = -6;
    }

    if(thisEvent->getnSelJet() > 7){
	jetPT8 = thisEvent->getSelJets()->at(7)->getp4()->Pt();
	jetEta8 = thisEvent->getSelJets()->at(7)->getp4()->Eta();
	jetBTagDisc8 = getJetCSV(thisEvent).at(7);
    } else{
	jetPT8 = -6;
	jetEta8 = -6;
	jetBTagDisc8 = -6;
    }

    if(thisEvent->getnbJet() > 2){
	bjetPT3 = thisEvent->getSelbJets()->at(2)->getp4()->Pt();
	bjetEta3 = thisEvent->getSelbJets()->at(2)->getp4()->Eta();
	bjetBTagDisc3 = getbJetCSV(thisEvent).at(2);
	bbjetHiggsMatched3 = thisEvent->getSelbJets()->at(2)->matchedtoHiggs;
	bbjetHiggsMatcheddR3 = thisEvent->getSelbJets()->at(2)->matchedtoHiggsdR;
	bbjetMinChiHiggsIndex3 = thisEvent->getSelbJets()->at(2)->minChiHiggsIndex;
    } else {
	bjetPT3 = -6;
	bjetEta3 = -6;
	bjetBTagDisc3 = -6;
	bbjetHiggsMatched3 = 0;
	bbjetHiggsMatcheddR3 = -6;
	bbjetMinChiHiggsIndex3 = -6;
	} 


    if(thisEvent->getnbJet() > 3){
	bjetPT4 = thisEvent->getSelbJets()->at(3)->getp4()->Pt();
	bjetEta4 = thisEvent->getSelbJets()->at(3)->getp4()->Eta();
	bjetBTagDisc4 = getbJetCSV(thisEvent).at(3);
	bbjetHiggsMatched4 = thisEvent->getSelbJets()->at(3)->matchedtoHiggs;
	bbjetHiggsMatcheddR4 = thisEvent->getSelbJets()->at(3)->matchedtoHiggsdR;
	bbjetMinChiHiggsIndex4 = thisEvent->getSelbJets()->at(3)->minChiHiggsIndex;
    } else {
	bjetPT4 = -6;
	bjetEta4 = -6;
	bjetBTagDisc4 = -6;
	bbjetHiggsMatched4 = 0;
	bbjetHiggsMatcheddR4 = -6;
	bbjetMinChiHiggsIndex4 = -6;
    }
    
    if(thisEvent->getnbJet() > 4){
	bjetPT5 = thisEvent->getSelbJets()->at(4)->getp4()->Pt();
	bjetEta5 = thisEvent->getSelbJets()->at(4)->getp4()->Eta();
	bjetBTagDisc5 = getbJetCSV(thisEvent).at(4);
	bbjetHiggsMatched5 = thisEvent->getSelbJets()->at(4)->matchedtoHiggs;
	bbjetHiggsMatcheddR5 = thisEvent->getSelbJets()->at(4)->matchedtoHiggsdR;
	bbjetMinChiHiggsIndex5 = thisEvent->getSelbJets()->at(4)->minChiHiggsIndex;
    } else {
	bjetPT5 = -6;
	bjetEta5 = -6;
	bjetBTagDisc5 = -6;
	bbjetHiggsMatched5 = 0;
	bbjetHiggsMatcheddR5 = -6;
	bbjetMinChiHiggsIndex5 = -6;
    }

    if(thisEvent->getnbJet() > 5){
	bjetPT6 = thisEvent->getSelbJets()->at(5)->getp4()->Pt();
	bjetEta6 = thisEvent->getSelbJets()->at(5)->getp4()->Eta();
	bjetBTagDisc6 = getbJetCSV(thisEvent).at(5);
	bbjetHiggsMatched6 = thisEvent->getSelbJets()->at(5)->matchedtoHiggs;
	bbjetHiggsMatcheddR6 = thisEvent->getSelbJets()->at(5)->matchedtoHiggsdR;
	bbjetMinChiHiggsIndex6 = thisEvent->getSelbJets()->at(5)->minChiHiggsIndex;
    } else {
	bjetPT6 = -6;
	bjetEta6 = -6;
	bjetBTagDisc6 = -6;
	bbjetHiggsMatched6 = 0;
	bbjetHiggsMatcheddR6 = -6;
	bbjetMinChiHiggsIndex6 = -6;
    }

    if(thisEvent->getnbJet() > 6){
	bjetPT7 = thisEvent->getSelbJets()->at(6)->getp4()->Pt();
	bjetEta7 = thisEvent->getSelbJets()->at(6)->getp4()->Eta();
	bjetBTagDisc7 = getbJetCSV(thisEvent).at(6);
	//bbjetHiggsMatched7 = thisEvent->getSelbJets()->at(6)->matchedtoHiggs;
	//bbjetHiggsMatcheddR7 = thisEvent->getSelbJets()->at(6)->matchedtoHiggsdR;
    } else {
	bjetPT7 = -6;
	bjetEta7 = -6;
	bjetBTagDisc7 = -6;
	//bbjetHiggsMatched7 = 0;
	//bbjetHiggsMatcheddR7 = -6;
    }


    if(thisEvent->getnbJet() > 7){
	bjetPT8 = thisEvent->getSelbJets()->at(7)->getp4()->Pt();
	bjetEta8 = thisEvent->getSelbJets()->at(7)->getp4()->Eta();
	bjetBTagDisc8 = getbJetCSV(thisEvent).at(7);
	//bbjetHiggsMatched8 = thisEvent->getSelbJets()->at(7)->matchedtoHiggs;
	//bbjetHiggsMatcheddR8 = thisEvent->getSelbJets()->at(7)->matchedtoHiggsdR;
    } else {
	bjetPT8 = -6;
	bjetEta8 = -6;
	bjetBTagDisc8 = -6;
	//bbjetHiggsMatched8 = 0;
	//bbjetHiggsMatcheddR8 = 0;
	} 


    if(thisEvent->getnLightJet() > 0){
    	lightjetPT1 = thisEvent->getSelLightJets()->at(0)->getp4()->Pt();
	lightjetEta1 = thisEvent->getSelLightJets()->at(0)->getp4()->Eta();
	lightjetBTagDisc1 = getlightJetCSV(thisEvent).at(0);
    } else{
	lightjetPT1 = -6;
	lightjetEta1 = -6;
	lightjetBTagDisc1 = -6;
    }

    if(thisEvent->getnLightJet() > 1){
    	lightjetPT2 = thisEvent->getSelLightJets()->at(1)->getp4()->Pt();
	lightjetEta2 = thisEvent->getSelLightJets()->at(1)->getp4()->Eta();
	lightjetBTagDisc2 = getlightJetCSV(thisEvent).at(1);
    } else{
	lightjetPT2 = -6;
	lightjetEta2 = -6;
	lightjetBTagDisc2 = -6;
    }

    if(thisEvent->getnLightJet() > 2){
    	lightjetPT3 = thisEvent->getSelLightJets()->at(2)->getp4()->Pt();
	lightjetEta3 = thisEvent->getSelLightJets()->at(2)->getp4()->Eta();
	lightjetBTagDisc3 = getlightJetCSV(thisEvent).at(2);
    } else{
	lightjetPT3 = -6;
	lightjetEta3 = -6;
	lightjetBTagDisc3 = -6;
    }


    weight= _weight;
    jetAverageMass = thisEvent->getSumSelJetMass()/thisEvent->getnSelJet();
    bjetAverageMass = thisEvent->getSumSelbJetMass()/thisEvent->getnbJet();
    lightJetAverageMass = thisEvent->getSumSelLightJetMass()/thisEvent->getnLightJet();
    bjetAverageMassSqr = (thisEvent->getSumSelbJetMass()*thisEvent->getSumSelbJetMass())/thisEvent->getnbJet();
    met = thisEvent->getMET()->getp4()->Pt();
    averageDeltaRjj = jetStat.meandR;
    averageDeltaRbb = bjetStat.meandR;
    averageDeltaRbj = bjStat.meandR;
    averageDeltaEtajj = jetStat.meandEta;
    averageDeltaEtabb = bjetStat.meandEta;
    averageDeltaEtabj = bjStat.meandEta;
    minDeltaRjj = jetStat.mindR;
    minDeltaRbb = bjetStat.mindR;
    minDeltaRbj = bjStat.mindR;
    maxDeltaEtabb = bjetStat.maxdEta;
    maxDeltaEtajj = jetStat.maxdEta;
    maxDeltaEtabj = bjStat.maxdEta;
    minDeltaRMassjj = jetStat.mindRMass;
    minDeltaRMassbb = bjetStat.mindRMass;
    minDeltaRMassbj = bjStat.mindRMass;
    minDeltaRpTjj = jetStat.mindRpT;
    minDeltaRpTbb = bjetStat.mindRpT;
    minDeltaRpTbj = bjStat.mindRpT;
    maxPTmassjjj = jjjMaxs.maxPTmass;
    maxPTmassjbb = jbbMaxs.maxPTmass;
    H0 = jetFoxWolfMom.h0;
    H1 = jetFoxWolfMom.h1;
    H2 = jetFoxWolfMom.h2;
    H3 = jetFoxWolfMom.h3;
    H4 = jetFoxWolfMom.h4;
    bH0 = bjetFoxWolfMom.h0;
    bH1 = bjetFoxWolfMom.h1;
    bH2 = bjetFoxWolfMom.h2;
    bH3 = bjetFoxWolfMom.h3;
    bH4 = bjetFoxWolfMom.h4;
    R1 = jetFoxWolfMom.r1;
    R2 = jetFoxWolfMom.r2;
    R3 = jetFoxWolfMom.r3;
    R4 = jetFoxWolfMom.r4;
    bR1 = bjetFoxWolfMom.r1;
    bR2 = bjetFoxWolfMom.r2;
    bR3 = bjetFoxWolfMom.r3;
    bR4 = bjetFoxWolfMom.r4;

    jetHT = thisEvent->getSumSelJetScalarpT();  
    bjetHT = thisEvent->getSumSelbJetScalarpT();
    lightjetHT = thisEvent->getSumSelLightJetScalarpT();
    jetNumber = thisEvent->getnSelJet();
    bjetNumber = thisEvent->getnbJet();
    lightjetNumber = thisEvent->getnLightJet();
    invMassZ1 = _bbMassMin1Z; //ZZ
    invMassZ2 = _bbMassMin2Z;
    chi2Z = _minChi2Z;
    invMassH1 = _bbMassMin1Higgs; //HH
    invMassH2 = _bbMassMin2Higgs;
    chi2Higgs = _minChi2Higgs;
    chi2HiggsZ = _minChi2HiggsZ; //ZH
    invMassHiggsZ1 = _bbMassMin1HiggsZ;
    invMassHiggsZ2 = _bbMassMin2HiggsZ;
    PTH1 = _bpTHiggs1;
    PTH2 = _bpTHiggs2;


    centralityjb = jbjetCent.centrality; 
    centralityjl = jlepCent.centrality; 
    aplanarity = thisEvent->eventShapeJet->getAplanarity();
    sphericity = thisEvent->eventShapeJet->getSphericity();
    transSphericity = thisEvent->eventShapeJet->getTransSphericity();
    cValue = thisEvent->eventShapeJet->getC();
    dValue = thisEvent->eventShapeJet->getD();
    baplanarity = thisEvent->eventShapeBjet->getAplanarity();
    bsphericity = thisEvent->eventShapeBjet->getSphericity();
    btransSphericity = thisEvent->eventShapeBjet->getTransSphericity();
    bcValue = thisEvent->eventShapeBjet->getC();
    bdValue = thisEvent->eventShapeBjet->getD();

    leptonPT1 = thisEvent->getSelLeptons()->at(0)->getp4()->Pt();
  //  leptonPT2 = thisEvent->getSelLeptons()->at(1)->getp4()->Pt();
    leptonEta1 = thisEvent->getSelLeptons()->at(0)->getp4()->Eta();
    //leptonEta2 = thisEvent->getSelLeptons()->at(1)->getp4()->Eta();
    leptonCharge1 = thisEvent->getSelLeptons()->at(0)->charge;
  //  leptonCharge2 = thisEvent->getSelLeptons()->at(1)->charge;
    leptonHT = thisEvent->getSelLeptonHT();
    ST = thisEvent->getSelLeptonST();


    if(thisEvent->getnSelMuon() > 0){
	muonPT1 = thisEvent->getSelMuons()->at(0)->getp4()->Pt();
	muonEta1 = thisEvent->getSelMuons()->at(0)->getp4()->Eta();
    } else {
	muonPT1 = -6;
	muonEta1 = -6;
    }
/*
    if(thisEvent->getnSelMuon() > 1){
	muonPT2 = thisEvent->getSelMuons()->at(1)->getp4()->Pt();
	muonEta2 = thisEvent->getSelMuons()->at(1)->getp4()->Eta();
	diMuonMass = thisEvent->getSelMuonsMass();
    } else {
	muonPT2 = -6;
	muonEta2 = -6;
	diMuonMass = -6;
    } 
*/
    if(thisEvent->getnSelElectron() > 0){
	elePT1 = thisEvent->getSelElectrons()->at(0)->getp4()->Pt();
	eleEta1 = thisEvent->getSelElectrons()->at(0)->getp4()->Eta();
    } else {
	elePT1 = -6;
	eleEta1 = -6;
    }
/*
    if(thisEvent->getnSelElectron() > 1){
	elePT2 = thisEvent->getSelElectrons()->at(1)->getp4()->Pt();
	eleEta2 = thisEvent->getSelElectrons()->at(1)->getp4()->Eta();
	diElectronMass = thisEvent->getSelElectronsMass();
    } else {
	elePT2 = -6;
	eleEta2 = -6;
	diElectronMass = -6;
    } 

  */  
    _inputTree->Fill();
}

void ttHHanalyzer::writeTree(){
    _of->file->cd();
    _treeDirs.at(0)->cd();
    _inputTree->Write();
    //    _inputTree->Delete();
    
}
//changed from here
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



int main(int argc, char** argv){
    commandLine cl(argc, argv);
    vector<string> filenames = fileNames(cl.filelist);
    double weight = cl.externalweight;   // Get global weight 
 
    // Create tree reader
    itreestream stream(filenames, "Events");
    if ( !stream.good() ) error("can't read root input files");

    eventBuffer ev(stream);
    std::cout << " Output filename: " << cl.outputfilename << std::endl;
    ////ttHHanalyzer analysis(cl.outputfilename, &ev, weight, true)
  
    // If you want to check or modify arguments,
    // Please check the [ src/tnm.cc ]
    // Arguments structure --> filelist, outputDirName, weight, Year, Data or MC, sampleName
    ttHHanalyzer analysis(cl.outputfilename, &ev, weight, true, cl.year, cl.isData, cl.sampleName);
    analysis.performAnalysis();

    ev.close();
    //    of.close();
    return 0;
}
