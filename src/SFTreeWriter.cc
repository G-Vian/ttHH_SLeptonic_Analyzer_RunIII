// Em src/SFTreeWriter.cc

#include "SFTreeWriter.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

// Método estático que cria e retorna a única instância
SFTreeWriter& SFTreeWriter::getInstance() {
    static SFTreeWriter instance; // Criada na primeira chamada, destruída no final do programa
    return instance;
}

// Construtor: apenas inicializa os ponteiros
SFTreeWriter::SFTreeWriter() 
    : _sf_output_file(nullptr), _sf_tree(nullptr), _isInitialized(false) {}

// Destrutor: chamado AUTOMATICAMENTE no final do programa para salvar tudo
SFTreeWriter::~SFTreeWriter() {
    if (_sf_output_file && _sf_tree) {
        std::cout << "[INFO] Salvando TTree de SFs em " << _sf_output_file->GetName() << std::endl;
        _sf_output_file->cd();
        _sf_tree->Write();
        _sf_output_file->Close();
        delete _sf_output_file;
    }
}

// Guarda o nome da amostra
void SFTreeWriter::setSampleName(const std::string& name) {
    if (_sampleName.empty()) {
        _sampleName = name;
    }
}

// Função para preencher a TTree
void SFTreeWriter::fill(float pt, float eta, float sf_trigger, float sf_reco, float sf_id, float sf_iso, int is_ele) {
    // Inicialização "preguiçosa": cria o arquivo e a TTree na primeira vez que 'fill' é chamado
    if (!_isInitialized) {
        if (_sampleName.empty()) {
            std::cerr << "ERRO: O nome da amostra (sampleName) não foi definido para o SFTreeWriter!" << std::endl;
            _isInitialized = true; // Evita repetir o erro
            return;
        }
        std::string filename = _sampleName + "_sf_ntuple.root";
        _sf_output_file = new TFile(filename.c_str(), "RECREATE");
        _sf_tree = new TTree("sf_tree", "TTree with lepton SFs");

        _sf_tree->Branch("lep_pt", &_lep_pt);
        _sf_tree->Branch("lep_eta", &_lep_eta);
        _sf_tree->Branch("sf_trigger", &_sf_trigger);
        _sf_tree->Branch("sf_reco", &_sf_reco);
        _sf_tree->Branch("sf_id", &_sf_id);
        _sf_tree->Branch("sf_iso", &_sf_iso);
        _sf_tree->Branch("lep_is_ele", &_lep_is_ele);
        
        _isInitialized = true;
    }
    
    // Se a árvore existe, preenche os valores
    if (_sf_tree) {
        _lep_pt = pt;
        _lep_eta = eta;
        _sf_trigger = sf_trigger;
        _sf_reco = sf_reco;
        _sf_id = id;
        _sf_iso = sf_iso;
        _lep_is_ele = is_ele;
        _sf_tree->Fill();
    }
}
