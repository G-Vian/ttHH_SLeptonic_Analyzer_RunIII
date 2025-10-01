// Em src/SFTreeWriter.h

#ifndef SFTreeWriter_h
#define SFTreeWriter_h

#include <string>

// Forward declarations para não precisar incluir os headers do ROOT aqui
class TFile;
class TTree;

class SFTreeWriter {
public:
    // Método para obter a única instância da classe
    static SFTreeWriter& getInstance();

    // Define o nome da amostra para o arquivo de saída
    void setSampleName(const std::string& name);

    // Preenche a TTree com os dados de um lépton
    void fill(float pt, float eta, float sf_trigger, float sf_reco, float sf_id, float sf_iso, int is_ele);

    // Deleta o construtor de cópia e o operador de atribuição para garantir que seja único
    SFTreeWriter(const SFTreeWriter&) = delete;
    SFTreeWriter& operator=(const SFTreeWriter&) = delete;

private:
    // Construtor e destrutor privados
    SFTreeWriter();
    ~SFTreeWriter();

    // Variáveis membro
    TFile* _sf_output_file;
    TTree* _sf_tree;
    std::string _sampleName;
    bool _isInitialized;

    // Variáveis para os branches
    float _lep_pt;
    float _lep_eta;
    float _sf_trigger;
    float _sf_reco;
    float _sf_id;
    float _sf_iso;
    int   _lep_is_ele;
};

#endif // SFTreeWriter_h
