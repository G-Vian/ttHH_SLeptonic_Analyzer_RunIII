#include <iostream>
#include <cstdlib> // para system()
#include <string>
#include <vector>

int main() {
    // Diretório de destino principal
    std::string baseDir = "/afs/cern.ch/user/g/gvian/JME_corrections";
    
    // 1. Cria a pasta principal se ela não existir
    std::string mkdirCmd = "mkdir -p " + baseDir;
    std::cout << "[INFO] Verificando/Criando diretório base: " << baseDir << std::endl;
    int ret = system(mkdirCmd.c_str());
    if (ret != 0) {
        std::cerr << "[ERRO] Falha ao criar o diretório base: " << baseDir << std::endl;
        return 1;
    }

    // 2. Lista de repositórios para clonar
    std::vector<std::string> repoUrls = {
        "https://gitlab.cern.ch/cms-analysis-corrections/JME/Run3-22CDSep23-Summer22-NanoAODv12.git",
        "https://gitlab.cern.ch/cms-analysis-corrections/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15.git",
        "https://gitlab.cern.ch/cms-analysis-corrections/JME/Run3-23DSep23-Summer23BPix-NanoAODv12.git",
        "https://gitlab.cern.ch/cms-analysis-corrections/JME/Run3-23CSep23-Summer23-NanoAODv12.git",
        "https://gitlab.cern.ch/cms-analysis-corrections/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12.git"
    };

    // Nomes das pastas locais para cada repositório
    std::vector<std::string> dirNames = {
        "Run3-22CDSep23-Summer22-NanoAODv12",
        "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
        "Run3-23DSep23-Summer23BPix-NanoAODv12",
        "Run3-23CSep23-Summer23-NanoAODv12",
        "Run3-22EFGSep23-Summer22EE-NanoAODv12"
    };

    // 3. Loop para clonar cada repositório
    for (size_t i = 0; i < repoUrls.size(); ++i) {
        std::string cloneDir = baseDir + "/" + dirNames[i];
        std::string gitCmd = "git clone " + repoUrls[i] + " " + cloneDir;

        std::cout << "\n-----------------------------------------------------" << std::endl;
        std::cout << "[INFO] Clonando repositório " << (i + 1) << "/" << repoUrls.size() << " para " << cloneDir << std::endl;
        
        ret = system(gitCmd.c_str());
        if (ret != 0) {
            std::cerr << "[AVISO] Falha ao clonar " << repoUrls[i] << "." << std::endl;
            std::cerr << "        Isso pode ter acontecido porque o diretório já existe. Verifique a pasta." << std::endl;
        } else {
            std::cout << "[INFO] Repositório clonado com sucesso." << std::endl;
        }
    }

    std::cout << "\n-----------------------------------------------------" << std::endl;
    std::cout << "[INFO] Processo de download finalizado." << std::endl;
    
    return 0;
}
