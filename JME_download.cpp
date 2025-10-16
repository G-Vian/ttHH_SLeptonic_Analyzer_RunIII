//// Nome do arquivo: JME_download.cpp
////
//// Compile com:
////    g++ JME_download.cpp -o JME_download
////
//// Execute com:
////    ./JME_download

#include <iostream>
#include <cstdlib> // para system()
#include <string>

int main() {
    // Diretório de destino onde o repositório será clonado
    // Vamos criar uma pasta específica para manter as coisas organizadas
    std::string baseDir = "/afs/cern.ch/user/g/gvian/JME_corrections";
    
    // 1. Cria a pasta principal se ela não existir
    std::string mkdirCmd = "mkdir -p " + baseDir;
    int ret = system(mkdirCmd.c_str());
    if (ret != 0) {
        std::cerr << "[ERRO] Falha ao criar o diretório: " << baseDir << std::endl;
        return 1;
    }

    // 2. Comando git clone com o novo link do repositório JME
    // O repositório será clonado para dentro da pasta JME_corrections/JME
    std::string gitUrl = "https://gitlab.cern.ch/cms-analysis-corrections/JME.git";
    std::string cloneDir = baseDir + "/JME";
    std::string gitCmd = "git clone " + gitUrl + " " + cloneDir;

    std::cout << "[INFO] Clonando o repositório para " << cloneDir << std::endl;
    ret = system(gitCmd.c_str());
    if (ret != 0) {
        std::cerr << "[ERRO] Falha ao clonar o repositório." << std::endl;
        std::cerr << "[INFO] Talvez o repositório já exista em " << cloneDir << std::endl;
        return 1;
    }

    std::cout << "[INFO] Repositório JME clonado com sucesso." << std::endl;
    return 0;
}
