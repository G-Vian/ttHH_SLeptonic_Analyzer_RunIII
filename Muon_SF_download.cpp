////compile : 
////   g++ Muon_SF_download.cpp -o Muon_SF_download
////execute 
////   ./Muon_SF_download

#include <iostream>
#include <cstdlib> // para system()
#include <string>

int main() {
    // Diretório de destino
    std::string baseDir = "/afs/cern.ch/user/g/gvian/muon_SF";
    
    // Cria a pasta principal se não existir
    std::string mkdirCmd = "mkdir -p " + baseDir;
    int ret = system(mkdirCmd.c_str());
    if (ret != 0) {
        std::cerr << "[ERROR] Failed to create directory: " << baseDir << std::endl;
        return 1;
    }

    // Comando git clone
    std::string gitCmd = "git clone https://gitlab.cern.ch/cms-muonPOG/muonefficiencies.git " + baseDir + "/muonefficiencies";

    std::cout << "[INFO] Cloning repository into " << baseDir << std::endl;
    ret = system(gitCmd.c_str());
    if (ret != 0) {
        std::cerr << "[ERROR] Failed to clone repository." << std::endl;
        return 1;
    }

    std::cout << "[INFO] Repository successfully cloned." << std::endl;
    return 0;
}
