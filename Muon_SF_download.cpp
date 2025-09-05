////compile : 
////   g++ Muon_SF_download.cpp -o Muon_SF_download -std=c++17
////execute 
////   ./Muon_SF_download


#include <iostream>
#include <filesystem>
#include <cstdlib>
#include <vector>
#include <string>

namespace fs = std::filesystem;

int main() {
    std::string targetDir = "/afs/cern.ch/user/g/gvian/muon_SF/muon_trigger_SF";

    // Cria a pasta se não existir
    if (!fs::exists(targetDir)) {
        if (!fs::create_directories(targetDir)) {
            std::cerr << "[ERROR] Failed to create directory: " << targetDir << std::endl;
            return 1;
        }
    }

    std::vector<std::string> urls = {
        "https://github.com/cms-muon/MuonEfficiencies/blob/main/Run3/2022/2022_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2022_eta_pt_schemaV2.json?raw=true",
        "https://github.com/cms-muon/MuonEfficiencies/blob/main/Run3/2022_EE/2022_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2022_EE_eta_pt_schemaV2.json?raw=true",
        "https://github.com/cms-muon/MuonEfficiencies/blob/main/Run3/2023/2023_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2023_eta_pt_schemaV2.json?raw=true",
        "https://github.com/cms-muon/MuonEfficiencies/blob/main/Run3/2023_BPix/2023_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2023_BPix_eta_pt_schemaV2.json?raw=true",
        "https://github.com/cms-muon/MuonEfficiencies/blob/main/Run3/2024/2024_Z/HLT/json/ScaleFactors_Muon_Z_HLT_2024_eta_pt_schemaV2.json?raw=true"
    };

    for (const auto& url : urls) {
        // Extrai o nome do arquivo
        auto pos = url.find_last_of('/');
        std::string fileName = url.substr(pos + 1);
        std::string fullPath = targetDir + "/" + fileName;

        // Baixa o arquivo usando wget
        std::string cmd = "wget -q -O " + fullPath + " " + url;
        std::cout << "[INFO] Downloading " << fileName << "..." << std::endl;
        int ret = system(cmd.c_str());
        if (ret != 0) {
            std::cerr << "[ERROR] Failed to download: " << url << std::endl;
        } else {
            std::cout << "[OK] Saved to: " << fullPath << std::endl;
        }
    }

    std::cout << "[INFO] All files downloaded to " << targetDir << std::endl;
    return 0;
}
