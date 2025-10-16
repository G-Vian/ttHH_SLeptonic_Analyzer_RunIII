#ifndef JET_CORRECTOR_H
#define JET_CORRECTOR_H

#include <string>
#include <vector>
#include <memory>
#include "tnm.h" // Para ter acesso às definições de 'pat::Jet' e 'objectJet'
#include "correction.h"
#include <TLorentzVector.h>

class JetCorrector {
public:
    // Construtor: carrega os arquivos de correção corretos com base no ano e tipo de amostra
    JetCorrector(const std::string& year, const std::string& dataOrMC);

    // Função principal: retorna o quadrivetor corrigido para um dado jato
    // systematic: "nominal", "JES_up:FlavorQCD", "JES_down:Total", etc.
    TLorentzVector getCorrectedP4(const pat::Jet& rawJet, double rho, const std::string& systematic = "nominal") const;

    // Função auxiliar para obter a lista de fontes de incerteza disponíveis
    const std::vector<std::string>& getUncertaintySources() const;


private:
    std::unique_ptr<correction::CorrectionSet> cset_;
    std::vector<std::string> jec_stack_names_;
    std::vector<std::string> uncertainty_sources_;
    bool isMC_;
};

#endif // JET_CORRECTOR_H
