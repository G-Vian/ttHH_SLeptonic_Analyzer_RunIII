#ifndef ttHHanalyzer_h
#define ttHHanalyzer_h
#include "tnm.h"
#include <cmath> 
#include <algorithm>
#include <TString.h>
#include <vector>
#include "TVector3.h"
#include "TH1F.h"
#include "TH2.h"//Trigger SF for electron (TSFel)
#include "TFile.h"
#include <iterator>
#include <string>  
//#include "EventShape/Class/src/EventShape.cc"
#include "EventShape/Class/interface/EventShape.h"
#include <TLorentzVector.h>
#include "TRandom3.h"
#include <unordered_map>
//#include "thhHypothesisCombinatorics.h"
//#include "HypothesisCombinatorics.h"
#include "include/tthHypothesisCombinatorics.h"
#include "include/HypothesisCombinatorics.h"
#include "fifo_map.hpp"

#include "json.hpp" // muon trigger SF (TSFmu)
using json = nlohmann::json;  /// this is for MUON trigger SF (TSFmu)
extern json muonTrigSFJson;

//using namespace ROOT::Math;
using nlohmann::fifo_map;

const float cLargeValue = 99999999999.;
const float cEps = 0.000000001; 
const float cHiggsMass = 125.38;
const float cZMass = 91.;

extern std::map<std::string, float> cut;







class objectPhysics { 
 public: 
    enum lFlavor{kNA, kEle, kMuon}; 
    explicit objectPhysics(const float pT, const float eta, const float phi, const float mass = 0){ 
	_p4.SetPtEtaPhiM(pT, eta, phi, mass); 
    }
    
    TLorentzVector * getp4(){
	return &_p4;
    }
    objectPhysics(){};
    void scale(float JES, bool up = true){
	_pxOffset = JES * _p4.Px();
	_pyOffset = JES * _p4.Py();
	_pzOffset = JES * _p4.Pz();
	_EOffset  = JES * _p4.E();
	if(up){
	    _p4.SetPxPyPzE(_p4.Px()+_pxOffset,_p4.Py()+_pyOffset,_p4.Pz()+_pzOffset, _p4.E()+_EOffset);
	} else {
	    _p4.SetPxPyPzE(_p4.Px()-_pxOffset,_p4.Py()-_pyOffset,_p4.Pz()-_pzOffset, _p4.E()-_EOffset);
	}
    }
    std::vector<float> getOffset(){
	std::vector<float> offset = {_pxOffset,_pyOffset,_pzOffset,_EOffset};
	return offset;
    }
    void subtractp4(const std::vector<float>& offset){
	_p4.SetPxPyPzE(_p4.Px()-offset[0],_p4.Py()-offset[1],_p4.Pz()-offset[2], _p4.E()-offset[3]);
    }
    void addp4(const std::vector<float>& offset){
	_p4.SetPxPyPzE(_p4.Px()+offset[0],_p4.Py()+offset[1],_p4.Pz()+offset[2], _p4.E()+offset[3]);
    }

//    void setYear(const int year){
//	_year = year;
 //   }

 private:
    TLorentzVector _p4;
    float _pxOffset = 0., _pyOffset = 0., _pzOffset = 0., _EOffset = 0.;
    std::string _year;
};


class objectGenPart:public objectPhysics { 
 public:    
    using objectPhysics::objectPhysics;
    bool hasHiggsMother = false;
    bool hasTopMother = false;
    bool matched = false; 
    float dRmatched = 10000.;
};


class objectJet:public objectPhysics {
 public:
    
    using objectPhysics::objectPhysics;
    bool matchedtoHiggs = false;
    float minChiHiggs = 0.;
    int minChiHiggsIndex = 0;
    float matchedtoHiggsdR = 0.;
    //int bTag;
    float bTagCSV, jetID, jetPUid;
    static constexpr float valbTagTight2017  = 0.7476;
    static constexpr float valbTagMedium2017 = 0.3040;
    static constexpr float valbTagLoose2017  = 0.0532;
    static constexpr float valbTagTight2018  = 0.7100; 
    static constexpr float valbTagMedium2018 = 0.2783; 
    static constexpr float valbTagLoose2018  = 0.0490; 
//robustParticleTransformer
/*
    static constexpr float valbTagTight2022  = 0.8482; //22
    static constexpr float valbTagMedium2022 = 0.4319; 
    static constexpr float valbTagLoose2022  = 0.0849; 
    static constexpr float valbTagTight2022E  = 0.8604; //22EE
    static constexpr float valbTagMedium2022E = 0.451; 
    static constexpr float valbTagLoose2022E  = 0.0897; 
    static constexpr float valbTagTight2023  = 0.7969; //23
    static constexpr float valbTagMedium2023 = 0.3487; 
    static constexpr float valbTagLoose2023  = 0.0681; 
    static constexpr float valbTagTight2023B  = 0.7994; //23Bpix
    static constexpr float valbTagMedium2023B = 0.3494; 
    static constexpr float valbTagLoose2023B  = 0.0683; 
    static constexpr float valbTagTight2024  = 0.8482; //24 (wrong values, update them)
    static constexpr float valbTagMedium2024 = 0.4319; 
    static constexpr float valbTagLoose2024  = 0.0849; 
    float valbTagTight;
    float valbTagMedium;
    float valbTagLoose;
    */

//PaticleNet
/*
    static constexpr float valbTagTight2022  =0.6734; //22
    static constexpr float valbTagMedium2022 = 0.245; 
    static constexpr float valbTagLoose2022  = 0.047; 
    static constexpr float valbTagTight2022E  = 0.6915; //22EE
    static constexpr float valbTagMedium2022E = 0.2605; 
    static constexpr float valbTagLoose2022E  = 0.0499; 
    static constexpr float valbTagTight2023  = 0.6172; //23
    static constexpr float valbTagMedium2023 =0.1917; 
    static constexpr float valbTagLoose2023  = 	0.0358; 
    static constexpr float valbTagTight2023B  = 0.6133; //23Bpix
    static constexpr float valbTagMedium2023B = 0.1919; 
    static constexpr float valbTagLoose2023B  = 0.0359; 
    static constexpr float valbTagTight2024  = 0.6133; //24 (wrong values, update them)
    static constexpr float valbTagMedium2024 =  0.1919; 
    static constexpr float valbTagLoose2024  =  0.0359; 
    float valbTagTight;
    float valbTagMedium;
    float valbTagLoose;
    
float getValbTagTight(const std::string& year) {
    if (year == "2017") {
        return valbTagTight2017;
    } else if (year == "2018") {
        return valbTagTight2018;
    } else if (year == "2022EE") {
        return valbTagTight2022E;
    } else if (year == "2022") {
        return valbTagTight2022;
    } else if (year == "2023B") {
        return valbTagTight2023B;
    } else if (year == "2023") {
        return valbTagTight2023;
    }  else if (year == "2024") {
        return valbTagTight2024;
    }
    return -1.0f;  // Valor padrão se o ano não for reconhecido
}

float getValbTagMedium(const std::string& year) {
    if (year == "2017") {
        return valbTagMedium2017;
    } else if (year == "2018") {
        return valbTagMedium2018;
    } else if (year == "2022EE") {
        return valbTagMedium2022E;
    } else if (year == "2022") {
        return valbTagMedium2022;
    } else if (year == "2023B") {
        return valbTagMedium2023B;
    } else if (year == "2023") {
        return valbTagMedium2023;
    } else if (year == "2024") {
        return valbTagMedium2024;
    }
    return -1.0f;  // Valor padrão se o ano não for reconhecido
}

float getValbTagLoose(const std::string& year) {
    if (year == "2017") {
        return valbTagLoose2017;
    } else if (year == "2018") {
        return valbTagLoose2018;
    } else if (year == "2022EE") {
        return valbTagLoose2022E;
    } else if (year == "2022") {
        return valbTagLoose2022;
    } else if (year == "2023B") {
        return valbTagLoose2023B;
    } else if (year == "2023") {
        return valbTagLoose2023;
    }else if (year == "2024") {
        return valbTagLoose2024;
    }
    return -1.0f;  // Valor padrão se o ano não for reconhecido
}

*/
//UParTAK4 --> 22 and 23 WP not available yet! I'm using PNet Values 

    static constexpr float valbTagTight2022  =0.6734; //22
    static constexpr float valbTagMedium2022 = 0.245; 
    static constexpr float valbTagLoose2022  = 0.047; 
    static constexpr float valbTagTight2022E  = 0.6915; //22EE
    static constexpr float valbTagMedium2022E = 0.2605; 
    static constexpr float valbTagLoose2022E  = 0.0499; 
    static constexpr float valbTagTight2023  = 0.6172; //23
    static constexpr float valbTagMedium2023 =0.1917; 
    static constexpr float valbTagLoose2023  = 	0.0358; 
    static constexpr float valbTagTight2023B  = 0.6133; //23Bpix
    static constexpr float valbTagMedium2023B = 0.1919; 
    static constexpr float valbTagLoose2023B  = 0.0359; 
    static constexpr float valbTagTight2024  = 0.4648; //24 (Correct values)
    static constexpr float valbTagMedium2024 =  0.1272; 
    static constexpr float valbTagLoose2024  = 0.0246; 
    float valbTagTight;
    float valbTagMedium;
    float valbTagLoose;
    
float getValbTagTight(const std::string& year) {
    if (year == "2017") {
        return valbTagTight2017;
    } else if (year == "2018") {
        return valbTagTight2018;
    } else if (year == "2022EE") {
        return valbTagTight2022E;
    } else if (year == "2022") {
        return valbTagTight2022;
    } else if (year == "2023B") {
        return valbTagTight2023B;
    } else if (year == "2023") {
        return valbTagTight2023;
    }  else if (year == "2024") {
        return valbTagTight2024;
    }
    return -1.0f;  // Valor padrão se o ano não for reconhecido
}

float getValbTagMedium(const std::string& year) {
    if (year == "2017") {
        return valbTagMedium2017;
    } else if (year == "2018") {
        return valbTagMedium2018;
    } else if (year == "2022EE") {
        return valbTagMedium2022E;
    } else if (year == "2022") {
        return valbTagMedium2022;
    } else if (year == "2023B") {
        return valbTagMedium2023B;
    } else if (year == "2023") {
        return valbTagMedium2023;
    } else if (year == "2024") {
        return valbTagMedium2024;
    }
    return -1.0f;  // Valor padrão se o ano não for reconhecido
}

float getValbTagLoose(const std::string& year) {
    if (year == "2017") {
        return valbTagLoose2017;
    } else if (year == "2018") {
        return valbTagLoose2018;
    } else if (year == "2022EE") {
        return valbTagLoose2022E;
    } else if (year == "2022") {
        return valbTagLoose2022;
    } else if (year == "2023B") {
        return valbTagLoose2023B;
    } else if (year == "2023") {
        return valbTagLoose2023;
    }else if (year == "2024") {
        return valbTagLoose2024;
    }
    return -1.0f;  // Valor padrão se o ano não for reconhecido
}


};

class objectMET:public objectPhysics {
 public:
    using objectPhysics::objectPhysics;
};

class objectBoostedJet:public objectPhysics {
 public:
    using objectPhysics::objectPhysics;
    float softDropMass;
};


class objectbJet:public objectPhysics {
 public:
    using objectPhysics::objectPhysics;
};

class objectLightJet:public objectPhysics {
 public:
    using objectPhysics::objectPhysics;
    };

class objectLep:public objectPhysics {
 public:
    using objectPhysics::objectPhysics;
    int charge;
    float miniPFRelIso;
    float pfRelIso03;
    float pfRelIso04;
    lFlavor flavor;
};

class event{
 public:
    event(){

    }

    struct evShapes{
	float objectP;
	float objectPT;
	float centrality;
    };

    struct maxObjects{
	float maxPT;
	float maxPTmass;
    };

    struct statObjects{
	float dR;
	float meandR;
	float mindR;
	float maxdR;
	float meandEta;
	float mindEta;
	float maxdEta;
	float meandPhi;
	float mindPhi;
	float maxdPhi;
	float mindRMass;
	float mindRpT;
    };

    struct foxWolframObjects{
	float h0;
	float h1;
	float h2;
	float h3;
	float h4;
	float r1;
	float r2;
	float r3;
	float r4;
    };

    EventShape * eventShapeJet, * eventShapeBjet;

    
    void addJet(objectJet * jet){
	_sumJetScalarpT+=fabs(jet->getp4()->Pt());
	_sumJetp4+= * jet->getp4();
	_jets.push_back(jet);
    }
    
    void selectJet(objectJet * jet){
	_sumSelJetScalarpT+=fabs(jet->getp4()->Pt());
	_sumSelJetMass+=jet->getp4()->M(); //
	_sumSelJetp4 += * jet->getp4();
	_selectJets.push_back(jet);
    }
    
    void selectLepton(objectLep * lepton){
	//	_sumLeptonScalarpT+=fabs(lepton->getp4()->Pt());  
	_selectLeptons.push_back(lepton);
    }

    void selectEle(objectLep * ele){
	ele->flavor = objectLep::kEle;
	_sumSelElectronScalarpT+=fabs(ele->getp4()->Pt());  
	_sumSelElectronp4 += *ele->getp4();  
	_selectLeptons.push_back(ele);
	_selectElectrons.push_back(ele);
    }

    void selectMuon(objectLep * muon){
	muon->flavor = objectLep::kMuon;
	_sumSelMuonScalarpT+=fabs(muon->getp4()->Pt());  
	_sumSelMuonp4 += *muon->getp4();  
	_selectLeptons.push_back(muon);
	_selectMuons.push_back(muon);
    }
    

    void selectHadronicHiggs(objectBoostedJet * boostedJet){
	//	_sumSelBoostedJetSoftDropMass+=fabs(boostedJet->softDropMass);
	_sumSelHadronicHiggsScalarpT+=fabs(boostedJet->getp4()->Pt());  
	_sumSelHadronicHiggsMass+=boostedJet->getp4()->M();  
	_sumHadronicHiggsp4 += *boostedJet->getp4();
	_selectHadronicHiggses.push_back(boostedJet);
    }

    void selectBoostedJet(objectBoostedJet * boostedJet){
	//	_sumSelBoostedJetSoftDropMass+=fabs(boostedJet->softDropMass);
	//	_sumSelHadronicHiggsScalarpT+=fabs(boostedJet->getp4()->Pt());  
	//_sumSelHadronicHiggsScalarMass+=fabs(boostedJet->getp4()->M());  
	//_sumHadronicHiggsp4 += *boostedJet->getp4();
	_selectBoostedJets.push_back(boostedJet);
    }


    void selectbJet(objectJet * jet){
	_sumSelbJetScalarpT+=fabs(jet->getp4()->Pt());  
	_sumSelbJetMass+=jet->getp4()->M();  
	_sumSelbJetp4 += *jet->getp4();
	_selectbJets.push_back(jet);
    }

    void selectLightJet(objectJet * jet){
	_sumSelLightJetScalarpT+=fabs(jet->getp4()->Pt());
	_sumSelLightJetMass+=jet->getp4()->M(); 
	_sumLightJetp4 += *jet->getp4();
	_selectLightJets.push_back(jet);
    }

    void selectLoosebJet(objectJet * jet){
	_loosebJets.push_back(jet);
    }

    void selectGenPart(objectGenPart * genPart){   
	_selectGenParts.push_back(genPart);
    } 


    void setMET(objectMET * met){
	_MET = met;
    }

    void setnVetoLepton(int nVeto){
	_nVetoLepton = nVeto;
    }

    float getSumJetScalarpT(){
	return _sumJetScalarpT;
    }

    float getSumSelJetScalarpT(){
	return _sumSelJetScalarpT;
    }

    float getSumSelJetpT(){
	return _sumSelJetp4.Pt();
    }

    float getSumSelHadronicHiggsScalarpT(){   
	return _sumSelHadronicHiggsScalarpT;
    }

    float getSumSelbJetScalarpT(){   
	return _sumSelbJetScalarpT;
    }

    float getSumSelbJetpT(){   
	return _sumSelbJetp4.Pt();
    }

    float getSumSelLightJetScalarpT(){   
	return _sumSelLightJetScalarpT;
    }

    float getSumSelLeptonScalarpT(){   
	return _sumSelMuonScalarpT + _sumSelElectronScalarpT;
    }

    float getSumSelLeptonpT(){   
	return _sumSelMuonp4.Pt() + _sumSelElectronp4.Pt();
    }

    float getSelLeptonHT(){   
	return _sumSelJetScalarpT + _sumSelMuonScalarpT + _sumSelElectronScalarpT;
    }

    float getSelLeptonST(){   
	return _sumSelJetScalarpT + _sumSelMuonScalarpT + _sumSelElectronScalarpT + getMET()->getp4()->Pt();
    }

    float getSumSelJetMass(){ 
	return _sumSelJetMass;
    }
    
    float getSumSelHadronicHiggsMass(){  
	return _sumSelHadronicHiggsMass;
    }

    float getSumSelbJetMass(){  
	return _sumSelbJetMass;
    }

    float getSumSelLightJetMass(){  
    	return _sumSelLightJetMass;
    }

    //    float getSelHadronicHiggsSoftDropMass(){
    //	return _selectHadronicHiggses.msoftdrop();
    // }


    float getSelLeptonsMass(){  
	return _sumSelMuonp4.M() + _sumSelElectronp4.M();
    }

    float getSelMuonsMass(){  
	return _sumSelMuonp4.M();
    }

    float getSelMuonsPT(){  
	return _sumSelMuonp4.Pt();
    }

    float getSelMuonsEta(){  
	return _sumSelMuonp4.Eta();
    }

    float getSelElectronsMass(){  
	return _sumSelElectronp4.M();
    } 

    float getSelElectronsPT(){  
	return _sumSelElectronp4.Pt();
    } 

    float getSelElectronsEta(){  
	return _sumSelElectronp4.Eta();
    } 

    int getnGenPart(){ 
	return _selectGenParts.size(); 
    }

    int getnbJet(){
	return _selectbJets.size(); 
    }

    int getnLightJet(){
    	return _selectLightJets.size(); 
    }

    int getnbLooseJet(){
	return _loosebJets.size(); 
    }

    int getnHadronicHiggs(){
	return _selectHadronicHiggses.size();
    }

    int getnSelJet(){
	return _selectJets.size(); 
    }

    int getnJet(){
	return _jets.size(); 
    }

    int getnSelElectron(){
	return _selectElectrons.size();
    }

    int getnSelMuon(){
	return _selectMuons.size(); 
    }

    bool orderLeptons(){
	if(_selectLeptons.size() != 2)
	    return false;
	else {
	    objectLep* tmp; 
	    if(_selectLeptons.at(0)->getp4()->Pt() < _selectLeptons.at(1)->getp4()->Pt()){
		//std::cout << "before " << _selectLeptons.at(0)->getp4()->Pt() << " " << _selectLeptons.at(1)->getp4()->Pt();
		tmp = _selectLeptons.at(0);
		_selectLeptons.at(0) = _selectLeptons.at(1);
		_selectLeptons.at(1) = tmp;
		//std::cout << "after " <<  _selectLeptons.at(0)->getp4()->Pt() << " " << _selectLeptons.at(1)->getp4()->Pt();
	    }
	    return true;
	}
    }

    int getnSelLepton(){
	return _selectLeptons.size();
    }

    int getnVetoLepton(){
	return _nVetoLepton;
    }

    objectMET* getMET(){
	return _MET;
    }

    std::vector<objectGenPart*>* getGenParts(){ 
	return &_selectGenParts;
    }

    std::vector<objectJet*>* getJets(){
	return &_jets;
    }


    std::vector<objectJet*>* getSelJets(){
	return &_selectJets;
    }

    std::vector<objectBoostedJet*>* getSelHadronicHiggses(){
	return &_selectHadronicHiggses;
    }

    std::vector<objectBoostedJet*>* getSelBoostedJets(){
	return &_selectBoostedJets;
    }

    std::vector<objectJet*>* getSelbJets(){
	return &_selectbJets;
    }

    std::vector<objectJet*>* getSelLightJets(){
    	return &_selectLightJets;
    }

    std::vector<objectJet*>* getLoosebJets(){
	return &_loosebJets;
    }

    std::vector<objectLep*>* getSelElectrons(){
	return &_selectElectrons;
    }

    std::vector<objectLep*>* getSelMuons(){
	return &_selectMuons;
    }

    std::vector<objectLep*>* getSelLeptons(){
	return &_selectLeptons;
    }

    // Centrality calculation //
    template <class object1, class object2>
	void getCentrality(std::vector<object1*>* cont1, std::vector<object2*>* cont2, evShapes& cent) {
	TLorentzVector obj1P4, obj2P4;
	float sumPT = 0., sumP = 0., centrality = 0.;

	if(cont1->size()  == 0 || cont2->size() == 0){
	    std::cout << "WTF!!!!" << std::endl;
	}
	for(int oindex=0; oindex < cont1->size(); oindex++){
	    obj1P4 = (*cont1->at(oindex)->getp4());
	    sumPT += obj1P4.Pt();
	    sumP  += obj1P4.P();
	}
	for(int iindex = 0; iindex < cont2->size(); iindex++){
	    obj2P4 = (*cont2->at(iindex)->getp4());
	    sumPT += obj2P4.Pt();
	    sumP  += obj2P4.P();
	}
	centrality = sumPT/sumP;
	cent.objectPT = sumPT;
	cent.objectP  = sumP;
	cent.centrality = centrality;
    }


    // Centrality calculation //                                                                                                                                                
    template <class object1, class object2>
	void getCentralityV2(std::vector<object1*>* cont1, std::vector<object2*>* cont2, evShapes& cent) {
        TLorentzVector obj1P4, obj2P4, sumP4;
        float  centrality = 0.;

        int nObject = 0.;
	if(cont1->size()  == 0 || cont2->size() == 0){
	    std::cout << "WTF!!!!" << std::endl;
        }
        for(int oindex=0; oindex < cont1->size(); oindex++){
            obj1P4 = (*cont1->at(oindex)->getp4());
            sumP4 += obj1P4;
        }
        for(int iindex = 0; iindex < cont2->size(); iindex++){
            obj2P4 = (*cont2->at(iindex)->getp4());
            sumP4 += obj2P4;
        }
        centrality = sumP4.Pt()/sumP4.P();
        cent.objectPT = sumP4.Pt();
        cent.objectP  = sumP4.P();
	cent.centrality = centrality;
    }


    template <class object1>
	void getMaxPTSame(std::vector<object1*>* cont1,  maxObjects& xxxMaxs) {
	float maxPT = 0., maxPTmass = 0., tmpPT = 0., tmpMass;
	TLorentzVector tmpP4; 
	int nObject = 0.;
	if(cont1->size()  == 0){
	    std::cout << "WTF!!!!" << std::endl;
	}
	for(int oindex=0; oindex < cont1->size(); oindex++){
	    for(int iindex = oindex+1; iindex < cont1->size(); iindex++){
		for(int mindex = oindex+2; mindex < cont1->size(); mindex++){
		    tmpP4 = (*cont1->at(oindex)->getp4() + *cont1->at(iindex)->getp4() + *cont1->at(mindex)->getp4());
		    tmpPT = tmpP4.Pt();
			//(*cont1->at(oindex)->getp4() + *cont1->at(iindex)->getp4() + *cont1->at(mindex)->getp4()).Pt();
		    tmpMass = tmpP4.M();
			//(*cont1->at(oindex)->getp4() + *cont1->at(iindex)->getp4() + *cont1->at(mindex)->getp4()).M();
		    if(maxPT < tmpPT){
			maxPT= tmpPT;
			maxPTmass = tmpMass;
		    }
		    nObject++;
		}
	    }
	}
	xxxMaxs.maxPT = maxPT;
	xxxMaxs.maxPTmass = maxPTmass;
    }


    template <class object1, class object2>
	void getMaxPTComb(std::vector<object1*>* cont1, std::vector<object2*>* cont2, maxObjects& xyyMaxs) {
	float maxPT = 0., maxPTmass = 0., tmpPT = 0., tmpMass;
	TLorentzVector tmpP4; 
	int nObject = 0.;
	if(cont1->size()  == 0 || cont2->size() == 0){
	    std::cout << "WTF!!!!" << std::endl;
	}
	for(int oindex=0; oindex < cont1->size(); oindex++){
	    for(int iindex = 0; iindex < cont2->size(); iindex++){
		for(int mindex = iindex+1; mindex < cont2->size(); mindex++){
		    tmpP4 = (*cont1->at(oindex)->getp4() + *cont2->at(iindex)->getp4() + *cont2->at(mindex)->getp4());
		    tmpPT = tmpP4.Pt();
		    tmpMass = tmpP4.M();
		    //  tmpPT = (*cont1->at(oindex)->getp4() + *cont2->at(iindex)->getp4() + *cont2->at(mindex)->getp4()).Pt();
		    //  tmpMass = (*cont1->at(oindex)->getp4() + *cont2->at(iindex)->getp4() + *cont2->at(mindex)->getp4()).M();
		    if(maxPT < tmpPT){
			maxPT= tmpPT;
			maxPTmass = tmpMass;
		    }
		    nObject++;
		}
	    }
	}
	xyyMaxs.maxPT = maxPT;
	xyyMaxs.maxPTmass = maxPTmass;
    } 

    template <class object1, class object2>
	void getStatsComb (std::vector<object1*>* cont1, std::vector<object2*>* cont2, statObjects& statsComb) {
	float mindR = 99999999999., mindEta = 99999999999., mindPhi = 99999999999., tmpdR = 0., tmpdPhi = 0., tmpdEta = 0., tmpMass = 0., tmpPT = 0, sumdEta = 0., sumdPhi = 0., sumdR = 0., maxdR = 0. , maxdPhi = 0., maxdEta = 0., mindRMass = 0., mindRpT = 0.;
	int nObject = 0;

	if(cont1->size()  == 0 || cont2->size() == 0){
	    std::cout << "WTF!!!!" << std::endl;
	}
	for(int oindex=0; oindex < cont1->size(); oindex++){
	    for(int iindex = 0; iindex < cont2->size(); iindex++){

		tmpdPhi = fabs(cont1->at(oindex)->getp4()->Phi() - cont2->at(iindex)->getp4()->Phi());
		tmpdEta = fabs(cont1->at(oindex)->getp4()->Eta() - cont2->at(iindex)->getp4()->Eta());
		tmpdR   = TMath::Sqrt(tmpdPhi*tmpdPhi + tmpdEta*tmpdEta);
		tmpMass = (*cont1->at(oindex)->getp4() + *cont2->at(iindex)->getp4()).M(); 
		tmpPT   = (*cont1->at(oindex)->getp4() + *cont2->at(iindex)->getp4()).Pt(); 

		sumdR   += tmpdR;
		sumdEta += tmpdEta;
		sumdPhi += tmpdPhi;
		if(mindR > tmpdR){
		    mindR = tmpdR;
		    mindRMass = tmpMass;
 		    mindRpT = tmpPT;
		}
		if(maxdR < tmpdR){
		    maxdR = tmpdR; 
		} 
		if(mindEta > tmpdEta){
		    mindEta = tmpdEta;
		}
		if(maxdEta < tmpdEta){
		    maxdEta = tmpdEta; 
		} 
		if(mindPhi > tmpdPhi){
		    mindPhi = tmpdPhi;
		}
		if(maxdPhi < tmpdPhi){
		    maxdPhi = tmpdPhi; 
		} 
		nObject++;
	    }
	}
	statsComb.dR        = tmpdR;
	statsComb.meandR    = sumdR   / (float) nObject;
	statsComb.meandEta  = sumdEta / (float) nObject;
	statsComb.meandPhi  = sumdPhi / (float) nObject;
	statsComb.mindR     = mindR;
	statsComb.mindEta   = mindEta;
	statsComb.mindPhi   = mindPhi;
	statsComb.maxdR     = maxdR;
	statsComb.maxdEta   = maxdEta;
	statsComb.maxdPhi   = maxdPhi;
	statsComb.mindRpT   = mindRpT;
	statsComb.mindRMass = mindRMass;
    }


    ///Starting to calculate Fox Wolfram moments///
    template <class object>
	void getFoxWolfram (std::vector<object*>* cont, foxWolframObjects& foxwolf){ //, double &h0, double &h1, double &h2, double &h3, double &h4){
	double jetEnergy = 0.0;
	//double costh;
	double h0 = 0.0, h1 = 0.0, h2 = 0.0, h3 = 0.0, h4 = 0.0;	
	double r1 = 0.0, r2 = 0.0, r3 = 0.0, r4 = 0.0;
	//	std::cout << "debug: " << cont->size() << std::endl;

	for (int oindex = 0; oindex < cont->size(); oindex++) {
	    jetEnergy += cont->at(oindex)->getp4()->E();
        }

	for(int oindex = 0; oindex < cont->size()-1; oindex++){
	    for(int iindex = oindex+1; iindex < cont->size(); iindex++){
		double costh = cos(cont->at(oindex)->getp4()->Angle(cont->at(iindex)->getp4()->Vect()));
		double p0 = 1.0;
		double p1 = costh;
		double p2 = 0.5*(3.0*costh*costh - 1.0);
		double p3 = 0.5*(5.0*costh*costh*costh - 3.0*costh);
		double p4 = 0.125*(35.0*costh*costh*costh*costh - 30.0*costh*costh + 3.0);
		double pipj = cont->at(oindex)->getp4()->P() * cont->at(iindex)->getp4()->P();
		h0 += (pipj/(jetEnergy*jetEnergy))*p0;
		h1 += (pipj/(jetEnergy*jetEnergy))*p1;
		h2 += (pipj/(jetEnergy*jetEnergy))*p2;
		h3 += (pipj/(jetEnergy*jetEnergy))*p3;
		h4 += (pipj/(jetEnergy*jetEnergy))*p4;
	    }
	}

	r1 = h1/h0;
	r2 = h2/h0;
	r3 = h3/h0;
	r4 = h4/h0;

	foxwolf.h0 = h0;
	foxwolf.h1 = h1;
	foxwolf.h2 = h2;
	foxwolf.h3 = h3;
	foxwolf.h4 = h4;
	foxwolf.r1 = r1;
	foxwolf.r2 = r2;
	foxwolf.r3 = r3;
	foxwolf.r4 = r4;
    }


    template <class object>
	void getStats (std::vector<object*>* cont, statObjects& stats){
	float mindR = 9999999999, mindEta = 9999999999, mindPhi = 99999999999., tmpdR = 0., tmpdPhi = 0., tmpdEta = 0., tmpMass = 0., tmpPT =0., sumdEta = 0., sumdPhi = 0., sumdR = 0., maxdR = 0. , maxdPhi = 0., maxdEta = 0., mindRMass = 0., mindRpT = 0.;
	int iindex = 0, nObject = 0;
	for(int oindex=0; oindex < cont->size(); oindex++){
	    for(iindex = oindex+1; iindex < cont->size(); iindex++){
		tmpdPhi = fabs(cont->at(iindex)->getp4()->Phi() - cont->at(oindex)->getp4()->Phi());
		tmpdEta = fabs(cont->at(iindex)->getp4()->Eta() - cont->at(oindex)->getp4()->Eta());
		tmpdR = TMath::Sqrt(tmpdPhi*tmpdPhi + tmpdEta*tmpdEta);
		tmpMass = (*cont->at(oindex)->getp4() + *cont->at(iindex)->getp4()).M(); 
		tmpPT   = (*cont->at(oindex)->getp4() + *cont->at(iindex)->getp4()).Pt(); 
		
		sumdR   += tmpdR;
		sumdEta += tmpdEta;
		sumdPhi += tmpdPhi;
		if(mindR > tmpdR){
		    mindR = tmpdR;
		    mindRMass = tmpMass;
		    mindRpT = tmpPT;
		}
		if(maxdR < tmpdR){
		    maxdR = tmpdR; 
		} 
		if(mindEta > tmpdEta){
		    mindEta = tmpdEta;
		}
		if(maxdEta < tmpdEta){
		    maxdEta = tmpdEta; 
		} 
		if(mindPhi > tmpdPhi){
		    mindPhi = tmpdPhi;
		} 
		if(maxdPhi < tmpdPhi){
		    maxdPhi = tmpdPhi; 
		} 
		nObject++;
	    }
	}
	stats.meandR    = sumdR   / (float) nObject;
	stats.meandEta  = sumdEta / (float) nObject;
	stats.meandPhi  = sumdPhi / (float) nObject;
	stats.mindR     = mindR;
	stats.mindEta   = mindEta;
	stats.mindPhi   = mindPhi;
	stats.maxdR     = maxdR;
	stats.maxdEta   = maxdEta;
	stats.maxdPhi   = maxdPhi;
	stats.mindRpT   = mindRpT;
	stats.mindRMass = mindRMass;
    }

    void summarize(){

	//	std::cout << "nJet: " << getnJet() << std::endl;//" jet scalar sum: " << _sumJetScalarpT << std::endl;
	//std::cout << "nSelectedJet: " << getnSelJet() << " jet selected scalar sum: " << _sumSelJetScalarpT << std::endl;
	//std::cout << "nbJet: " << getnbJet() << std::endl;
    //	std::cout << "nElectron: "<< getnElectron() << " nMuon: " << getnMuon() << " nLepton: " << getnLepton() << std::endl;
    //	statObjects jetStat;
    //	getStats(getSelJets(), jetStat);
	//std::cout << "min Jet dr: " << jetStat.mindR << std::endl;
    }

    float getbTagSys(){
	return _bTagSysW;
    }

    void setbTagSys(float bTagSysWeight){
	_bTagSysW = bTagSysWeight;
    }

    void setTrigger(bool accept){
	_trigger = accept;
    }
    
    bool getTriggerAccept(){
	return _trigger;
    }

    void setFilter(bool clean){
	_filter = clean;
    }

    bool getMETFilter(){
	return _filter;
    }

    void setPV(int passPV){
	_pv = passPV;
    }

    float getPVvalue(){
	return _pv;
    }

 private:
    std::vector<objectGenPart*>  _selectGenParts; 
    std::vector<objectJet*>                _jets;
    std::vector<objectJet*>               _bjets;
    objectMET*                              _MET;
    std::vector<objectLep*>               _muons;
    std::vector<objectLep*>           _electrons; 
    std::vector<objectJet*>          _selectJets;
    std::vector<objectJet*>         _selectbJets;
    std::vector<objectBoostedJet*>   _selectHadronicHiggses;
    std::vector<objectBoostedJet*>   _selectBoostedJets;
    std::vector<objectJet*>     _selectLightJets;
    std::vector<objectJet*>          _loosebJets;
    std::vector<objectLep*>     _selectElectrons; 
    std::vector<objectLep*>         _selectMuons; 
    std::vector<objectLep*>       _selectLeptons;

    bool _trigger = false, _filter = false;    
    int _pv = -1;

    float  _sumJetScalarpT=0., _sumSelJetScalarpT=0., _sumSelbJetScalarpT=0.,_sumSelHadronicHiggsScalarpT=0., _sumSelLightJetScalarpT=0., _sumSelMuonScalarpT=0., _sumSelElectronScalarpT=0., _sumSelJetMass=0., _sumSelbJetMass=0., _sumSelHadronicHiggsMass=0., _sumSelLightJetMass=0., _bTagSysW = 1. , _sumSelHadronicHiggsSoftDropMass=0;

    //    int nJet = 0, nbJet = 0, nSelJet = 0, nSelbJet = 0;
    int _nVetoLepton = 0;
    TLorentzVector _sumJetp4, _sumSelJetp4, _sumSelbJetp4, _sumHadronicHiggsp4, _sumLightJetp4, _sumSelMuonp4, _sumSelElectronp4; 
};

class ttHHanalyzer {
 public:
    enum sysName { kJES, kJER, kbTag, noSys };
    ttHHanalyzer(const std::string & cl, eventBuffer * ev, float weight = 1., bool systematics = false,
 		  std::string year= "nothing", std::string DataOrMC = "nothing", std::string sampleName = "nothing"){
	_weight = weight;
	_initialWeight = weight;   //Trigger SF for electron (TSFel)
	_ev = ev;
	_cl = cl;
	_sys = systematics;
	_of = new outputFile(_cl);
	_year= year;
	_DataOrMC = DataOrMC;
	_sampleName = sampleName;
    static const int LOG_INTERVAL = 1000;

	initHistograms();	
	initTree();
	initSys();
       	std::string dummy = "";
	HypoComb = new tthHypothesisCombinatorics(std::string("data/blrbdtweights_80X_V4/weights_64.xml"), std::string(""));
    }
    void createObjects(event*,sysName,bool);
    bool selectObjects(event*);
    void analyze(event*);
    void process(event*, sysName, bool);
    void loop(sysName, bool);
    void performAnalysis();
    void fillHistos(event * thisevent);
    void writeHistos();
    void fillTree(event * thisevent);
    void writeTree();
    TH1F * hmet,* hmetPhi, *hmetEta, *hAvgDeltaRjj, *hAvgDeltaRbb,*hAvgDeltaRbj, *hAvgDeltaEtajj, *hAvgDeltaEtabb, *hAvgDeltaEtabj, *hminDeltaRjj, *hminDeltaRbb, *hminDeltaRbj,  *hminDeltaRpTjj, *hminDeltaRpTbb, *hminDeltaRpTbj, *hminDeltaRMassjj, *hminDeltaRMassbb,*hminDeltaRMassbj, *hmaxDeltaEtajj, *hmaxDeltaEtabb, *hmaxDeltaEtabj, *hmaxPTmassjbb, *hmaxPTmassjjj, *hjetAverageMass, *hBjetAverageMass, *hHadronicHiggsAverageMass, *hLightJetAverageMass, *hBjetAverageMassSqr, *hHadronicHiggsSoftDropMass1, *hHadronicHiggsSoftDropMass2, *hjetHT, *hBjetHT, *hHadronicHiggsHT, *hLightJetHT, *hjetNumber, *hBjetNumber, *hHadronicHiggsNumber, *hLightJetNumber, *hInvMassHadW, *hInvMassZ1, *hInvMassZ2,*hInvMassHSingleMatched,*hInvMassHSingleNotMatched ,*hChi2HiggsSingleNotMatched, *hChi2HiggsSingleMatched , *hInvMassH1, *hInvMassH2,*hInvMassHZ1, *hInvMassHZ2, *hInvMassH1mChi, *hInvMassH2mChi,*hPTH1, *hPTH2, *hChi2Higgs, *hChi2HiggsZ, *hChi2HadW, *hChi2Z, *hAplanarity, *hSphericity, *hTransSphericity, *hCvalue, *hDvalue, *hBjetAplanarity, *hBjetSphericity, *hBjetTransSphericity ,*hBjetCvalue, *hBjetDvalue, *hCentralityjl, *hCentralityjb, *hleptonNumber,  *hElecNumber, *hMuonNumber, *hLeptonPT1, *hMuonPT1, *hElePT1, *hLeptonPhi1, *hMuonPhi1, *hElePhi1, *hLeptonEta1, *hMuonEta1, *hEleEta1, *hLeptonPT2, *hMuonPT2, *hElePT2, *hLeptonPhi2, *hMuonPhi2, *hElePhi2, *hLeptonEta2, *hMuonEta2, *hEleEta2, *hLepCharge1, *hLepCharge2, *hleptonHT, *hST, *hDiMuonMass, *hDiElectronMass, *hDiMuonPT, *hDiElectronPT, *hDiMuonEta, *hDiElectronEta, *hH0, *hH1, *hH2, *hH3, *hH4, *hR1, *hR2, * hR3, *hR4, *hBjetH0, *hBjetH1, *hBjetH2, *hBjetH3, *hBjetH4, *hBjetR1, *hBjetR2, * hBjetR3, *hBjetR4, *hCutFlow, *hCutFlow_w,
	*hInvMassHH1Matched,
	*hInvMassHH1NotMatched,
	*hInvMassHH2Matched,
        *hInvMassHH2NotMatched,
        *hChi2HHNotMatched,
        *hChi2HHMatched;


    tthHypothesisCombinatorics * HypoComb; 


    //fifo_map<std::string,int> cutflow{{"noCut", 0}, {"nlepton==2", 0}, {"nOpositeChargedLep", 0}, {"nMassCut", 0}, {"MET>40", 0}};
    //fifo_map<std::string,float> cutflow_w{{"noCut", 0}, {"nlepton==2", 0}, {"nOpositeChargedLep", 0}, {"nMassCut", 0}, {"MET>40", 0}};
    fifo_map<std::string,int> cutflow{ {"noCut", 0}, {"nHLTrigger", 0}, {"Muon_Trigger", 0}, {"Elec_Trigger", 0}, {"nFilter", 0}, {"nPV", 0}, {"njets>5", 0}, {"nbjets>4", 0}, {"nlepton==1", 0}, {"MET>20", 0}};
    fifo_map<std::string,int> cutflow_w{ {"noCut", 0}, {"nHLTrigger", 0}, {"Muon_Trigger", 0}, {"Elec_Trigger", 0}, {"nFilter", 0}, {"nPV", 0}, {"njets>5", 0}, {"nbjets>4", 0}, {"nlepton==1", 0}, {"MET>20", 0}};




    //    std::unordered_map<std::string, int> cutflow {{"noCut", 0}, {"njets>3", 0}, {"nbjets>2", 0}, {"nlepton==2", 0}, {"nOpositeChargedLep", 0}, {"nMassCut", 0}, {"nTotal", 0}};

 private: 
    std::unique_ptr<std::ofstream> sf_summary_log_file;  // Log para summary SF
    std::unique_ptr<std::ofstream> event_log_file; // logs
    std::unique_ptr<std::ofstream> sf_log_file; // logs
    long long total_electrons_processed = 0;
    long long total_muons_processed = 0;
    bool _sys;
    float _weight;
    float _initialWeight;  //Trigger SF for electron (TSFel)
//    int _year;
//    std::string _data;
//    std::string _sample;
    std::string _DataOrMC, _year, _sampleName; 
    TH1D * _hJES, * _hbJES, *_hbJetEff, *_hJetEff, *_hSysbTagM ;
    TString _pathJES = "HL_YR_JEC.root";
    TString _nameJES = "TOTAL_DIJET_AntiKt4EMTopo_YR2018";
    TString _namebJES = "TOTAL_BJES_AntiKt4EMTopo_YR2018";
    static const int nHistsJets = 8;
    static const int nHistsbJets = 6;
    static const int nHistsLightJets = 6;
    std::vector<TH1F*> hjetsPTs, hjetsEtas, hbjetsPTs, hbjetsEtas, hLightJetsPTs, hLightJetsEtas, hjetsBTagDisc, hbjetsBTagDisc, hLightJetsBTagDisc;
    event::evShapes jlepCent, jbjetCent;
    event::maxObjects jbbMaxs, jjjMaxs;
    event::statObjects jetStat, bjetStat, bjStat, ljetStat, lbjetStat, genPbjetStat; 
    event::foxWolframObjects jetFoxWolfMom, bjetFoxWolfMom;
    std::string _cl;
    eventBuffer * _ev;
    std::vector<event*> events;
    outputFile * _of;
    float _bbMassMinSHiggsNotMatched, _bbMassMinSHiggsMatched, _minChi2SHiggsNotMatched = 999999999. , _minChi2SHiggsMatched = 999999999.; 
    float _bbMassMinHH1NotMatched, _bbMassMinHH1Matched,_bbMassMinHH2NotMatched, _bbMassMinHH2Matched, _minChi2HHNotMatched = 999999999. , _minChi2HHMatched = 999999999.; 

    float _bbMassMin1Higgs, _bbMassMin2Higgs, _minChi2Higgs = 999999999.;
    float _bpTHiggs1, _bpTHiggs2;
    float _bbMassMin1HiggsZ, _bbMassMin2HiggsZ, _minChi2HiggsZ = 999999999.;
    float _bbMassMin1Z, _bbMassMin2Z, _minChi2Z = 999999999.;
    TRandom3 _rand;

	///////////////////////Electron ID (MVA ISO WP90) (IDSFel)
	void initEleIDSF();
	float getEleIDSF(float eta, float phi, float pt, float &unc);
	
	TH2F *h2_eleIDSF_others = nullptr;
	TH2F *h2_eleIDSF_2023B_Hole = nullptr;
	TH2F *h2_eleIDSF_2023B_NoHole = nullptr;
	
	float eleIDSFUncertainty = 0.0;

	///////////////////////Reco electron (RSFel)
    float getEleRecoSF(float eta, float pt, float &unc);
    void initRecoSF();

    TH2F *h2_eleRecoSF_low  = nullptr;
    TH2F *h2_eleRecoSF_mid  = nullptr;
    TH2F *h2_eleRecoSF_high = nullptr;

    float recoSFUncertainty = 0.0;

	///////////////////////Trigger electron (TSFel)
	// Eletron trigger SF
	TFile* eleTrigSFFile = nullptr;
	TH2F* h2_eleTrigSF = nullptr;
	TH2F* h2_eleTrigSF_unc = nullptr;
	TH2F* h2_effMC = nullptr;
	
	void initTriggerSF();
	float getEleTrigSF(float eta, float pt, float& sf_unc);
	float triggerSFUncertainty = 0.0;
	float eleTriggerSF = 1.0;  // para armazenar o SF do evento atual
	int _entryInLoop = -1;
	
	// Histogramas para elétrons
	TH1F* h_sf_vs_pt         = nullptr;
	TH1F* h_sf_vs_eta        = nullptr;
	TH1F* h_effMC_vs_pt      = nullptr;
	TH1F* h_effMC_vs_eta     = nullptr;
	TH1F* h_sf_vs_pt_sum     = nullptr;
	TH1F* h_sf_vs_pt_count   = nullptr;
	TH1F* h_sf_vs_eta_sum    = nullptr;
	TH1F* h_sf_vs_eta_count  = nullptr;
	TH1F* h_effMC_vs_pt_sum  = nullptr;
	TH1F* h_effMC_vs_pt_count= nullptr;
	TH1F* h_effMC_vs_eta_sum = nullptr;
	TH1F* h_effMC_vs_eta_count = nullptr;
	
	TH1F* h_sf_vs_pt_avg     = nullptr;
	TH1F* h_sf_vs_eta_avg    = nullptr;
	TH1F* h_effMC_vs_pt_avg  = nullptr;
	TH1F* h_effMC_vs_eta_avg = nullptr;
	
	// Muon trigger SF
	TH1F* h_sf_muon_vs_pt         = nullptr;
	TH1F* h_sf_muon_vs_eta        = nullptr;
	TH1F* h_sf_muon_vs_pt_sum     = nullptr;
	TH1F* h_sf_muon_vs_pt_count   = nullptr;
	TH1F* h_sf_muon_vs_eta_sum    = nullptr;
	TH1F* h_sf_muon_vs_eta_count  = nullptr;
	TH1F* h_sf_muon_vs_pt_avg     = nullptr;
	TH1F* h_sf_muon_vs_eta_avg    = nullptr;
	
	void initMuonHLTriggerSF();
	float getMuonTrigSF(float eta, float pt);

///////////////////////////


    void diMotherReco(const TLorentzVector & dPar1p4,const TLorentzVector & dPar2p4,const TLorentzVector & dPar3p4,const TLorentzVector & dPar4p4, const float mother1mass, const float  mother2mass, float & _minChi2,float & _bbMassMin1, float & _bbMassMin2);
    void motherReco(const TLorentzVector & dPar1p4,const TLorentzVector & dPar2p4, const float mother1mass, float & _minChi2,float & _bbMassMin1);

    /*    std::vector<double> getJetCutFlow(event *thisevent){
	int jetCounter = 0;
	std::vector<double> jetCutFlow;
        for(int i = 0; i < _ev->size(); i++){
	    jetCounter += 1;
            jetCutFlow.push_back(jetCounter);
	}
        return jetCutFlow;
	}*/


    std::vector<double> getJetCSV(event *thisevent){
	std::vector<double> jetCSV;
        for(int i = 0; i < thisevent->getnSelJet(); i++){
            jetCSV.push_back(thisevent->getSelJets()->at(i)->bTagCSV);
	}
        return jetCSV;
    }

    std::vector<double> getbJetCSV(event *thisevent){
	std::vector<double> bjetCSV;
        for(int i = 0; i < thisevent->getnbJet(); i++){
            bjetCSV.push_back(thisevent->getSelbJets()->at(i)->bTagCSV);
	}
        return bjetCSV;
    }

    std::vector<double> getlightJetCSV(event *thisevent){
	std::vector<double> lightjetCSV;
        for(int i = 0; i < thisevent->getnLightJet(); i++){
            lightjetCSV.push_back(thisevent->getSelLightJets()->at(i)->bTagCSV);
	}
        return lightjetCSV;
    }

    std::vector<TLorentzVector> getJetP4(event *thisevent){
	std::vector<TLorentzVector> jetP4;
	jetP4.reserve(thisevent->getnSelJet());
        for(const auto thisJet: *thisevent->getSelJets()){
            jetP4.push_back(*thisJet->getp4());
	}
        return jetP4;
    }

    std::vector<TLorentzVector> getLepP4(event *thisevent){
	std::vector<TLorentzVector> lepP4;
	lepP4.reserve(thisevent->getnSelLepton());
        for(const auto thisLep: *thisevent->getSelLeptons()){
            lepP4.push_back(*thisLep->getp4());
	}
        return lepP4;
    }


    /*    float getBTagValue(event *thisevent, int jetNo){
	float bTagValue;
	if(thisevent->getnSelJet() <= jetNo-1){
	    return -1;
	}

	if(bool(thisevent->getSelJets()->at(jetNo-1)->bTag & (1<<2)) == true){
	    bTagValue = 0.99;
	} else if(bool(thisevent->getSelJets()->at(jetNo-1)->bTag & (1<<1)) == true){
	    bTagValue = 0.66;
	} else if(bool(thisevent->getSelJets()->at(jetNo-1)->bTag & (1<<0)) == true){
	    bTagValue = 0.33;
	} else {
	    bTagValue = 0;
	}
	return bTagValue;
	}*/


    float getSysJES(TH1* hSys, float pT){
	return hSys->GetBinContent(hSys->FindBin(pT));
    } 

    float getSysJER(float sigma){
	return _rand.Gaus(sigma/2., sigma);
    } 

    void initSys(){
	TFile *_fJES = TFile::Open(_pathJES);
	_hJES = (TH1D*)_fJES->Get(_nameJES);
	_hbJES = (TH1D*)_fJES->Get(_namebJES);
	int npTbin = 9;
	float sysbTagM[] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.016, 0.018, 0.023, 0.046 };
	float pTBinEdges[] = { 30, 50, 70, 100, 140, 200, 300, 600, 1000, 3000 };
	_hbJetEff = new TH1D("bjeteff","bjet efficiency", npTbin, pTBinEdges);
	_hJetEff = new TH1D("jeteff","jet efficiency", npTbin, pTBinEdges);
	_hSysbTagM = new TH1D("bTagMSys","btag medium systematics", npTbin, pTBinEdges);
	int nbinsx = _hSysbTagM->GetXaxis()->GetNbins();
	for(int bind = 1; bind < nbinsx+1; bind++){
	    _hSysbTagM->SetBinContent(bind, sysbTagM[bind-1]);
	}
    }

    void getbJetEffMap(){
	_hbJetEff->Divide(_hJetEff);
    }

    std::vector<TDirectory *> _histoDirs; 
    std::vector<TDirectory *> _treeDirs; 

    void initHistograms(sysName sysType = noSys, bool up = false){
	hjetsPTs.resize(nHistsJets); hjetsEtas.resize(nHistsJets); hbjetsPTs.resize(nHistsbJets); hbjetsEtas.resize(nHistsbJets); hLightJetsPTs.resize(nHistsLightJets), hLightJetsEtas.resize(nHistsLightJets), hjetsBTagDisc.resize(nHistsJets), hbjetsBTagDisc.resize(nHistsbJets), hLightJetsBTagDisc.resize(nHistsLightJets);



	hCutFlow = new TH1F("cutflow", "N_{cutFlow}", cutflow.size(), 0, cutflow.size()-1);
	hCutFlow_w = new TH1F("cutflow_w", "N_{weighted}", cutflow.size(), 0, cutflow.size()-1);

	    
	TString trail = "";
	if(sysType == kbTag){
	    if(up) trail += "btag_up";
	    else trail += "btag_down";
	} else if(sysType == kJES){
	    if(up) trail += "JES_up";
	    else trail += "JES_down";
	}else if(sysType == kJER){
	    if(up) trail += "JER_up";
	    else trail += "JER_down";
	}

	_of->file->cd();
	std::vector<TDirectory *> tmpDirs; 	
	TDirectory *jet = _of->file->mkdir("jet"+trail);
	tmpDirs.push_back(jet); 
	jet->cd();

	hmet = new TH1F("met"+trail, "MET"+trail, 50, 0, 1000);
	hmetPhi = new TH1F("metPhi"+trail, "MET #phi"+trail, 50, -5, 5);
	hmetEta = new TH1F("metEta"+trail, "MET #eta"+trail, 50, -5, 5);


	for(int i=0; i < nHistsJets; i++){
	    if(i < 3){
		hjetsPTs.at(i)  = new TH1F(TString::Format("jetPT%d",(i+1))+trail, TString::Format("jet%d p_{T} [GeV]",i+1)+trail, 50, 0, 1500);
		hjetsEtas.at(i) = new TH1F(TString::Format("jetEta%d",(i+1))+trail, TString::Format("jet%d #eta",i+1)+trail, 50, -4, 4);
	    } else  {
		hjetsPTs.at(i)  = new TH1F(TString::Format("jetPT%d",(i+1))+trail, TString::Format("jet%d p_{T} [GeV]",i+1)+trail, 50, 0, 350);
		hjetsEtas.at(i) = new TH1F(TString::Format("jetEta%d",(i+1))+trail, TString::Format("jet%d #eta",i+1)+trail, 50, -4, 4);

	    }
	    hjetsBTagDisc.at(i)  = new TH1F(TString::Format("jetBTagDisc%d",(i+1))+trail, TString::Format("jet%d btagDisc" ,i+1)+trail, 50, 0, 1);

	}

	for(int i=0; i < nHistsbJets; i++){
	    if(i < 3){
		hbjetsPTs.at(i) = new TH1F(TString::Format("bjetPT%d",(i+1))+trail, TString::Format("bjet%d p_{T} [GeV]",i+1)+trail, 50, 0, 1500);
		hbjetsEtas.at(i) = new TH1F(TString::Format("bjetEta%d",(i+1))+trail, TString::Format("bjet%d #eta",i+1)+trail, 50, -4, 4);
	    } else {
	    hbjetsPTs.at(i) = new TH1F(TString::Format("bjetPT%d",(i+1))+trail, TString::Format("bjet%d p_{T} [GeV]",i+1)+trail, 50, 0, 350);
	    hbjetsEtas.at(i) = new TH1F(TString::Format("bjetEta%d",(i+1))+trail, TString::Format("bjet%d #eta",i+1)+trail, 50, -4, 4);
	    }
	    hbjetsBTagDisc.at(i)  = new TH1F(TString::Format("bjetBTagDisc%d",(i+1))+trail, TString::Format("bjet%d btagDisc" ,i+1)+trail, 50, 0, 1);
	}


	for(int i=0; i < nHistsLightJets; i++){
	    if(i < 3){
		hLightJetsPTs.at(i) = new TH1F(TString::Format("lightJetPT%d",(i+1))+trail, TString::Format("lightJet%d p_{T} [GeV]",i+1)+trail, 50, 0, 400);
		hLightJetsEtas.at(i) = new TH1F(TString::Format("lightJetEta%d",(i+1))+trail, TString::Format("lightJet%d #eta",i+1)+trail, 50, -4, 4);
	    } else {hLightJetsPTs.at(i) = new TH1F(TString::Format("lightJetPT%d",(i+1))+trail, TString::Format("lightJet%d p_{T} [GeV]",i+1)+trail, 50, 0, 150);
		hLightJetsEtas.at(i) = new TH1F(TString::Format("lightJetEta%d",(i+1))+trail, TString::Format("lightJet%d #eta",i+1)+trail, 50, -4, 4);
	    }
	    hLightJetsBTagDisc.at(i)  = new TH1F(TString::Format("lightJetBTagDisc%d",(i+1))+trail, TString::Format("lightJet%d btagDisc" ,i+1)+trail, 50, 0, 1);
	}
   

	hAvgDeltaRjj = new TH1F("deltaRavgjj"+trail, "#DeltaR_{jj}^{avg}"+trail, 50, 0, 5);
	hAvgDeltaRbb = new TH1F("deltaRavgbb"+trail, "#DeltaR_{bb}^{avg}"+trail, 50, 0, 5);
	hAvgDeltaRbj = new TH1F("deltaRavgbj"+trail, "#DeltaR_{bj}^{avg}"+trail, 50, 0, 5);
	hAvgDeltaEtajj = new TH1F("deltaEtaavgjj"+trail, "#Delta#eta_{jj}^{avg}"+trail, 50, 0, 3);
	hAvgDeltaEtabb = new TH1F("deltaEtaavgbb"+trail, "#Delta#eta_{bb}^{avg}"+trail, 50, 0, 3);
	hAvgDeltaEtabj = new TH1F("deltaEtaavgbj"+trail, "#Delta#eta_{bj}^{avg}"+trail, 50, 0, 3);
	hminDeltaRjj = new TH1F("deltaRminjj"+trail, "#DeltaR_{jj}^{min}"+trail, 50, 0, 3);
	hminDeltaRbb = new TH1F("deltaRminbb"+trail, "#DeltaR_{bb}^{min}"+trail, 50, 0, 3);
	hminDeltaRbj = new TH1F("deltaRminbj"+trail, "#DeltaR_{bj}^{min}"+trail, 50, 0, 3);
	hminDeltaRpTjj = new TH1F("pTdeltaRminjj"+trail, "#DeltaR_{jj, p_{T}}^{min}"+trail, 50, 0, 600);
	hminDeltaRpTbb = new TH1F("pTdeltaRminbb"+trail, "#DeltaR_{bb, p_{T}}^{min}"+trail, 50, 0, 600);
	hminDeltaRpTbj = new TH1F("pTdeltaRminbj"+trail, "#DeltaR_{bj, p_{T}}^{min}"+trail, 50, 0, 600);
	hminDeltaRMassjj = new TH1F("massDeltaRminjj"+trail, "#DeltaR_{jj, mass}^{min}"+trail, 50, 0, 100);
	hminDeltaRMassbb = new TH1F("massDeltaRminbb"+trail, "#DeltaR_{bb, mass}^{min}"+trail, 50, 0, 100);
	hminDeltaRMassbj = new TH1F("massDeltaRminbj"+trail, "#DeltaR_{bj, mass}^{min}"+trail, 50, 0, 100);
	hmaxDeltaEtabb = new TH1F("deltaEtamaxbb"+trail, "#Delta#eta_{bb}^{max}"+trail, 50, 0, 5);
	hmaxDeltaEtajj = new TH1F("deltaEtamaxjj"+trail, "#Delta#eta_{jj}^{max}"+trail, 50, 0, 5);
	hmaxDeltaEtabj = new TH1F("deltaEtamaxbj"+trail, "#Delta#eta_{bj}^{max}"+trail, 50, 0, 5);
	hmaxPTmassjbb = new TH1F("maxPTmassjbb"+trail, "m_{jbb}^{max p_{T}}"+trail, 50, 0, 100); 
	hmaxPTmassjjj = new TH1F("maxPTmassjjj"+trail, "m_{jjj}^{max p_{T}}"+trail, 50, 0, 100); 
	hjetAverageMass = new TH1F("jetAvgMass"+trail, "m_{j}^{avg}"+trail, 50, 0, 60);
	hBjetAverageMass = new TH1F("jetBAvgMass"+trail, "m_{b}^{avg}"+trail, 50, 0, 1000);
	hHadronicHiggsAverageMass = new TH1F("higgsHadAvgMass"+trail, "m_{H_{had}}^{avg}"+trail, 50, 0, 60);
	hLightJetAverageMass = new TH1F("jetLightAvgMass"+trail, "m_{light}^{avg}"+trail, 50, 0, 60);
	hBjetAverageMassSqr = new TH1F("jetBAvgMassSqr"+trail, "(m^{2})_{b}^{avg}"+trail, 50, 0, 2500);
	hHadronicHiggsSoftDropMass1 = new TH1F("higgsHadSoftDropMass1"+trail, "msoftdrop_{H_{had}}"+trail, 50, 0, 400);
	hHadronicHiggsSoftDropMass2 = new TH1F("higgsHadSoftDropMass2"+trail, "msoftdrop_{H_{had}}"+trail, 50, 0, 300);
	hjetHT = new TH1F("jetHT"+trail, "H_{T} [GeV]"+trail, 50, 0, 3500);
	hBjetHT = new TH1F("jetBHT"+trail, "H_{T}^{b} [GeV]"+trail, 50, 0, 3500); 
	hHadronicHiggsHT = new TH1F("jetHadronicHiggsHT"+trail, "H_{T}^{H_{had}} [GeV]"+trail, 50, 0, 4000); 
	hLightJetHT = new TH1F("jetLightHT"+trail, "H_{T}^{light} [GeV]"+trail, 50, 0, 1200); 
	hjetNumber = new TH1F("jetNumber"+trail, "N_{jet}"+trail, 13, 2, 15);
	hBjetNumber = new TH1F("jetBNumber"+trail, "N_{bjet}"+trail, 8, 2, 10); 
	hHadronicHiggsNumber = new TH1F("jetHadronicHiggsNumber"+trail, "N_{H_{had}}"+trail, 8, 2, 10); 
	hLightJetNumber = new TH1F("jetLightNumber"+trail, "N_{lightJet}"+trail, 10, 0, 10); 
	hInvMassHadW = new TH1F("invMass_hadW"+trail, "m_{W,had}"+trail, 50, 0, 400);
	hInvMassZ1 = new TH1F("invMass_Z1"+trail, "m_{Z,1} [GeV]"+trail, 50, 0, 400);
	hInvMassZ2 = new TH1F("invMass_Z2"+trail, "m_{Z,2} [GeV]"+trail, 50, 0, 400);
	hInvMassH1 = new TH1F("invMass_Higgs1"+trail, "m_{H,1} [GeV]"+trail, 50, 0, 500);
	hInvMassH2 = new TH1F("invMass_Higgs2"+trail, "m_{H,2} [GeV]"+trail, 50, 0, 500);
	hInvMassH1mChi = new TH1F("invMass_Higgs1_mChi"+trail, "m_{H,1} min(#chi^{2})"+trail, 50, 0, 500);
	hInvMassH2mChi = new TH1F("invMass_Higgs2_mChi"+trail, "m_{H,2} min(#chi^{2})"+trail, 50, 0, 500);

	hPTH1 = new TH1F("pT_Higgs1"+trail, "p_{T(H,1)} [GeV]"+trail, 50, 0, 500);
	hPTH2 = new TH1F("pT_Higgs2"+trail, "p_{T(H,2)} [GeV]"+trail, 50, 0, 500);
	hInvMassHZ1 = new TH1F("invMass_HiggsZ1"+trail, "m^{Z}_{H,1} [GeV]"+trail, 50, 0, 500);
	hInvMassHZ2 = new TH1F("invMass_HiggsZ2"+trail, "m^{Z}_{H,2} [GeV]"+trail, 50, 0, 500);
	hChi2Higgs = new TH1F("chi2Higgs"+trail, "#chi^{2}_{HH}"+trail, 50, 0, 5000);
	hChi2Z = new TH1F("chi2Z"+trail, "#chi^{2}_{ZZ}"+trail, 50, 0, 5000);
	hChi2HiggsZ = new TH1F("chi2HiggsZ"+trail, "#chi^{2}_{HZ}"+trail, 50, 0, 5000);
	//hChi2HadW = new TH1F("chi2HadW"+trail, "#chi^{2}_{W,had}"+trail, 50, 0, 400);  
	hInvMassHSingleMatched = new TH1F("invMass_HiggsMatched"+trail, "m_{H,matched} [GeV]"+trail, 50, 0, 500); 
        hInvMassHSingleNotMatched  = new TH1F("invMass_HiggsNotMatched"+trail, "m_{H,unmatched} [GeV]"+trail, 50, 0, 500); 
        hChi2HiggsSingleNotMatched = new TH1F("chi2HiggsNotMatched"+trail, "#chi^{2}_{H,unmatched}"+trail, 50, 0, 10);
        hChi2HiggsSingleMatched =new TH1F("chi2HiggsMatched"+trail, "#chi^{2}_{H,matched}"+trail, 50, 0, 10);

	hInvMassHH1Matched = new TH1F("invMass_HH1Matched"+trail, "m_{H1,matched} [GeV]"+trail, 50, 0, 500); 
        hInvMassHH1NotMatched  = new TH1F("invMass_HH1NotMatched"+trail, "m_{H1,unmatched} [GeV]"+trail, 50, 0, 500); 
	hInvMassHH2Matched = new TH1F("invMass_HH2Matched"+trail, "m_{H2,matched} [GeV]"+trail, 50, 0, 500); 
        hInvMassHH2NotMatched  = new TH1F("invMass_HH2NotMatched"+trail, "m_{H2,unmatched} [GeV]"+trail, 50, 0, 500); 
        hChi2HHNotMatched = new TH1F("chi2HHNotMatched"+trail, "#chi^{2}_{H,unmatched}"+trail, 50, 0, 10);
        hChi2HHMatched =new TH1F("chi2HHMatched"+trail, "#chi^{2}_{H,matched}"+trail, 50, 0, 10);


	hAplanarity = new TH1F("aplanarity"+trail, "A"+trail, 50, 0, 0.5);  
	hSphericity = new TH1F("sphericity"+trail, "S"+trail, 50, 0, 1);  
	hTransSphericity = new TH1F("transSphericity"+trail, "S_{#perp}"+trail, 50, 0, 1);  
	hCvalue = new TH1F("C"+trail, "C value"+trail, 50, 0, 1);
	hDvalue = new TH1F("D"+trail, "D value"+trail, 50, 0, 1); 
	hCentralityjb = new TH1F("centralityjb"+trail, "centrality_{jb}"+trail, 50, 0, 1); 
	hCentralityjl = new TH1F("centralityjl"+trail, "centrality_{jl}"+trail, 50, 0, 1); 
	
	hH0 = new TH1F("H0"+trail, "H_{0}"+trail, 50, 0.2, 0.45); 
	hH1 = new TH1F("H1"+trail, "H_{1}"+trail, 50, -0.2, 0.45); 
	hH2 = new TH1F("H2"+trail, "H_{2}"+trail, 50, -0.2, 0.3); 
	hH3 = new TH1F("H3"+trail, "H_{3}"+trail, 50, -0.2, 0.3); 
	hH4 = new TH1F("H4"+trail, "H_{4}"+trail, 50, -0.2, 0.3); 

	hR1 = new TH1F("R1"+trail, "R_{1}"+trail, 50, 0, 1); 
	hR2 = new TH1F("R2"+trail, "R_{2}"+trail, 50, 0, 1); 
	hR3 = new TH1F("R3"+trail, "R_{3}"+trail, 50, 0, 1); 
	hR4 = new TH1F("R4"+trail, "R_{4}"+trail, 50, 0, 1); 

	hBjetH0 = new TH1F("H0_bjet"+trail, "H_{0,bjet}"+trail, 50, -0.2, 0.45); 
	hBjetH1 = new TH1F("H1_bjet"+trail, "H_{1,bjet}"+trail, 50, -0.2, 0.45); 
	hBjetH2 = new TH1F("H2_bjet"+trail, "H_{2,bjet}"+trail, 50, -0.2, 0.3); 
	hBjetH3 = new TH1F("H3_bjet"+trail, "H_{3,bjet}"+trail, 50, -0.2, 0.3); 
	hBjetH4 = new TH1F("H4_bjet"+trail, "H_{4,bjet}"+trail, 50, -0.2, 0.3); 

	hBjetR1 = new TH1F("R1_bjet"+trail, "R_{1,bjet}"+trail, 50, 0, 1); 
	hBjetR2 = new TH1F("R2_bjet"+trail, "R_{2,bjet}"+trail, 50, 0, 1); 
	hBjetR3 = new TH1F("R3_bjet"+trail, "R_{3,bjet}"+trail, 50, 0, 1); 
	hBjetR4 = new TH1F("R4_bjet"+trail, "R_{4,bjet}"+trail, 50, 0, 1); 

	hBjetAplanarity = new TH1F("aplanarity_bjet"+trail, "A_{bjet}"+trail, 50, 0, 0.5);  
	hBjetSphericity = new TH1F("sphericity_bjet"+trail, "S_{bjet}"+trail, 50, 0, 1);  
	hBjetTransSphericity = new TH1F("transSphericity_bjet"+trail, "S_{#perp, bjet}"+trail, 50, 0, 1);  
	hBjetCvalue = new TH1F("C_bjet"+trail, "C value_{bjet}"+trail, 50, 0, 1);
	hBjetDvalue = new TH1F("D_bjet"+trail, "D value_{bjet}"+trail, 50, 0, 1); 


	    
	_of->file->cd();
	TDirectory *lepton = _of->file->mkdir("Lepton"+trail);
	tmpDirs.push_back(lepton);
	lepton->cd();

/////////////Electron Trigger SF //////////////////  (TSFel)
	h_sf_vs_pt     = new TH1F("h_sf_vs_pt",     "SF vs pT;Electron pT [GeV];SF",        100, 0, 700);
	h_sf_vs_eta    = new TH1F("h_sf_vs_eta",    "SF vs Eta;Electron #eta;SF",           100, -5, 5);
	h_effMC_vs_pt  = new TH1F("h_effMC_vs_pt",  "EffMC vs pT;Electron pT [GeV];Eff.",   100, 0, 700);
	h_effMC_vs_eta = new TH1F("h_effMC_vs_eta", "EffMC vs Eta;Electron #eta;Eff.",      100, -5, 5);
	h_sf_vs_pt_sum    = new TH1F("h_sf_vs_pt_sum",    "Sum SF vs pT",        100, 0, 700);
	h_sf_vs_pt_count  = new TH1F("h_sf_vs_pt_count",  "Count SF vs pT",      100, 0, 700);
	h_sf_vs_eta_sum   = new TH1F("h_sf_vs_eta_sum",   "Sum SF vs eta",       100, -5, 5);
	h_sf_vs_eta_count = new TH1F("h_sf_vs_eta_count", "Count SF vs eta",     100, -5, 5);
	
	h_effMC_vs_pt_sum    = new TH1F("h_effMC_vs_pt_sum",    "Sum Eff vs pT",      100, 0, 700);
	h_effMC_vs_pt_count  = new TH1F("h_effMC_vs_pt_count",  "Count Eff vs pT",    100, 0, 700);
	h_effMC_vs_eta_sum   = new TH1F("h_effMC_vs_eta_sum",   "Sum Eff vs eta",     100, -5, 5);
	h_effMC_vs_eta_count = new TH1F("h_effMC_vs_eta_count", "Count Eff vs eta",   100, -5, 5);
	
	//////////////////////////////////////////////////
	 /////////////Muon Trigger SF //////////////////  (TSFmu)
	h_sf_muon_vs_pt        = new TH1F("h_sf_muon_vs_pt", "Muon SF vs pT;Muon pT [GeV];SF", 100, 0, 700);
	h_sf_muon_vs_eta       = new TH1F("h_sf_muon_vs_eta", "Muon SF vs Eta;Muon #eta;SF", 100, -5, 5);
	h_sf_muon_vs_pt_sum    = new TH1F("h_sf_muon_vs_pt_sum", "Sum SF vs pT", 100, 0, 700);
	h_sf_muon_vs_pt_count  = new TH1F("h_sf_muon_vs_pt_count", "Count SF vs pT", 100, 0, 700);
	h_sf_muon_vs_eta_sum   = new TH1F("h_sf_muon_vs_eta_sum", "Sum SF vs eta", 100, -5, 5);
	h_sf_muon_vs_eta_count = new TH1F("h_sf_muon_vs_eta_count", "Count SF vs eta", 100, -5, 5);
	h_sf_muon_vs_pt_avg    = new TH1F("h_sf_muon_vs_pt_avg", "Avg SF vs pT", 100, 0, 700);
	h_sf_muon_vs_eta_avg   = new TH1F("h_sf_muon_vs_eta_avg", "Avg SF vs eta", 100, -5, 5);
///////////////////////////////////////////////
	hLepCharge1 = new TH1F("lepCharge1"+trail, "lepCh1"+trail, 4, -2, 2);
	hLepCharge2 = new TH1F("lepCharge2"+trail, "lepCh2"+trail, 4, -2, 2);
	hLepCharge1 = new TH1F("lepCharge1"+trail, "lepCh1"+trail, 4, -2, 2);
	hLepCharge2 = new TH1F("lepCharge2"+trail, "lepCh2"+trail, 4, -2, 2);
	    
	hleptonNumber = new TH1F("lepNumber"+trail, "N_{lep}"+trail, 4, 0, 4);
	hElecNumber = new TH1F("ElecNumber"+trail, "N_{Elec}"+trail, 4, 0, 4);
	hMuonNumber = new TH1F("MuonNumber"+trail, "N_{Muon}"+trail, 4, 0, 4);
	    
	hleptonNumber = new TH1F("lepNumber"+trail, "N_{lep}"+trail, 4, 0, 4);
	hleptonHT = new TH1F("leptonHT"+trail, "H_{T}^{lep} [GeV]"+trail, 50, 0, 3500);
	hST = new TH1F("ST"+trail, "S_{T} [GeV]"+trail, 50, 0, 2000);
	hDiMuonMass = new TH1F("diMuonMass"+trail, "m_{#mu#mu} [GeV]"+trail, 50, 0, 200);
	hDiElectronMass = new TH1F("diEleMass"+trail, "m_{ee} [GeV]"+trail, 50, 0, 200);
	hDiMuonPT = new TH1F("diMuonPT"+trail, "p_{T, #mu#mu} [GeV]"+trail, 50, 0, 400);
	hDiElectronPT = new TH1F("diElePT"+trail, "p_{T, ee} [GeV]"+trail, 50, 0, 400);
	hDiMuonEta = new TH1F("diMuonEta"+trail, "#eta_{#mu#mu}"+trail, 50, -3, 3);
	hDiElectronEta = new TH1F("diEleEta"+trail, "#eta_{ee}"+trail, 50, -3, 3);

	hLeptonPT1 = new TH1F("leptonPT1"+trail, "lepton p_{T,1}"+trail, 50, 0, 700);
	hMuonPT1 = new TH1F("muonPT1"+trail, "muon p_{T,1}"+trail, 50, 0, 700);
	hElePT1 = new TH1F("elePT1"+trail, "ele p_{T,1}"+trail, 50, 0, 700);

	hLeptonPhi1 = new TH1F("leptonPhi1"+trail, "lepton #phi_{1}"+trail, 50, -4, 4);
	hMuonPhi1 = new TH1F("muonPhi1"+trail, "muon #phi_{1}"+trail, 50, -4, 4);
	hElePhi1 = new TH1F("elePhi1"+trail, "ele #phi_{1}"+trail, 50, -4, 4);

	hLeptonEta1 = new TH1F("leptonEta1"+trail, "lepton #eta_{1}"+trail, 50, -3, 3);
	hMuonEta1 = new TH1F("muonEta1"+trail, "muon #eta_{1}"+trail, 50, -3, 3);
	hEleEta1 = new TH1F("eleEta1"+trail, "ele #eta_{1}"+trail, 50, -3, 3); 

	hLeptonPT2 = new TH1F("leptonPT2"+trail, "lepton p_{T,2}"+trail, 50, 0, 500);
	hMuonPT2 = new TH1F("muonPT2"+trail, "muon p_{T,2}"+trail, 50, 0, 500);
	hElePT2 = new TH1F("elePT2"+trail, "ele p_{T,2}"+trail, 50, 0, 500);
	hLeptonPhi2 = new TH1F("leptonPhi2"+trail, "lepton #phi_{2}"+trail, 50, -4, 4);
	hMuonPhi2 = new TH1F("muonPhi2"+trail, "muon #phi_{2}"+trail, 50, -4, 4);
	hElePhi2 = new TH1F("elePhi2"+trail, "ele #phi_{2}"+trail, 50, -4, 4);
	hLeptonEta2 = new TH1F("leptonEta2"+trail, "lepton #eta_{2}"+trail, 50, -3, 3);
	hMuonEta2 = new TH1F("muonEta2"+trail, "muon #eta_{2}"+trail, 50, -3, 3);
	hEleEta2 = new TH1F("eleEta2"+trail, "ele #eta_{2}"+trail, 50, -3, 3); 
	
	_histoDirs = tmpDirs;
    }

    
    TTree * _inputTree;
    float jetPT1, jetPT2, jetPT3, jetPT4, jetPT5, jetPT6, jetPT7, jetPT8, bjetPT1, bjetPT2, bjetPT3, bjetPT4, bjetPT5, bjetPT6, bjetPT7, bjetPT8, jetEta1, jetEta2, jetEta3, jetEta4, jetEta5, jetEta6, jetEta7, jetEta8, bjetEta1, bjetEta2, bjetEta3, bjetEta4, bjetEta5, bjetEta6, bjetEta7, bjetEta8, lightjetPT1, lightjetPT2, lightjetPT3, lightjetEta1, lightjetEta2, lightjetEta3, jetBTagDisc1, jetBTagDisc2, jetBTagDisc3, jetBTagDisc4, jetBTagDisc5, jetBTagDisc6, jetBTagDisc7, jetBTagDisc8, bjetBTagDisc1,bjetBTagDisc2, bjetBTagDisc3, bjetBTagDisc4, bjetBTagDisc5, bjetBTagDisc6, bjetBTagDisc7, bjetBTagDisc8, lightjetBTagDisc1, lightjetBTagDisc2, lightjetBTagDisc3, met, metPhi, metEta, averageDeltaRjj, averageDeltaRbb, averageDeltaEtajj, averageDeltaEtabb, minDeltaRjj, minDeltaRbb, maxDeltaEtajj, maxDeltaEtabb, jetAverageMass, bjetAverageMass, lightJetAverageMass, bjetAverageMassSqr, jetHT, bjetHT, lightjetHT, invMassHadW, invMassZ1, invMassZ2, invMassH1, invMassH2, chi2Higgs, chi2HiggsZ, chi2HadW, chi2Z, invMassHiggsZ1, invMassHiggsZ2, PTH1, PTH2, weight, aplanarity, sphericity, transSphericity, cValue, dValue, baplanarity, centralityjb, centralityjl, bsphericity, btransSphericity, bcValue, bdValue, leptonEta1, muonEta1, eleEta1, leptonPT1, muonPT1, elePT1, leptonEta2, muonEta2, eleEta2, leptonPT2, muonPT2, elePT2, diElectronMass, diMuonMass, leptonHT, ST, leptonCharge1, leptonCharge2, H0, H1, H2, H3, H4, bH0, bH1, bH2, bH3, bH4, R1, R2, R3, R4, bR1, bR2, bR3, bR4, maxPTmassjbb, maxPTmassjjj, minDeltaRpTbb, minDeltaRpTjj, minDeltaRpTbj, minDeltaRMassjj, minDeltaRMassbj, minDeltaRMassbb, averageDeltaRbj,  averageDeltaEtabj, minDeltaRbj, maxDeltaEtabj, bbjetHiggsMatched1, bbjetHiggsMatched2, bbjetHiggsMatched3, bbjetHiggsMatched4, bbjetHiggsMatched5, bbjetHiggsMatched6,bbjetHiggsMatcheddR1, bbjetHiggsMatcheddR2, bbjetHiggsMatcheddR3, bbjetHiggsMatcheddR4, bbjetHiggsMatcheddR5, bbjetHiggsMatcheddR6, bbjetMinChiHiggsIndex1, bbjetMinChiHiggsIndex2, bbjetMinChiHiggsIndex3, bbjetMinChiHiggsIndex4, bbjetMinChiHiggsIndex5, bbjetMinChiHiggsIndex6;

    int jetNumber, bjetNumber, lightjetNumber; 
    void initTree(sysName sysType = noSys, bool up = false){
	_of->file->cd();
	std::vector<TDirectory*> tmpDirs;
	TString trail = "";
	if(sysType == kbTag){
	    if(up) trail += "_btag_up";
	    else trail += "_btag_down";
	} else if(sysType == kJES){
	    if(up) trail += "_JES_up";
	    else trail += "_JES_down";
	} else if(sysType == kJER){
	    if(up) trail += "_JER_up";
	    else trail += "_JER_down";
	}
	
	TDirectory *tree = _of->file->mkdir("Tree"+trail);
	tree->cd();
	tmpDirs.push_back(tree);
	
        _inputTree = new  TTree("Tree","tree for dnn inputs");

	_inputTree->Branch("jetPT1", &jetPT1, "jetPT1/f");
	_inputTree->Branch("jetPT2", &jetPT2, "jetPT2/f");	  
	_inputTree->Branch("jetPT3", &jetPT3, "jetPT3/f");
	_inputTree->Branch("jetPT4", &jetPT4, "jetPT4/f");	 
	_inputTree->Branch("jetPT5", &jetPT5, "jetPT5/f");
	_inputTree->Branch("jetPT6", &jetPT6, "jetPT6/f");	
	_inputTree->Branch("jetPT7", &jetPT7, "jetPT7/f");	
	_inputTree->Branch("jetPT8", &jetPT8, "jetPT8/f");	
	_inputTree->Branch("bjetPT1", &bjetPT1, "bjetPT1/f"); 
	_inputTree->Branch("bjetPT2", &bjetPT2, "bjetPT2/f");
	_inputTree->Branch("bjetPT3", &bjetPT3, "bjetPT3/f");	
	_inputTree->Branch("bjetPT4", &bjetPT4, "bjetPT4/f");	 
	_inputTree->Branch("bjetPT5", &bjetPT5, "bjetPT5/f");
	_inputTree->Branch("bjetPT6", &bjetPT6, "bjetPT6/f");
	_inputTree->Branch("bjetPT7", &bjetPT7, "bjetPT7/f");
	_inputTree->Branch("bjetPT8", &bjetPT8, "bjetPT8/f");
	_inputTree->Branch("lightjetPT1", &lightjetPT1, "lightjetPT1/f");
	_inputTree->Branch("lightjetPT2", &lightjetPT2, "lightjetPT2/f");	  
	_inputTree->Branch("lightjetPT3", &lightjetPT3, "lightjetPT3/f");	
	_inputTree->Branch("jetEta1", &jetEta1, "jetEta1/f");  
	_inputTree->Branch("jetEta2", &jetEta2, "jetEta2/f");
	_inputTree->Branch("jetEta3", &jetEta3, "jetEta3/f");	
	_inputTree->Branch("jetEta4", &jetEta4, "jetEta4/f");	  
	_inputTree->Branch("jetEta5", &jetEta5, "jetEta5/f");	
	_inputTree->Branch("jetEta6", &jetEta6, "jetEta6/f");
	_inputTree->Branch("jetEta7", &jetEta7, "jetEta7/f");
	_inputTree->Branch("jetEta8", &jetEta8, "jetEta8/f");
	_inputTree->Branch("bjetEta1", &bjetEta1, "bjetEta1/f");	 
	_inputTree->Branch("bjetEta2", &bjetEta2, "bjetEta2/f");	
	_inputTree->Branch("bjetEta3", &bjetEta3, "bjetEta3/f");	  
	_inputTree->Branch("bjetEta4", &bjetEta4, "bjetEta4/f");	  
	_inputTree->Branch("bjetEta5", &bjetEta5, "bjetEta5/f");	
	_inputTree->Branch("bjetEta6", &bjetEta6, "bjetEta6/f");  
	_inputTree->Branch("bjetEta7", &bjetEta7, "bjetEta7/f");  
	_inputTree->Branch("bjetEta8", &bjetEta8, "bjetEta8/f");  
	_inputTree->Branch("lightjetEta1", &lightjetEta1, "lightjetEta1/f");
	_inputTree->Branch("lightjetEta2", &lightjetEta2, "lightjetEta2/f");	  
	_inputTree->Branch("lightjetEta3", &lightjetEta3, "lightjetEta3/f");	
	_inputTree->Branch("jetBTagDisc1", &jetBTagDisc1, "jetBTagDisc1/f");  
	_inputTree->Branch("jetBTagDisc2", &jetBTagDisc2, "jetBTagDisc2/f");  
	_inputTree->Branch("jetBTagDisc3", &jetBTagDisc3, "jetBTagDisc3/f");  
	_inputTree->Branch("jetBTagDisc4", &jetBTagDisc4, "jetBTagDisc4/f");  
	_inputTree->Branch("jetBTagDisc5", &jetBTagDisc5, "jetBTagDisc5/f");  
	_inputTree->Branch("jetBTagDisc6", &jetBTagDisc6, "jetBTagDisc6/f");  
	_inputTree->Branch("jetBTagDisc7", &jetBTagDisc7, "jetBTagDisc7/f");  
	_inputTree->Branch("jetBTagDisc8", &jetBTagDisc8, "jetBTagDisc8/f");  
	_inputTree->Branch("bjetBTagDisc1", &bjetBTagDisc1, "bjetBTagDisc1/f");  
	_inputTree->Branch("bjetBTagDisc2", &bjetBTagDisc2, "bjetBTagDisc2/f");  
	_inputTree->Branch("bjetBTagDisc3", &bjetBTagDisc3, "bjetBTagDisc3/f");  
	_inputTree->Branch("bjetBTagDisc4", &bjetBTagDisc4, "bjetBTagDisc4/f");  
	_inputTree->Branch("bjetBTagDisc5", &bjetBTagDisc5, "bjetBTagDisc5/f");  
	_inputTree->Branch("bjetBTagDisc6", &bjetBTagDisc6, "bjetBTagDisc6/f");  
	_inputTree->Branch("bjetBTagDisc7", &bjetBTagDisc7, "bjetBTagDisc7/f");  
	_inputTree->Branch("bjetBTagDisc8", &bjetBTagDisc8, "bjetBTagDisc8/f");  
	_inputTree->Branch("lightjetBTagDisc1", &lightjetBTagDisc1, "lightjetBTagDisc1/f");  
	_inputTree->Branch("lightjetBTagDisc2", &lightjetBTagDisc2, "lightjetBTagDisc2/f");  
	_inputTree->Branch("lightjetBTagDisc3", &lightjetBTagDisc3, "lightjetBTagDisc3/f");  
	_inputTree->Branch("met", &met, "met/f");    
	_inputTree->Branch("metPhi", &metPhi, "metPhi/f");    
	_inputTree->Branch("averageDeltaRjj", &averageDeltaRjj, "averageDeltaRjj/f"); 
	_inputTree->Branch("averageDeltaRbb", &averageDeltaRbb, "averageDeltaRbb/f"); 
	_inputTree->Branch("averageDeltaEtajj", &averageDeltaEtajj, "averageDeltaEtajj/f"); 
	_inputTree->Branch("averageDeltaEtabb", &averageDeltaEtabb, "averageDeltaEtabb/f"); 
	_inputTree->Branch("minDeltaRjj", &minDeltaRjj, "minDeltaRjj/f"); 
	_inputTree->Branch("minDeltaRbb", &minDeltaRbb, "minDeltaRbb/f"); 
	_inputTree->Branch("maxDeltaEtabb", &maxDeltaEtabb, "maxDeltaEtabb/f"); 
	_inputTree->Branch("maxDeltaEtajj", &maxDeltaEtajj, "maxDeltaEtajj/f"); 
	_inputTree->Branch("maxDeltaEtabj", &maxDeltaEtabj, "maxDeltaEtabj/f"); 
	_inputTree->Branch("minDeltaRbj", &minDeltaRbj, "minDeltaRbj/f"); 
	_inputTree->Branch("averageDeltaEtabj", &averageDeltaEtabj, "averageDeltaEtabj/f"); 
	_inputTree->Branch("averageDeltaRbj", &averageDeltaRbj, "averageDeltaRbj/f"); 
	_inputTree->Branch("minDeltaRMassjj", &minDeltaRMassjj, "minDeltaRMassjj/f"); 
	_inputTree->Branch("minDeltaRMassbb", &minDeltaRMassbb, "minDeltaRMassbb/f"); 
	_inputTree->Branch("minDeltaRMassbj", &minDeltaRMassbj, "minDeltaRMassbj/f"); 
	_inputTree->Branch("minDeltaRpTjj", &minDeltaRpTjj, "minDeltaRpTjj/f"); 
	_inputTree->Branch("minDeltaRpTbb", &minDeltaRpTbb, "minDeltaRpTbb/f"); 
	_inputTree->Branch("minDeltaRpTbj", &minDeltaRpTbj, "minDeltaRpTbj/f"); 
	_inputTree->Branch("maxPTmassjjj", &maxPTmassjjj, "maxPTmassjjj/f"); 
	_inputTree->Branch("maxPTmassjbb", &maxPTmassjbb, "maxPTmassjbb/f"); 
	_inputTree->Branch("H0", &H0, "H0/f");
	_inputTree->Branch("H1", &H1, "H1/f");
	_inputTree->Branch("H2", &H2, "H2/f");
	_inputTree->Branch("H3", &H3, "H3/f");
	_inputTree->Branch("H4", &H4, "H4/f");
	_inputTree->Branch("bH0", &bH0, "bH0/f");
	_inputTree->Branch("bH1", &bH1, "bH1/f");
	_inputTree->Branch("bH2", &bH2, "bH2/f");
	_inputTree->Branch("bH3", &bH3, "bH3/f");
	_inputTree->Branch("bH4", &bH4, "bH4/f");
	_inputTree->Branch("R1", &R1, "R1/f");
	_inputTree->Branch("R2", &R2, "R2/f");
	_inputTree->Branch("R3", &R3, "R3/f");
	_inputTree->Branch("R4", &R4, "R4/f");
	_inputTree->Branch("bR1", &bR1, "bR1/f");
	_inputTree->Branch("bR2", &bR2, "bR2/f");
	_inputTree->Branch("bR3", &bR3, "bR3/f");
	_inputTree->Branch("bR4", &bR4, "bR4/f");
	_inputTree->Branch("jetAverageMass", &jetAverageMass, "jetAverageMass/f");
	_inputTree->Branch("bjetAverageMass", &bjetAverageMass, "bjetAverageMass/f");
	_inputTree->Branch("bjetAverageMassSqr", &bjetAverageMassSqr, "bjetAverageMassSqr/f");
	_inputTree->Branch("jetHT", &jetHT, "jetHT/f");
	_inputTree->Branch("bjetHT", &bjetHT, "bjetHT/f");
	_inputTree->Branch("lightjetHT", &lightjetHT, "lightjetHT/f");
	_inputTree->Branch("jetNumber", &jetNumber, "jetNumber/i");
	_inputTree->Branch("bjetNumber", &bjetNumber, "bjetNumber/i");
	_inputTree->Branch("lightjetNumber", &lightjetNumber, "lightjetNumber/i");
	_inputTree->Branch("invMassZ1", &invMassZ1, "invMassZ1/f");
	_inputTree->Branch("invMassZ2", &invMassZ2, "invMassZ2/f");
	_inputTree->Branch("invMassH1", &invMassH1, "invMassH1/f");
	_inputTree->Branch("invMassH2", &invMassH2, "invMassH2/f");
	_inputTree->Branch("chi2Higgs", &chi2Higgs, "chi2Higgs/f");
	_inputTree->Branch("chi2HadW", &chi2HadW, "chi2HadW/f");
	_inputTree->Branch("chi2Z", &chi2Z, "chi2Z/f");
	_inputTree->Branch("chi2HiggsZ", &chi2HiggsZ, "chi2HiggsZ/f");
	_inputTree->Branch("invMassHiggsZ1", &invMassHiggsZ1, "invMassHiggsZ1/f");
	_inputTree->Branch("invMassHiggsZ2", &invMassHiggsZ2, "invMassHiggsZ2/f");
	_inputTree->Branch("PTH1", &PTH1, "PTH1/f");
	_inputTree->Branch("PTH2", &PTH2, "PTH2/f");


	_inputTree->Branch("centralityjl", &centralityjl, "centralityjl/f");
	_inputTree->Branch("centralityjb", &centralityjb, "centralityjb/f");
	_inputTree->Branch("aplanarity", &aplanarity, "aplanarity/f");
	_inputTree->Branch("sphericity", &sphericity, "sphericity/f");
	_inputTree->Branch("transSphericity", &transSphericity, "transSphericity/f");
	_inputTree->Branch("cValue", &cValue, "cValue/f");
	_inputTree->Branch("dValue", &dValue, "dValue/f");
	_inputTree->Branch("baplanarity", &baplanarity, "baplanarity/f");
	_inputTree->Branch("bsphericity", &bsphericity, "bsphericity/f");
	_inputTree->Branch("btransSphericity", &btransSphericity, "btransSphericity/f");
	_inputTree->Branch("bcValue", &bcValue, "bcValue/f");
	_inputTree->Branch("bdValue", &bdValue, "bdValue/f");

	_inputTree->Branch("weight", &weight, "weight/f");

	_inputTree->Branch("leptonPT1", &leptonPT1, "leptonPT1/f");
	_inputTree->Branch("muonPT1", &muonPT1, "muonPT1/f");
	_inputTree->Branch("elePT1", &elePT1, "elePT1/f");
	_inputTree->Branch("leptonEta1", &leptonEta1, "leptonEta1/f");
	_inputTree->Branch("muonEta1", &muonEta1, "muonEta1/f");
	_inputTree->Branch("eleEta1", &eleEta1, "eleEta1/f");
	_inputTree->Branch("leptonPT2", &leptonPT2, "leptonPT2/f");
	_inputTree->Branch("muonPT2", &muonPT2, "muonPT2/f");
	_inputTree->Branch("elePT2", &elePT2, "elePT2/f");
	_inputTree->Branch("leptonEta2", &leptonEta2, "leptonEta2/f");
	_inputTree->Branch("muonEta2", &muonEta2, "muonEta2/f");
	_inputTree->Branch("eleEta2", &eleEta2, "eleEta2/f");
	_inputTree->Branch("diElectronMass", &diElectronMass, "diElectronMass/f");
	_inputTree->Branch("diMuonMass", &diMuonMass, "diMuonMass/f");
	_inputTree->Branch("leptonHT", &leptonHT, "leptonHT/f");
	_inputTree->Branch("ST", &ST, "ST/f");
	_inputTree->Branch("leptonCharge1", &leptonCharge1, "leptonCharge1/f");
	_inputTree->Branch("leptonCharge2", &leptonCharge2, "leptonCharge2/f");

	_inputTree->Branch("bbjetHiggsMatched1", &bbjetHiggsMatched1, "bbjetHiggsMatched1/f");
	_inputTree->Branch("bbjetHiggsMatched2", &bbjetHiggsMatched2, "bbjetHiggsMatched2/f");
	_inputTree->Branch("bbjetHiggsMatched3", &bbjetHiggsMatched3, "bbjetHiggsMatched3/f");
	_inputTree->Branch("bbjetHiggsMatched4", &bbjetHiggsMatched4, "bbjetHiggsMatched4/f");
	_inputTree->Branch("bbjetHiggsMatched5", &bbjetHiggsMatched5, "bbjetHiggsMatched5/f");
	_inputTree->Branch("bbjetHiggsMatched6", &bbjetHiggsMatched6, "bbjetHiggsMatched6/f");

	_inputTree->Branch("bbjetHiggsMatcheddR1", &bbjetHiggsMatcheddR1, "bbjetHiggsMatcheddR1/f");
	_inputTree->Branch("bbjetHiggsMatcheddR2", &bbjetHiggsMatcheddR2, "bbjetHiggsMatcheddR2/f");
	_inputTree->Branch("bbjetHiggsMatcheddR3", &bbjetHiggsMatcheddR3, "bbjetHiggsMatcheddR3/f");
	_inputTree->Branch("bbjetHiggsMatcheddR4", &bbjetHiggsMatcheddR4, "bbjetHiggsMatcheddR4/f");
	_inputTree->Branch("bbjetHiggsMatcheddR5", &bbjetHiggsMatcheddR5, "bbjetHiggsMatcheddR5/f");
	_inputTree->Branch("bbjetHiggsMatcheddR6", &bbjetHiggsMatcheddR6, "bbjetHiggsMatcheddR6/f");

	_inputTree->Branch("bbjetMinChiHiggsIndex1", &bbjetMinChiHiggsIndex1, "bbjetMinChiHiggsIndex1/f");
	_inputTree->Branch("bbjetMinChiHiggsIndex2", &bbjetMinChiHiggsIndex2, "bbjetMinChiHiggsIndex2/f");
	_inputTree->Branch("bbjetMinChiHiggsIndex3", &bbjetMinChiHiggsIndex3, "bbjetMinChiHiggsIndex3/f");
	_inputTree->Branch("bbjetMinChiHiggsIndex4", &bbjetMinChiHiggsIndex4, "bbjetMinChiHiggsIndex4/f");
	_inputTree->Branch("bbjetMinChiHiggsIndex5", &bbjetMinChiHiggsIndex5, "bbjetMinChiHiggsIndex5/f");
	_inputTree->Branch("bbjetMinChiHiggsIndex6", &bbjetMinChiHiggsIndex6, "bbjetMinChiHiggsIndex6/f");

	_treeDirs = tmpDirs;
    }
};	
#endif
