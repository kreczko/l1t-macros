#ifndef TL1EVENTCLASS_H
#define TL1EVENTCLASS_H

#include <string>
#include <memory>

#include <TMath.h>
#include <TVector2.h>

#include "TL1DataClass.h"

//#include "../160510_l1t-integration-v46p0/CMSSW_8_0_7/src/L1Trigger/L1TNtuples/interface/L1AnalysisL1CaloTowerDataFormat.h"
//#include "../160510_l1t-integration-v46p0/CMSSW_8_0_7/src/L1Trigger/L1TNtuples/interface/L1AnalysisRecoJetDataFormat.h"
//#include "../160510_l1t-integration-v46p0/CMSSW_8_0_7/src/L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetDataFormat.h"
//#include "../160510_l1t-integration-v46p0/CMSSW_8_0_7/src/L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
//#include "../160510_l1t-integration-v46p0/CMSSW_8_0_7/src/L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
//#include "../160510_l1t-integration-v46p0/CMSSW_8_0_7/src/L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetFilterDataFormat.h"
//#include "../160510_l1t-integration-v46p0/CMSSW_8_0_7/src/L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"

class TL1EventClass
{
    public:
        TL1EventClass(std::string inDir);
        bool Next();
        void GetEntry(int i);

        bool MuonFilter() const;
        bool JetFilter(const int & iJet) const;
        bool SumsFilter() const;

        void RecalculateVariables();
        void GetL1Sums();
        void GetL1Jets();
        void GetFPL1Mht();

        void RecalculatedL1Ett();
        void RecalculateRecoEtt();

        double fRecalcMht, fRecalcMhtPhi, fRecalcHtt;
        double fL1Met, fL1Mht, fL1Ett, fL1Htt, fL1MetPhi, fL1MhtPhi;
        std::vector<double> fL1JetEt, fL1JetPhi, fL1JetEta;

        bool fMhtPassFlag;

        double fRecalcL1Ett, fRecalcRecoEtt;
        
        double fFPL1Mht, fFPL1MhtPhi;

        std::shared_ptr<L1Analysis::L1AnalysisL1CaloTowerDataFormat>   fCaloTowers;
        std::shared_ptr<L1Analysis::L1AnalysisRecoJetDataFormat>       fJets;
        std::shared_ptr<L1Analysis::L1AnalysisRecoMetDataFormat>       fSums;
        std::shared_ptr<L1Analysis::L1AnalysisL1UpgradeDataFormat>     fUpgrade;
        std::shared_ptr<L1Analysis::L1AnalysisRecoVertexDataFormat>    fVertex;
        std::shared_ptr<L1Analysis::L1AnalysisRecoMetFilterDataFormat> fMetFilters;
        std::shared_ptr<L1Analysis::L1AnalysisRecoMuon2DataFormat>     fMuons;     

    private:
        std::shared_ptr<TL1DataClass> caloTower, jetReco, metFilterReco, muonReco, recoTree, upgrade;
        Long64_t fPos;
};

TL1EventClass::TL1EventClass(std::string inDir) :
    caloTower(new TL1DataClass("l1CaloTowerTree/L1CaloTowerTree",inDir)),
    jetReco(new TL1DataClass("l1JetRecoTree/JetRecoTree",inDir)),
    metFilterReco(new TL1DataClass("l1MetFilterRecoTree/MetFilterRecoTree",inDir)),
    muonReco(new TL1DataClass("l1MuonRecoTree/Muon2RecoTree",inDir)),
    recoTree(new TL1DataClass("l1RecoTree/RecoTree",inDir)),
    upgrade(new TL1DataClass("l1UpgradeTree/L1UpgradeTree",inDir)),
    fCaloTowers(new L1Analysis::L1AnalysisL1CaloTowerDataFormat()),
    fJets      (new L1Analysis::L1AnalysisRecoJetDataFormat()),
    fSums      (new L1Analysis::L1AnalysisRecoMetDataFormat()),
    fUpgrade   (new L1Analysis::L1AnalysisL1UpgradeDataFormat()),
    fVertex    (new L1Analysis::L1AnalysisRecoVertexDataFormat()),
    fMetFilters(new L1Analysis::L1AnalysisRecoMetFilterDataFormat()),
    fMuons     (new L1Analysis::L1AnalysisRecoMuon2DataFormat()),
    fPos(0)
{
    caloTower->SetBranchAddress("L1CaloTower", &fCaloTowers);
    jetReco->SetBranchAddress("Jet", &fJets);
    jetReco->SetBranchAddress("Sums", &fSums);
    metFilterReco->SetBranchAddress("MetFilters", &fMetFilters);
    muonReco->SetBranchAddress("Muon", &fMuons);
    recoTree->SetBranchAddress("Vertex", &fVertex);
    upgrade->SetBranchAddress("L1Upgrade", &fUpgrade);

    fMhtPassFlag = true;
}

bool TL1EventClass::Next()
{
    if( fPos >= caloTower->GetEntries()-1 ) return false;
    caloTower->LoadTree(fPos);
    jetReco->LoadTree(fPos);
    metFilterReco->LoadTree(fPos);
    muonReco->LoadTree(fPos);
    recoTree->LoadTree(fPos);
    upgrade->LoadTree(fPos);

    caloTower->GetEntry(fPos);
    jetReco->GetEntry(fPos);
    metFilterReco->GetEntry(fPos);
    muonReco->GetEntry(fPos);
    recoTree->GetEntry(fPos);
    upgrade->GetEntry(fPos);
    
    ++fPos;

    return true;
}

void TL1EventClass::GetEntry(int i)
{
    caloTower->LoadTree(i);
    jetReco->LoadTree(i);
    metFilterReco->LoadTree(i);
    muonReco->LoadTree(i);
    recoTree->LoadTree(i);
    upgrade->LoadTree(i);

    caloTower->GetEntry(i);
    jetReco->GetEntry(i);
    metFilterReco->GetEntry(i);
    muonReco->GetEntry(i);
    recoTree->GetEntry(i);
    upgrade->GetEntry(i);
}

bool TL1EventClass::MuonFilter() const
{
    bool pass(false);
    for(unsigned jMuon=0; jMuon<fMuons->nMuons; ++jMuon)
    {
        double pt = fMuons->pt[jMuon];
        double iso = fMuons->iso[jMuon];
        double isLoose = fMuons->isLooseMuon[jMuon];
        if( pt>=20.0 && iso<=0.1 && isLoose==1.0 )
            pass = true;
    }
    return pass;
}

bool TL1EventClass::JetFilter(const int & iJet) const
{
    double eta = fJets->eta[iJet];
    double nhef = fJets->nhef[iJet];
    double pef = fJets->pef[iJet];
    double mef = fJets->mef[iJet];
    double chef = fJets->chef[iJet];
    double eef = fJets->eef[iJet];
    int chMult = fJets->chMult[iJet];
    int nhMult = fJets->nhMult[iJet];
    int phMult = fJets->phMult[iJet];
    int elMult = fJets->elMult[iJet];
    int muMult = fJets->muMult[iJet];
    bool jetPass;
    if (abs(eta)<=3.00 && nhef<0.90 && pef<0.90 && (chMult+nhMult+phMult+elMult+muMult)>1 && mef<0.8 
            && (abs(eta)>2.4||(abs(eta)<=2.4 && chef>0 && (chMult+elMult+muMult)>0 && eef<0.90))){
        jetPass = true;}//this jet has passed
    else{jetPass = false;}//this jet has failed, or is outside the central zone           
    return (jetPass && muMult==0);
}

bool TL1EventClass::SumsFilter() const
{
    return fMetFilters->hbheNoiseFilter;
}

void TL1EventClass::RecalculateVariables()
{
    std::shared_ptr<TVector2> mht(new TVector2(0.,0.));
    double jetEt(0.0);
    unsigned jetCount(0);
    for(int iJet=0; iJet<fJets->nJets; ++iJet)
    {
        bool passJetFilter = this->JetFilter(iJet);
        if( fJets->etCorr[iJet] < 30.0 ) continue;
        if( !passJetFilter )
        {
            fMhtPassFlag = false;
            break;
        }
        std::shared_ptr<TVector2> jet(new TVector2(0.,0.));
        jet->SetMagPhi(fJets->etCorr[iJet], fJets->phi[iJet]);
        jetEt += fJets->etCorr[iJet];
        ++jetCount;

        *(mht.get()) -= *(jet.get());
    }
    fRecalcMht = mht->Mod();
    fRecalcMhtPhi = mht->Phi();
    fRecalcHtt = jetEt;
}

void TL1EventClass::RecalculatedL1Ett()
{
    int iEtaMax = 28;
    Double_t sumEt = 0.0;
    for(unsigned jTower=0; jTower<fCaloTowers->nTower; ++jTower)
    {
        int ieta = fCaloTowers->ieta[jTower];
        int iet  = fCaloTowers->iet[jTower];
        if( abs(ieta) <= iEtaMax )
            sumEt += 0.5 * (double)iet;
    }
    fRecalcL1Ett = sumEt;
}

void TL1EventClass::RecalculateRecoEtt()
{
    double jetEt(0.0);
    for(int iJet=0; iJet<fJets->nJets; ++iJet)
        if( fJets->etCorr[iJet] >= 30.0 )
            jetEt += fJets->etCorr[iJet];
    fRecalcRecoEtt = jetEt;
}

void TL1EventClass::GetL1Sums()
{
    for(int iter=0; iter<fUpgrade->nSums; ++iter)
    {
        double et = fUpgrade->sumEt[iter];
        double phi = fUpgrade->sumPhi[iter];
        if( fUpgrade->sumBx[iter] == 0 )
        {
            if(fUpgrade->sumType[iter] == L1Analysis::kTotalEt)   fL1Ett = et;
            if(fUpgrade->sumType[iter] == L1Analysis::kTotalHt)   fL1Htt = et;
            if(fUpgrade->sumType[iter] == L1Analysis::kMissingEt)
            {
                fL1Met = et;
                fL1MetPhi = phi;
            }
            if(fUpgrade->sumType[iter] == L1Analysis::kMissingHt)
            {
                fL1Mht = et;
                fL1MhtPhi = phi;
            }
        }
    }
}

void TL1EventClass::GetL1Jets()
{
    fL1JetEt.clear();
    fL1JetEta.clear();
    fL1JetPhi.clear();
    for(unsigned iJet=0; iJet<fUpgrade->nJets; ++iJet)
    {
        if( fUpgrade->jetBx[iJet] == 0 )
        {
            fL1JetEt.push_back(fUpgrade->jetEt[iJet]);
            fL1JetEta.push_back(fUpgrade->jetEta[iJet]);
            fL1JetPhi.push_back(fUpgrade->jetPhi[iJet]);
        }
    }
}

void TL1EventClass::GetFPL1Mht()
{
    double sumEx(0.0), sumEy(0.0);
    for(unsigned iJet=0; iJet<fL1JetEt.size(); ++iJet)
    {
        if( abs(fL1JetEta[iJet]) >= 3.0 ) continue;
        if( fL1JetEt[iJet] < 30.0 ) continue;
        sumEx += fL1JetEt[iJet]*TMath::Cos(fL1JetPhi[iJet]);
        sumEy += fL1JetEt[iJet]*TMath::Sin(fL1JetPhi[iJet]);
    }
    fFPL1Mht = TMath::Sqrt(sumEx*sumEx + sumEy*sumEy);
    fFPL1MhtPhi = TMath::Pi()+TMath::ATan(sumEy/sumEx);
    if( sumEx > 0.0 ) fFPL1MhtPhi = TMath::Pi() - fFPL1MhtPhi;
}

#endif
