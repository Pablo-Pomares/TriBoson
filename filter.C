// **ROOT**
// Version: ROOT 6.32.04
// Description: Filter for events containing less than 4 muons and 
//              all jets that passes a btag
//
// Created: 30 Nov 2024
// Last Modified: 12 Jan 2025
// Author: Pablo Pomares
// Email: pablo.pomaresv@alumno.buap.mx
//

#include <string>
#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

void filter(const std::string& file_dir, const std::string& output_name){
  TFile *f_old = TFile::Open(file_dir.c_str());
  TTree* t_old;
  f_old->GetObject("Events", t_old);

  const int nentries = t_old->GetEntries();

  // Deactivate all branches
  t_old->SetBranchStatus("*", 0);

  // Activate desired branches
  for (auto actBranchName : {"run", "event", "Muon_charge", "Muon_dxy", "Muon_dxyErr", "Muon_isGlobal", 
                             "Muon_isTracker", "Muon_pt", "Muon_phi", "Muon_eta", "nMuon", "MET_phi",
                             "MET_pt", "MET_significance", "nJet", "Jet_btagCSVV2", "Jet_btagDeepB"}){
    t_old->SetBranchStatus(actBranchName, 1);
  };

  // Entry selection
  UInt_t nMuon, nJet;
  Float_t Jet_btagDeepB[20], MET_significance;
  Bool_t Muon_isGlobal[10];
  Int_t Muon_charge[10];
  t_old->SetBranchAddress("nMuon", &nMuon);
  t_old->SetBranchAddress("nJet", &nJet);
  t_old->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB);
  t_old->SetBranchAddress("Muon_isGlobal", &Muon_isGlobal);
  t_old->SetBranchAddress("Muon_charge", &Muon_charge);
  t_old->SetBranchAddress("MET_significance", &MET_significance);

  TFile newfile(output_name.c_str(), "recreate"); // Creates new file
  auto t_new = t_old->CloneTree(0);

  for (UInt_t i=0; i<nentries; i++){
    t_old->GetEntry(i);

    int sumCharge = 0;
    for (UInt_t j=0; j<nMuon; j++){
      sumCharge += Muon_charge[j];
    };
    bool chargeViolation = sumCharge;

    bool passbTag = true;
    for (UInt_t j=0; j<nJet; j++){
      if (Jet_btagDeepB[j] > 0.5847){
        passbTag = false;
      };
    };

    bool passnMuon = false;
    if (nMuon == 4){
      passnMuon = true;
    };

    bool passAllGlobal = true;
    for (UInt_t j=0; j<nMuon; j++){
      if (Muon_isGlobal[j] == 0){
        passAllGlobal = false;
      };
    };

    bool passMETsignificance = false;
    if (MET_significance > 0.1){
      passMETsignificance = true;
    };

    if (passAllGlobal && passbTag && passnMuon && !chargeViolation && passMETsignificance) {
      t_new->Fill();
    };

  };

  t_new->Print();
  newfile.Write();
}
