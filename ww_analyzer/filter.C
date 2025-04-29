// **ROOT**
// Version: ROOT 6.32.04
// Description: Searches for events containing 2 muons and that pass
//              other search parameters.
//
// Created: 25 Feb 2025
// Last Modified: 28 Feb 2025
// Author: Pablo Pomares
// Email: pablo.pomaresv@alumno.buap.mx
//

#include <string>
#include <iostream>
#include <vector>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

// Function to calculate the radial distance of a muon to a jet
Float_t deltaR(Float_t eta1, Float_t eta2, Float_t phi1, Float_t phi2){
  Float_t Deta2 = (eta1 - eta2)*(eta1 - eta2);
  Float_t Dphi2 = (phi1 - phi2)*(phi1 - phi2);
  Float_t DR = TMath::Sqrt(Deta2 + Dphi2);
  
  return DR;
}

void filter(const std::string& file_dir, const std::string& output_name){
  TFile *f_old = TFile::Open(file_dir.c_str());
  TTree* t_old;
  f_old->GetObject("Events", t_old);

  const int nentries = t_old->GetEntries();

  // Deactivate all branches
  t_old->SetBranchStatus("*", 0);

  // Activate desired branches
  for (auto actBranchName : {"run", "event", "Muon_charge", "Muon_dxy", "Muon_dxyErr", "Muon_isGlobal",
                             "Muon_pt", "Muon_phi", "Muon_eta", "nMuon", "MET_phi",
                             "MET_pt", "MET_significance", "nJet", "Jet_btagCSVV2", "Jet_btagDeepB",
                             "Jet_pt", "Jet_phi", "Jet_eta"}){
    t_old->SetBranchStatus(actBranchName, 1);
  };

  // Entry selection
  UInt_t nMuon, nJet;
  Float_t Jet_btagDeepB[40], Jet_pt[40], Jet_phi[40], Jet_eta[40], MET_pt,
          Muon_pt[25], Muon_eta[25], Muon_phi[25];
  Bool_t Muon_isGlobal[25];
  Int_t Muon_charge[25];
//  Float_t* Jet_btagDeepB = nullptr;
//  Float_t* Jet_pt = nullptr;
//  Float_t* Jet_phi = nullptr;
//  Float_t* Jet_eta = nullptr;
//  Float_t MET_pt;
//  Float_t* Muon_pt = nullptr;
//  Float_t* Muon_eta = nullptr;
//  Float_t* Muon_phi = nullptr;
//  Bool_t* Muon_isGlobal = nullptr;

  t_old->SetBranchAddress("nMuon", &nMuon);
  t_old->SetBranchAddress("Muon_pt", &Muon_pt);
  t_old->SetBranchAddress("Muon_eta", &Muon_eta);
  t_old->SetBranchAddress("Muon_phi", &Muon_phi);
  t_old->SetBranchAddress("nJet", &nJet);
  t_old->SetBranchAddress("Jet_pt", &Jet_pt);
  t_old->SetBranchAddress("Jet_eta", &Jet_eta);
  t_old->SetBranchAddress("Jet_phi", &Jet_phi);
  t_old->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB);
  t_old->SetBranchAddress("Muon_isGlobal", &Muon_isGlobal);
  t_old->SetBranchAddress("Muon_charge", &Muon_charge);
  t_old->SetBranchAddress("MET_pt", &MET_pt);

  int num_passnMuon = 0;
  int num_passMETpt = 0;
  int num_passnJet = 0;
  int num_passAllGlobal = 0;
  int num_passtotal_charge = 0;
  int num_passdeltaR = 0;
  int num_passLeadingMuonpt = 0;
  int num_passbtag = 0;

  TFile newfile(output_name.c_str(), "recreate"); // Creates new file
  auto t_new = t_old->CloneTree(0);

  for (int i=0; i<nentries; i++){
    t_old->GetEntry(i);
    
    // Checks that only two muons in the event. If it fails the event is skipped.
    bool passnMuon = (nMuon == 2);
    if (!passnMuon){continue;} else {num_passnMuon += 1;};
    
    // Due to the WW production, a minimum missing transverse energy is requiered due to the
    // prescence of neutrinos.
    bool passMETpt = false;
    if (MET_pt > 55.0){
      passMETpt = true;
    };
    if (!passMETpt) {continue;} else {num_passMETpt += 1;};

    // For a better event selection, only global muons are used.
    bool passAllGlobal = true;
    for (UInt_t j=0; j<nMuon; j++){
      if (Muon_isGlobal[j] == 0){
        passAllGlobal = false;
      };
    };
    if (!passAllGlobal) {continue;} else {num_passAllGlobal += 1;};

    // Opposite charge of the muons is needed 
    int sumCharge = 0;
    for (UInt_t j=0; j<nMuon; j++){
      sumCharge += Muon_charge[j];
    };
    bool total_charge = sumCharge;
    if (total_charge) {continue;} else {num_passtotal_charge += 1;};

    // We select the events were the leading muon is greater than 25 GeV and the subleading muon is 
    // greater than 20 GeV
    std::sort(&Muon_pt[0], &Muon_pt[2], std::greater<Float_t>());
    float subLeadingMuons = Muon_pt[1];
    float leadingMuons = Muon_pt[0];
    bool passLeadingMuonpt = ((leadingMuons > 25.0) && (subLeadingMuons > 20.0));
    if (!passLeadingMuonpt){continue;} else {
      num_passLeadingMuonpt += 1;
    };

    // A medium b-tag is implemented to reduce tt background
    bool passbTag = true;
    if (nJet > 0){
      for (UInt_t j=0; j<nJet; j++){
        if (Jet_btagDeepB[j] > 0.5847){
          passbTag = false;
        };
      };
    };
    if (!passbTag) {continue;} else {num_passbtag += 1;};
    
    // If radial distance lesser than 0.4, jet is rejected
    if (nJet > 0){
      for (UInt_t j=0; j<nMuon; j++){
        for (UInt_t k=0; k<nJet; k++){
          Float_t dR = deltaR(Muon_eta[j], Jet_eta[k], Muon_phi[j], Jet_phi[k]);
          if (dR < 0.4){
            nJet -= 1;
          }
        }
      };
    };

    // Checks for the number of jets. Less or equal than one is requiered
    bool passnJet = (nJet < 2);
    if (!passnJet){continue;} else {num_passnJet += 1;};

    t_new->Fill();
  };

  t_new->Print();
  newfile.Write();

  cout << "Pass nMuon criteria: " << num_passnMuon << "\n";
  cout << "Pass MET_pt criteria: " << num_passMETpt << "\n";
  cout << "Pass AllGlobal criteria: " << num_passAllGlobal << "\n";
  cout << "Pass total charge criteria: " << num_passtotal_charge << "\n";
  cout << "Pass Muon pt criteria: " << num_passLeadingMuonpt << "\n";
  cout << "Pass b tag criteria: " << num_passbtag << "\n";
  cout << "Pass nJet criteria: " << num_passnJet << "\n";
}
