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

  TFile newfile(output_name.c_str(), "recreate"); // Creates new file
  auto t_new = t_old->CloneTree(0);

  for (UInt_t i=0; i<nentries; i++){
    t_old->GetEntry(i);

    // Checks that only two muons in the event. If it fails the event is skipped.
    bool passnMuon = (nMuon == 4);
    if (!passnMuon){continue;};

    // For a better event selection, only global muons are used.
    bool passAllGlobal = true;
    for (UInt_t j=0; j<nMuon; j++){
      if (Muon_isGlobal[j] == 0){
        passAllGlobal = false;
      };
    };
    if (!passAllGlobal) {continue;};

    // Opposite charge of the muons is needed 
    int charge_sum = 0;
    for (UInt_t j=0; j<nMuon; j++){
      charge_sum += Muon_charge[j];
    };
    bool total_charge = (charge_sum == 0);
    if (!total_charge) {continue;};

    // We select the events were the leading muon is greater than 25 GeV and the subleading muon is 
    // greater than 20 GeV
    std::sort(&Muon_pt[0], &Muon_pt[5], std::greater<Float_t>());
    float subLeadingMuonsSum = Muon_pt[0] + Muon_pt[1];
    float leadingMuonsSum = Muon_pt[2] + Muon_pt[3];
    bool passLeadingMuonpt = true;
    if ((leadingMuonsSum < 50.0) && (subLeadingMuonsSum < 20.0)){
      passLeadingMuonpt = false;
      //break;
    };
    if (!passLeadingMuonpt){continue;};

    // A medium b-tag is implemented to reduce tt background
    bool passbTag = true;
    if (nJet > 0){
      for (UInt_t j=0; j<nJet; j++){
        if (Jet_btagDeepB[j] > 0.5847){
          passbTag = false;
          //break;
        };
      };
    };
    if (!passbTag) {continue;};

    // A radial distance greater than 0.4 is desired
    bool passdeltaR = true;
    if (nJet > 0){
      for (UInt_t j=0; j<nMuon; j++){
        for (UInt_t k=0; k<nJet; k++){
          Float_t dR = deltaR(Muon_eta[j], Jet_eta[k], Muon_phi[j], Jet_phi[k]);
          if (dR < 0.4){
            nJet--;
          }
        }
      };
    };

    // Due to the WW production, a minimum missing transverse energy is requiered due to the
    // prescence of neutrinos.
    if (MET_pt < 20.0){
      continue;
    };

    t_new->Fill();
  };

  t_new->Print();
  newfile.Write();
}
