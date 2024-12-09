// **ROOT**
// Version: ROOT 6.32.04
// Description: Filter for events containing less than 4 muons and 
//              all jets that passes a btag
//
// Created: 30 Nov 2024
// Last Modified: 01 Dec 2024
// Author: Pablo Pomares
// Email: pablo.pomaresv@alumno.buap.mx
//

void filter(){
  TFile *f_old = TFile::Open("root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/127C2975-1B1C-A046-AABF-62B77E757A86.root");
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
  Float_t Jet_btagDeepB[20];
  Bool_t Muon_isGlobal[10];
  Int_t Muon_charge[10];
  t_old->SetBranchAddress("nMuon", &nMuon);
  t_old->SetBranchAddress("nJet", &nJet);
  t_old->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB);
  t_old->SetBranchAddress("Muon_isGlobal", &Muon_isGlobal);
  t_old->SetBranchAddress("Muon_charge", &Muon_charge);

  TFile newfile("prueba1.root", "recreate");
  auto t_new = t_old->CloneTree(0);

  for (int i=0; i<nentries; i++){
    t_old->GetEntry(i);

    int sumCharge = 0;
    for (int j=0; j<nMuon; j++){
      sumCharge += Muon_charge[j];
    };
    bool chargeViolation = sumCharge;

    bool passbTag = true;
    for (int j=0; j<nJet; j++){
      if (Jet_btagDeepB[j] > 0.5847){
        passbTag = false;
      };
    };

    bool passnMuon = false;
    if (nMuon == 4){
      passnMuon = true;
    };

    bool passAllGlobal = true;
    for (int j=0; j<nMuon; j++){
      if (Muon_isGlobal[j] == 0){
        passAllGlobal = false;
      };
    };
    if (passAllGlobal && passbTag && passnMuon && !chargeViolation) {
      t_new->Fill();
    };

  };

  t_new->Print();
  newfile.Write();
}
