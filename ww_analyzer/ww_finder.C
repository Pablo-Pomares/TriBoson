// **ROOT**
// Version: ROOT 6.32.04
// Description: Searches for events containing 2 muons and that pass
//              other search parameters.
//
// Created: 28 Feb 2025
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
#include "../utils.h"

using namespace std;

void ww_finder(const std::string& file_dir) {
  TFile *f = TFile::Open(file_dir.c_str(), "READ");
  TTree* t;
  f->GetObject("Events", t);
  
  const int nentries = t->GetEntries();

  // Entry selection
  UInt_t nMuon, run;
  ULong64_t event;
  Int_t Muon_charge[2];
  Float_t Muon_pt[2], Muon_phi[2], Muon_dxy[2], 
          Muon_eta[2];
  Float_t MET_pt, MET_phi;

  t->SetBranchAddress("run", &run);
  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("nMuon", &nMuon);
  t->SetBranchAddress("Muon_charge", &Muon_charge);
  t->SetBranchAddress("Muon_pt", &Muon_pt);
  t->SetBranchAddress("Muon_phi", &Muon_phi);
  t->SetBranchAddress("Muon_dxy", &Muon_dxy);
  t->SetBranchAddress("Muon_eta", &Muon_eta);
  t->SetBranchAddress("MET_pt", &MET_pt);
  t->SetBranchAddress("MET_phi", &MET_phi);

  int ww = 0;

  for (int i=0; i<nentries; i++) {
    t->GetEntry(i);
    WW_Analyzer Info(Muon_pt, Muon_phi, Muon_eta, MET_pt, MET_phi);
    if (Info.pass_all) {
      ww += 1;
    }
  }
  double ww_cross_section = 12*ww/luminosity;
  cout << "Del total de eventos: " << nentries << "\n";
  cout << "Se tienen " << ww << " candidatos a WW" << "\n";
  cout << "Lo que nos da una luminosidad de: " << ww_cross_section << "\n";
}
