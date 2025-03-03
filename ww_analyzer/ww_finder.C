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

const double pi = TMath::Pi();
const double z_mass = 91.1880; // GeV/c2 pm 0.0020 per PDG(2025) 

// Checks if the perpendicular component of p^miss_T with respect to a given muon is greater or 
// equal to 20 GeV if it is requiered.
bool pt_miss_proj_test(Float_t MET_pt, Float_t MET_phi, Float_t Muon_phi){
  float pt_miss_proj = 0.0;
  double min_dphi = pi/2;

  // Measures difference in angles between a given muon and the missing transverse energy.
  float dphi = TMath::Abs(MET_phi - Muon_phi);
  if (dphi > TMath::Pi()) {
      dphi = 2.0 * TMath::Pi() - dphi;
    }
  if (dphi < min_dphi){
    pt_miss_proj = MET_pt*TMath::Sin(dphi);
  };

  bool pass_pt_miss_proj = (pt_miss_proj > 20.0);

  return pass_pt_miss_proj;
}

class Analyzer{
  private:
  Float_t diff_to_Z = 0.0; // to filter Drell Yan background

  public:
  // Conditions WW events must follow
  bool pass_mll = true;
  bool pass_pt_miss_proj = true;
  bool pass_notdy = true;

  // All in one, for convinience
  bool pass_all = (pass_mll && pass_pt_miss_proj && pass_notdy); 

  Analyzer(Float_t Muon_pt[2], Float_t Muon_phi[2], Float_t Muon_eta[2], Float_t MET_pt, Float_t MET_phi){
    Float_t mll = inv_mass(Muon_pt[0], Muon_pt[1], Muon_phi[0], Muon_phi[1], Muon_eta[0], Muon_phi[1]);
    if (mll > 40) {pass_mll = false;};

    diff_to_Z = TMath::Abs(mll - z_mass);
    if (diff_to_Z < 15) {pass_notdy = false;};

    for (int i=0; i<2; i++){
      if (!pass_pt_miss_proj) {continue;};
      pass_pt_miss_proj = pt_miss_proj_test(MET_pt, MET_phi, Muon_phi[i]);
    }
  }
};

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
  int no_ww = 0;
  int hell_yeah = 0;

  for (int i=0; i<nentries; i++) {
    t->GetEntry(i);
    Analyzer Info(Muon_pt, Muon_phi, Muon_eta, MET_pt, MET_phi);

    if (Info.pass_all) {
      hell_yeah += 1;
    }
  }
  cout << hell_yeah << "\n";
}
