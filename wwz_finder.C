// **ROOT**
// Version: ROOT 6.32.04
// Description: For now, it identifies events with 1, 2 or no Z bosons
//
// Created: 08 Dec 2024
// Last Modified: 10 Jan 2025
// Author: Pablo Pomares
// Email: pablo.pomaresv@alumno.buap.mx

// Function for the invariant mass

#include "lester_mt2_bisect.h"

const Float_t muon_mass = 0.1056583; // MeV/c2

Double_t inv_mass(Float_t pt1, Float_t pt2, Float_t phi1, Float_t phi2, Float_t eta1, Float_t eta2) {

  Double_t eta_diff = eta1 - eta2;
  Double_t phi_diff = phi1 - phi2;
  Double_t pt_prod = 2*pt1*pt2;

  Double_t m2 = pt_prod*(TMath::CosH(eta_diff) - TMath::Cos(phi_diff));
  Double_t m = TMath::Sqrt(m2);
  
  return m;
}


// Find the number of Z bosons on a event and its mass.
// If >1 it returns (for now) a single mass. However, this is unimportant because event is discarded.
class Z_finder{
  public:
    UInt_t num_z;
    Double_t masses[2];
    std::array<int, 2> z_muon_index;

    Z_finder(Float_t muon_pt[4], Float_t muon_phi[4], Float_t muon_eta[4], Int_t muon_charge[4])
      : num_z(0), masses{0.0, 0.0}, z_muon_index{0, 0}
    {
      for (int i=0; i<4; i++){
        int local_z = 0; // Number of Z candidates for that share a muon
        int z_local_index[3] {-1, -1, -1};
        Double_t local_masses[2] = {0., 0.};

        for (int j=i+1; j<4; j++){
          // Muon properties
          Float_t pt1 = muon_pt[i];
          Float_t pt2 = muon_pt[j];
          Float_t phi1 = muon_phi[i];
          Float_t phi2 = muon_phi[j];
          Float_t eta1 = muon_eta[i];
          Float_t eta2 = muon_eta[j];
          Int_t q1 = muon_charge[i];
          Int_t q2 = muon_charge[j];

          // Check for different charge
          bool same_charge = q1 + q2;
          z_local_index[0] = i;
          z_muon_index[0] = i;

          if (!same_charge) {
            Double_t m = inv_mass(pt1, pt2, phi1, phi2, eta1, eta2);
            if (m > 81.2 && m < 101.2){
              local_masses[local_z] = m;
              z_local_index[local_z+1] = j;
              local_z++;
            };
          };
        };

        // If more there is more than one local Z choose the one with greater
        // phi difference
        if (local_z == 2){
          float phi1 = TMath::Abs(muon_phi[z_local_index[0]]);
          float phi2 = TMath::Abs(muon_phi[z_local_index[1]]);
          float phi3 = TMath::Abs(muon_phi[z_local_index[2]]);

          float diff1_2 = phi1 - phi2;
          float diff1_3 = phi1 - phi3;

          if (diff1_2 > diff1_3){
            masses[num_z] = local_masses[0];
            z_muon_index[1] = z_local_index[1];
          }
          else {
            masses[num_z] = local_masses[1];
            z_muon_index[1] = z_local_index[2];
          };
        };

        if (local_z) {num_z++;};
      }
    }
};

class Vec_comp{
  public:
    Double_t px;
    Double_t py;

    Vec_comp(Double_t pt, Double_t phi){
      px = pt*TMath::Cos(phi);
      py = pt*TMath::Sin(phi);
    }
};

//Muons that are not part of the Z boson
std::tuple<int, int> free_muons(std::array<int, 2> arr){
  int free_index[2];
  int nfree = 0;
  bool isfree;

  for (int i=0; i<4; i++){
    isfree = true;
    for (int j=0; j<2; j++){
      if (i==arr[j]) {isfree = false;}
    }

    if (isfree){
      free_index[nfree] = i;
      nfree++;
    }
  }
  std::tuple<int, int> free_muons = {free_index[0], free_index[1]};

  return free_muons;
};

void wwz_finder(const std::string& file_name){
  TFile *f = TFile::Open(file_name.c_str(), "READ");
  TTree* t;
  f->GetObject("Events", t);
  
  const int nentries = t->GetEntries();

  UInt_t nMuon, run;
  ULong64_t event;
  Int_t Muon_charge[4];
  Float_t Muon_pt[4], Muon_phi[4], Muon_dxy[4], 
          Muon_eta[4];
  Float_t MET_pt, MET_phi;

  t->SetBranchAddress("run", &run);
  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("nMuon", &nMuon);
  t->SetBranchAddress("Muon_charge", &Muon_charge);
  t->SetBranchAddress("Muon_pt", &Muon_pt);
  t->SetBranchAddress("Muon_phi", &Muon_phi);
  t->SetBranchAddress("Muon_dxy", &Muon_dxy);
  t->SetBranchAddress("MET_pt", &MET_pt);
  t->SetBranchAddress("MET_phi", &MET_phi);

  int single_z = 0;
  int double_z = 0;
  int no_z = 0;

  for (int i=0; i<nentries; i++){
    t->GetEntry(i);
    Z_finder Info(Muon_pt, Muon_phi, Muon_eta, Muon_charge);
    int num_z = Info.num_z;
    std::array<int, 2> index = Info.z_muon_index;
    if (num_z == 1){
      single_z++;

      auto freemuons = free_muons(index); // free muons
      int f1 = std::get<0>(freemuons);
      int f2 = std::get<1>(freemuons);
      float pt_visA = Muon_pt[f1];
      float phi_visA = Muon_phi[f1];
      float pt_visB = Muon_pt[f2];
      float phi_visB = Muon_phi[f2];

      Vec_comp p_visA(pt_visA, phi_visA);
      Vec_comp p_visB(pt_visB, phi_visB);
      Vec_comp p_miss(MET_pt, MET_phi);
      
      asymm_mt2_lester_bisect::disableCopyrightMessage();

      double mT2 = asymm_mt2_lester_bisect::get_mT2(muon_mass, p_visA.px, p_visA.py,
                                                    muon_mass, p_visB.px, p_visB.py,
                                                    p_miss.px, p_miss.py, 0, 0, 0);
      std::cout << "Con valor mT2: " << mT2 << std::endl;
    }
    else if (num_z == 2){double_z++;}
    else {no_z++;}
  }
  std::cout << "Se tienen " << single_z << " eventos con 1 Z." << std::endl;
  std::cout << "Se tienen " << double_z << " eventos con 2 Z." << std::endl;
  std::cout << "Se tienen " << no_z << " eventos con ningún Z." << std::endl;
}
