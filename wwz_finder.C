// **ROOT**
// Version: ROOT 6.32.04
// Description: For now, it identifies events with 1, 2 or no Z bosons
//
// Created: 08 Dec 2024
// Last Modified: 08 Dec 2024
// Author: Pablo Pomares
// Email: pablo.pomaresv@alumno.buap.mx

// Function for the invariant mass

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
std::tuple<UInt_t, Double_t> z_finder(Float_t muon_pt[4], Float_t muon_phi[4], Float_t muon_eta[4], Int_t muon_charge[4]){
  UInt_t num_z = 0;
  Double_t masses[2];
  for (int i=0; i<4; i++){
    int local_z = 0; // if >1 check for the greatest eta diff
    int z_index[3];
    Double_t local_masses[2] = {0., 0.};

    for (int j=i+1; j<4; j++){
      Float_t pt1 = muon_pt[i];
      Float_t pt2 = muon_pt[j];
      Float_t phi1 = muon_phi[i];
      Float_t phi2 = muon_phi[j];
      Float_t eta1 = muon_eta[i];
      Float_t eta2 = muon_eta[j];
      Int_t q1 = muon_charge[i];
      Int_t q2 = muon_charge[j];

      bool same_charge = q1 + q2;
      z_index[0] = i;

      if (!same_charge) {
        Double_t m = inv_mass(pt1, pt2, phi1, phi2, eta1, eta2);
        if (m > 81.2 && m < 101.2){
          local_masses[local_z] = m;
          z_index[local_z+1] = j;
          local_z++;
          };
        }
    };

    if (local_z == 2){
      float phi1 = TMath::Abs(muon_phi[z_index[0]]);
      float phi2 = TMath::Abs(muon_phi[z_index[1]]);
      float phi3 = TMath::Abs(muon_phi[z_index[2]]);

      float diff1_2 = phi1 - phi2;
      float diff1_3 = phi1 - phi3;

      if (diff1_2 > diff1_3){
        masses[num_z] = local_masses[0];
      }
      else {
        masses[num_z] = local_masses[1];
      };
    };

    if (local_z) {num_z++;};
  }

  std::tuple<UInt_t, Double_t> result = {num_z, masses[0]};

  return result;
}

void wwz_finder(){
  TFile *f = TFile::Open("prueba1.root", "READ");
  TTree* t;
  f->GetObject("Events", t);
  
  const int nentries = t->GetEntries();

  UInt_t nMuon, run;
  ULong64_t event;
  Int_t Muon_charge[4];
  Float_t Muon_pt[4], Muon_phi[4], Muon_dxy[4], 
          Muon_eta[4], MET_pt[4], MET_phi[4];

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
    auto [num_z, mass] = z_finder(Muon_pt, Muon_phi, Muon_eta, Muon_charge);
    if (num_z == 1){single_z++;}
    else if (num_z == 2){double_z++;}
    else {no_z++;}
  }
  std::cout << "Se tienen " << single_z << " eventos con 1 Z." << std::endl;
  std::cout << "Se tienen " << double_z << " eventos con 2 Z." << std::endl;
  std::cout << "Se tienen " << no_z << " eventos con ningún Z." << std::endl;
}
