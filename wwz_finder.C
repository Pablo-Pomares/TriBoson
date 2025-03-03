//***ROOT*******************************************************************
//*Version: ROOT 6.32.04                                                   *
//*Description: Finds events that have 1 Z and calculates mT2 and mll from *
//*             the remainig 2 muons                                       *
//*                                                                        *
//*Created: 08 Dec 2024                                                    *
//*Last Modified: 22 Feb 2025                                              *
//*Author: Pablo Pomares                                                   *
//*Email: pablo.pomaresv@alumno.buap.mx                                    *
//**************************************************************************


#include <vector>
#include <iostream>
#include "lester_mt2_bisect.h"
#include "utils.h"

using namespace std;

void wwz_finder(const std::string& file_dir){
  TFile *f = TFile::Open(file_dir.c_str(), "READ");
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
  t->SetBranchAddress("Muon_eta", &Muon_eta);
  t->SetBranchAddress("MET_pt", &MET_pt);
  t->SetBranchAddress("MET_phi", &MET_phi);

  int single_z_noWW = 0;
  int single_z_posWW = 0;
  int no_z = 0;
  int posWWZ = 0;

  vector<float> pt_visAVec;
  vector<float> phi_visAVec;
  vector<float> eta_visAVec;
  vector<float> pt_visBVec;
  vector<float> phi_visBVec;
  vector<float> eta_visBVec;
  vector<float> pt_METVec;
  vector<float> phi_METVec;
  vector<float> mT2Vec;
  vector<float> mllVec;
  vector<int> runVec;
  vector<long> eventVec;

  vector<int> wwzrunVec;
  vector<long> wwzeventVec;

  for (int i=0; i<nentries; i++){
    t->GetEntry(i);
    Z_finder Info(Muon_pt, Muon_phi, Muon_eta, Muon_charge);
    bool hasZ = Info.hasZ;
    std::array<int, 2> index = Info.z_muon_index;
    if (hasZ){

      auto freemuons = free_muons(index); // free muons
      int f1 = std::get<0>(freemuons);
      int f2 = std::get<1>(freemuons);
      float pt_visA = Muon_pt[f1];
      float phi_visA = Muon_phi[f1];
      float eta_visA = Muon_eta[f1];
      float pt_visB = Muon_pt[f2];
      float phi_visB = Muon_phi[f2];
      float eta_visB = Muon_eta[f2];

      if ((pt_visA > 25.0 || pt_visB > 25.0) && (pt_visA > 20.0 && pt_visB > 20.0)){

        single_z_posWW += 1;
        Vec_comp p_visA(pt_visA, phi_visA);
        Vec_comp p_visB(pt_visB, phi_visB);
        Vec_comp p_miss(MET_pt, MET_phi);
        
        asymm_mt2_lester_bisect::disableCopyrightMessage();

        double mT2 = asymm_mt2_lester_bisect::get_mT2(muon_mass, p_visA.px, p_visA.py,
                                                      muon_mass, p_visB.px, p_visB.py,
                                                      p_miss.px, p_miss.py, 0, 0, 0);
        Double_t mll = inv_mass(pt_visA, pt_visB, phi_visA, phi_visB, eta_visA, eta_visB);
        pt_visAVec.push_back(pt_visA);
        phi_visAVec.push_back(phi_visA);
        eta_visAVec.push_back(eta_visA);
        pt_visBVec.push_back(pt_visB);
        phi_visBVec.push_back(phi_visB);
        eta_visBVec.push_back(eta_visB);
        pt_METVec.push_back(MET_pt);
        phi_METVec.push_back(MET_phi);
        mT2Vec.push_back(mT2);
        mllVec.push_back(mll);
        runVec.push_back(run);
        eventVec.push_back(event);

        bool passmT2 = (mT2 >= 60) && (mT2 <= 80);
        bool passMET_pt = (MET_pt >= 60) && (MET_pt <= 80);
        if ((mll > 100.0) && passMET_pt && passMET_pt) {
          posWWZ += 1;
          wwzrunVec.push_back(run);
          wwzeventVec.push_back(event);
        }
      }
      else {
        single_z_noWW += 1;
      }
    }
    else {no_z++;}
  }
  std::cout << "Se tiene(n) " << posWWZ << " evento(s) con WWZ" << std::endl;
  if (posWWZ){
    for (int i=0; i<wwzeventVec.size(); i++){
      cout << "   En " << "run: " << wwzrunVec[i] << " event: " << wwzeventVec[i] << "\n"; 
    }
  }
  std::cout << "Se tiene(n) " << single_z_posWW << " evento(s) con 1 Z y posibilidad de WW." << std::endl;
  //std::cout << "Se tiene(n) " << single_z_noWW << " evento(s) con 1 Z sin WW." << std::endl;
  //std::cout << "Se tienen " << no_z << " eventos con ningún Z." << std::endl;
  ofstream outfile("posWWZdata.csv", ofstream::out);
  save_file(outfile, runVec, eventVec, pt_visAVec, phi_visAVec, eta_visAVec, pt_visBVec, phi_visBVec, eta_visBVec, pt_METVec, phi_METVec, mT2Vec, mllVec);
}
