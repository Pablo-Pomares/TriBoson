//***ROOT******************************************************************
//*Version: ROOT 6.32.04                                                  *
//*Description: Finds events that have 1 Z and calculates mT2 and mll from*
//*             the remainig 2 muons                                      *
//*                                                                       *
//*Created: 08 Dec 2024                                                   *
//*Last Modified: 07 Feb 2025                                             *
//*Author: Pablo Pomares                                                  *
//*Email: pablo.pomaresv@alumno.buap.mx                                   *
//*************************************************************************


#include <vector>
#include <iostream>
#include "lester_mt2_bisect.h"
using namespace std;

const Float_t muon_mass = 0.1056583; // MeV/c2

//Saves relevant data to csv file
void save_file(ofstream& out, vector<int> run, vector<long> event, vector<float> pt_visA, vector<float> phi_visA, vector<float> eta_visA,
               vector<float> pt_visB, vector<float> phi_visB, vector<float> eta_visB,
               vector<float> pt_MET, vector<float> phi_MET, vector<float> mT2, vector<float> mll){
  int vsize = pt_visA.size();
  if (out.is_open()){
    out << "run" << "," << "event" << "," << "pt_visA" << "," << "phi_visA" << "," << "eta_visA" << "," << "pt_visB" << "," << "phi_visB" << "," << "eta_visB" << "," << "pt_MET" << "," << "phi_MET" << "," << "mT2" << ","<< "mll" << "\n";  
    for (int i = 0; i < vsize; i++){
      out << run[i] << "," << event[i] << "," << pt_visA[i] << "," << phi_visA[i] << "," << eta_visA[i] << "," << pt_visB[i] << "," <<  phi_visB[i] << "," << eta_visB[i] << "," << pt_MET[i] << "," << phi_MET[i] << "," <<  mT2[i] << "," <<  mll[i] << "\n";
    }
    out.close();
    cout << "File written succesfully!\n";
  }
  else {
    cout << "Error opening file\n";
  }
}


// Function for the invariant mass
Double_t inv_mass(Float_t pt1, Float_t pt2, Float_t phi1, Float_t phi2, Float_t eta1, Float_t eta2) {

  Double_t eta_diff = eta1 - eta2;
  Double_t phi_diff = phi1 - phi2;
  Double_t pt_prod = 2*pt1*pt2;

  Double_t m2 = pt_prod*(TMath::CosH(eta_diff) - TMath::Cos(phi_diff));
  Double_t m = TMath::Sqrt(m2);
  
  return m;
}


// Find the Z bosons on a event and its mass.
// If it finds more than one candidate it selects the more probable one. 
// i. e., the one with the greater absolute difference of phi.
class Z_finder{
  private:
    int z_temp_index[4] = {-1, -1, -1, -1}; // If more than one particle found in Z range the first pair                                             //
                                            // would correspond to the index of first founded Z and the last pair for the second one.
    int z_temp_count = 0;
    UInt_t num_z;

  public:
    bool hasZ = false;
    Double_t massZ = 0.0;
    std::array<int, 2> z_muon_index;
    Double_t pi = TMath::Pi();

    Z_finder(Float_t muon_pt[4], Float_t muon_phi[4], Float_t muon_eta[4], Int_t muon_charge[4])
      : num_z(0), z_muon_index{-1, -1}
    {
      for (int i=0; i<4; i++){
        int local_z = 0; // Number of Z candidates for that share a muon
        int z_local_index[3] = {-1, -1, -1};
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

          if (!same_charge) {
            Double_t m = inv_mass(pt1, pt2, phi1, phi2, eta1, eta2);
            if (m > 81.2 && m < 101.2){
              local_masses[local_z] = m;
              z_local_index[local_z+1] = j;
              local_z++;
            };
          };
        };

        if (local_z){
          z_temp_index[z_temp_count] = z_local_index[0];
          num_z += 1;
          z_temp_count += 1;
        }

        // If more there is more than one local Z choose the one whose
        // phi difference is closest to pi
        if (local_z == 2){
          float phi1 = TMath::Abs(muon_phi[z_local_index[0]]);
          float phi2 = TMath::Abs(muon_phi[z_local_index[1]]);
          float phi3 = TMath::Abs(muon_phi[z_local_index[2]]);

          float diff1_2 = TMath::Abs(phi1 + phi2 - pi);
          float diff1_3 = TMath::Abs(phi1 + phi3 - pi);

          if (diff1_2 < diff1_3){
            massZ = local_masses[0];
            z_temp_index[z_temp_count] = z_local_index[1];
          }
          else {
            massZ = local_masses[1];
            z_temp_index[z_temp_count] = z_local_index[2];
          };
        }
        else if (local_z == 1) {
         z_temp_index[z_temp_count] = z_local_index[1]; 
        } 
      }
      
      // Due to the prescence of MET, finding 2 Z is not possible.
      // Therefore, we choose the one whose muons are in a more opposite in direction
      if (num_z == 2) {
        float phi1 = TMath::Abs(muon_phi[z_temp_index[0]]);
        float phi2 = TMath::Abs(muon_phi[z_temp_index[1]]);
        float phi3 = TMath::Abs(muon_phi[z_temp_index[2]]);
        float phi4 = TMath::Abs(muon_phi[z_temp_index[3]]);

        float diff1_2 = TMath::Abs(phi1 + phi2 - pi);
        float diff3_4 = TMath::Abs(phi3 + phi4 - pi);

        if (diff1_2 < diff3_4){
          z_muon_index[0] = z_temp_index[0];
          z_muon_index[1] = z_temp_index[1];
        }
        else {
          z_muon_index[0] = z_temp_index[2];
          z_muon_index[1] = z_temp_index[3];
        }
      }
      if (num_z) {hasZ = true;}
    }
};

//Gets xy components from pT and phi
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
  ofstream outfile("data2.csv", ofstream::out);
  save_file(outfile, runVec, eventVec, pt_visAVec, phi_visAVec, eta_visAVec, pt_visBVec, phi_visBVec, eta_visBVec, pt_METVec, phi_METVec, mT2Vec, mllVec);
}
