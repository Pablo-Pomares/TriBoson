//***ROOT*******************************************************************
//*Version: ROOT 6.32.04                                                   *
//*Description: Classes and functions used throughout the different codes  *
//*                                                                        *
//*Created: 08 Dec 2024                                                    *
//*Last Modified: 22 Feb 2025                                              *
//*Author: Pablo Pomares                                                   *
//*Email: pablo.pomaresv@alumno.buap.mx                                    *
//**************************************************************************

#include <vector>
#include <iostream>

using namespace std;

const Float_t muon_mass = 0.1056583; // MeV/c2
const double z_mass = 91.1880; // GeV/c2 pm 0.0020 per PDG(2025) 
const double pi = TMath::Pi();
const double luminosity = 16.290713420;
const double brTo4mu = 0.0004138;

//------------------General------------------

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

//------------------WWZ------------------

// Find the Z bosons on a event and its mass.
// If it finds more than one candidate it selects the more probable one. 
// i. e., the one with the greater absolute difference of phi.
class Z_finder{
  private:
    int z_temp_index[4] = {-1, -1, -1, -1}; // If more than one particle found in Z range the first pair
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

// Returns invariant mass of the higher and the lower pair
class Four_muon_mll{
	public:
    Double_t high_mass = 0.0; // mass of the highest pair
    Double_t low_mass = numeric_limits<Double_t>::max(); // mass of the lowest pair
    array<int, 2> high_index;
    array<int, 2> low_index;

    Four_muon_mll(Float_t muon_pt[4], Float_t muon_phi[4], Float_t muon_eta[4], Int_t muon_charge[4])
      : high_index{-1, -1}, low_index{-1, -1}
    {
			for (int i=0; i<4; i++){
				for (int j=i+1; j<4; j++){
          Float_t pt1 = muon_pt[i];
          Float_t pt2 = muon_pt[j];
          Float_t phi1 = muon_phi[i];
          Float_t phi2 = muon_phi[j];
          Float_t eta1 = muon_eta[i];
          Float_t eta2 = muon_eta[j];
          Int_t q1 = muon_charge[i];
          Int_t q2 = muon_charge[j];

          // Check for different charge
          bool opposite_charge = q1 + q2;

          if (opposite_charge) {
            Double_t m = inv_mass(pt1, pt2, phi1, phi2, eta1, eta2);
						if (m > high_mass){
							high_index = {i, j};
							high_mass = m;
						}
						if (m < low_mass){
							low_index = {i, j};
							low_mass = m;
						} 
					}
				}
      }
		}
};

float get_pt_4l (Float_t Muon_pt[4], Float_t Muon_phi[4]){

  Vec_comp mu1(Muon_pt[0], Muon_phi[0]);
  Vec_comp mu2(Muon_pt[1], Muon_phi[1]);
  Vec_comp mu3(Muon_pt[2], Muon_phi[2]);
  Vec_comp mu4(Muon_pt[3], Muon_phi[3]);
  float fourl_pt = TMath::Sqrt((mu1.px + mu2.px + mu3.px + mu4.px)*(mu1.px + mu2.px + mu3.px + mu4.px) + (mu1.py + mu2.py + mu3.py + mu4.py)*(mu1.py + mu2.py + mu3.py + mu4.py));
  bool pass_4l_pt = (fourl_pt > 40); 

  return fourl_pt;
}

//------------------WW------------------

// Searches for WW candidates.
class WW_Analyzer{
  public:
  // Conditions WW events must follow
  Float_t mll = 0;
  Float_t diff_to_Z = 0.0; // to filter Drell Yan background
  bool pass_mll = false;
  bool pass_pt_miss_proj = false;
  bool pass_notdy = false;
  bool pass_pairpT = false;
  float pt_miss_proj = 0.0;
  // All in one, for convinience
  bool pass_all;

  WW_Analyzer(Float_t Muon_pt[2], Float_t Muon_phi[2], Float_t Muon_eta[2], Float_t MET_pt, Float_t MET_phi){
    mll = inv_mass(Muon_pt[0], Muon_pt[1], Muon_phi[0], Muon_phi[1], Muon_eta[0], Muon_phi[1]);
    if (mll > 40) {pass_mll = true;};

    diff_to_Z = TMath::Abs(mll - z_mass);
    if (diff_to_Z > 15) {pass_notdy = true;};

    for (int i=0; i<2; i++){
      if (pass_pt_miss_proj) {continue;};
      pass_pt_miss_proj = WW_Analyzer::pt_miss_proj_test(MET_pt, MET_phi, Muon_phi[i]);
    }

    Vec_comp mu1(Muon_pt[0], Muon_phi[0]);
    Vec_comp mu2(Muon_pt[1], Muon_phi[1]);
    float pair_pt = TMath::Sqrt((mu1.px + mu2.px)*(mu1.px + mu2.px) + (mu1.py + mu2.py)*(mu1.py + mu2.py));
    pass_pairpT = (pair_pt > 30.0);

    pass_all = (pass_mll && pass_pt_miss_proj && pass_notdy && pass_pairpT); 
  }

  // Checks if the perpendicular component of p^miss_T with respect to a given muon is greater or 
  // equal to 20 GeV if it is requiered.
  bool pt_miss_proj_test(Float_t MET_pt, Float_t MET_phi, Float_t Muon_phi){
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
};
