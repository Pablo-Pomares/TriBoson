//***ROOT******************************************************************
//*Version: ROOT 6.32.04                                                  *
//*Description: Graphs                                                    *                                      
//*                                                                       *
//*Created: 22 Dec 2024                                                   *
//*Last Modified: 09 Mar 2025                                             *
//*Author: Pablo Pomares                                                  *
//*Email: pablo.pomaresv@alumno.buap.mx                                   *
//*************************************************************************

#include <TROOT.h>
#include <TStyle.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLatex.h"

#include "utils.h"
#include "lester_mt2_bisect.h"

using namespace std;

void graph(const string& file_name){

  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineStyle(0);

  TFile *f = TFile::Open(file_name.c_str(), "READ");
  TTree *t;
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

// Events of interest
  TCanvas* c1 = new TCanvas("c1", "", 800, 600); //Canvas for mT2
  TCanvas* c2 = new TCanvas("c2", "", 800, 600); //Canvas for high mll
  TCanvas* c3 = new TCanvas("c3", "", 800, 600); //Canvas for low mll
  TH1F* h1 = new TH1F("h1", "sel pT4l", 30, 0, 90);
  TH1F* h2 = new TH1F("h2", "sel High mll", 100, 0, 200);
  TH1F* h3 = new TH1F("h3", "sel Low mll", 100, 0, 200);

// All events
  TCanvas* c4 = new TCanvas("c4", "", 800, 600); //Canvas for mT2
  TCanvas* c5 = new TCanvas("c5", "", 800, 600); //Canvas for high mll
  TCanvas* c6 = new TCanvas("c6", "", 800, 600); //Canvas for low mll
  TCanvas* c7 = new TCanvas("c7", "", 800, 600); //Canvas for METpT
  TH1F* h4 = new TH1F("h4", "all pT4l", 30, 0, 90);
  TH1F* h5 = new TH1F("h5", "all High mll", 100, 0, 200);
  TH1F* h6 = new TH1F("h6", "all Low mll", 100, 0, 200);
  TH1F* h7 = new TH1F("h7", "all METpT", 100, 0, 200);

  for (int i = 0; i < nentries; i++){
    t->GetEntry(i);
    float pT_4l = get_pt_4l(Muon_pt, Muon_phi);
    Z_finder Info1(Muon_pt, Muon_phi, Muon_eta, Muon_charge);
    std::array<int, 2> index = Info1.z_muon_index;
    Four_muon_mll Info2(Muon_pt, Muon_phi, Muon_eta, Muon_charge);

		int f1 = std::get<0>(Info2.low_index);
		int f2 = std::get<1>(Info2.low_index);
		float pt_visA = Muon_pt[f1];
		float phi_visA = Muon_phi[f1];
		float eta_visA = Muon_eta[f1];
		float pt_visB = Muon_pt[f2];
		float phi_visB = Muon_phi[f2];
		float eta_visB = Muon_eta[f2];

    float free_muon_pt[2] = {pt_visA, pt_visB};
    float free_muon_phi[2] = {phi_visA, phi_visB};
    float free_muon_eta[2] = {eta_visA, eta_visB};

		Vec_comp p_visA(pt_visA, phi_visA);
		Vec_comp p_visB(pt_visB, phi_visB);
		Vec_comp p_miss(MET_pt, MET_phi);

    WW_Analyzer Info_ww(free_muon_pt, free_muon_phi, free_muon_eta, MET_pt, MET_phi);
    if (Info1.hasZ) {


      auto freemuons = free_muons(Info1.z_muon_index); // free muons
      int f1 = std::get<0>(freemuons);
      int f2 = std::get<1>(freemuons);
      float pt_visA_int = Muon_pt[f1];
      float phi_visA_int = Muon_phi[f1];
      float eta_visA_int = Muon_eta[f1];
      float pt_visB_int = Muon_pt[f2];
      float phi_visB_int = Muon_phi[f2];
      float eta_visB_int = Muon_eta[f2];

      bool passMET_pt = (MET_pt > 70);
      bool min_pt_z = (((Muon_pt[index[0]] > 20) && (Muon_pt[index[1]] > 20)) && ((Muon_pt[index[0]] > 25) || (Muon_pt[index[1]] > 25)));
      bool min_pt_ww = (((pt_visA > 20) && (pt_visB > 20)) && ((pt_visA > 25) || (pt_visB > 25)));
      bool min_muon_pt = (min_pt_z && min_pt_ww);
      if (passMET_pt && (pT_4l > 40) && Info_ww.pass_all && min_muon_pt) {
        h1->Fill(pT_4l);
        h2->Fill(Info2.high_mass);
        h3->Fill(Info2.low_mass);
      }
    }
    h4->Fill(pT_4l);
    h5->Fill(Info2.high_mass);
    h6->Fill(Info2.low_mass);
    h7->Fill(MET_pt);
  }
  
  c1->cd();
  h1->Draw();

  c2->cd();
  h2->Draw();

  c3->cd();
  h3->Draw();
  
  c4->cd();
  h4->Draw();

  c5->cd();
  h5->Draw();

  c6->cd();
  h6->Draw();

  c7->cd();
  h7->Draw();
}
