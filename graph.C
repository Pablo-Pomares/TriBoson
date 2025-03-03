//***ROOT******************************************************************
//*Version: ROOT 6.32.04                                                  *
//*Description: Graphs                                                    *                                      
//*                                                                       *
//*Created: 22 Dec 2024                                                   *
//*Last Modified: 22 Feb 2025                                             *
//*Author: Pablo Pomares                                                  *
//*Email: pablo.pomaresv@alumno.buap.mx                                   *
//*************************************************************************

#include <TROOT.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLatex.h"

#include "utils.h"
#include "lester_mt2_bisect.h"

using namespace std;

void graph(const string& file_name){
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


  TCanvas* c1 = new TCanvas("c1", "", 800, 600); //Canvas for mT2
  TCanvas* c2 = new TCanvas("c2", "", 800, 600); //Canvas for high mll
  TCanvas* c3 = new TCanvas("c3", "", 800, 600); //Canvas for low mll
  TH1F* h1 = new TH1F("h1", "mT2", 30, 0, 90);
  TH1F* h2 = new TH1F("h2", "High mll", 100, 0, 200);
  TH1F* h3 = new TH1F("h3", "Low mll", 100, 0, 200);

  for (int i = 0; i < nentries; i++){
    t->GetEntry(i);
    Four_muon_mll Info(Muon_pt, Muon_phi, Muon_eta, Muon_charge);

		int f1 = std::get<0>(Info.low_index);
		int f2 = std::get<1>(Info.low_index);
		float pt_visA = Muon_pt[f1];
		float phi_visA = Muon_phi[f1];
		float eta_visA = Muon_eta[f1];
		float pt_visB = Muon_pt[f2];
		float phi_visB = Muon_phi[f2];
		float eta_visB = Muon_eta[f2];

		Vec_comp p_visA(pt_visA, phi_visA);
		Vec_comp p_visB(pt_visB, phi_visB);
		Vec_comp p_miss(MET_pt, MET_phi);

		asymm_mt2_lester_bisect::disableCopyrightMessage();

		double mT2 = asymm_mt2_lester_bisect::get_mT2(muon_mass, p_visA.px, p_visA.py,
																									muon_mass, p_visB.px, p_visB.py,
																									p_miss.px, p_miss.py, 0, 0, 0);

    h1->Fill(mT2);
    h2->Fill(Info.high_mass);
    h3->Fill(Info.low_mass);
  }
  
  c1->cd();
  h1->Draw();

  c2->cd();
  h2->Draw();

  c3->cd();
  h3->Draw();
}
