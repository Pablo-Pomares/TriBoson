//***ROOT******************************************************************
//*Version: ROOT 6.32.04                                                  *
//*Description: Graphs                                                    *                                      
//*                                                                       *
//*Created:  05 Mar 2025                                                   *
//*Last Modified: 09 Mar 2025                                             *
//*Author: Pablo Pomares                                                  *
//*Email: pablo.pomaresv@alumno.buap.mx                                   *
//*************************************************************************

#include <TROOT.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLatex.h"

#include "../utils.h"
#include "../lester_mt2_bisect.h"

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
  TH1F* h1 = new TH1F("h1", "mT2", 30, 5, 120);
  TH1F* h2 = new TH1F("h2", "mll", 100, 0, 200);

  TCanvas* c3 = new TCanvas("c3", "", 800, 600); //Canvas for mT2
  TCanvas* c4 = new TCanvas("c4", "", 800, 600); //Canvas for high mll
  TH1F* h3 = new TH1F("h3", "passing mT2", 30, 5, 120);
  TH1F* h4 = new TH1F("h4", "passing mll", 100, 0, 200);

  for (int i = 0; i < nentries; i++){
    t->GetEntry(i);
    WW_Analyzer Info(Muon_pt, Muon_phi, Muon_eta, MET_pt, MET_phi);
		float pt_visA = Muon_pt[0];
		float phi_visA = Muon_phi[0];
		float eta_visA = Muon_eta[0];
		float pt_visB = Muon_pt[1];
		float phi_visB = Muon_phi[1];
		float eta_visB = Muon_eta[1];

		Vec_comp p_visA(pt_visA, phi_visA);
		Vec_comp p_visB(pt_visB, phi_visB);
		Vec_comp p_miss(MET_pt, MET_phi);

		asymm_mt2_lester_bisect::disableCopyrightMessage();

		double mT2 = asymm_mt2_lester_bisect::get_mT2(muon_mass, p_visA.px, p_visA.py,
																									muon_mass, p_visB.px, p_visB.py,
																									p_miss.px, p_miss.py, 0, 0, 0);

    h1->Fill(mT2);
    h2->Fill(Info.mll);

    if (Info.pass_all){
      h3->Fill(mT2);
      h4->Fill(Info.mll);
    }
  }
  
  c1->cd();
  h1->Draw();
  c1->Print("img/all_mt2.png");

  c2->cd();
  h2->Draw();
  c1->Print("img/all_mll.png");

  c3->cd();
  h3->Draw();
  c1->Print("img/passing_mt2.png");

  c4->cd();
  h4->Draw();
  c1->Print("img/passing_mll.png");
}
