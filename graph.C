#include <TROOT.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLatex.h"

using namespace std;

void showoff(const string& file_name){
  TTree *t = new TTree("t", "t");
  t->ReadFile(file_name.c_str(), "pt_visA/F:phi_visA/F:eta_visA/F:pt_visB/F:phi_visB/F:eta_visB/F:pt_MET/F:phi_MET/F:mT2/F:mll/F");
  
  const int nentries = t->GetEntries();
  Float_t pt_visA, pt_visB, phi_visA, phi_visB, mT2, mll;

  t->SetBranchAddress("pt_visA", &pt_visA);
  t->SetBranchAddress("phi_visA", &phi_visA);
  t->SetBranchAddress("pt_visB", &pt_visB);
  t->SetBranchAddress("phi_visB", &phi_visB);
  t->SetBranchAddress("mT2", &mT2);
  t->SetBranchAddress("mll", &mll);

  TCanvas* c1 = new TCanvas("c1", "", 800, 600); //Canvas for mT2
  TCanvas* c2 = new TCanvas("c2", "", 800, 600); //Canvas for mll
  TH1F* h1 = new TH1F("h1", "mT2", 30, 0, 90);
  TH1F* h2 = new TH1F("h2", "mll", 100, 0, 200);

  for (int i = 0; i < nentries; i++){
    t->GetEntry(i);
    h1->Fill(mT2);
    h2->Fill(mll);
  }
  
  c1->cd();
  h1->Draw();

  c2->cd();
  h2->Draw();
}
