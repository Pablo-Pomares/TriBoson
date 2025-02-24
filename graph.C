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

using namespace std;

void graph(const string& file_name){
  TFile *f = TFile::Open(file_name.c_str(), "READ")
  TTree *t;
  
  const int nentries = t->GetEntries();

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
