void prueba_filtro1(){
  TFile *f_old = TFile::Open("root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/127C2975-1B1C-A046-AABF-62B77E757A86.root");
  TTree* t_old;
  f_old->GetObject("Events", t_old);

  // Número de entradas
  const int nentries = t_old->GetEntries();

  // Desactivo las todas las branches
  t_old->SetBranchStatus("*", 0);

  // Activo aquellas que quiero
  for (auto actBranchName : {"Muon_charge", "Muon_dxy", "Muon_dxyErr", "Muon_isGlobal", 
                             "Muon_isTracker", "Muon_pt", "Muon_phi", "Muon_eta", "nMuon", "MET_phi",
                             "MET_pt", "MET_significance", "nJet", "Jet_btagCSVV2", "Jet_btagDeepB"}){
    t_old->SetBranchStatus(actBranchName, 1);
  };

  TFile newfile("prueba1.root", "recreate");
  auto t_new = t_old->CloneTree(0);

  UInt_t nMuon;
  t_old->SetBranchAddress("nMuon", &nMuon);

  for (int i=0; i<nentries; i++){
    t_old->GetEntry(i);
    if (nMuon == 4){
      t_new->Fill();
    };
  };

  t_new->Print();
  newfile.Write();
}
