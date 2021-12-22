#include <vector>

void CondenseEvent(TString fileName) {

  TFile *f = TFile::Open(Form("%s",fileName.Data()),"READ");
  std::cout << " condensing file: " << fileName.Data() << std::endl;
  if (!f) { std::cout << " exiting no file " << std::endl; return; }

  TTree *t_rpd; f->GetObject("RPD1tree",t_rpd);
  TTree *t_zdc1; f->GetObject("ZDC1tree",t_zdc1);
  TTree *t_zdc2; f->GetObject("ZDC2tree",t_zdc2);
  TTree *t_zdc3; f->GetObject("ZDC3tree",t_zdc3);
  TTree *t_zdc4; f->GetObject("ZDC4tree",t_zdc4);
  TTree *t_eventData; f->GetObject("EventData",t_eventData);
  TTree *t_eventGen; f->GetObject("EventGen",t_eventGen);

  //RPD
  std::vector<double> *energy_rpd = 0;
  std::vector<int> *nCherenkovs_rpd = 0;
  std::vector<int> *channel_rpd = 0;
  std::vector<double> *time_rpd = 0;
  std::vector<double> energy_rpd_out;
  std::vector<int> nCherenkovs_rpd_out;
  std::vector<int> channel_rpd_out;
  std::vector<double> time_rpd_out;
  t_rpd->SetBranchAddress("energy",&energy_rpd);
  t_rpd->SetBranchAddress("nCherenkovs",&nCherenkovs_rpd);
  t_rpd->SetBranchAddress("channel",&channel_rpd);
  t_rpd->SetBranchAddress("time",&time_rpd);
  int nentries_rpd = (int) t_rpd->GetEntries();
  //ZDC1
  std::vector<int> *nCherenkovs_zdc1 = 0;
  std::vector<int> nCherenkovs_zdc1_out;
  t_zdc1->SetBranchAddress("nCherenkovs",&nCherenkovs_zdc1);
  int nentries_zdc1 = (int) t_zdc1->GetEntries();
  //ZDC2
  std::vector<int> *nCherenkovs_zdc2 = 0;
  std::vector<int> nCherenkovs_zdc2_out;
  t_zdc2->SetBranchAddress("nCherenkovs",&nCherenkovs_zdc2);
  int nentries_zdc2 = (int) t_zdc2->GetEntries();
  //ZDC3
  std::vector<int> *nCherenkovs_zdc3 = 0;
  std::vector<int> nCherenkovs_zdc3_out;
  t_zdc3->SetBranchAddress("nCherenkovs",&nCherenkovs_zdc3);
  int nentries_zdc3 = (int) t_zdc3->GetEntries();
  //ZDC4
  std::vector<int> *nCherenkovs_zdc4 = 0;
  std::vector<int> nCherenkovs_zdc4_out;
  t_zdc4->SetBranchAddress("nCherenkovs",&nCherenkovs_zdc4);
  int nentries_zdc4 = (int) t_zdc4->GetEntries();
  //EventData
  double gunPosX;
  double gunPosY;
  double gunPosZ;
  std::vector<double> gunPosX_out;
  std::vector<double> gunPosY_out;
  std::vector<double> gunPosZ_out;
  t_eventData->SetBranchAddress("gunPosX",&gunPosX);
  t_eventData->SetBranchAddress("gunPosY",&gunPosY);
  t_eventData->SetBranchAddress("gunPosZ",&gunPosZ);
  int nentries_eventData = (int) t_eventData->GetEntries();
  //EventGen
  int nPart;
  std::vector<int> *pdgid = 0;
  std::vector<double> *px = 0;
  std::vector<double> *py = 0;
  std::vector<double> *pz = 0;
  int nPart_out = 0;
  int pdgid_out;
  std::vector<double> px_out;
  std::vector<double> py_out;
  std::vector<double> pz_out;
  t_eventGen->SetBranchAddress("nPart",&nPart);
  t_eventGen->SetBranchAddress("pdgid",&pdgid);
  t_eventGen->SetBranchAddress("px",&px);
  t_eventGen->SetBranchAddress("py",&py);
  t_eventGen->SetBranchAddress("pz",&pz);
  int nentries_eventGen = (int) t_eventGen->GetEntries();

  // RPD
  for (int i = 0; i < nentries_rpd; ++i) {
    t_rpd->GetEntry(i);
    for (int j = 0; j < energy_rpd->size(); ++j) {
      energy_rpd_out.push_back(energy_rpd->at(j));
    }
    for (int j = 0; j < nCherenkovs_rpd->size(); ++j) {
      nCherenkovs_rpd_out.push_back(nCherenkovs_rpd->at(j));
    }
    for (int j = 0; j < channel_rpd->size(); ++j) {
      channel_rpd_out.push_back(channel_rpd->at(j));
    }
    for (int j = 0; j < time_rpd->size(); ++j) {
      time_rpd_out.push_back(time_rpd->at(j));
    }
  }
  // ZDC1
  for (int i = 0; i < nentries_zdc1; ++i) {
    t_zdc1->GetEntry(i);
    for (int j = 0; j < nCherenkovs_zdc1->size(); ++j) {
      nCherenkovs_zdc1_out.push_back(nCherenkovs_zdc1->at(j));
    }
  }
  // ZDC2
  for (int i = 0; i < nentries_zdc2; ++i) {
    t_zdc2->GetEntry(i);
    for (int j = 0; j < nCherenkovs_zdc2->size(); ++j) {
      nCherenkovs_zdc2_out.push_back(nCherenkovs_zdc2->at(j));
    }
  }
  // ZDC3
  for (int i = 0; i < nentries_zdc3; ++i) {
    t_zdc3->GetEntry(i);
    for (int j = 0; j < nCherenkovs_zdc3->size(); ++j) {
      nCherenkovs_zdc3_out.push_back(nCherenkovs_zdc3->at(j));
    }
  }
  // ZDC4
  for (int i = 0; i < nentries_zdc4; ++i) {
    t_zdc4->GetEntry(i);
    for (int j = 0; j < nCherenkovs_zdc4->size(); ++j) {
      nCherenkovs_zdc4_out.push_back(nCherenkovs_zdc4->at(j));
    }
  }
  // EventData
  for (int i = 0; i < nentries_eventData; ++i) {
    t_eventData->GetEntry(i);
    gunPosX_out.push_back(gunPosX);
    gunPosY_out.push_back(gunPosY);
    gunPosZ_out.push_back(gunPosZ);

  }
  // EventGen
  for (int i = 0; i < nentries_eventGen; ++i) {
    t_eventGen->GetEntry(i);
    nPart_out += nPart;
    pdgid_out = pdgid->at(0);
    for (int j = 0; j < px->size(); ++j) {
      px_out.push_back(px->at(j));
      py_out.push_back(py->at(j));
      pz_out.push_back(pz->at(j));
    }
  }
  f->Close();
  delete f;
  
  TFile *f_out = TFile::Open(Form("%s",fileName.Data()),"RECREATE");
  //EventData
  TTree *t_eventData_out = new TTree("EventData","EventData");
  t_eventData_out->Branch("gunPosX",&gunPosX_out);
  t_eventData_out->Branch("gunPosY",&gunPosY_out);
  t_eventData_out->Branch("gunPosZ",&gunPosZ_out);
  t_eventData_out->Fill();
  //ZDC1
  TTree *t_zdc1_out = new TTree("ZDC1tree","ZDC1tree");
  t_zdc1_out->Branch("nCherenkovs",&nCherenkovs_zdc1_out);
  t_zdc1_out->Fill();
  //ZDC2
  TTree *t_zdc2_out = new TTree("ZDC2tree","ZDC2tree");
  t_zdc2_out->Branch("nCherenkovs",&nCherenkovs_zdc2_out);
  t_zdc2_out->Fill();
  //ZDC3
  TTree *t_zdc3_out = new TTree("ZDC3tree","ZDC3tree");
  t_zdc3_out->Branch("nCherenkovs",&nCherenkovs_zdc3_out);
  t_zdc3_out->Fill();
  //ZDC4
  TTree *t_zdc4_out = new TTree("ZDC4tree","ZDC4tree");
  t_zdc4_out->Branch("nCherenkovs",&nCherenkovs_zdc4_out);
  t_zdc4_out->Fill();
  //RPD
  TTree *t_rpd_out = new TTree("RPD1tree","RPD1tree");
  t_rpd_out->Branch("energy",&energy_rpd_out);
  t_rpd_out->Branch("nCherenkovs",&nCherenkovs_rpd_out);
  t_rpd_out->Branch("channel",&channel_rpd_out);
  t_rpd_out->Branch("time",&time_rpd_out);
  t_rpd_out->Fill();
  //EventGen
  TTree *t_eventGen_out = new TTree("EventGen","EventGen");
  t_eventGen_out->Branch("nPart",&nPart_out);
  t_eventGen_out->Branch("px",&px_out);
  t_eventGen_out->Branch("py",&py_out);
  t_eventGen_out->Branch("pz",&pz_out);
  t_eventGen_out->Branch("pdgid",&pdgid_out);
  t_eventGen_out->Fill();

  std::cout << " Writing new condensed file " << std::endl;  
  f_out->Write();
  delete f_out;
}
