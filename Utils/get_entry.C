#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"


void get_entry(bool isA){
  //RPD file
  TFile *input;
  if (isA) input = new TFile("A_myGeneration.root");
  else input = new TFile("B_myGeneration.root");
  //Variables to read-out the tree
  TTree* particle = (TTree*)input->Get("Particle");
  std::cout << particle->GetEntries() <<  endl;

}
