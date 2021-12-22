//Toy Model for v1 and the Reaction Plane Angle using v1 = coefficient * momentum
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TLegend.h"

int Global_pTChoice;
int Spectators_Number;

using namespace std;

void Draw2DPlot(TH2D* h2, string draw_options, bool logz, string x_title,
                string label, string out_file, string y_axis){
  TCanvas* c1 = new TCanvas(x_title.c_str(),"c2",500,500);
  c1->cd();
  h2->SetMarkerSize(0.5);
  h2->SetMarkerStyle(20);
  if(x_title != "keep") h2->GetXaxis()->SetTitle(x_title.c_str());
  h2->GetXaxis()->SetTitleSize(0.05);
  if(y_axis != "keep") h2->GetYaxis()->SetTitle(y_axis.c_str());
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->Draw(draw_options.c_str());
  if(logz) {
    gPad->SetLogz();
    std::cout << h2->GetMinimum(1e-8) << std::endl;
    h2->GetZaxis()->SetRangeUser( h2->GetMinimum(1e-8),h2->GetMaximum()*2);
  }
  gPad->SetTopMargin(0.12);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.12);
  label = "";
  gStyle->SetOptStat(0);
  c1->Print(out_file.c_str());
  gStyle->SetOptStat("rme");
  delete c1;
}

void DrawPlot(vector <TH1D*> h1, bool logy, string x_title, string y_axis, string label, string out_file){

  TCanvas* c1 = new TCanvas(x_title.c_str(),"c1",500,500);
  c1->cd();

  TLegend*leg = new TLegend(0.4, 0.85,0.95, 0.98 );
  label = "";
  //h1->SetFillColorAlpha(color,0.3);
  h1.at(0)->GetXaxis()->SetTitle(x_title.c_str());
  h1.at(0)->GetXaxis()->SetTitleSize(0.05);
  h1.at(0)->GetYaxis()->SetTitle(y_axis.c_str());
  h1.at(0)->GetYaxis()->SetTitleSize(0.05);
  h1.at(0)->Draw();
  leg->AddEntry(h1.at(0),h1.at(0)->GetTitle(),"l");

  if(logy){
    gPad->SetLogy();
     h1.at(0)->GetYaxis()->SetRangeUser(1e-6, h1.at(0)->GetMaximum()*100.);
  }
  if(h1.size() > 1){
    for(int i = 1; i < (int)h1.size(); i++){
      h1.at(i)->Draw("SAME");
      leg->AddEntry(h1.at(i),h1.at(i)->GetTitle(),"l");
    }
  }
  gPad->SetTopMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.12);

  TLatex *lat = new TLatex();
  lat->SetTextFont(72);
  lat->SetTextSize(0.04);
  lat->SetTextFont(42);
  lat->DrawLatexNDC(.2,.92,label.c_str());
  leg->SetTextSize(0.035
  );
  leg->Draw();
  gPad->Update();
  TPaveStats *st = (TPaveStats*)h1.at(0)->FindObject("stats");
  st->SetX1NDC(0.13);
  st->SetX2NDC(0.4);
  st->SetY1NDC(0.85);
  st->SetY2NDC(0.98);
  gPad->Update();
  gStyle->SetOptStat("rme");
  c1->Print(out_file.c_str());
  delete c1;

}

//Version 3.0 - Both ZDC included in the simulation + Fermi momentum toy generation
//pT distribution now coming from the Fermi momentum simulation - therefore the Fermi Momentum generation is being enabled using

void V1Gen( int events = 1, bool debug = false ){

  int spectatorRange[2];

  spectatorRange[0] = 10;
  spectatorRange[1] = 50;

  //Per nucleon energy in the HI collision, in GeV
  double Ptot = 2760.;
  double neutron_mass = 0.939565; // GeV
  //  double z_distance = 127000. // distance from IP in mm for run 4
  
  double z_distance = 141000.; // distance from IP in mm for run 4
  
  TGraph *grpos = new TGraph(); //The TGraph was used to plot momentum space of the particles
  grpos->SetMarkerColor(2);
  grpos->SetTitle("ZDC+ Position");

  TGraph *grmom = new TGraph(); //The TGraph was used to plot momentum space of the particles
  grmom->SetMarkerColor(2);
  grmom->SetTitle("ZDC+ P_{x} vs P_{y}");

  TGraph *grposB = new TGraph(); //The TGraph was used to plot momentum space of the particles
  grposB->SetMarkerColor(1);
  grposB->SetTitle("ZDC- Position");

  TGraph *grmomB = new TGraph(); //The TGraph was used to plot momentum space of the particles
  grmomB->SetMarkerColor(1);
  grmomB->SetTitle("ZDC- P_{x} vs P_{y}");

  double FermiMomentum, FermiMomentumMax = 0, Pref = 0, gamma = 0, beta = 0;
  
  //===================================
  //Global pt selection
  FermiMomentumMax = 0.265; //in GeV  -- provided by Brian, 80s paper reference
  //Computing relativistic factors
  gamma = sqrt(Ptot*Ptot+neutron_mass*neutron_mass)/neutron_mass;
  beta = sqrt(1 - pow(1/gamma, 2));
  cout << "You've picked the Fermi Momentum generation! " << endl;
  cout << "Gamma: " << gamma << endl;
  cout << "Beta: " << beta << endl;


  //PARAMETER FOR GENERATION
  double pTnuclearModule;
  double pTNuclearComponents[2];
  double pTNuclearComponentsB[2];
  double reactionPlaneAngle; //Defining the reaction plane angle along which the pT of v1 will be applied

  int particles = -1;   //# of neutrons going to ZDC + (or A)
  int particlesB = -1;   //# of neutrons going to ZDC - (or C)
  

  double bufferPtGen;
  //Random number generator engine -- with seed set to 0, guaranteed to be unique to time and space
  TRandom3* myRand = new TRandom3(0);
  //Containers
  vector < double > pTplus;         //Vector containing pT of neutrons going towards ZDC + (or A)
  vector < double > azimuth_plus;   //Vector containing azimuthal angle of neutrons going towards ZDC + (or A)
  vector < double > pTminus;        //Vector containing pT of neutrons going towards ZDC - (or C)
  vector < double > azimuth_minus;  //Vector containing azimuthal angle of neutrons going towards ZDC - (or C)
  vector < TLorentzVector* > neutronPlus;
  vector < TLorentzVector* > neutronMinus;

  //Debug tools
  vector < TLorentzVector* > orig_neutronPlus;
  vector < TLorentzVector* > orig_neutronMinus;
  
  vector < int > vPdgid, vPdgid_b, vStatus, vStatus_b, vNPart, vNPart_b;
  vector < double > vNeutronMass;
  vector < double > vNeutronMass_b;
  double ImpactParameter = 0;
  int pdgid = 2112;
  int status = 1;
  int nPart = -1;
  int nPart_b = -1;

    
  //Buffer variables for kinematic computation
  vector < double > vPx, vPy, vPz, vE;
  vector <double >  vPx_b, vPy_b, vPz_b, vE_b;
  //  TFile* myFileA = new TFile(Form("A_myGeneration_%d.root",processNumber),"RECREATE");
  TFile* myFileA = new TFile("A_myGeneration.root","RECREATE");
  double avgX = 0;
  double avgY = 0;
  double avgX_b = 0;
  double avgY_b = 0;
  vector < double > vX, vY, vX_b, vY_b;
  double reactionPlaneGeneratedAngle, reactionPlaneGeneratedAngleB;

  TTree* tree = new TTree("Particle","myGeneration - ZDC +");

  tree->Branch("nNeutrons", &particles, "nNeutrons/I");
  tree->Branch("RP_true_value", &reactionPlaneAngle, "RP_true_value/D");
  tree->Branch("RP_gen_value", &reactionPlaneGeneratedAngle, "RP_gen_value/D");
  tree->Branch("avgX", &avgX, "avgX/D");
  tree->Branch("avgY", &avgY, "avgY/D");   
  tree->Branch("nPart", &nPart);
  tree->Branch("x", &vX);
  tree->Branch("y", &vY);
  tree->Branch("px", &vPx);
  tree->Branch("py", &vPy);
  tree->Branch("pz", &vPz);
  tree->Branch("pdgid", &vPdgid);
  tree->Branch("status", &vStatus);
  tree->Branch("ImpactParameter", &ImpactParameter);
  tree->Branch("m", &vNeutronMass);
  tree->Branch("E", &vE);
  tree->Branch("pt_nuclear", &pTnuclearModule);



  //  TFile* myFileB = new TFile(Form("B_myGeneration_%d.root",processNumber),"RECREATE");
  TFile* myFileB = new TFile("B_myGeneration.root","RECREATE");
  //
  TTree* treeB = new TTree("Particle","myGeneration - ZDC -");

  treeB->Branch("nNeutrons", &particlesB, "nNeutrons/I");
  treeB->Branch("RP_true_value", &reactionPlaneAngle, "RP_true_value/D");
  treeB->Branch("RP_gen_value", &reactionPlaneGeneratedAngleB, "RP_gen_value/D");
  treeB->Branch("nPart", &nPart_b);
  treeB->Branch("avgX", &avgX_b, "avgX/D");
  treeB->Branch("avgY", &avgY_b, "avgY/D");   
  treeB->Branch("x", &vX_b);
  treeB->Branch("y", &vY_b);
  treeB->Branch("px", &vPx_b);
  treeB->Branch("py", &vPy_b);
  treeB->Branch("pz", &vPz_b);
  treeB->Branch("pdgid", &vPdgid_b);
  treeB->Branch("status", &vStatus_b);
  treeB->Branch("ImpactParameter", &ImpactParameter);
  treeB->Branch("m", &vNeutronMass_b);
  treeB->Branch("E", &vE_b);
  treeB->Branch("pt_nuclear", &pTnuclearModule);


  double diffBgen;
  double diffAgen;

  for(int i = 0; i< events; i++)
  {
    std::cout << " event: " << i << std::endl;
    //===================================
        //randomize the nulcear pt
        pTnuclearModule = 0.005 + myRand->Rndm()*0.045;
    
        //Clearing containers

        neutronPlus.clear();
        orig_neutronPlus.clear();
        neutronMinus.clear();
        orig_neutronMinus.clear();
        vX_b.clear();
        vY_b.clear();
        vX.clear();
        vY.clear();
	vPx_b.clear();
        vPy_b.clear();
        vPz_b.clear();
        vE_b.clear();
        vPx.clear();
        vPy.clear();
	vPy_b.clear();
        vE.clear();
        vNeutronMass.clear();
        vNeutronMass_b.clear();
	vPdgid.clear();
	vStatus.clear();
	vPdgid_b.clear();
	vStatus_b.clear();
	
        reactionPlaneAngle = -666.;
        reactionPlaneGeneratedAngle = -666.;
        reactionPlaneGeneratedAngleB = -666.;
	
        //===================================
        //Spectators Block - Here is defined how many spectators we have for the event
	// uniform distribution between 0 and 1
	particles = spectatorRange[0] + myRand->Uniform() * (spectatorRange[1] - spectatorRange[0]);
	particlesB = spectatorRange[0] + myRand->Uniform() * (spectatorRange[1] - spectatorRange[0]);
	//	std::cout << " nParticles " << particles << " nParticlesB " << particlesB << std::endl;
	// if using this to simulate multiple neutrons at same time. set nPart to particles, otherwise nPart = 1
	//	nPart = particles;
	nPart = 1;
	//	nPart_b = particlesB;
	nPart_b = 1;

        //===================================
        //pT nuclear block - Here we extract a direction for the pT nuclear and we compute components -- defined b/w [-Pi,Pi]
        reactionPlaneAngle = -TMath::Pi()+(myRand->Rndm()*TMath::TwoPi());
        pTNuclearComponents[0] = pTnuclearModule*TMath::Cos(reactionPlaneAngle);
        pTNuclearComponents[1] = pTnuclearModule*TMath::Sin(reactionPlaneAngle);
	// directed flow assumes the deflection angle in the arm is equal and opposite in azimuth
        pTNuclearComponentsB[0] = -pTNuclearComponents[0];
        pTNuclearComponentsB[1] = -pTNuclearComponents[1];
	
	//Storing event info into vector
        TLorentzVector sumNeut = TLorentzVector(0,0,0,0);
	
        if(debug) cout << "Event " << i << endl;
        for(int j = 0; j < particles; j++)
      	{
	  if (debug) std::cout << " Particle: " << j << std::endl;
          //Extraction of pT for particle i

	  double Px, Py, Pz, E;

	  //Fermi momentum simulation
	  // WHY +1?????
	  Pref = FermiMomentumMax+1.;	    
	  while(Pref > FermiMomentumMax ){
	    // fermi momentum randomly selected between [-FMM, FMM]

	    Px = -FermiMomentumMax+myRand->Rndm()*2*FermiMomentumMax;
	    Py = -FermiMomentumMax+myRand->Rndm()*2*FermiMomentumMax;
	    Pz = -FermiMomentumMax+myRand->Rndm()*2*FermiMomentumMax;
	    E = sqrt(Px*Px+Py*Py+Pz*Pz+neutron_mass*neutron_mass);
	    Pref = sqrt(Px*Px+Py*Py+Pz*Pz);
	   
	  }
	  pTplus.push_back(sqrt(Px*Px+Py*Py));

          orig_neutronPlus.push_back(new TLorentzVector(Px, Py, Pz, E));
	  orig_neutronPlus.back()->Boost(0.,0.,beta);
	  Px = orig_neutronPlus.back()->Px();
	  Py = orig_neutronPlus.back()->Py();
	  Pz = orig_neutronPlus.back()->Pz();
	  E  = orig_neutronPlus.back()->E();

          //Now adding the contribution of pT nuclear
          Px += pTNuclearComponents[0];
          Py += pTNuclearComponents[1];
          neutronPlus.push_back(new TLorentzVector(Px, Py, Pz, E));
          if(neutronPlus.back()->Pt() > (FermiMomentumMax + pTnuclearModule)) cout << " WARNING - momentum conservation violated " << endl;
	  // 	  std::cout << "Final Px " << Px << " E " << E << std::endl;
	  vPx.push_back(Px);
	  vPy.push_back(Py);
	  vPz.push_back(Pz);
	  vE.push_back(E);
	  vPdgid.push_back(pdgid);
	  vStatus.push_back(status);
	  vNeutronMass.push_back(neutron_mass);
	  if (j == particles - 1) {
	    for(int k = 0; k < (int)neutronPlus.size(); k++){
	      //	      std::cout << " vX " << (neutronPlus.at(k)->Px()/neutronPlus.at(k)->Pz())*z_distance << std::endl;
	      double x = (neutronPlus.at(k)->Px()/neutronPlus.at(k)->Pz())*z_distance;
	      double y = (neutronPlus.at(k)->Py()/neutronPlus.at(k)->Pz())*z_distance;
	      vX.push_back(x);
	      vY.push_back(y);
	      avgX += x;
	      avgY += y;
	      sumNeut += (*neutronPlus.at(k));

	    }
	    reactionPlaneGeneratedAngle = sumNeut.Phi();
	    avgX /= particles;
	    avgY /= particles;
	    //	    std::cout << " RXN PLANE ANGLE !!!!!!!!!! " << reactionPlaneGeneratedAngle << std::endl;
	  }
	  tree->Fill();

      	}//End of loop on particles of ZDC + (A)

	//Storing event info into vector
        TLorentzVector sumNeutB = TLorentzVector(0,0,0,0);

        for(int j = 0; j < particlesB; j++)
      	{
	  double Px_b, Py_b, Pz_b, E_b;
	  // Fermi momentum simulation
	  Pref = FermiMomentumMax+1.;
	  while(Pref > FermiMomentumMax ){
	    // fermi momentum different for each arm
	    Px_b = -FermiMomentumMax+myRand->Rndm()*2*FermiMomentumMax;
	    Py_b = -FermiMomentumMax+myRand->Rndm()*2*FermiMomentumMax;
	    Pz_b = -FermiMomentumMax+myRand->Rndm()*2*FermiMomentumMax;
	    E_b = sqrt(Px_b*Px_b+Py_b*Py_b+Pz_b*Pz_b+neutron_mass*neutron_mass);
	    Pref = sqrt(Px_b*Px_b+Py_b*Py_b+Pz_b*Pz_b);
	  }
	  pTminus.push_back(sqrt(Px_b*Px_b+Py_b*Py_b));

	  // Apply boost
          orig_neutronMinus.push_back(new TLorentzVector(Px_b, Py_b, Pz_b, E_b));
	  orig_neutronMinus.back()->Boost(0.,0.,beta);
	  Px_b = orig_neutronMinus.back()->Px();
	  Py_b = orig_neutronMinus.back()->Py();
	  Pz_b = orig_neutronMinus.back()->Pz();
	  E_b  = orig_neutronMinus.back()->E();

          // Now adding the contribution of pT nuclear
          Px_b += pTNuclearComponentsB[0];
          Py_b += pTNuclearComponentsB[1];
          neutronMinus.push_back(new TLorentzVector(Px_b, Py_b, Pz_b, E_b));
          if(neutronMinus.back()->Pt() > (FermiMomentumMax + pTnuclearModule)) cout << " WARNING - momentum conservation violated " << endl;
	  vPx_b.push_back(Px_b);
	  vPy_b.push_back(Py_b);
	  vPz_b.push_back(Pz_b);
	  vE_b.push_back(E_b);
	  vPdgid_b.push_back(pdgid);
	  vStatus_b.push_back(status);
	  vNeutronMass_b.push_back(neutron_mass);

	  if (j == particlesB - 1) {
	    for(int k = 0; k < (int)neutronMinus.size(); k++){
	      double x = (neutronMinus.at(k)->Px()/neutronMinus.at(k)->Pz())*(-1*z_distance);
	      double y = (neutronMinus.at(k)->Py()/neutronMinus.at(k)->Pz())*(-1*z_distance);
	      vX_b.push_back(x);
	      vY_b.push_back(y);
	      avgX_b += x;
	      avgY_b += y;
	      sumNeutB += (*neutronMinus.at(k));
	    }
	    reactionPlaneGeneratedAngleB = sumNeutB.Phi();
	    avgX_b /= particlesB;
	    avgY_b /= particlesB;

	    //	    std::cout << " RXN PLANE ANGLE B !!!!!!!!!! " << reactionPlaneGeneratedAngleB << std::endl;
	  }
	  treeB->Fill();
      	}//End of loop on particles of ZDC + (A)



	// uncomment these (and comment out previous fills) to have multi-neutron events
	//        tree->Fill();
	//        treeB->Fill();
  }//End of event loop


//DrawPlot(vector <TH1D*> h1, bool logy, string x_title, string y_axis, string label, string out_file)



  myFileB->cd();
  treeB->Write();

  myFileA->cd();
  tree->Write();


  myFileA->Close();
  myFileB->Close();


  return;
}
