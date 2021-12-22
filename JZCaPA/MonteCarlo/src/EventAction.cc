//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \ingroup mc
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
/// \author Chad Lantz
/// \date 16 April 2020


#include "EventAction.hh"
#include "FiberHit.hh"
#include "FiberSD.hh"
#include "AnalysisManager.hh"

#include "G4UnitsTable.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction( )
: G4UserEventAction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  G4cout << ">>> Begin event " << evt->GetEventID() << G4endl;
  fEventNo = evt->GetEventID();
  m_analysisManager = AnalysisManager::getInstance();

  m_RPDdblVec = m_analysisManager->GetRPDdoubleVectors( );
  m_RPDintVec = m_analysisManager->GetRPDintVectors   ( );
  m_ZDCdblVec = m_analysisManager->GetZDCdoubleVectors( );
  m_ZDCintVec = m_analysisManager->GetZDCintVectors   ( );

  fEventSeed1 = CLHEP::HepRandom::getTheSeeds()[0];
  fEventSeed2 = CLHEP::HepRandom::getTheSeeds()[1];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt){
  G4cout << "Event " << evt->GetEventID() << " ended. Processing hits..." << G4endl;
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  FiberHitsCollection* HC = 0;
  G4int nCollections =  HCE->GetNumberOfCollections();
  
  //Bail if there are no hit collections
  if(!HCE) return;

  //Loop over the collections we do have
  for(int hitsCollID = 0; hitsCollID < nCollections; hitsCollID++) {
    HC = (FiberHitsCollection*)(HCE->GetHC(hitsCollID));
    G4String name = HC->GetSDname();

    FiberSD* sd = (FiberSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector( name );
    if( sd->OpticalIsOn() ) ProcessOpticalHitCollection( HC );
    else                    ProcessHitCollection( HC );

  }// end hit collection loop
  //Give event information to analysisManager and fill the Ntuples
  m_analysisManager->SetEventNo( fEventNo );
  G4PrimaryVertex* pVert = evt->GetPrimaryVertex();
  m_analysisManager->SetGunPosition( pVert->GetX0(), pVert->GetY0(), pVert->GetZ0() );
  m_analysisManager->SetEventSeeds( fEventSeed1, fEventSeed2 );
  m_analysisManager->FillNtuples();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Optical off + Nonreduced
void EventAction::ProcessHitCollection( FiberHitsCollection* HC ){

  G4int nCherenkovsSum = 0;
  G4double eDepSum = 0.0;

  G4String name = HC->GetSDname();
  FiberSD* sd = (FiberSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector( name );
  int modNo = sd->GetModNum();

  // If you're looking for reduced tree stuff, go to FiberSD and AnalysisManager.
  // FiberSD doesn't add hits to the collection in reduced tree mode, so there
  // is nothing to process here
  if( sd->IsReduced() || sd->IsMLReduced() ){
    G4cout << name << " nHits = " << sd->GetNhits() << G4endl;
    return;
  }

  int n_hit = HC->entries();
  G4cout << name << " nHits = " << n_hit << G4endl;

  for ( G4int i = 0 ; i < n_hit; i++){

    G4double      eDep          = (*HC)[i]->getEdep();

    G4int         rodNo         = (*HC)[i]->getRodNb();
    G4int         nCherenkovs   = (*HC)[i]->getNCherenkovs(); // This is the number of cherenkovs in a single step within the SD
    G4int         trackID       = (*HC)[i]->getTrackID();


    G4ThreeVector position      = (*HC)[i]->getPos();
    G4ThreeVector origin        = (*HC)[i]->getOrigin();
    G4ThreeVector momentum      = (*HC)[i]->getMomentum();
    G4double      energy        = (*HC)[i]->getEnergy();
    G4double      velocity      = (*HC)[i]->getVelocity();
    G4double      beta          = (*HC)[i]->getBeta();
    G4double      charge        = (*HC)[i]->getCharge();
    G4double      time          = (*HC)[i]->getTime();

    G4int         pid           = (*HC)[i]->getParticle()->GetPDGEncoding();
    
    // Goal of this block of code: only save hit once per charged particle per SD
    // Requires: summing data from each step per track per SD
    // Since geant tracks a single particle all the way through to its death, we can achieve this by summing contiguous steps from the same track in the same volume
    // Only save the hit on the track's final step in the SD. Then reset the Cherenkov and EDep counters to zero
    nCherenkovsSum += nCherenkovs;
    eDepSum += eDep;

    G4int nextTrack = i+1;

    if (nextTrack != n_hit) {
      G4int nextRodNo    = (*HC)[nextTrack]->getRodNb();
      G4int nextTrackID  = (*HC)[nextTrack]->getTrackID();
      
      if (trackID == nextTrackID && rodNo == nextRodNo)  {
	continue;
      }
    }        

    if( sd->IsRPD() ){
      //doubles
      m_RPDdblVec->at(modNo-1).at(0). push_back( position.x() );
      m_RPDdblVec->at(modNo-1).at(1). push_back( position.y() );
      m_RPDdblVec->at(modNo-1).at(2). push_back( position.z() );
      m_RPDdblVec->at(modNo-1).at(3). push_back( momentum.x() );
      m_RPDdblVec->at(modNo-1).at(4). push_back( momentum.y() );
      m_RPDdblVec->at(modNo-1).at(5). push_back( momentum.z() );
      m_RPDdblVec->at(modNo-1).at(6). push_back( energy       );
      m_RPDdblVec->at(modNo-1).at(7). push_back( velocity     );
      m_RPDdblVec->at(modNo-1).at(8). push_back( beta         );
      m_RPDdblVec->at(modNo-1).at(9). push_back( eDepSum      );
      m_RPDdblVec->at(modNo-1).at(10).push_back( charge       );
      m_RPDdblVec->at(modNo-1).at(11).push_back( time         );

      //ints
      m_RPDintVec->at(modNo-1).at(0).push_back( rodNo          );
      m_RPDintVec->at(modNo-1).at(1).push_back( modNo          );
      m_RPDintVec->at(modNo-1).at(2).push_back( nCherenkovsSum );
      m_RPDintVec->at(modNo-1).at(3).push_back( trackID        );
      m_RPDintVec->at(modNo-1).at(4).push_back( pid            );
    // end if RPD
    } else{
      if( sd->IsZDC() ){
        //doubles
        m_ZDCdblVec->at(modNo-1).at(0). push_back( position.x() );
        m_ZDCdblVec->at(modNo-1).at(1). push_back( position.y() );
        m_ZDCdblVec->at(modNo-1).at(2). push_back( position.z() );
        m_ZDCdblVec->at(modNo-1).at(3). push_back( momentum.x() );
        m_ZDCdblVec->at(modNo-1).at(4). push_back( momentum.y() );
        m_ZDCdblVec->at(modNo-1).at(5). push_back( momentum.z() );
        m_ZDCdblVec->at(modNo-1).at(6). push_back( energy       );
        m_ZDCdblVec->at(modNo-1).at(7). push_back( velocity     );
        m_ZDCdblVec->at(modNo-1).at(8). push_back( beta         );
	m_ZDCdblVec->at(modNo-1).at(9). push_back( eDepSum      );
        m_ZDCdblVec->at(modNo-1).at(10).push_back( charge       );
        m_ZDCdblVec->at(modNo-1).at(11).push_back( time         );

        //ints
        m_ZDCintVec->at(modNo-1).at(0).push_back( rodNo          );
        m_ZDCintVec->at(modNo-1).at(1).push_back( modNo          );
	m_ZDCintVec->at(modNo-1).at(2).push_back( nCherenkovsSum );
        m_ZDCintVec->at(modNo-1).at(3).push_back( trackID        );
        m_ZDCintVec->at(modNo-1).at(4).push_back( pid            );
      }// end if ZDC
    }// end else (!RPD)

    // Reset Cherenkov and EDep counters to zero after every saved hit
    nCherenkovsSum = 0;
    eDepSum = 0.;
  }//end of hit loop
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Optical on + Nonreduced
void EventAction::ProcessOpticalHitCollection ( FiberHitsCollection* HC ){
  G4String name = HC->GetSDname();
  FiberSD* sd = (FiberSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector( name );
  int modNo = sd->GetModNum();

  // If you're looking for reduced tree stuff, go to FiberSD and AnalysisManager.
  // FiberSD doesn't add hits to the collection in reduced tree mode, so there
  // is nothing to process here
  if( sd->IsReduced() || sd->IsMLReduced() ){
    G4cout << name << " nOpticalHits = " << sd->GetNhits() << G4endl;
    return;
  }


  int n_hit = HC->entries();
  G4cout << name << " nOpticalHits = " << n_hit << G4endl;
  for ( G4int i = 0 ; i < n_hit; i++){

    G4ThreeVector origin   = (*HC)[i]->getOrigin();
    G4ThreeVector momentum = (*HC)[i]->getMomentum();
    G4double      time     = (*HC)[i]->getTime();
    G4int         rodNo    = (*HC)[i]->getRodNb();
    G4double      energy   = (*HC)[i]->getEnergy();
    
    if( sd->IsZDC() ){//ZDC hitsCollID, check to be sure/symmetric
      //double
      m_ZDCdblVec->at(modNo-1).at(0). push_back( origin.x()   );
      m_ZDCdblVec->at(modNo-1).at(1). push_back( origin.y()   );
      m_ZDCdblVec->at(modNo-1).at(2). push_back( origin.z()   );
      m_ZDCdblVec->at(modNo-1).at(3). push_back( momentum.x() );
      m_ZDCdblVec->at(modNo-1).at(4). push_back( momentum.y() );
      m_ZDCdblVec->at(modNo-1).at(5). push_back( momentum.z() );
      m_ZDCdblVec->at(modNo-1).at(6). push_back( energy       );
      m_ZDCdblVec->at(modNo-1).at(7). push_back( time         );

      //int
      m_ZDCintVec->at(modNo-1).at(0).push_back( rodNo );
    }//end fill ZDC vectors
    if( sd->IsRPD() ){//RPD hitsCollID, check to be sure/symmetric
      //double
      m_RPDdblVec->at(modNo-1).at(0). push_back( origin.x()   );
      m_RPDdblVec->at(modNo-1).at(1). push_back( origin.y()   );
      m_RPDdblVec->at(modNo-1).at(2). push_back( origin.z()   );
      m_RPDdblVec->at(modNo-1).at(3). push_back( momentum.x() );
      m_RPDdblVec->at(modNo-1).at(4). push_back( momentum.y() );
      m_RPDdblVec->at(modNo-1).at(5). push_back( momentum.z() );
      m_RPDdblVec->at(modNo-1).at(6). push_back( time         );

      //int
      m_RPDintVec->at(modNo-1).at(0). push_back( rodNo );
    }//end fill RPD vectors
  }// end hit loop

  //Fill the total number of cherenkovs produced in the event
  if( sd->IsZDC() ) m_analysisManager->FillZDCnCherenkovs( modNo, sd->GetNCherenkovs() );
  if( sd->IsRPD() ) m_analysisManager->FillRPDnCherenkovs( modNo, sd->GetNCherenkovs() );
}
