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
/// \file AnalysisManager.hh
/// \brief Definition of the AnalysisManager class
/// \author Chad Lantz
/// \date 16 April 2020

// Note: Relies on use of G4AnalysisManager. See: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Analysis/managers.html#analysis-manager


#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "SteppingAction.hh"
#include "ModTypeZDC.hh"
#include "ModTypeRPD.hh"

#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisManager* AnalysisManager::analysisManager = NULL;

AnalysisManager* AnalysisManager::getInstance(void)
{
    if (analysisManager == NULL) {
        analysisManager = new AnalysisManager();
    }
    return analysisManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisManager::AnalysisManager()
 : m_FactoryOn(false)
{
  m_analysisManager = G4AnalysisManager::Instance();
  m_lastStepVec = new std::vector< G4ThreeVector >;
  m_gunPos = new G4ThreeVector;

  m_Pi0Vert = new std::vector< G4ThreeVector >;
  m_Pi0Mom = new std::vector< G4ThreeVector >;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisManager::~AnalysisManager()
{
  delete m_ZDCdblVec;
  delete m_RPDdblVec;
  delete m_ZDCintVec;
  delete m_RPDintVec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::Book( G4String fileName )
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in AnalysisManager.hh
  G4cout << "Using " << m_analysisManager->GetType() << " type Analysis Manager" << G4endl;
  m_analysisManager->SetVerboseLevel(1);
  //  m_analysisManager->SetNtupleMerging(true);

  // Open an output file
  //
  G4bool fileOpen;
  if(fileName == ""){
    fileOpen = m_analysisManager->OpenFile("output");
  } else {
    fileOpen = m_analysisManager->OpenFile(fileName);
  }
  if (! fileOpen) {
    G4cerr << "\n---> AnalysisManager::Book(): cannot open "
           << m_analysisManager->GetFileName() << G4endl;
    return;
  }

  //Get information about the detector configuration from DetectorConstruction
  DetectorConstruction* detectorConstruction = (DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  std::vector< ModTypeZDC* > *ZDCvec = detectorConstruction->GetZDCvec();
  std::vector< ModTypeRPD* > *RPDvec = detectorConstruction->GetRPDvec();
  OPTICAL   = detectorConstruction->GetOpticalFlag();
  PI0   = detectorConstruction->GetPI0Flag();

  //Make vectors for the detectors we have
  //Indecies are [module#][dataType][dataPoint]
  m_ZDCdblVec = new std::vector< std::vector< std::vector<double> > >( ZDCvec->size() );
  m_ZDCintVec = new std::vector< std::vector< std::vector< int  > > >( ZDCvec->size() );
  m_RPDdblVec = new std::vector< std::vector< std::vector<double> > >( RPDvec->size() );
  m_RPDintVec = new std::vector< std::vector< std::vector< int  > > >( RPDvec->size() );

  SteppingAction* steppingAction = (SteppingAction*)G4RunManager::GetRunManager()->GetUserSteppingAction();
  steppingAction->SetLastStepVec( m_lastStepVec, &m_lastStepPidVec );
  steppingAction->SetPi0Mom( m_Pi0Mom );
  steppingAction->SetPi0Vertex( m_Pi0Vert );


  //Create the event data tree (Tree 0)
  MakeEventDataTree();

  //Create ZDC trees and branches (Tree 1 - nZDCs)
  m_ZDCfiberVec.resize( ZDCvec->size(), 0 );
  m_ZDCtimeVec.resize( ZDCvec->size(), 0 );
  for(uint i = 0; i < ZDCvec->size(); i++){

    // Create the vector that will store nCherenkovs for the SD if necessary.
    // Must be done here because the vector address is needed to create
    // the output tree, but the SD has not been created yet
    if( ZDCvec->at(i)->GetReducedTreeFlag() ){
      m_ZDCfiberVec[i] = new std::vector< G4int  >( ZDCvec->at(i)->GetnFibers(), 0 );
      m_ZDCtimeVec[i]  = new std::vector< G4double  >( 128, 0 );
    }
    if( ZDCvec->at(i)->GetMLReducedTreeFlag() ){
      m_ZDCfiberVec[i] = new std::vector< G4int  >( 1, 0 );
    }

    MakeZDCTree( ZDCvec->at(i)->GetModNum() - 1, m_ZDCfiberVec[i], m_ZDCtimeVec[i], ZDCvec->at(i)->GetReducedTreeFlag(), ZDCvec->at(i)->GetMLReducedTreeFlag(), ZDCvec->at(i)->GetOpticalFlag() );

  }//end ZDC loop

  //Create RPD trees (Tree nZDCs - nZDCs + nRPDs)
  // RPDvec->size() is just the number of modules
  m_RPDfiberVec.resize( RPDvec->size(), 0 );
  m_YOriginVec.resize( RPDvec->size(), 0 );
  m_EnergyVec.resize( RPDvec->size(), 0 );
  m_channelVec.resize( RPDvec->size(), 0 );
  m_timeHitVec.resize( RPDvec->size(), 0 );
  m_IncidenceAngleVec.resize( RPDvec->size(), 0 );
  m_RPDfiberGenVec.resize( RPDvec->size(), 0 );
  m_RPDtimeVec.resize( RPDvec->size(), 0 );
  for(uint i = 0; i < RPDvec->size(); i++){

    // Create the vector that will store nCherenkovs for the SD if necessary.
    // Must be done here because the vector address is needed to create
    // the output tree, but the SD has not been created yet
    if( RPDvec->at(i)->GetReducedTreeFlag() ){
      m_RPDfiberVec[i] = new std::vector< G4int  >( RPDvec->at(i)->GetnFibers(), 0 );
      m_RPDtimeVec[i]  = new std::vector< G4double  >( 128, 0 );
    }

    else if ( RPDvec->at(i)->GetMLReducedTreeFlag() ){
      m_IncidenceAngleVec[i] = new std::vector< G4double  >;
      m_EnergyVec[i] = new std::vector< G4double  >;
      m_channelVec[i] = new std::vector< G4int  >;
      m_timeHitVec[i] = new std::vector< G4double  >;
      m_YOriginVec[i] = new std::vector< G4double  >;
      m_RPDfiberGenVec[i] = new std::vector< G4int  >( RPDvec->at(i)->GetnChannels(), 0 );
      m_RPDfiberVec[i] = new std::vector< G4int  >( RPDvec->at(i)->GetnChannels(), 0 );
      m_RPDtimeVec[i]  = new std::vector< G4double  >( RPDvec->at(i)->GetnChannels(), 0 );
    }
    else {
      m_RPDfiberVec[i] = new std::vector< G4int  >( RPDvec->at(i)->GetnFibers(), 0 );
    }
    G4bool rpdOpticalFlag = RPDvec->at(i)->GetOpticalFlag();
    G4bool rpdIsOptical = 0;
    if (rpdOpticalFlag) rpdIsOptical = 1;
    MakeRPDTree( RPDvec->at(i)->GetModNum() - 1, i, RPDvec->at(i)->GetReducedTreeFlag(), RPDvec->at(i)->GetMLReducedTreeFlag(), rpdIsOptical);

  }//end RPD loop

  //Turn on the data factory I guess
  m_FactoryOn = true;

  G4cout << "\n----> Output file is open in "
         << m_analysisManager->GetFileName() << "."
         << m_analysisManager->GetFileType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::Save()
{
  if (! m_FactoryOn) return;

  m_analysisManager->Write();
  m_analysisManager->CloseFile();

  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;

  delete m_analysisManager;
  m_FactoryOn = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::FillNtuples(){

  // Fill the last step vectors

  //G4cout << "Before lastStep loop" << G4endl;
  for(uint i = 0; i < m_lastStepVec->size(); i++){
    m_lastStepXVec.push_back( m_lastStepVec->at(i).x() );
    m_lastStepYVec.push_back( m_lastStepVec->at(i).y() );
    m_lastStepZVec.push_back( m_lastStepVec->at(i).z() );
  }
  for(uint i = 0; i < m_Pi0Mom->size(); i++){
    m_Pi0MomX.push_back( m_Pi0Mom->at(i).x() );
    m_Pi0MomY.push_back( m_Pi0Mom->at(i).y() );
    m_Pi0MomZ.push_back( m_Pi0Mom->at(i).z() );
  }

  for(uint i = 0; i < m_Pi0Vert->size(); i++){
    m_Pi0VertX.push_back( m_Pi0Vert->at(i).x() );
    m_Pi0VertY.push_back( m_Pi0Vert->at(i).y() );
    m_Pi0VertZ.push_back( m_Pi0Vert->at(i).z() );
  }



  m_analysisManager->FillNtupleIColumn( m_eventDataNtupleNo, 0, m_eventNo );

  m_analysisManager->FillNtupleDColumn( m_eventDataNtupleNo, 1, m_gunPos->x() );
  m_analysisManager->FillNtupleDColumn( m_eventDataNtupleNo, 2, m_gunPos->y() );
  m_analysisManager->FillNtupleDColumn( m_eventDataNtupleNo, 3, m_gunPos->z() );
  // Note: Geant seeds are unsigned ints but we save it to the tree as a signed int to use the FillNtupleIColumn function. So if you want to look at a particular event, convert the signed int from the ntuple to unsigned int using:
  // int trueSeed = (int16_t) signedSeed; // defined in <stdint.h>
  int seed1 = (short) m_EventSeed1;
  m_analysisManager->FillNtupleIColumn( m_eventDataNtupleNo, 4, seed1 );
  int seed2 = (short) m_EventSeed2;
  m_analysisManager->FillNtupleIColumn( m_eventDataNtupleNo, 5, seed2 );

  // fill ntuples  //
  for(int i = 0; i < m_analysisManager->GetNofNtuples(); i++){
    m_analysisManager->AddNtupleRow(i);
  }

  // Clear ZDC vectors //
  for(uint i = 0; i < m_ZDCdblVec->size(); i++){
    for(uint j = 0; j < m_ZDCdblVec->at(i).size(); j++){
      m_ZDCdblVec->at(i).at(j).clear();
    }
    for(uint j = 0; j < m_ZDCintVec->at(i).size(); j++){
      m_ZDCintVec->at(i).at(j).clear();
    }
  }
  // Clear RPD vectors //
  for(uint i = 0; i < m_RPDdblVec->size(); i++){
    for(uint j = 0; j < m_RPDdblVec->at(i).size(); j++){
      m_RPDdblVec->at(i).at(j).clear();
    }
    for(uint j = 0; j < m_RPDintVec->at(i).size(); j++){
      m_RPDintVec->at(i).at(j).clear();
    }
  }
  m_lastStepVec->clear();
  m_lastStepXVec.clear();
  m_lastStepYVec.clear();
  m_lastStepZVec.clear();
  m_lastStepPidVec.clear();

  m_Pi0Mom->clear();
  m_Pi0MomX.clear();
  m_Pi0MomY.clear();
  m_Pi0MomZ.clear();

  m_Pi0Vert->clear();
  m_Pi0VertX.clear();
  m_Pi0VertY.clear();
  m_Pi0VertZ.clear();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::FillZDCnCherenkovs( int zdcNo, int nCherenkovs ){
  std::cout << " filling ZDC n Cherenkovs " << std::endl;
  m_analysisManager->FillNtupleIColumn( m_ZDCnTupleNo[zdcNo - 1], 0, nCherenkovs );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::FillRPDnCherenkovs( int rpdNo, int nCherenkovs ){
  std::cout << " filling RPD n Cherenkovs " << std::endl;
  m_analysisManager->FillNtupleIColumn( m_RPDnTupleNo[rpdNo - 1], 0, nCherenkovs );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::MakeEventDataTree( ){

  m_eventDataNtupleNo = m_analysisManager->GetNofNtuples();
  m_analysisManager->CreateNtuple( "EventData", "Event Data");

  //Integer
  m_analysisManager->CreateNtupleIColumn( m_eventDataNtupleNo, "EventNo"  );

  //Doubles
  m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "gunPosX"  );
  m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "gunPosY"  );
  m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "gunPosZ"  );

  // is really an unsigned int
  m_analysisManager->CreateNtupleIColumn( m_eventDataNtupleNo, "EventSeed1"  );
  m_analysisManager->CreateNtupleIColumn( m_eventDataNtupleNo, "EventSeed2"  );

  //std::vector< double >
  //  m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "lastStepX",   m_lastStepXVec   );
  //  m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "lastStepY",   m_lastStepYVec   );
  //  m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "lastStepZ",   m_lastStepZVec   );


  //std::vector< int >
  //  m_analysisManager->CreateNtupleIColumn( m_eventDataNtupleNo, "lastStepPID", m_lastStepPidVec );

  if(PI0){
      //std::vector< double >
      m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "pi0Px",   m_Pi0MomX );
      m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "pi0Py",   m_Pi0MomY );
      m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "pi0Pz",   m_Pi0MomZ );

      //std::vector< double >
      m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "pi0Vx",   m_Pi0VertX );
      m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "pi0Vy",   m_Pi0VertY );
      m_analysisManager->CreateNtupleDColumn( m_eventDataNtupleNo, "pi0Vz",   m_Pi0VertZ );
  }

  m_analysisManager->FinishNtuple( );

  G4cout << "Created EventData tree with ntuple number " << m_eventDataNtupleNo << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::MakeEventGenTree( std::vector< std::vector<int>* > &intVec , std::vector< std::vector<double>* > &dblVec, G4int type ){

  m_eventGenNtupleNo = m_analysisManager->GetNofNtuples();

  m_analysisManager->CreateNtuple( "EventGen", "Event Generator Data");

  if(type == 0){ // CRMC generated events
    //Integer
    m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo, "nPart"  );
    //    m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo, "nSpec"  );
    //    m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo, "model"  );

    //Doubles
    //    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "impactParameter"  );

    // vectors below need to be resized to single neutron and value moved to 0 position
    ////std::vector< double >
    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "px", *dblVec[0] );
    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "py", *dblVec[1] );
    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "pz", *dblVec[2] );
    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo,  "E", *dblVec[3] );
    //    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo,  "m", *dblVec[4] );


    //vector< int > branches
    m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo,      "pdgid", *intVec[0] );
    //    m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo, "CRMCstatus", *intVec[1] );
    //    m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo, "keptStatus", *intVec[2] );
  }
  else if(type == 1){ // Toy pt Generator
    //Integer
    m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo, "nSpectators"  );

    //Doubles
    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "ptCollision"  );
    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "ptBreakup"    );
    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "rpAngle"      );

    ////std::vector< double >
    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "px", *dblVec[0] );
    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "py", *dblVec[1] );
    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "pz", *dblVec[2] );
    m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo,  "E", *dblVec[3] );

  }

  m_analysisManager->FinishNtuple( );

  G4cout << "Created eventGen tree with ntuple number " << m_eventGenNtupleNo << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::FillEventGenTree( int nPart, int nSpec, int model, double impactParam ){
  m_analysisManager->FillNtupleIColumn( m_eventGenNtupleNo, 0, nPart );
  //  m_analysisManager->FillNtupleIColumn( m_eventGenNtupleNo, 1, nSpec );
  //  m_analysisManager->FillNtupleIColumn( m_eventGenNtupleNo, 2, model );

  //Doubles
  //  m_analysisManager->FillNtupleDColumn( m_eventGenNtupleNo, 3, impactParam );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::FillEventGenTree( int nSpectators, double ptCollision, double ptBreakup, double rpAngle ){
  m_analysisManager->FillNtupleIColumn( m_eventGenNtupleNo, 0, nSpectators );

  //Doubles
  m_analysisManager->FillNtupleDColumn( m_eventGenNtupleNo, 1, ptCollision );
  m_analysisManager->FillNtupleDColumn( m_eventGenNtupleNo, 2, ptBreakup   );
  m_analysisManager->FillNtupleDColumn( m_eventGenNtupleNo, 3, rpAngle     );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::MakeZDCTree( G4int zdcNo, std::vector< int >* nCherenkovVec, std::vector< double >* timeVec, G4bool reducedTree, G4bool MLReducedTree, G4bool thisIsOptical ){
  G4int nTupleNo = m_analysisManager->GetNofNtuples();
  m_ZDCnTupleNo.push_back( nTupleNo );
  char name[20];
  sprintf(name,"ZDC%dtree",zdcNo+1);
  m_analysisManager->CreateNtuple( name, "ZDC data");

  // If nCherenkovVec is defined for this SD the user has set the ReducedTree flag
  // Otherwise generate the full tree
  if( reducedTree ){

    // This branch is a vector nFibers long that contains the number of cherenkovs
    // produced in each fiber or photons that arrive at the top of each fiber during an event
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", *nCherenkovVec );

    // This branch acts as a histogram to store the timing of photon arrival or production
    // depending on optical settings
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "timeHist",    *timeVec       );

  }
  else if (MLReducedTree ){

    // This branch is a vector nFibers long that contains the number of cherenkovs
    // produced in each fiber or photons that arrive at the top of each fiber during an event
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", *nCherenkovVec );

  }
  else { //Full tree

    if(thisIsOptical){


      // This branch is a vector nFibers long that contains the number of cherenkovs
      // produced in each fiber during an event. Make this first so it's columnId is 0
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs"                               );

      //vector< int > branch
      m_ZDCintVec->at(zdcNo).resize(1);
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_ZDCintVec->at(zdcNo).at(0) );

      //Resize the vector for the number of optical branches storing double vectors
      m_ZDCdblVec->at(zdcNo).resize(7);
      //vector< double > branches
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_ZDCdblVec->at(zdcNo).at(0) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_ZDCdblVec->at(zdcNo).at(1) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_ZDCdblVec->at(zdcNo).at(2) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Px",          m_ZDCdblVec->at(zdcNo).at(3) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Py",          m_ZDCdblVec->at(zdcNo).at(4) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Pz",          m_ZDCdblVec->at(zdcNo).at(5) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "time",        m_ZDCdblVec->at(zdcNo).at(6) );

    //The next branches are only created if the optical flag is off
    } else {
      m_ZDCdblVec->at(zdcNo).resize(12);
      //Non-optical vector< double > branches
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_ZDCdblVec->at(zdcNo).at(0) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_ZDCdblVec->at(zdcNo).at(1) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_ZDCdblVec->at(zdcNo).at(2) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Px",          m_ZDCdblVec->at(zdcNo).at(3) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Py",          m_ZDCdblVec->at(zdcNo).at(4) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Pz",          m_ZDCdblVec->at(zdcNo).at(5) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "energy",      m_ZDCdblVec->at(zdcNo).at(6) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "velocity",    m_ZDCdblVec->at(zdcNo).at(7) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "beta",        m_ZDCdblVec->at(zdcNo).at(8) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "eDep",        m_ZDCdblVec->at(zdcNo).at(9) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "charge",      m_ZDCdblVec->at(zdcNo).at(10) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "time",        m_ZDCdblVec->at(zdcNo).at(11) );

      //Non-optical vector< int > branches
      m_ZDCintVec->at(zdcNo).resize(5);
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_ZDCintVec->at(zdcNo).at(0) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "modNo",       m_ZDCintVec->at(zdcNo).at(1) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_ZDCintVec->at(zdcNo).at(2) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "trackID",     m_ZDCintVec->at(zdcNo).at(3) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "pid",         m_ZDCintVec->at(zdcNo).at(4) );
    }
  }//end if nCherenkovVec
  m_analysisManager->FinishNtuple( );

  G4cout << "Created " << name << " tree with ntuple number " << nTupleNo << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void AnalysisManager::MakeRPDTree( G4int rpdNo, std::vector< double >* YOriginVec, std::vector< int >* nGenCherenkovVec, std::vector< int >* nCherenkovVec, std::vector< double >* timeVec, G4bool reducedTree, G4bool MLReducedTree, G4bool thisIsOptical ){
void AnalysisManager::MakeRPDTree( G4int rpdNo, G4int modNum, G4bool reducedTree, G4bool MLReducedTree, G4bool thisIsOptical ){
  G4int nTupleNo = m_analysisManager->GetNofNtuples();
  m_RPDnTupleNo.push_back( nTupleNo );

  char name[20];
  sprintf(name,"RPD%dtree",rpdNo+1);
  m_analysisManager->CreateNtuple( name, "RPD data");

  // If nCherenkovVec is defined for this SD the user has set the ReducedTree flag
  // Otherwise generate the full tree
  if( reducedTree ){

      // This branch is a vector nFibers long that contains the number of cherenkovs
      // produced in each fiber or photons that arrive at the top of each fiber during an event
    //      m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", *nCherenkovVec );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", *(m_RPDfiberVec[modNum]) );

      // This branch acts as a histogram to store the timing of photon arrival or production
      // depending on optical settings
    //     m_analysisManager->CreateNtupleDColumn( nTupleNo, "timeHist",    *timeVec       );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "timeHist",    *(m_RPDtimeVec[modNum])       );

  }
  else if( MLReducedTree ){
    // Branch of vectors nHits long to keep y origins

    //    m_analysisManager->CreateNtupleDColumn( nTupleNo, "incidenceAngle", *(m_IncidenceAngleVec[modNum]) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "energy", *(m_EnergyVec[modNum]) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "channel", *(m_channelVec[modNum]) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "time", *(m_timeHitVec[modNum]) );
    //    m_analysisManager->CreateNtupleDColumn( nTupleNo, "yOrigin", *(m_YOriginVec[modNum]) );
    // Branch of vectors nChannels long that contains the number of cherenkovs per channel for each event
    //    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nGenCherenkovs", *(m_RPDfiberGenVec[modNum]) );
    // Branch of vectors nChannels long that contains the number of cherenkovs per channel for each event
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", *(m_RPDfiberVec[modNum]));
    // This branch acts as a histogram to store the timing of photon arrival or production
    // depending on optical settings
    //    m_analysisManager->CreateNtupleDColumn( nTupleNo, "timeHist",    *timeVec       );
    //    m_analysisManager->CreateNtupleDColumn( nTupleNo, "timeHist",    *(m_RPDtimeVec[modNum])       );
  }
  else { //Full tree


    //This is done first so it can have an columnId of 0 if optical is on
    if(thisIsOptical){
      // This branch is a vector nFibers long that contains the number of cherenkovs
      // produced in each fiber during an event. Make this first so it's columnId is 0
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs"                               );

      //vector< int > branch
      m_RPDintVec->at(rpdNo).resize(1);
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_RPDintVec->at(rpdNo).at(0) );

      //Resize the vector for the number of optical branches storing double vectors
      //      m_RPDdblVec->at(rpdNo).resize(7);
      m_RPDdblVec->at(rpdNo).resize(8);
      //vector< double > branches
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_RPDdblVec->at(rpdNo).at(0) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_RPDdblVec->at(rpdNo).at(1) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_RPDdblVec->at(rpdNo).at(2) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Px",          m_RPDdblVec->at(rpdNo).at(3) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Py",          m_RPDdblVec->at(rpdNo).at(4) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Pz",          m_RPDdblVec->at(rpdNo).at(5) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "energy",      m_RPDdblVec->at(rpdNo).at(6) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "time",        m_RPDdblVec->at(rpdNo).at(7) );

    //The next branches are only created if the optical flag is off
    } else {
      m_RPDdblVec->at(rpdNo).resize(12); //Increase the size to store more branches
      //Non-optical vector< double > branches
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_RPDdblVec->at(rpdNo).at(0) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_RPDdblVec->at(rpdNo).at(1) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_RPDdblVec->at(rpdNo).at(2) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Px",          m_RPDdblVec->at(rpdNo).at(3) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Py",          m_RPDdblVec->at(rpdNo).at(4) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Pz",          m_RPDdblVec->at(rpdNo).at(5) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "energy",      m_RPDdblVec->at(rpdNo).at(6) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "velocity",    m_RPDdblVec->at(rpdNo).at(7) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "beta",        m_RPDdblVec->at(rpdNo).at(8) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "eDep",        m_RPDdblVec->at(rpdNo).at(9) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "charge",      m_RPDdblVec->at(rpdNo).at(10));
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "time",        m_RPDdblVec->at(rpdNo).at(11));

      //Non-optical vector< int > branches
      m_RPDintVec->at(rpdNo).resize(5);
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_RPDintVec->at(rpdNo).at(0) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "modNo",       m_RPDintVec->at(rpdNo).at(1) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_RPDintVec->at(rpdNo).at(2) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "trackID",     m_RPDintVec->at(rpdNo).at(3) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "pid",         m_RPDintVec->at(rpdNo).at(4) );
    }
  }//end if nCherenkovVec
  m_analysisManager->FinishNtuple( nTupleNo );

  G4cout << "Created " << name << " tree with ntuple number " << nTupleNo << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector< G4int >* AnalysisManager::GetFiberVector( G4bool ZDC, G4bool RPD, G4int modNum  ){

       if(ZDC) return m_ZDCfiberVec[ modNum - 1 ];
  else if(RPD) return m_RPDfiberVec[ modNum - 1 ];

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector< G4int >* AnalysisManager::GetFiberGenVector( G4int modNum  ){

  return m_RPDfiberGenVec[ modNum - 1 ];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector< G4double > * AnalysisManager::GetYOriginVector( G4int modNum  ){

  return m_YOriginVec[ modNum - 1 ];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector< G4double > * AnalysisManager::GetIncidenceAngleVector( G4int modNum  ){

  return m_IncidenceAngleVec[ modNum - 1 ];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector< G4double > * AnalysisManager::GetEnergyVector( G4int modNum  ){

  return m_EnergyVec[ modNum - 1 ];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::vector< G4int > * AnalysisManager::GetChannelVector( G4int modNum  ){

  return m_channelVec[ modNum - 1 ];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::vector< G4double > * AnalysisManager::GetTimeHitVector( G4int modNum  ){

  return m_timeHitVec[ modNum - 1 ];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector< G4double >* AnalysisManager::GetTimeVector( G4bool ZDC, G4bool RPD, G4int modNum  ){

  if(ZDC) return m_ZDCtimeVec[ modNum - 1 ];
  else if(RPD) return m_RPDtimeVec[ modNum - 1 ];

  return 0;
}
