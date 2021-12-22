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
/// \defgroup mc MonteCarlo
/// \ingroup mc
/// \file zdc.cc


#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4VModularPhysicsList.hh"
#include "G4FastSimulationPhysics.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "G4OpticalPhysics.hh"

#include "TRandom3.h"

#include "Randomize.hh"
#include <climits>

/*
*/
namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " zdc [-m macro ] [-t nThreads] [-r seed] [-i inputFileName] [-o outputFileName] [-p physicsParameterization]"
           << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 13 ) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String output = "";
  G4String input = "";
  G4bool   gflash = false;
  G4bool   fastOptical = false;
  G4long   myseed = -1;
  G4long   myseed1 = -1;
  G4long   myseed2 = -1;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif

  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m"  ) macro      = argv[i+1];
    else if ( G4String(argv[i]) == "-o"  ) output     = argv[i+1];
    else if ( G4String(argv[i]) == "-i"  ) input      = argv[i+1];
    else if ( G4String(argv[i]) == "-r"  ) myseed     = atoi(argv[i+1]);
    else if ( G4String(argv[i]) == "-r1"  ) myseed1     = atoi(argv[i+1]);
    else if ( G4String(argv[i]) == "-r2"  ) myseed2     = atoi(argv[i+1]);    
    else if ( G4String(argv[i]) == "-f"  ){
      G4String buffer = argv[i+1];
      if( buffer.contains("t") ) fastOptical = true;
    }
    else if ( G4String(argv[i]) == "-p"  ){
      G4String buffer = argv[i+1];
      if( buffer.contains("t") ) gflash = true;
    }
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
       nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }

  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if ( macro.size() == 0 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) runManager->SetNumberOfThreads(nThreads);
  G4cout << "Using G4MULTITHREADED" << G4endl;
#else
  G4RunManager * runManager = new G4RunManager;
#endif
  //Set the seed
    if (myseed != -1) {
      G4Random::setTheSeed(myseed);
    }
    else if (myseed1 != -1 && myseed2 != -1) {
      long seeds[2] = {myseed1,myseed2};      
      G4Random::setTheSeeds(seeds);
    }
    // default method is to let ROOT set seeds based off TUUID object (guaranteed unique to time and space)
    else {

      long seeds[2];
      TRandom3* myRand1 = new TRandom3(0);
      long seed = myRand1->Integer(std::numeric_limits<unsigned int>::max());
      seeds[0] = seed;
      TRandom3* myRand2 = new TRandom3(0);
      seed = myRand2->Integer(std::numeric_limits<unsigned int>::max());
      seeds[1] = seed;
  
      G4Random::setTheSeeds(seeds);
      runManager->SetRandomNumberStore(false);
      delete myRand1; delete myRand2;
    }

  // Set mandatory initialization classes
  
  // Detector construction
  runManager->SetUserInitialization( new DetectorConstruction( gflash ) );

  // Default ATLAS physics list
  G4VModularPhysicsList* physicsList = new QGSP_BERT;
  // G4VModularPhysicsList* physicsList = new FTFP_BERT;
  
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);

  // currently only configured for the RPD
  if (fastOptical) { 
    G4FastSimulationPhysics* fastSimPhysics = new G4FastSimulationPhysics();
    fastSimPhysics->ActivateFastSimulation("opticalphoton");
    physicsList->RegisterPhysics(fastSimPhysics);
  }
  
  runManager->SetUserInitialization(physicsList);

  // User action initialization
  runManager->SetUserInitialization( new ActionInitialization( output ) );
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( macro.size() ) {
    // batch mode
    if(input != ""){
      UImanager->ApplyCommand("/beam/input " + input);
    }
    UImanager->ExecuteMacroFile(macro);
  }
  else {
    // interactive mode
    if(input != ""){
      G4String command = "/beam/input " + input;
      UImanager->ApplyCommand(command);
    }
    UImanager->ExecuteMacroFile("init_vis.mac");
    ui->SessionStart();
    delete ui;
  }


  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
