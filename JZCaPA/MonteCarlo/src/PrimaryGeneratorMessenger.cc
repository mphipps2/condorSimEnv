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
/// \ingroup mc
/// \file PrimaryGeneratorMessenger.cc
/// \brief Implementation of the PrimaryGeneratorMessenger class
//
//

#include "PrimaryGeneratorMessenger.hh"

#include <sstream>
#include <iostream>

#include "G4OpticalSurface.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"


PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction * Gen)
:G4UImessenger(),fGenerator(Gen)
{
  fGeneratorDir = new G4UIdirectory("/beam/");
  fGeneratorDir->SetGuidance("");

  fLHCDir = new G4UIdirectory("/beam/LHC/");
  fLHCDir->SetGuidance("");

  fSPSDir = new G4UIdirectory("/beam/SPS/");
  fSPSDir->SetGuidance("");

  fFNALDir = new G4UIdirectory("/beam/FNAL/");
  fFNALDir->SetGuidance("");

  fBeamTypeCmd = new G4UIcmdWithAString("/beam/type", this);
  fBeamTypeCmd->SetGuidance("Set the beam type. LHC, SPS, or FNAL");
  fBeamTypeCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fBeamTypeCmd->SetToBeBroadcasted(false);

  fCRMCmodelCmd = new G4UIcmdWithAString("/beam/GeneratorModel", this);
  fCRMCmodelCmd->SetGuidance("Generate events from CRMC using the selected model");
  fCRMCmodelCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fCRMCmodelCmd->SetParameterName("Event generator model",true);
  fCRMCmodelCmd->SetToBeBroadcasted(false);

  fInputFileCmd = new G4UIcmdWithAString("/beam/input", this);
  fInputFileCmd->SetGuidance("Set an input file with pre-generated events");
  fInputFileCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fInputFileCmd->SetParameterName("Input file name",true);
  fInputFileCmd->SetToBeBroadcasted(false);

  fPtDistCmd = new G4UIcmdWithAString("/beam/FragmentationPtDistribution", this);
  fPtDistCmd->SetGuidance("Set the distribution for pt due to fragmentation. Options are: gaus, histogram, and function");
  fPtDistCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fPtDistCmd->SetParameterName("Distribution type",true);
  fPtDistCmd->SetToBeBroadcasted(false);

  fMultiplicityDistCmd = new G4UIcmdWithAString("/beam/MultiplicityDistribution", this);
  fMultiplicityDistCmd->SetGuidance("Set the distribution for spectator multiplicity. Options are histogram and function");
  fMultiplicityDistCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fMultiplicityDistCmd->SetParameterName("Distribution type",true);
  fMultiplicityDistCmd->SetToBeBroadcasted(false);

  fBeamPosCmd = new G4UIcmdWith3VectorAndUnit("/beam/pos", this);
  fBeamPosCmd->SetGuidance("Set the origin position of the beam");
  fBeamPosCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fBeamPosCmd->SetToBeBroadcasted(false);
  fBeamPosCmd->SetParameterName("X","Y","Z",true);
  fBeamPosCmd->SetDefaultUnit("mm");

  fPseudorapitityCutCmd = new G4UIcmdWithADouble("/beam/SetMinimmPseudorapidity", this);
  fPseudorapitityCutCmd->SetGuidance("Minimum pseudorapidity to be considered by the event generator");
  fPseudorapitityCutCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fPseudorapitityCutCmd->SetToBeBroadcasted(false);

  fCollisionPtMeanCmd = new G4UIcmdWithADouble("/beam/SetMeanCollisionPt", this);
  fCollisionPtMeanCmd->SetGuidance("Set the mean pt kick due to fragmentation in the ToyV1 generator");
  fCollisionPtMeanCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fCollisionPtMeanCmd->SetToBeBroadcasted(false);

  fFragmentationPtMeanCmd = new G4UIcmdWithADouble("/beam/SetMeanFragmentationPt", this);
  fFragmentationPtMeanCmd->SetGuidance("Set the mean pt kick due to the collision in the ToyV1 generator");
  fFragmentationPtMeanCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fFragmentationPtMeanCmd->SetToBeBroadcasted(false);

  fProjectBeamCmd = new G4UIcmdWithADoubleAndUnit("/beam/projectBeam", this);
  fProjectBeamCmd->SetGuidance("Set the plane in Z which the beam will be projected to");
  fProjectBeamCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fProjectBeamCmd->SetToBeBroadcasted(false);
  fProjectBeamCmd->SetDefaultUnit("mm");

  fnPrimariesCmd = new G4UIcmdWithAnInteger("/beam/nPrimaries", this);
  fnPrimariesCmd->SetGuidance("Set the plane in Z which the beam will be projected to");
  fnPrimariesCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fnPrimariesCmd->SetToBeBroadcasted(false);

  fMinSpectatorsCmd = new G4UIcmdWithAnInteger("/beam/MinSpectators", this);
  fMinSpectatorsCmd->SetGuidance("Set the minimum number of spectators per side for the ToyV1 generator");
  fMinSpectatorsCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fMinSpectatorsCmd->SetToBeBroadcasted(false);

  fMaxSpectatorsCmd = new G4UIcmdWithAnInteger("/beam/MaxSpectators", this);
  fMaxSpectatorsCmd->SetGuidance("Set the maximum number of spectators per side for the ToyV1 generator");
  fMaxSpectatorsCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fMaxSpectatorsCmd->SetToBeBroadcasted(false);

  fRandomizeRPCmd = new G4UIcmdWithABool("/beam/randomizeReactionPlane", this);
  fRandomizeRPCmd->SetGuidance("Randomize the reaction plane angle each event");
  fRandomizeRPCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRandomizeRPCmd->SetDefaultValue(true);
  fRandomizeRPCmd->SetToBeBroadcasted(false);

  //LHC beam commands
  fVerticalCrossingCmd = new G4UIcmdWithADoubleAndUnit("/beam/LHC/verticalCrossingAngle", this);
  fVerticalCrossingCmd->SetGuidance("Set vertical crossing angle for LHC beam");
  fVerticalCrossingCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fVerticalCrossingCmd->SetToBeBroadcasted(false);
  fVerticalCrossingCmd->SetDefaultUnit("mrad");

  fHorizontalCrossingCmd = new G4UIcmdWithADoubleAndUnit("/beam/LHC/horizontalCrossingAngle", this);
  fHorizontalCrossingCmd->SetGuidance("Set horizontal crossing angle for LHC beam");
  fHorizontalCrossingCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fHorizontalCrossingCmd->SetToBeBroadcasted(false);
  fHorizontalCrossingCmd->SetDefaultUnit("mrad");



}

/*
 *
 */
PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger(){
  delete fBeamTypeCmd;
  delete fCRMCmodelCmd;
  delete fInputFileCmd;
  delete fPtDistCmd;
  delete fMultiplicityDistCmd;
  delete fPseudorapitityCutCmd;
  delete fCollisionPtMeanCmd;
  delete fFragmentationPtMeanCmd;
  delete fVerticalCrossingCmd;
  delete fHorizontalCrossingCmd;
  delete fProjectBeamCmd;
  delete fRandomizeRPCmd;
  delete fnPrimariesCmd;
  delete fMinSpectatorsCmd;
  delete fMaxSpectatorsCmd;
  delete fBeamPosCmd;
}

/*
 *
 */
void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == fBeamTypeCmd){
    newValue.toLower();
    fGenerator->SetBeamType(newValue);
  }
  if(command == fCRMCmodelCmd){
    newValue.toLower();
    fGenerator->SetGeneratorModel(newValue);
  }
  if(command == fInputFileCmd){
    fGenerator->SetInputFileName(newValue);
  }
  if(command == fPtDistCmd){
    newValue.toLower();
    fGenerator->SetFragmentationPtDist(newValue);
  }
  if(command == fMultiplicityDistCmd){
    newValue.toLower();
    fGenerator->SetMultiplicityDist(newValue);
  }
  else if(command == fPseudorapitityCutCmd){
    fGenerator->SetPseudorapidityCut( fPseudorapitityCutCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fCollisionPtMeanCmd){
    fGenerator->SetCollisionPtMean( fCollisionPtMeanCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fFragmentationPtMeanCmd){
    fGenerator->SetFragmentationPtMean( fFragmentationPtMeanCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fVerticalCrossingCmd){
    fGenerator->SetVerticalCrossingAngle( fVerticalCrossingCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fHorizontalCrossingCmd){
    fGenerator->SetHorizontalCrossingAngle( fHorizontalCrossingCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fProjectBeamCmd){
    fGenerator->SetProjectionPlane( fProjectBeamCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fRandomizeRPCmd){
    fGenerator->SetRandomizeRP( fRandomizeRPCmd->GetNewBoolValue(newValue) );
  }
  else if(command == fBeamPosCmd){
    fGenerator->SetBeamPos( new G4ThreeVector( fBeamPosCmd->GetNew3VectorValue(newValue) ) );
  }
  else if(command == fnPrimariesCmd){
    fGenerator->SetnPrimaries( fnPrimariesCmd->GetNewIntValue(newValue) );
  }
  else if(command == fMinSpectatorsCmd){
    fGenerator->SetMinNspectators( fMinSpectatorsCmd->GetNewIntValue(newValue) );
  }
  else if(command == fMaxSpectatorsCmd){
    fGenerator->SetMaxNspectators( fMaxSpectatorsCmd->GetNewIntValue(newValue) );
  }
}
