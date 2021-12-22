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
// $Id: SteppingAction.hh 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file SteppingAction.hh
/// \brief Definition of the SteppingAction class

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <vector>

/// Stepping action class
///

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction();
  virtual ~SteppingAction();

  inline  void SetOpticalFlag( G4bool arg ){ OPTICAL = arg; }
  inline  void SetPI0Flag( G4bool arg ){ PI0 = arg; }

  inline  void SetLastStepVec( std::vector< G4ThreeVector >* vec, std::vector< int >* _lastStepPidVec){ m_lastStepVec = vec;  m_lastStepPidVec = _lastStepPidVec; }
  inline  void SetPi0Mom( std::vector< G4ThreeVector >* vec ){ m_Pi0Mom = vec; }
  inline  void SetPi0Vertex( std::vector< G4ThreeVector >* vec ){ m_Pi0Vert = vec; }

  // method from the base class
  virtual void UserSteppingAction(const G4Step*);

private:
  std::vector< G4ThreeVector >* m_lastStepVec;
  std::vector< G4ThreeVector >* m_Pi0Mom;
  std::vector< G4ThreeVector >* m_Pi0Vert;
  std::vector< int >* m_lastStepPidVec;

  G4int prevTrackID;
  G4int TIR_count;
  G4bool   OPTICAL;
  G4bool   PI0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
