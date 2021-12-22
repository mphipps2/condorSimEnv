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
/// \file include/PrimaryGeneratorMessenger_h.hh
/// \brief Definition of the PrimaryGeneratorMessenger_h class
//
//
//

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3VectorAndUnit;


class PrimaryGeneratorMessenger: public G4UImessenger{
  public:

    PrimaryGeneratorMessenger(PrimaryGeneratorAction* );
   ~PrimaryGeneratorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:

    PrimaryGeneratorAction*    fGenerator;

    G4UIdirectory*             fGeneratorDir;
    G4UIdirectory*             fLHCDir;
    G4UIdirectory*             fSPSDir;
    G4UIdirectory*             fFNALDir;

    G4UIcmdWithAString*        fBeamTypeCmd;
    G4UIcmdWithAString*        fCRMCmodelCmd;
    G4UIcmdWithAString*        fInputFileCmd;
    G4UIcmdWithAString*        fPtDistCmd;
    G4UIcmdWithAString*        fMultiplicityDistCmd;
    G4UIcmdWithADouble*        fPseudorapitityCutCmd;
    G4UIcmdWithADouble*        fCollisionPtMeanCmd;
    G4UIcmdWithADouble*        fFragmentationPtMeanCmd;
    G4UIcmdWithADoubleAndUnit* fVerticalCrossingCmd;
    G4UIcmdWithADoubleAndUnit* fHorizontalCrossingCmd;
    G4UIcmdWithADoubleAndUnit* fProjectBeamCmd;
    G4UIcmdWithAnInteger*      fnPrimariesCmd;
    G4UIcmdWithAnInteger*      fMinSpectatorsCmd;
    G4UIcmdWithAnInteger*      fMaxSpectatorsCmd;
    G4UIcmdWithABool*          fRandomizeRPCmd;
    G4UIcmdWith3VectorAndUnit* fBeamPosCmd;


};
#endif
