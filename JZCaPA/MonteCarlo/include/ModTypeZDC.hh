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
// Author: Michael Phipps

#ifndef ModTypeZDC_h
#define ModTypeZDC_h 1


#include "globals.hh"
#include "G4PVPlacement.hh"
#include "Materials.hh"

#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4Material;

/// Detector construction class to define materials and geometry.

class ModTypeZDC
{
public:
  ModTypeZDC(const int, ModTypeZDC*);
  ModTypeZDC(const int, G4LogicalVolume*, G4ThreeVector* );
  ~ModTypeZDC();

  virtual void  Construct();
  virtual void  ConstructSDandField();

  virtual void  ConstructDetector();


  inline  void  SetPosition            ( G4ThreeVector* vec ){ delete m_pos;       m_pos       = vec; }
  inline  void  SetFiberDiameters      ( G4ThreeVector* vec ){ delete m_fiberDiam; m_fiberDiam = vec; }
  inline  void  SetAbsorberDimensions  ( G4ThreeVector* vec ){ delete m_absDim;    m_absDim    = vec; }
  inline  void  SetnAbsorbers          ( G4int          arg ){ m_nAbsorbers        = arg; }
  inline  void  SetSteelAbsHeight      ( G4double       arg ){ STEEL_ABSORBER = true; m_SteelAbsHeight = arg; }
  inline  void  SetHousingThickness    ( G4double       arg ){ m_HousingThickness  = arg; }
  inline  void  SetGapThickness        ( G4double       arg ){ m_GapThickness      = arg; }
  inline  void  SetOpticalFlag         ( G4bool         arg ){ OPTICAL             = arg; }
  inline  void  SetOverlapsFlag        ( G4bool         arg ){ CHECK_OVERLAPS      = arg; }
  inline  void  SetReducedTreeFlag     ( G4bool         arg ){ REDUCED_TREE        = arg; }
  inline  void  SetMLReducedTreeFlag   ( G4bool         arg ){ ML_REDUCED_TREE     = arg; }
  virtual void  SetHousingMaterial     ( G4String  material );
  virtual void  SetAbsorberMaterial    ( G4String  material );



  inline  G4LogicalVolume* GetAbsorberLogicalVolume( ){ return m_WLogical;    }
  inline  G4Material*      GetAbsorberMaterial     ( ){ return m_matAbsorber; }
  inline  G4ThreeVector*   GetPosition             ( ){ return m_pos;         }
  inline  G4int            GetModNum               ( ){ return m_modNum;      }
  inline  G4int            GetnFibers              ( ){ return m_nFibers;     }
  inline  G4bool           GetOpticalFlag          ( ){ return OPTICAL;       }
  inline  G4bool           GetReducedTreeFlag      ( ){ return REDUCED_TREE;  }
  inline  G4bool           GetMLReducedTreeFlag      ( ){ return ML_REDUCED_TREE;  }

protected:
  const G4int      m_modNum;
  G4int            m_nAbsorbers;
  G4int            m_nFibers;
  G4ThreeVector*   m_pos;
  G4ThreeVector*   m_fiberDiam;
  G4ThreeVector*   m_absDim;
  G4double         m_HousingThickness;
  G4double         m_GapThickness;
  G4double         m_SteelAbsHeight;
  G4double         m_topOfVolume;
  G4double         m_bottomOfVolume;
  G4bool           OPTICAL;
  G4bool           CHECK_OVERLAPS;
  G4bool           STEEL_ABSORBER;
  G4bool           REDUCED_TREE;
  G4bool           ML_REDUCED_TREE;
  Materials*       m_materials;
  G4Material*      m_matAbsorber;
  G4Material*      m_matHousing;
  G4LogicalVolume* m_logicMother;

protected:

  G4VSolid*          m_ModuleBox;
  G4LogicalVolume*   m_ModuleLogical;
  G4VPhysicalVolume* m_ModulePhysical;

  G4VSolid*          m_AirScoringBox;
  G4LogicalVolume*   m_AirScoringLogical;
  G4VPhysicalVolume* m_AirScoringPhysical;

  G4VSolid*          m_HousingBox;
  G4LogicalVolume*   m_HousingLogical;
  G4VPhysicalVolume* m_HousingPhysical;

  G4VSolid*          m_W;
  G4LogicalVolume*   m_WLogical;
  std::vector < G4VPhysicalVolume* > m_WPhysical;

  G4VSolid*          m_Steel;
  G4LogicalVolume*   m_SteelLogical;
  std::vector < G4VPhysicalVolume* > m_SteelPhysical;

  // Vertical quartz strips (these strips are full -- no partial segments)
  G4VSolid*          m_FiberCoreTube;
  G4LogicalVolume*   m_FiberCoreLogical;
  std::vector< std::vector < G4VPhysicalVolume* > > m_FiberCorePhysical;

  G4VSolid*          m_CladdingTube;
  G4LogicalVolume*   m_CladdingLogical;
  std::vector< std::vector < G4VPhysicalVolume* > > m_FiberCladdingPhysical;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
