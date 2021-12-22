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
// $Id: DetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "ModTypeZDC.hh"
#include "ModTypeRPD.hh"
#include "XMLSettingsReader.hh"
#include "Materials.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4Cache.hh"
#include "G4MagneticField.hh"

#include <vector>


class G4Box;
class G4Para;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

class G4UniformMagField;
class PurgMagTabulatedField3D;

class G4Region;

class GFlashHomoShowerParameterisation;
class GFlashShowerModel;
class GFlashHitMaker;
class GFlashParticleBounds;


class Survey {

 public :

	/** Type of detector - ZDC or RPD **/
    G4String detector;
    /** x_pos for this survey */
    double x_pos;
	/** y_pos for this survey */
    double y_pos;
	/** z_pos for this survey */
    double z_pos;
	/** cos_x for this survey */
    double cos_x;
	/** cos_y for this survey */
    double cos_y;
	/** cos_z for this survey */
    double cos_z;

};

class Alignment {

 public:

    /** X position of the Desy Table **/
    double x_table;
    /** Y position of the Desy Table **/
    double y_table;
    /** First detector met by the beam **/
    std::string upstream_Det;
    /** Second detector met by the beam **/
    std::string mid_Det;
    /** Third detector met by the beam **/
    std::string downstream_Det;
    /** GOLIATH magnet status **/
    bool magnet_On;
    /** Target in **/
    bool target_In;
    /** Lead absorber in **/
    bool lead_In;

};

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  DetectorConstruction(G4bool);
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual G4VPhysicalVolume* ConstructWorldVolume(G4double x, G4double y, G4double z);
  virtual G4VPhysicalVolume* ConstructSPSTestBeam();
  virtual G4VPhysicalVolume* ManualConstruction();
  virtual void ConstructSDandField();

  virtual void    LoadConfigurationFile ( G4String _inFile = "" );
  virtual void    LoadAlignmentFile     ( G4String _inFile = "" );
  virtual void    AddZDC                ( G4ThreeVector* position = NULL );
  virtual void    AddRPD                ( G4ThreeVector* position = NULL );
  virtual void    SetConfigFileName     ( G4String _name ){ m_configFileName = _name; }
  virtual void    SetRunNumber          ( G4int _run ){ m_runNumber = _run; }
  inline  void    SetOverlapsFlag       ( G4bool arg ){CHECK_OVERLAPS = arg;}
  inline  void    SetOpticalFlag        ( G4bool arg ){OPTICAL = arg;}
  inline  void    ForcePosition         ( G4bool arg ){ForceDetectorPosition = arg;}
  inline  void    SetPI0Flag            ( G4bool arg ){PI0 = arg;}

  virtual Survey* GetSurvey             ( G4String name );
  inline  G4bool  GetOverlapsFlag       ( ){ return CHECK_OVERLAPS;  }
  inline  G4bool  GetOpticalFlag        ( ){ return OPTICAL;         }
  inline  G4bool  GetPI0Flag            ( ){ return PI0;             }
  inline  std::vector< ModTypeZDC* >* GetZDCvec ( ){ return &m_ZDCvec; }
  inline  std::vector< ModTypeRPD* >* GetRPDvec ( ){ return &m_RPDvec; }

  //Manual World Volume
  inline  void SetWorldDimensions      ( G4ThreeVector* vec ){ ConstructWorldVolume( vec->x(), vec->y(), vec->z() ); }

  //For manual ZDCs
  inline  void SetZDCPosition          ( G4ThreeVector* vec ){ m_ZDCvec.at(currentZDC-1)->SetPosition(vec);           }
  inline  void SetZDCFiberDiameters    ( G4ThreeVector* vec ){ m_ZDCvec.at(currentZDC-1)->SetFiberDiameters(vec);     }
  inline  void SetZDCAbsorberDimensions( G4ThreeVector* vec ){ m_ZDCvec.at(currentZDC-1)->SetAbsorberDimensions(vec); }
  inline  void SetZDCnAbsorbers        ( G4int          arg ){ m_ZDCvec.at(currentZDC-1)->SetnAbsorbers(arg);         }
  inline  void SetZDCHousingThickness  ( G4double       arg ){ m_ZDCvec.at(currentZDC-1)->SetHousingThickness(arg);   }
  inline  void SetZDCGapThickness      ( G4double       arg ){ m_ZDCvec.at(currentZDC-1)->SetGapThickness(arg);       }
  inline  void SetZDCSteelAbsHeight    ( G4double       arg ){ m_ZDCvec.at(currentZDC-1)->SetSteelAbsHeight(arg);     }
  inline  void SetZDCOpticalFlag       ( G4bool         arg ){ m_ZDCvec.at(currentZDC-1)->SetOpticalFlag(arg);        }
  inline  void SetZDCOverlapsFlag      ( G4bool         arg ){ m_ZDCvec.at(currentZDC-1)->SetOverlapsFlag(arg);       }
  inline  void SetZDCReducedTreeFlag   ( G4bool         arg ){ m_ZDCvec.at(currentZDC-1)->SetReducedTreeFlag(arg);    }
  inline  void SetZDCMLReducedTreeFlag   ( G4bool         arg ){ m_ZDCvec.at(currentZDC-1)->SetMLReducedTreeFlag(arg);    }
  inline  void SetZDCHousingMaterial   ( G4String       arg ){ m_ZDCvec.at(currentZDC-1)->SetHousingMaterial(arg);    }
  inline  void SetZDCAbsorberMaterial  ( G4String       arg ){ m_ZDCvec.at(currentZDC-1)->SetAbsorberMaterial(arg);   }
  inline  void SetCurrentZDC           ( G4int          arg ){ currentZDC = arg; }
  virtual void DuplicateZDC            ( G4int       module );


  //For manual RPDs
  inline  void SetRPDPosition          ( G4ThreeVector* vec ){ m_RPDvec.at(currentRPD-1)->SetPosition(vec);        }
  inline  void SetRPDFiberDimensions   ( G4ThreeVector* vec ){ m_RPDvec.at(currentRPD-1)->SetFiberDiameters(vec);  }
  inline  void SetRPDRotation          ( G4double       arg ){ m_RPDvec.at(currentRPD-1)->SetRPDRotation(arg);     }
  inline  void SetRPDHousingThickness  ( G4double       arg ){ m_RPDvec.at(currentRPD-1)->SetHousingThickness(arg);}
  inline  void SetRPDFiberPitchX       ( G4double       arg ){ m_RPDvec.at(currentRPD-1)->SetFiberPitchX(arg);     }
  inline  void SetRPDFiberPitchZ       ( G4double       arg ){ m_RPDvec.at(currentRPD-1)->SetFiberPitchZ(arg);     }
  inline  void SetRPDTileSize          ( G4double       arg ){ m_RPDvec.at(currentRPD-1)->SetTileSize(arg);        }
  inline  void SetRPDMinWallThickness  ( G4double       arg ){ m_RPDvec.at(currentRPD-1)->SetMinWallThickness(arg);}
  inline  void SetRPDReadoutDistance   ( G4double       arg ){ m_RPDvec.at(currentRPD-1)->SetReadoutDistance(arg); }
  inline  void SetRPDDetectorType      ( G4String       arg ){ m_RPDvec.at(currentRPD-1)->SetDetectorType(arg);    }
  inline  void SetRPDOpticalFlag       ( G4bool         arg ){ m_RPDvec.at(currentRPD-1)->SetOpticalFlag(arg);     }
  inline  void SetRPDOverlapsFlag      ( G4bool         arg ){ m_RPDvec.at(currentRPD-1)->SetOverlapsFlag(arg);    }
  inline  void SetRPDReducedTreeFlag   ( G4bool         arg ){ m_RPDvec.at(currentRPD-1)->SetReducedTreeFlag(arg);    }
  inline  void SetRPDMLReducedTreeFlag ( G4bool         arg ){ m_RPDvec.at(currentRPD-1)->SetMLReducedTreeFlag(arg);    }
  inline  void SetRPD_NRows            ( G4int          arg ){ m_RPDvec.at(currentRPD-1)->SetNRows(arg);           }
  inline  void SetRPD_NColumns         ( G4int          arg ){ m_RPDvec.at(currentRPD-1)->SetNColumns(arg);           }
  inline  void SetRPD_NCyclesPerTile   ( G4int          arg ){ m_RPDvec.at(currentRPD-1)->SetNCyclesPerTile(arg);           }
  inline  void SetCurrentRPD           ( G4int          arg ){ currentRPD = arg; }
  virtual void DuplicateRPD            ( G4int       module );


protected:
  G4ThreeVector*          m_WorldDimensions;
  G4String                m_configFileName;
  Materials*              m_materials;

  /* World objects */
  G4Box*                  m_solidWorld;
  G4LogicalVolume*        m_logicWorld;
  G4VPhysicalVolume*      m_physWorld;

  /* Lead target objects */
  G4LogicalVolume*        logic_leadTarget;
  G4LogicalVolume*        logic_leadBlock;

  /* Configuration parser objects and storage */
  XMLSettingsReader*      m_XMLparser;
  Alignment* 	            m_alignment;
  Survey*                 m_survey;
  G4int                   m_runNumber;
  std::vector < Survey* > m_surveyEntries;

  /* Number of each detector */
  std::vector< ModTypeZDC* > m_ZDCvec;
  std::vector< ModTypeRPD* > m_RPDvec;
  G4int currentZDC;
  G4int currentRPD;


private:
	G4Cache<G4MagneticField*> fField;  //pointer to the thread-local fields

  /* Run condition flags */
  G4bool m_GFlash;
  G4bool OPTICAL;
  G4bool ForceDetectorPosition;
  G4bool CHECK_OVERLAPS;
  G4bool PI0;

  //Fast simulation
  std::vector< G4Region* > m_Region;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
