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
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef FiberSD_h
#define FiberSD_h 1

#include "G4VSensitiveDetector.hh"
#include "FiberHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class FiberSD : public G4VSensitiveDetector
{
public:
  FiberSD(G4String, G4int,G4bool);
  ~FiberSD();

  void HistInitialize();

  void   Initialize     ( G4HCofThisEvent* );
  G4bool ProcessHits    ( G4Step*, G4TouchableHistory* );
  void   EndOfEvent     ( G4HCofThisEvent* );
  void   SetReducedTree ( G4int _nFibers, G4int _nChannels );
  void   SetMLReducedTree (G4int _nFibers, G4int _nChannels );
  void   FillTimeVector ( G4int fiberNo, G4double time, G4int weight = 1 );
  void   FillChannelTimeVector ( G4int channelNo, G4double time, G4int weight = 1 );

  inline G4bool   OpticalIsOn            ( ){ return OPTICAL;        }
  inline G4bool   IsZDC                  ( ){ return ZDC;            }
  inline G4bool   IsRPD                  ( ){ return RPD;            }
  inline G4bool   IsReduced              ( ){ return REDUCED_TREE;   }
  inline G4bool   IsMLReduced            ( ){ return ML_REDUCED_TREE;}
  inline G4double GetTopOfVolume         ( ){ return m_topOfVolume;  }
  inline G4int    GetNCherenkovs         ( ){ return m_nCherenkovs;  }
  inline G4int    GetModNum              ( ){ return m_modNum;       }
  inline G4int    GetNhits               ( ){ return m_nHits;        }
  inline void     SetTopOfVolume         ( G4double _top  ){ m_topOfVolume = _top; }
  inline void     SetBottomOfVolume         ( G4double _bottom  ){ m_bottomOfVolume = _bottom; }
  inline void     SetnFibers             ( G4int _nFibers ){ m_nFibers = _nFibers;}
  inline G4int    GetRPDChannelMapping (G4int _fiberNum) {return _fiberNum / m_nFiberPerChannel;} 
  inline G4int    GetZDCChannelMapping (G4int _fiberNum) {return 0;} 
  G4double GetIncidenceAngle(const G4Step *aStep);
  G4double GetPostStepIncidenceAngle(const G4Step *aStep);
  G4double GetTrackIncidenceAngle(const G4Step *aStep);
  inline std::vector< G4double > * GetIncidenceAngleVec() {return m_incidenceAngleVec;}
  
private:
  int HCID;
  G4double m_modCoreIndexRefraction;
  FiberHitsCollection* fiberCollection;
  G4bool is_downward_light;
  G4bool is_upward_light;
  G4int m_modNum;
  G4int m_nCherenkovs;
  G4int m_nFibers;
  G4int m_nChannels;
  G4int m_nFiberPerChannel;
  G4int m_nHits;
  G4bool OPTICAL;
  G4bool REDUCED_TREE;
  G4bool ML_REDUCED_TREE;
  G4bool ZDC, RPD;
  G4double m_topOfVolume;
  G4double m_bottomOfVolume;
  std::vector< G4int >* m_genCherenkovVec;
  std::vector< G4int >* m_cherenkovVec;
  std::vector< G4double > * m_yOriginVec;
  std::vector< G4double > * m_energyVec;
  std::vector< G4int > * m_channelVec;
  std::vector< G4double >* m_timeVec;
  std::vector< G4double >* m_timeHitVec;
  
  std::vector< G4double >* m_incidenceAngleVec;
  // Number of TIR reflections photon has taken in fiber
  G4int    m_TIR_count;
  G4int    m_prevTrackID;
  G4bool   m_lastFiveCore;
  G4bool   m_isTIR;
  G4double m_yOrigin;
  G4double m_firstStepTime;
  G4double m_avgStepTime;
  G4double m_avgTimeLength;
  G4double m_avgStepLength;
  // Average vertical step length
  G4double m_avgYStepLength;
  // Number of upward reflections before checking capture condition 
  //  static const int m_captureThreshold = 7;
  static const int m_captureThreshold = 5;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
