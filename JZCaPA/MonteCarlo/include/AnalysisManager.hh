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
/// \file AnalysisManager.hh
/// \brief Selection of the analysis technology
/// \author Chad Lantz
/// \date 16 April 2020

#ifndef AnalysisManager_h
#define AnalysisManager_h 1

#include "g4root.hh"
//#include "g4cvs.hh"
//#include "g4xml.hh"

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4ThreeVector.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class AnalysisManager
{
public:
  static AnalysisManager* getInstance(void);
  ~AnalysisManager();

  void Book( G4String fileName );
  void Save();

  void FillNtuples( );
  void FillZDCnCherenkovs( int zdcNo, int nCherenkovs );
  void FillRPDnCherenkovs( int rpdNo, int nCherenkovs );
  void FillEventGenTree  ( int nPart, int nSpec, int model, double impactParam );
  void FillEventGenTree  ( int nSpectators, double ptCollision, double ptBreakup, double rpAngle );

  void MakeEventDataTree ( );
  void MakeZDCTree       ( G4int zdcNo, std::vector< int >* nCherenkovVec, std::vector< double >* timeVec, G4bool reducedTree, G4bool MLReducedTree, G4bool thisIsOptical );
  void MakeRPDTree       ( G4int rpdNo, G4int modNum, G4bool reducedTree, G4bool MLReducedTree, G4bool thisIsOptical );
  void MakeEventGenTree  ( std::vector< std::vector<int>* > &intVec , std::vector< std::vector<double>* > &dblVec, G4int type );

  inline G4bool         GetOpticalFlag( ){ return OPTICAL; }
  inline G4bool         GetPI0Flag    ( ){ return PI0;     }
  std::vector< G4int >* GetFiberVector( G4bool ZDC, G4bool RPD, G4int modNum );
  std::vector< G4int >* GetFiberGenVector( G4int modNum );
  std::vector< G4double >* GetIncidenceAngleVector( G4int modNum );
  std::vector< G4double >* GetYOriginVector( G4int modNum );
  std::vector< G4double >* GetEnergyVector( G4int modNum );
  std::vector< G4int >*    GetChannelVector( G4int modNum );
  std::vector< G4double >* GetTimeVector ( G4bool ZDC, G4bool RPD, G4int modNum );
  std::vector< G4double >* GetTimeHitVector ( G4int modNum );

  inline void SetGunPosition ( G4double x, G4double y, G4double z ){ m_gunPos->set(x,y,z); }
  inline void SetEventSeeds ( unsigned int seed1, unsigned int seed2){ m_EventSeed1 = seed1; m_EventSeed2 = seed2; }
  inline void SetEventNo     ( G4int _eventNo ){ m_eventNo = _eventNo; }

  inline  std::vector< std::vector< std::vector<double> > >*  GetRPDdoubleVectors( ){return m_RPDdblVec;}
  inline  std::vector< std::vector< std::vector< int  > > >*  GetRPDintVectors   ( ){return m_RPDintVec;}
  inline  std::vector< std::vector< std::vector<double> > >*  GetZDCdoubleVectors( ){return m_ZDCdblVec;}
  inline  std::vector< std::vector< std::vector< int  > > >*  GetZDCintVectors   ( ){return m_ZDCintVec;}


private:
  AnalysisManager();
  static AnalysisManager* analysisManager;
  G4int m_eventNo;
  G4int m_eventDataNtupleNo;
  G4int m_eventGenNtupleNo;
  G4bool m_FactoryOn;
  G4bool OPTICAL;
  G4bool PI0;
  unsigned int m_EventSeed1;
  unsigned int m_EventSeed2;
  G4ThreeVector* m_gunPos;
  G4AnalysisManager* m_analysisManager;
  DetectorConstruction* m_detectorConstruction;
  std::vector< std::vector< std::vector<double> > > *m_ZDCdblVec, *m_RPDdblVec;
  std::vector< std::vector< std::vector< int  > > > *m_ZDCintVec, *m_RPDintVec;
  std::vector< std::vector< G4int  >* > m_ZDCfiberVec, m_RPDfiberVec, m_RPDfiberGenVec, m_channelVec ;
  std::vector< std::vector< G4double >* >   m_ZDCtimeVec, m_RPDtimeVec, m_YOriginVec, m_EnergyVec, m_IncidenceAngleVec, m_timeHitVec;
  std::vector< G4ThreeVector >* m_lastStepVec, *m_Pi0Mom, *m_Pi0Vert;
  std::vector< int >  m_lastStepPidVec, m_ZDCnTupleNo, m_RPDnTupleNo;
  std::vector< double > m_lastStepXVec, m_lastStepYVec, m_lastStepZVec;
  std::vector< double > m_Pi0MomX, m_Pi0MomY, m_Pi0MomZ;
  std::vector< double > m_Pi0VertX, m_Pi0VertY, m_Pi0VertZ;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
