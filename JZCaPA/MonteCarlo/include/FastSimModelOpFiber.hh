#ifndef FastSimModelOpFiber_h
#define FastSimModelOpFiber_h 1

#include "G4VFastSimulationModel.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4GenericMessenger.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4Material.hh"


struct FastFiberData {
public:
  FastFiberData(G4int, G4double, G4ThreeVector);
  ~FastFiberData() {}

  FastFiberData& operator=(const FastFiberData &right);
  void reset();
  inline void setDeltaPosition(G4double delta) {deltaPosition = delta;}

  G4int trackID;
  G4double globalTime;
  G4ThreeVector globalPosition;
  G4double deltaPosition;
private:

};


class FastSimModelOpFiber : public G4VFastSimulationModel {
public:
  FastSimModelOpFiber(G4String, G4Region*, G4double, G4int, G4int);
  ~FastSimModelOpFiber();

  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  virtual G4bool ModelTrigger(const G4FastTrack&);
  virtual void DoIt(const G4FastTrack&, G4FastStep&);

  void SetCoreMaterial(G4Material* mat) { fCoreMaterial = mat; }
  G4bool Compare_doubles(double x, double y, double epsilon = 0.0000001f);
  G4double GetTrackIncidenceAngle(const G4Track *track);
  G4double GetPrestepIncidenceAngle(const G4Track *track);
  G4double GetPoststepIncidenceAngle(const G4Track *track);
  
private:
  void DefineCommands();

  bool checkTotalInternalReflection(const G4Track* track);
  G4double CalculateVelocityForOpticalPhoton(const G4Track* track);
  void setOpBoundaryProc(const G4Track* track);
  void getCoreMaterial(const G4Track* track);
  void reset();

  FastFiberData mStepPrevious;
  FastFiberData mStepCurrent;
  
  G4GenericMessenger* fMessenger;
  G4OpBoundaryProcess* fOpBoundaryProc;
  G4Material* fCoreMaterial;
  G4bool fProcAssigned;
  G4int fSafetySteps;
  G4int fSafetyReflections;
  G4double fTrkLength;
  G4double fNtransport;
  G4double fMaxTransport;
  G4double fTransportUnit;
  G4ThreeVector fFiberAxis;
  G4ThreeVector fDownwardFiberAxis;
  G4bool fActiveArea;
  G4bool fReadoutArea;
  G4bool fIsDownwardLight;
  G4double fReadoutFiberLen;
  G4ThreeVector fFiberEnd;
  G4bool fKill;
  G4bool fTransported;
  G4bool f_is_run3;
  G4bool f_is_run4;
  G4int fNtotIntRefl;
  G4int fTrackId;
  G4int fNFibers;
  G4int fNChannels;
  G4int fNFibersPerChannel;
};

#endif




