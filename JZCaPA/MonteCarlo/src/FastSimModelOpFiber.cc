// Note: this fast physics class is called in between a step, before an action is taken. So current step length is 0 and we don't have a postStepPoint yet (by default it's equivalent to the preStepPoint)

#include "FastSimModelOpFiber.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4OpProcessSubType.hh"
#include "G4Tubs.hh"

#include "TMath.h"
#include "TRandom3.h"

#include <cmath>

FastFiberData::FastFiberData(G4int id, G4double globTime, G4ThreeVector pos) {
  trackID = id;
  globalTime = globTime;
  globalPosition = pos;
}

FastFiberData& FastFiberData::operator=(const FastFiberData &right) {
  trackID = right.trackID;
  globalTime = right.globalTime;
  globalPosition = right.globalPosition;
  deltaPosition = right.deltaPosition;

  return *this;
}

void FastFiberData::reset() {
  trackID = 0;
  globalTime = 0;
  globalPosition.set(0,0,0);
}


FastSimModelOpFiber::FastSimModelOpFiber(G4String name, G4Region* envelope, G4double readoutFiberLen, G4int nFibers, G4int nChannels)
  : G4VFastSimulationModel(name,envelope),
    mStepPrevious(0,0.,G4ThreeVector(0)),
    mStepCurrent(0,0.,G4ThreeVector(0))
{
  f_is_run3 = true;
  f_is_run4 = false;
  fOpBoundaryProc = NULL;
  fCoreMaterial = NULL;
  fProcAssigned = false;
  fSafetySteps = 2;
  //  fSafetyReflections = 2;
  //  fSafetyReflections = 3;
  fSafetyReflections = 4;
  fTrkLength = 0.;
  fNtransport = 0.;
  fTransportUnit = 0.;
  fFiberAxis = G4ThreeVector(0);
  fDownwardFiberAxis = G4ThreeVector(0);
  fKill = false;
  fNtotIntRefl = 0;
  fTransported = false;
  fTrackId = 0;
  fActiveArea = false;
  fReadoutFiberLen = readoutFiberLen;
  fNFibers = nFibers;
  fNChannels = nChannels;
  fNFibersPerChannel = fNFibers / fNChannels;
  fIsDownwardLight = false;
  DefineCommands();
}

FastSimModelOpFiber::~FastSimModelOpFiber() {}

G4bool FastSimModelOpFiber::IsApplicable(const G4ParticleDefinition& type) {
  return &type == G4OpticalPhoton::OpticalPhotonDefinition();
}

G4bool FastSimModelOpFiber::ModelTrigger(const G4FastTrack& fasttrack) {
  const G4Track* track = fasttrack.GetPrimaryTrack();

  if (fTrackId != track->GetTrackID()) {
    reset();
    fTrackId = track->GetTrackID();
  }

  getCoreMaterial(track);

  auto matPropTable = fCoreMaterial->GetMaterialPropertiesTable();

  if ( !matPropTable ) {
    std::cout << " !matPropTable " << std::endl;
    return false;
  }


  G4TouchableHandle theTouchable = track->GetTouchableHandle();
  auto fiberPos = theTouchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0.,0.,0.));
  fFiberAxis = theTouchable->GetHistory()->GetTopTransform().Inverse().TransformAxis(G4ThreeVector(0.,0.,1.));
  fDownwardFiberAxis = theTouchable->GetHistory()->GetTopTransform().Inverse().TransformAxis(G4ThreeVector(0.,0.,1.));
  fDownwardFiberAxis.rotateZ(M_PI);
  // total accummulated track length (across all prior steps)
  fTrkLength = track->GetTrackLength();
  G4Tubs* tubs = static_cast<G4Tubs*>(theTouchable->GetSolid());
  G4double fiberLen = 2.*tubs->GetZHalfLength();

  G4ThreeVector deltaPos = mStepCurrent.globalPosition - mStepPrevious.globalPosition;

  // fTransportUnit is the straightline vertical magnitude of the step size
  fTransportUnit = deltaPos.dot(fFiberAxis);


  // flag direction of light
  if ( fTransportUnit < 0. ) {
    fIsDownwardLight = true;
    // flip axis for fiber unit vector
    fFiberAxis = fDownwardFiberAxis;
  }
  else {
    fIsDownwardLight = false;
  }

  G4int channelNum = track->GetTouchableHandle()->GetCopyNumber(0) / fNFibersPerChannel;
  G4double activeFiberLen;
  // Warning: This part is hard coded for now. If fiber lengths change this needs manually changed
  if (f_is_run3) {
    if (channelNum > 11) activeFiberLen = 45.6;
    else if (channelNum <= 11 && channelNum > 7) activeFiberLen = 34.2;
    else if (channelNum <= 7 && channelNum > 3) activeFiberLen = 22.8;
    else activeFiberLen = 11.4;
  }
  else if (f_is_run4) {
    if (channelNum > 11) activeFiberLen = 38.4;
    else if (channelNum <= 11 && channelNum > 7) activeFiberLen = 28.8;
    else if (channelNum <= 7 && channelNum > 3) activeFiberLen = 19.2;
    else activeFiberLen = 9.6;
  }

  // fiberPos is vector giving fiber center in global coordinates
  // fiberAxis is a unit y vector (0,1,0)
  // fiberLen is 507 or 11.4, 22.8, 34.2, 45.6 (depending whether its a readout fiber or from active area)
  // fiberEnd is just a y transform taking you to the center of fiber end. x and z from fiberPosVec preserved in fiberEndVec
  if (f_is_run3) {
    if (std::floor(fiberLen) == 11 || std::floor(fiberLen) == 22 || std::floor(fiberLen) == 34 || std::floor(fiberLen) == 45) {   
      fActiveArea = true;
    }
    else {
      fActiveArea = false;
    }
  }
  else if (f_is_run4) {
    if (std::floor(fiberLen) == 9 || std::floor(fiberLen) == 19 || std::floor(fiberLen) == 28 || std::floor(fiberLen) == 38) {   
      fActiveArea = true;
    }
    else {
      fActiveArea = false;
    }

  }

  fFiberEnd.set(0,0,0);
  // if it's upward take light to top where it's recorded. else it goes to the bottom where it has a chance to be recaptured
  // depending which way we're going and which fiber we're in we may need to splice on either the active fiber length or the readout fiber length

  // readout area/downward light
  if (!fActiveArea && fIsDownwardLight) fFiberEnd = fiberPos + fFiberAxis*fiberLen/2 + fFiberAxis*activeFiberLen;
  // active area/downward light
  else if (fActiveArea && fIsDownwardLight) fFiberEnd = fiberPos + fFiberAxis*fiberLen/2;
  // active area/upward light
  else if (fActiveArea && !fIsDownwardLight) fFiberEnd = fiberPos + fFiberAxis*fiberLen/2 + fFiberAxis*fReadoutFiberLen;
  // active area/upward light
  else if (!fActiveArea && !fIsDownwardLight)  fFiberEnd = fiberPos + fFiberAxis*fiberLen/2 ;

  if ( !checkTotalInternalReflection(track) ) {
    return false; // nothing to do if the previous status is not total internal reflection
  }

  // toEnd is a vector giving the straight line distance from current position to center of fiber end
  auto toEnd = fFiberEnd - track->GetPosition();
  // toEndAxis gives scalar magnitude of that straight line distance to center of fiber end from current position
  double toEndAxis = toEnd.dot(fFiberAxis);
  // fMaxTransport gives the number of steps needed to reach the end
  fMaxTransport = std::floor(toEndAxis/TMath::Abs(fTransportUnit));
  fNtransport = fMaxTransport - fSafetySteps;

  // fNtransport is estimate of how many total internal reflections the track will take to reach fiber end
  // require at least n = fSafetySteps of total internal reflections at the end. ie) if we're very close to the end let geant handle things, no transport
  if ( fNtransport < 0. ) {
    reset();

    return false;
  }

  if ( matPropTable->GetProperty(kABSLENGTH) ) {
    double attLength = matPropTable->GetProperty(kABSLENGTH)->Value( track->GetDynamicParticle()->GetTotalMomentum() );
    double totalPathLength =  deltaPos.mag()*fMaxTransport;
    double survivalProb = exp((-1 * totalPathLength) / attLength);
    // initializing to 0 invokes random seed
    TRandom3 *r3 = new TRandom3(0);
    G4bool survived = r3->Binomial(1,survivalProb);
    if (!survived) {
      fKill = true;
      return true;
    }
  }

  fKill = false;

  return true;
}

void FastSimModelOpFiber::DoIt(const G4FastTrack& fasttrack, G4FastStep& faststep) {
  auto track = fasttrack.GetPrimaryTrack();

  // either a stopped track or one that would be absorbed during transport
  if (fKill) {
    faststep.KillPrimaryTrack();
    reset();
    return;
  }

  //total time traveled -- fTrkLength will go to
  // vertical magnitude of step size * number of steps remaining * fiber unit vector

  G4ThreeVector deltaPos = mStepCurrent.globalPosition - mStepPrevious.globalPosition;
  double deltaTime = mStepCurrent.globalTime - mStepPrevious.globalTime;
  //fTransportUnit: vertical magnitude of the step size;
  //fMaxTransport: # of steps to bypass
  //fFiberAxis: unit vector pointing down fiber axis (in direction of light)
  G4ThreeVector posShift = TMath::Abs(fTransportUnit)*fMaxTransport*fFiberAxis;
  double timeShift = deltaTime*fMaxTransport;
  G4ThreeVector finalPos = track->GetPosition() + posShift;
  G4double finalTime = track->GetGlobalTime() + timeShift;

  auto fiberPos = track->GetTouchableHandle()->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0.,0.,0.));

  /*
  G4double incidentAngle = fOpBoundaryProc->GetIncidentAngle();
  G4double incidentDegrees = (incidentAngle / (2*M_PI)) * 360;
  if (incidentDegrees > 90) incidentDegrees = 180 - incidentDegrees;
  */
  G4double incidentAngle2 = GetTrackIncidenceAngle(track);
  if (incidentAngle2 > 90) incidentAngle2 = 180 - incidentAngle2;

  faststep.ProposePrimaryTrackFinalPosition( finalPos, false );
  faststep.ProposePrimaryTrackFinalTime( finalTime );
  faststep.ProposePrimaryTrackFinalKineticEnergy( track->GetKineticEnergy() );
  faststep.ProposePrimaryTrackFinalMomentumDirection( track->GetMomentumDirection(), false );
  faststep.ProposePrimaryTrackFinalPolarization( track->GetPolarization(), false );
  //  fTransported = true;
  reset();

  return;
}

bool FastSimModelOpFiber::checkTotalInternalReflection(const G4Track* track) {
  if (!fProcAssigned) { // locate OpBoundaryProcess only once
    setOpBoundaryProc(track);
  }

  if ( track->GetTrackStatus()==fStopButAlive || track->GetTrackStatus()==fStopAndKill ) {
    std::cerr << " trackstatus: stopped in checkTIR " << std::endl;
    return false;
  }
  /*
  G4int theStatus = fOpBoundaryProc->GetStatus();
  // matches calculated incident angle as long as we are in core/cladding
  G4double incidentAngle = fOpBoundaryProc->GetIncidentAngle();
  G4double incidentDegrees = (incidentAngle / (2*M_PI)) * 360;
  if (incidentDegrees > 90) incidentDegrees = 180 - incidentDegrees;
  */
  G4double incidentAngle2 = GetTrackIncidenceAngle(track);
  if (incidentAngle2 > 90) incidentAngle2 = 180 - incidentAngle2;

  // 22.8 is the transition b/w active and readout fibers. The normal is relative to the center of the fiber now rather than the side. This causes this a captured track to appear uncaptured during this step
  if (f_is_run3) {
    if (fNtotIntRefl > 0) {
      if (Compare_doubles(track->GetPosition().y(),22.8)) {   
	return false;
      }
      // if light reaches end of fiber and we're not going to transport, we need to reset the counters to deal with recaptured light correctly
      else if (Compare_doubles(track->GetPosition().y(),fFiberEnd.y()) && !Compare_doubles(fFiberEnd.y(),22.8)) {
	reset();
	return false;
      }
    }
  }
  else if (f_is_run4) {
    if (fNtotIntRefl > 0) {
      if (Compare_doubles(track->GetPosition().y(),19.2)) {
	return false;
      }
      // if light reaches end of fiber and we're not going to transport, we need to reset the counters to deal with recaptured light correctly
      else if (Compare_doubles(track->GetPosition().y(),fFiberEnd.y()) && !Compare_doubles(fFiberEnd.y(),19.2)) {
	reset();
	return false;
      }
    }

  }

  // UPDATE THIS TO CALCULATE CRITICAL ANGLE FROM MATPROPERTIESTABLE
  //  if (incidentDegrees >= 82.148) {
  if (incidentAngle2 >= 82.148) {
    // this doesnt work b/c for some reason the step length isnt updated yet at this point and so it shows stepLength = 0, which it always is in the cladding to core step which registers as a StepTooSmall
    //  if ( fOpBoundaryProc->GetStatus()==G4OpBoundaryProcessStatus::TotalInternalReflection ) {
    if ( fTrackId != track->GetTrackID() ) { // reset everything if when encountered a different track
      reset();
    }

    fTrackId = track->GetTrackID();
    fNtotIntRefl++;
    G4int trackID = track->GetTrackID();
    G4double globalTime = track->GetGlobalTime();
    G4ThreeVector globalPosition = track->GetPosition();
    mStepPrevious = mStepCurrent;
    mStepCurrent = FastFiberData(trackID,globalTime,globalPosition);
    auto deltaPos = mStepCurrent.globalPosition - mStepPrevious.globalPosition;
    mStepCurrent.setDeltaPosition(deltaPos.mag());
    if ( fNtotIntRefl > fSafetyReflections ) { // require at least n = fSafety of total internal reflections at the beginning
      // the 22.8 check above should get rid of any cases of the step length changing within a TIR
      if (!Compare_doubles(mStepCurrent.deltaPosition,mStepPrevious.deltaPosition)) {
	std::cout << "Warning: step length has changed for this TIR track!!!! Current step " << mStepCurrent.deltaPosition << " Previous Step: " << mStepPrevious.deltaPosition << " y: " << globalPosition.y() << " if fiberLength changed make sure you change hardcoded values in FastSimModelOpFiber! " << std::endl;
	return false;
      }
      return true;
    }
    else {
      //      std::cout << " exiting tir due to less than safety. nTotalIntRefl " << fNtotIntRefl << " safety " << fSafety << std::endl;
    }
  }
  else {
    // to be safe, reset whenever we dont have TIR
    reset();
  }

  return false;
}


void FastSimModelOpFiber::setOpBoundaryProc(const G4Track* track) {
  G4ProcessManager* pm = track->GetDefinition()->GetProcessManager();
  auto postStepProcessVector = pm->GetPostStepProcessVector();

  for (int np = 0; np < (int) postStepProcessVector->entries(); np++) {
    auto theProcess = (*postStepProcessVector)[np];
    if ( theProcess->GetProcessType()!=fOptical || theProcess->GetProcessSubType()!=G4OpProcessSubType::fOpBoundary ) continue;

    fOpBoundaryProc = (G4OpBoundaryProcess*)theProcess;
    fProcAssigned = true;

    break;
  }

  return;
}

G4double FastSimModelOpFiber::CalculateVelocityForOpticalPhoton(const G4Track* track) {
  G4double velocity = CLHEP::c_light;
  G4bool update_groupvel = false;
  G4MaterialPropertyVector* groupvel = nullptr;

  // check if previous step is in the same volume
  // and get new GROUPVELOCITY table if necessary
  if ( (fCoreMaterial!=nullptr) && (groupvel==nullptr) ) {
    if ( fCoreMaterial->GetMaterialPropertiesTable() != nullptr ) {
      groupvel = fCoreMaterial->GetMaterialPropertiesTable()->GetProperty("GROUPVEL");
    }
    update_groupvel = true;
  }

  if ( groupvel != nullptr ) {
    // light velocity = c/(rindex+d(rindex)/d(log(E_phot)))
    // values stored in GROUPVEL material properties vector
    G4double current_momentum = track->GetDynamicParticle()->GetTotalMomentum();
    if ( update_groupvel ) {
      velocity = groupvel->Value(current_momentum);
    }
  }

  return velocity;
}

void FastSimModelOpFiber::getCoreMaterial(const G4Track* track) {
  auto physVol = track->GetVolume();
  auto logVol = physVol->GetLogicalVolume();

  if ( logVol->GetNoDaughters()==0 ) {
    fCoreMaterial = logVol->GetMaterial();
  } else {
    auto corePhys = logVol->GetDaughter(0);
    fCoreMaterial = corePhys->GetLogicalVolume()->GetMaterial();
  }
}

void FastSimModelOpFiber::reset() {
  fTrkLength = 0.;
  fNtransport = 0.;
  fTransportUnit = 0.;
  //  fTransported = false;
  fFiberAxis = G4ThreeVector(0);
  fDownwardFiberAxis = G4ThreeVector(0);
  fIsDownwardLight = false;
  fKill = false;
  fNtotIntRefl = 0;
  fTrackId = 0;
  fActiveArea = false;
  mStepPrevious.reset();
  mStepCurrent.reset();
}

void FastSimModelOpFiber::DefineCommands() {
  fMessenger = new G4GenericMessenger(this, "/DRsim/fastOp/", "fast Op transport model control");
  G4GenericMessenger::Command& safetyCmd = fMessenger->DeclareProperty("safety",fSafetyReflections,"min number of total internal reflection");
  safetyCmd.SetParameterName("safety",true);
  safetyCmd.SetDefaultValue("2.");
}

G4bool FastSimModelOpFiber::Compare_doubles(double x, double y, double epsilon) {
  //  std::cout << " current deltaP " << x << " previous deltaP " << y << " difference " << fabs(x - y) << " epsilon " << epsilon << std::endl;
  if(fabs(x - y) < epsilon)
    return true; //they are same
  return false; //they are not same
}


G4double FastSimModelOpFiber::GetTrackIncidenceAngle(const G4Track *track) {
  //  G4ThreeVector photonDirection = postStep->GetMomentum() / postStep->GetMomentum().mag();
  G4ThreeVector photonDirection = track->GetMomentumDirection();
  G4ThreeVector stepPos = track->GetPosition();

  const G4VTouchable *touchable = track->GetTouchable();

  const G4RotationMatrix *rotation = touchable->GetRotation();
  G4RotationMatrix rotation_inv = rotation->inverse();
  G4ThreeVector translation = touchable->GetTranslation();
  G4VSolid *sector = touchable->GetSolid();

  G4ThreeVector posLocal = *rotation * (stepPos - translation);
  G4ThreeVector normal =  - sector->SurfaceNormal(posLocal);

  G4ThreeVector photonDirectionLocal = *rotation * photonDirection;

  G4double incidenceAngle = acos( normal.dot(photonDirectionLocal) );
  G4double incidenceDegrees = (incidenceAngle / (2*M_PI)) * 360;

  return incidenceDegrees;
}

G4double FastSimModelOpFiber::GetPrestepIncidenceAngle(const G4Track *track)
{
  G4StepPoint *preStep = track->GetStep()->GetPreStepPoint();
  G4ThreeVector photonDirection = preStep->GetMomentum() / preStep->GetMomentum().mag();
  G4ThreeVector stepPos = preStep->GetPosition();

  const G4VTouchable *touchable = preStep->GetTouchable();

  const G4RotationMatrix *rotation = touchable->GetRotation();
  G4RotationMatrix rotation_inv = rotation->inverse();
  G4ThreeVector translation = touchable->GetTranslation();
  G4VSolid *sector = touchable->GetSolid();

  G4ThreeVector posLocal = *rotation * (stepPos - translation);
  G4ThreeVector normal =  - sector->SurfaceNormal(posLocal);

  G4ThreeVector photonDirectionLocal = *rotation * photonDirection;

  G4double incidenceAngle = acos( normal.dot(photonDirectionLocal) );
  G4double incidenceDegrees = (incidenceAngle / (2*M_PI)) * 360;

  return incidenceDegrees;
}

G4double FastSimModelOpFiber::GetPoststepIncidenceAngle(const G4Track *track)
{
  G4StepPoint *postStep = track->GetStep()->GetPostStepPoint();
  G4ThreeVector photonDirection = postStep->GetMomentum() / postStep->GetMomentum().mag();
  G4ThreeVector stepPos = postStep->GetPosition();

  const G4VTouchable *touchable = postStep->GetTouchable();

  const G4RotationMatrix *rotation = touchable->GetRotation();
  G4RotationMatrix rotation_inv = rotation->inverse();
  G4ThreeVector translation = touchable->GetTranslation();
  G4VSolid *sector = touchable->GetSolid();

  G4ThreeVector posLocal = *rotation * (stepPos - translation);
  G4ThreeVector normal =  - sector->SurfaceNormal(posLocal);

  G4ThreeVector photonDirectionLocal = *rotation * photonDirection;

  G4double incidenceAngle = acos( normal.dot(photonDirectionLocal) );
  G4double incidenceDegrees = (incidenceAngle / (2*M_PI)) * 360;

  return incidenceDegrees;
}
