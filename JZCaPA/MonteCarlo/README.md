## ZDC R&D MonteCarlo

 Instructions from Mike Phipps (taken from his repo, michael.william.phipps@cern.ch). Updated by Chad Lantz in April 2020 for MC2.0 update

 This simulation is intended to allow standalone (non-ATHENA) simulation of the current ZDC geometry, as well as configurable geometries that could be used for prototype studies. The structure is intended to closely resemble examples/extended/analysis/AnaEx01 with the addition of several simulated environments and configurations stored in XML format.

 #### DetectorConstruction

   - DetectorConstruction configures the simulation geometry based on the macro commands executed in geometry.mac. Multiple environments are/will be available for selection with the command ```/Detector/SelectEnvionment```.

   - For test beam environments, the particular environment must be set via ```/Detector/SelectEnvionment```, as well as the run number via ```/Detector/RunNumber```

   - Manual construction can also be accomplished by using the macro command ```/Detector/ForcePosition true```. geometry.mac currently contains all detector construction commands available to the user. These commands include things like size and composition of absorbers, fibers and housing.

   - **Please note:** When positioning manually it is up to the user to ensure all detector geometry as well as the beam source is located inside the mother volume.


#### PhysicsList

   - This package defaults to the standard ATLAS hadronic physics list which is is FTFP_Bert (Fritiof string model (>\~5 GeV); Bertini-style cascade (<~10 GeV)). This can be changed via macro command ```/Physics/SelectList```. The other hadronic list currently configured is QGSP_Bert.

   - By default, FTFP_Bert comes with the following EM and low energy lists also turned on: G4EmStandard, G4GAmmaLeptoNuclearPhys, Decay, hElasticWEL_CHIPS, hInelastic FTFP_BERT, stopping, ionInelasticFTFP_BIC, neutronTrackingCut.

   - Pion and muon decay have also been turned on within the PhysicsList.cc module

#### ActionInitialization

   - Action initialization begins within ActionInitialization.cc, where PrimaryGeneratorAction, RunAction, EventAction and SteppingAction are all declared.

   - The RunAction module is just used to initialize the output via AnalysisManager and pass the optical flag to SteppingAction at the start of each run.

   - The EventAction module is used primarily at the end of each event process hit collections and fill the output nTuples.

   - The SteppingAction module is used to determine the location which the primary particles were killed, as well as kill optical photons in volumes where the optical flag is off.

#### PrimaryGeneratorAction

   - The PrimaryGeneratorAction module uses G4ParticleGun to initialize the default positioning and configuration of the particle generator. By default the gun is set at the edge of the first module, centered in the (x,y) transverse plane.

   - Expansion of this class will provide a generator configuration for each beam environment offered by DetectorConstruction.

#### FiberSD

   - The detector response is handled in the FiberSD module. Each quartz rod is set as a sensitive detector and hits within these volumes are parsed in the Fiber Sensitive Detector class (FiberSD.cc) and FiberHits are created and filled in Collection vectors. This allows retention of detailed information about the hit including the particular fiber it was registered in, the module, the radiator gap number, the trackID, the particleID, as well as kinematic information.

   - This is also the stage where individual hits generate Cherenkov radiation and the number of Cherenkov photons created is counted.

#### VISUALISATION

   - The visualization manager is set via the G4VisExecutive class
   in the main() function in exampleB1.cc.    
   The initialisation of the drawing is done via a set of /vis/ commands
   in the macro vis.mac. This macro is automatically read from
   the main function when the program is used in interactive running mode.

   - By default, vis.mac opens an OpenGL viewer (```/vis/open OGL```).
   The user can change the initial viewer by commenting out this line
   and instead uncommenting one of the other ```/vis/open``` statements, such as
   HepRepFile or DAWNFILE (which produce files that can be viewed with the
   HepRApp and DAWN viewers, respectively).  Note that one can always
   open new viewers at any time from the command line.  For example, if
   you already have a view in, say, an OpenGL window with a name
   "viewer-0", then
   ```/vis/open DAWNFILE```
   then to get the same view
    ```/vis/viewer/copyView viewer-0```
   or to get the same view *plus* scene-modifications
      ```/vis/viewer/set/all viewer-0```
   then to see the result
      ```/vis/viewer/flush```

   - The DAWNFILE, HepRepFile drivers are always available
   (since they require no external libraries), but the OGL driver requires
   that the Geant4 libraries have been built with the OpenGL option.

   - From Release 9.6 the vis.mac macro has additional commands
   that demonstrate additional functionality of the vis system, such as
   displaying text, axes, scales, date, logo and shows how to change
   viewpoint and style.  Consider copying these to other examples or
   your application.  To see even more commands use help or
   ls or browse the available UI commands in the Application
   Developers Guide, Section 7.1.

   - For more information on visualization, including information on how to
   install and run DAWN, OpenGL and HepRApp, see the visualization tutorials,
   for example,
   http://geant4.slac.stanford.edu/Presentations/vis/G4[VIS]Tutorial/G4[VIS]Tutorial.html
   (where [VIS] can be replaced by DAWN, OpenGL and HepRApp)

   - The tracks are automatically drawn at the end of each event, accumulated
   for all events and erased at the beginning of the next run.

#### USER INTERFACES

   - The user command interface is set via the G4UIExecutive class
   in the main() function in zdc.cc. This script is the driver for your simulation


 C- HOW TO RUN

  - Execute zdc in the 'interactive mode' with visualization (must be run from the $JZCaPA/bin directory):
  ```
  $ ./zdc
  ```
  and type in the commands from run1.mac line by line:  
  ```
  Idle> /run/beamOn 1
  Idle> ...
  Idle> exit
  ```
  or
  ```
  Idle> /control/execute run1.mac
  ....
  Idle> exit
  ```

  - Execute zdc in the 'batch' mode from macro files (without visualization)
  ```
  $ ./zdc -m run1.mac -o Output/resultsFile.root -i Input/myGeneration_A.root
  ```
