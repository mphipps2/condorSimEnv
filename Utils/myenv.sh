#!/bin/bash 
export PATH=/software/root-6.16.00/bin:/software/geant-4.10.5/bin:/software/geant-4.10.5/bin:/software/cmake-3.14.1/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/${USER}/.local/bin:/home/${USER}/bin
export JZCaPA=${PWD}
export ROOTSYS=/software/root-6.16.00
export LD_LIBRARY_PATH=/software/root-6.16.00/lib:/software/geant-4.10.5/lib64:/software/geant-4.10.5/lib64 
export GEANT4_DIR=/software/geant-4.10.5/  
export GEANT4_SHARE=${GEANT4_DIR}/share/Geant4-10.5.0/data/ 
export G4LEDATA=${GEANT4_SHARE}/G4EMLOW7.7
export G4SAIDXSDATA=${GEANT4_SHARE}/G4SAIDDATA2.0
export G4NEUTRONXSDATA=${GEANT4_SHARE}/G4PARTICLEXS1.1/neutron
export G4PARTICLEXSDATA=${GEANT4_SHARE}/G4PARTICLEXS1.1
export G4LEVELGAMMADATA=${GEANT4_SHARE}/PhotonEvaporation5.3
export G4ENSDFSTATEDATA=${GEANT4_SHARE}/G4ENSDFSTATE2.2
