#!/bin/bash
#hostname;

source /software/root-6.16.00/bin/thisroot.sh
source myenv.sh
#generate file
echo "root -l -q -b 'V1Gen.C'"
root -l -q -b 'V1Gen.C'
#create a new run file
ls
cp $2 A_$7_$2

echo $(root -l -q -b  'get_entry.C(1)' |tail -n 1)
#set geant run macro to independently simulate n neutrons
echo "/run/beamOn $(root -l -q -b  'get_entry.C(1)' |tail -n 1)" >> A_$7_$2
echo "./zdc $1 A_$7_$2 $3 A_$4 $5 $6 -i A_$8"
./zdc $1 A_$7_$2 $3 A_$4 $5 $6 "-i" A_$8
ls
echo "root -l -q -b 'CondenseEvent.C("'A_$4'")'"
root -l -q -b 'CondenseEvent.C("'A_$4'")'

# run digitizer
echo "root -l -q -b 'outputDigitization.cc("'A_$4'")'"
root -l -q -b 'outputDigitization.cc("'A_$4'")'

mkdir OutputA
# run jzcapa over the waveform
echo "./rpdMLTrainingAnalysis 1 WF_A_$4"
./rpdMLTrainingAnalysis 1 WF_A_$4 OutputA/
mv A_myGeneration.root OutputA/event$7A_gen.root
mv A_output_$7.root OutputA/event$7A_sim.root
mv OutputA/output1.root OutputA/event$7A_jzcapa.root
hadd OutputA/event$7A.root OutputA/event$7A_gen.root OutputA/event$7A_sim.root OutputA/event$7A_jzcapa.root

#create a new run file
cp $2 B_$7_$2
#run geant simulation side C w/ nEvents equal to number of neutrons generated in V1Gen.C
echo "/run/beamOn $(root -l -q -b  'get_entry.C(0)' |tail -n 1)" >> B_$7_$2
echo "./zdc $1 B_$7_$2 $3 B_$4 $5 $6 -i B_$8"
./zdc $1 B_$7_$2 $3 B_$4 $5 $6 "-i" B_$8
ls
# geant events come with each neutron as separate entry. This script makes it a true N neutron single event
echo "root -l -q -b 'CondenseEvent.C("'B_$4'")'"
root -l -q -b 'CondenseEvent.C("'B_$4'")'

# run digitizer
echo "root -l -q -b 'outputDigitization.cc("'B_$4'")'"
root -l -q -b 'outputDigitization.cc("'B_$4'")'

# run jzcapa over the waveform
mkdir OutputB
echo "./rpdMLTrainingAnalysis 1 WF_B_$4"
./rpdMLTrainingAnalysis 1 WF_B_$4 OutputB/
mv B_myGeneration.root OutputB/event$7B_gen.root
mv B_output_$7.root OutputB/event$7B_sim.root
mv OutputB/output1.root OutputB/event$7B_jzcapa.root
hadd OutputB/event$7B.root OutputB/event$7B_gen.root OutputB/event$7B_sim.root OutputB/event$7B_jzcapa.root
