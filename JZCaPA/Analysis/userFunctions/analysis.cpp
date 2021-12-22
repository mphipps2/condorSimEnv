/** @file AnalysisExample.cpp
 *  @brief Example of a simple analysis: read and plot all the w
 *
 *
 *  @author Yakov Kulinich, Riccardo Longo
 *  @bug No known bugs.
 */

#include "DataReader.h"
#include "WFAnalysis.h"

using namespace std;

int main(int argc, char *argv[]){

  // can put arguments here, for now I will just do defaults.

  int nCh    = 20;   // 5 DRS4 x 4 ch/board - 16 RPD channels
  int nSamp  = 1024; // Default number of samples?
  int runNum = 54;   // !! Change for your test !!

  string fNameIn = "TreeZDCBeamTestRun54.root"; // !! Change for your test !!

  // DataReader is the main class. It reads data and also
  // has analysis classes in it. User should only have to
  // modify the analysis classes and add output in them.
  // User has to add their analysis to DataReader.
  DataReader* r = new DataReader( nCh, nSamp, fNameIn, runNum );

  r->AddAnalysis( new WFAnalysis() );

  r->Run();
  
  delete r;
  
  return 0;
}
