/** @file AnalysisExample.cpp
 *  @brief Example of a simple analysis: read and plot all the w
 *
 *
 *  @author Yakov Kulinich, Riccardo Longo
 *  @bug No known bugs.
 */

#include "DataReader.h"
#include "WFAnalysis.h"
#include "ZDCAnalysis.h"
#include "RPDAnalysis.h"
#include "EventTimer.h"

using namespace std;


int main(int argc, char *argv[]){

  // can put arguments here, for now I will just do defaults.

  int nCh    = 20;   // 5 DRS4 x 4 ch/board - 16 RPD channels
  int nSamp  = 1024; // Default number of samples?
  int runNum = 99;   // !! Change for your test !!

  string fNameIn = Form("ZDCBeamTestRun%d.root", runNum); // !! Change for your test !!

  // DataReader is the main class. It reads data and also
  // has analysis classes in it. User should only have to
  // modify the analysis classes and add output in them.
  // User has to add their analysis to DataReader.
  DataReader* r = new DataReader( nCh, nSamp, fNameIn, runNum );
  r->SetVerbosity(1);

  r->SelectDetectorForAnalysis(true,true,false);
  r->AddPreAnalysis( new WFAnalysis() );
  r->AddDetectorAnalysis( new ZDCAnalysis() );
  //r->AddDetectorAnalysis( new RPDAnalysis() );
  r->LoadConfigurationFile();
  r->LoadAlignmentFile();
  r->LoadTimingFile();
  r->EnablePlotLabel();
  
  EventTimer timer(1000, r, kFALSE);
  timer.TurnOn();
  
  r->Run();
  
  timer.TurnOff();
  std::cout << std::endl << "Finished!" << std::endl;
  
  delete r;
  
  return 0;
}
