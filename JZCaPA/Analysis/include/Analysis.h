/** @file Analysis
 *  @brief Function prototypes for Analysis
 *
 *  This contains the prototypes and members 
 *  for Analysis
 *
 *  @author Yakov Kulinich
 *  @bug No known bugs.
 */

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <string>
#include <vector>
#include "TVirtualPad.h"

#include "Detector.h"
#include "Containers.h"
#include "Visualizer.h"

class TFile;
class TH1;
class TH2;
class TH3;
class TTree;

class Analysis{

 public:
  Analysis( ){};
  virtual ~Analysis( ){};

  virtual void Initialize     ( ) = 0;
  virtual void Initialize     ( std::vector < Detector* > ) = 0;
  virtual void SetupHistograms( ) = 0;
  virtual void AnalyzeEvent   ( ) = 0;
  virtual void AnalyzeEvent   ( const std::vector< TH1* >& ) = 0;
  virtual void AnalyzeEvent   ( const std::vector< std::vector< float > >& ) = 0;
  virtual void AnalyzeEvent   ( const std::vector< Channel* > ) = 0;
  virtual void SetVerbosity   ( int _level ){ m_verbose = _level; }
  virtual void SetBranches    ( TTree* _tree ){ m_AnalysisTree = _tree; }
  virtual void Finalize       ( ) = 0;
  virtual void AssignVisualizer ( Visualizer *_viz ){ m_viz = _viz; }
  
 protected:
  /** Verbosity level */
  int m_verbose = 0;
  /** Output tree */
  TTree *m_AnalysisTree = 0;
  /** Visualizer for plots **/
  Visualizer* m_viz = 0;
  /** Alignment information for the given run */
  Alignment* m_alignment;
};

#endif
