/** @file WFAnalysis
 *  @brief Function prototypes for WFAnalysis
 *
 *  This contains the prototypes and members
 *  for WFAnalysis
 *
 *  @author Yakov Kulinich
 *  @bug No known bugs.
 */

#ifndef WFANALYSIS_H
#define WFANALYSIS_H

#include "Analysis.h"
#include "Containers.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TLine.h"


class WFAnalysis : public Analysis{

 public :
  WFAnalysis( );
  WFAnalysis( int _sensitivity, double _threshold ){ m_diffSens = _sensitivity; m_Tmultiple = _threshold; }
  WFAnalysis( bool _DRS4, bool _FFT = false ){ m_DRS4 = _DRS4; m_FFT = _FFT; }
  virtual ~WFAnalysis( );

  virtual void   Initialize     ( );
  virtual void   Initialize     ( std::vector < Detector* > ){};
  virtual void   SetupHistograms( );
  virtual void   GetDifferential( Channel* Ch );
  virtual double GetRMS         ( Channel* Ch ) ;
  virtual void   GetPedestal    ( Channel* ch );
  virtual void   GetCharge      ( Channel* ch );
  virtual void   FindHitWindow  ( Channel* ch );
  virtual void   ZeroSuppress   ( Channel* ch );
  virtual void   LowPassFilter  ( Channel* ch, TH1D* hIn = 0 );
  virtual void   DRS4Cal        ( Channel* ch );
  virtual void   AnalyzeEvent   ( ){};
  virtual void   AnalyzeEvent   ( const std::vector< TH1* >& );
  virtual void   AnalyzeEvent   ( const std::vector< std::vector< float > >& );
  virtual void   AnalyzeEvent   ( const std::vector< Channel* > vCh );
  virtual void   Finalize       ( );

  /** Histogram used to get derivative RMS */
  TH1D *hRMS = 0;
  TH1D *hRMSrpd = 0;
  bool m_isRPD;
  bool m_DRS4;
  bool m_FFT;

 private :
  /** Sensitivity level to hits (differentiation window) */
  int m_diffSens;
  /** Hit threshold multiplier */
  double m_Tmultiple;
  /** Frequency threshold for low pass filter */
  int fCutoff;
  /** Sensitivity level to hits (differentiation window) for RPD channels*/
  int m_RPDdiffSens = 3;
  /** Hit threshold multiplier for RPD channels*/
  double m_RPDTmultiple = 3.5;
  /** Frequency threshold for low pass filter for RPD channels*/
  int m_RPDfCutoff = 50;
  /** Sensitivity level to hits (differentiation window) for ZDC channels*/
  int m_ZDCdiffSens = 7;
  /** Hit threshold multiplier for ZDC channels*/
  double m_ZDCTmultiple = 3.5;
  /** Frequency threshold for low pass filter for ZDC channels*/
  int m_ZDCfCutoff = 50;

  /** Histogram used to find pedestal */
  TH1D *hPed = 0;
  /** Fit used to get RMS value */
  TF1 *f = 0;
  /** TF1 used to correct for DRS4 non-linearity */
  TF1 *nlFit = 0;
  /** Input resistance of the DRS4 in ohms */
  double Rin = 50.0;
};

#endif
