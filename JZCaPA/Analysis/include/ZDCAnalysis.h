/** @file ZDCAnalysis
 *  @brief Function prototypes for ZDCAnalysis
 *
 *  This contains the prototypes and members 
 *  for ZDCAnalysis
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */

#ifndef ZDCANALYSIS_H
#define ZDCANALYSIS_H

#include "Analysis.h"
#include "TH2D.h"
#include "Containers.h"
#include "Detector.h"

class ZDCAnalysis : public Analysis{

 public :
  ZDCAnalysis( );
  virtual ~ZDCAnalysis( );
  
  virtual void   Initialize     ( );
  virtual void   Initialize     ( std::vector < Detector* > _vDet);
  virtual void   SetupHistograms( );
  virtual void   AnalyzeEvent   ( );
  virtual void   AnalyzeEvent   ( const std::vector< TH1* >& ){};
  virtual void   AnalyzeEvent   ( const std::vector< std::vector< float > >& ){};
  virtual void   AnalyzeEvent   ( const std::vector< Channel* > ){};
  virtual void   SetBranches    ( TTree* _tree );
  virtual void   Finalize       ( );


///// Energy histos

  /** Ratio between the Charge detected in the two ZDCs **/
  TH1D* hChargeRatio;
  /** Sum of ZDC charges **/
  TH1D* hChargeSum;
  /** Ratio between the peak heights of the two ZDCs **/
  TH1D* hPeakRatio;
  /** Charge correlation */
  TH2D *hCharge;
  /** ZDC1 Charge */
  TH1D *hCharge1;
  /** ZDC2 Charge */
  TH1D *hCharge2;
  /** Peak height correlation */
  TH2D *hPeak;
  /** ZDC1 peak height */
  TH1D *hPeak1;
  /** ZDC2 peak height */
  TH1D *hPeak2;
  /** Differential peak correlation */
  TH2D *hDpeak;
  /** ZDC1 peak height */
  TH1D *hDpeak1;
  /** ZDC2 peak height */
  TH1D *hDpeak2;
  /** Charge-Peak correlation, ZDC1 */
  TH2D* hChargePeakZDC1;
  /** Charge-Peak correlation, ZDC2 */
  TH2D* hChargePeakZDC2;
  
///// Timing histos

  /* ZDC1 peak center */
  TH1D* hArrival1;
  /* ZDC2 peak center */
  TH1D* hArrival2;
  /* Time of flight ZDC1 - ZDC2 */
  TH1D* hToF;
  

 private :
  /** Pointer to the first ZDC module */
  Detector *m_zdc1 = 0;
  /** Pointer to the second ZDC module */
  Detector *m_zdc2 = 0;
  /** Pointer to the first ZDC channel */
  Channel* zdc1 = 0;
  /** Pointer to the second ZDC channel */
  Channel* zdc2 = 0;
 
  
};

#endif
