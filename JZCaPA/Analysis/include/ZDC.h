/** @file ZDC
 *  @brief Function prototypes for ZDC
 *
 *  This contains the prototypes and members
 *  for ZDC
 *
 *  @author Chad Lantz, Riccardo Longo
 *  @bug No known bugs.
 */

#ifndef ZDC_H
#define ZDC_H

#include "Detector.h"

class ZDC : public Detector{

 public:
  ZDC( );
  ZDC(std::vector<Channel *> _readOut, int _runNumber, int _zdcNumber);
  virtual ~ZDC( );

  virtual void PrintMap  ( );
  void SetNumber ( int _number ) { m_Number = _number; }

 private:
  int m_Number;
};

#endif
