/** @file Detector
 *  @brief Function prototypes for Detector
 *
 *  This contains the prototypes and members 
 *  for Detector
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */

#ifndef DETECTOR_H
#define DETECTOR_H

#include <vector>
#include <string>

#include "Containers.h"
#include "TTree.h"

class Detector{

 public:
  Detector( );
  Detector( std::vector< Channel* > _element ){ m_Element = _element;}
  virtual ~Detector( );

  
  virtual Channel* GetElement  (int row, int column);
  virtual Channel* GetElement  (std::string _name);
  virtual std::vector < Channel* > GetChannelsVector () { return m_Element; }
  virtual double*    GetPosition ( ) { return m_Position; }
  virtual double*    GetAngle    ( ) { return m_Angle; }
  virtual Alignment* GetAlignment( ) { return m_Alignment; }

  virtual void     SetNSamples ( int _nSamples )  { m_nSamp = _nSamples; }
  virtual void     SetElement  ( Channel* _entry) { m_Element.push_back(_entry); }
  virtual void     SetPosition (double x, double y, double z) { m_Position[0] = x; m_Position[1] = y; m_Position[2] = z; }
  virtual void     SetAngle    (double _cosx = 0, double _cosy = 0, double _cosz = 0) { m_Angle[0] = _cosx; m_Angle[1] = _cosy; m_Angle[2] = _cosz; }
  virtual void     SetBranches ( TTree* _dataTree );
  virtual void     SetAlignment( Alignment* _alignment ){ m_Alignment = _alignment; }
  virtual void     DeclareHistograms ( );
  virtual void     FillHistograms ( );

  virtual void     PrintMap    ( ) = 0;
  
  private:
  /** Vector of channels associated to the dector **/
  std::vector< Channel* > m_Element;
  /** Three element array with x, y, and z of some pre-defined point on the detector **/
  double m_Position[3];
  /** Three element array of angle about the x, y, and z axis **/
  double m_Angle[3];
  /** Number of samples per channel **/
  int m_nSamp = 1024;
  /** Alignment of the Testbeam */
  Alignment* m_Alignment = 0;

  };

#endif
