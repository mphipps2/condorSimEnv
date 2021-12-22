/** @file RPD
 *  @brief Function prototypes for RPD
 *
 *  This contains the prototypes and members
 *  for RPD
 *
 *  @author Chad Lantz, Riccardo Longo
 *  @bug No known bugs.
 */

#ifndef RPD_H
#define RPD_H

#include "Containers.h"
#include "Detector.h"

class RPD : public Detector{

 public:
  RPD( );
  RPD( std::vector< Channel* > _readOut, int _runNumber );
  ~RPD( );

  Channel* GetElement(int row, int column);

  void SetnRows(int rows){nRows=rows; nElements = nRows * nColumns; ResizeSortedElements(); };
  void SetnCols(int cols){nColumns=cols; nElements = nRows * nColumns; ResizeSortedElements(); };
  void ResizeSortedElements();

  virtual void PrintMap( );

 private:
  /** Number of rows **/
  int nRows = 4;
  /** Number of columns **/
  int nColumns = 4;
  /** Translation of row elements for 2018 testbeam
      Converts the tile mapping (row,col)
            (1,4) (1,3) (1,2) (1,1)          (0,0) (0,1) (0,2) (0,3)
      from  (2,4) (2,3) (2,2) (2,1)   to     (1,0) (1,1) (1,2) (1,3)
            (3,4) (3,3) (3,2) (3,1)          (2,0) (2,1) (2,2) (2,3)
            (4,4) (4,3) (4,2) (4,1)          (3,0) (3,1) (3,2) (3,3)

      As viewed from the beam souce
  **/
  int rowTranslation [4];
  /** Translation of column elements for 2018 testbeam
      Converts the tile mapping (row,col)
            (1,4) (1,3) (1,2) (1,1)          (0,0) (0,1) (0,2) (0,3)
      from  (2,4) (2,3) (2,2) (2,1)   to     (1,0) (1,1) (1,2) (1,3)
            (3,4) (3,3) (3,2) (3,1)          (2,0) (2,1) (2,2) (2,3)
            (4,4) (4,3) (4,2) (4,1)          (3,0) (3,1) (3,2) (3,3)

  As viewed from the beam souce**/
  int colTranslation [4];
  /** Total elements **/
  int nElements = nRows * nColumns;
  /** 2D vector of channels sorted in a [row][column] format **/
  std::vector< std::vector< Channel* > > m_SortedElements;

};

#endif
