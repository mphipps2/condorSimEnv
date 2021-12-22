/** @file EventTimer
 *  @brief Reimplementation of TTimer
 *
 *  This contains the prototypes and members 
 *  for EventTimer
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */
 
#ifndef EVENTTIMER_H
#define EVENTTIMER_H

#include "TTimer.h"
#include "DataReader.h"
 
class EventTimer : public TTimer {
  public:
    EventTimer( );
    EventTimer(Long_t milliSec = 0, DataReader* obj = 0, Bool_t mode = kTRUE);
    ~EventTimer();
    
    virtual Bool_t Notify();

    DataReader* m_object;
    Long_t m_rate;
};

#endif