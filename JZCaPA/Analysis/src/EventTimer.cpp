/** @ingroup ana
 *  @file EventTimer.cpp
 *  @brief Implementation of EventTimer.
 *
 *  Function definitions for EventTimer are provided.
 *  This class is a reimplementation of TTimer with
 *  a new notify method for JZCaPA
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */

#include "EventTimer.h"

/** @brief Default constructor for EventTimer
 *
 */
EventTimer::EventTimer( ) : TTimer( ){

}

/** @brief Custom constructor for EventTimer
 *
 *  Makes a TTimer and sets the object member of EventTimer
 */
EventTimer::EventTimer(Long_t milliSec, DataReader* obj, Bool_t mode)
    : TTimer(milliSec, mode){
    m_object = obj;
    m_rate = milliSec;
}

/** @brief Default destructor for EventTimer
 *
 */
EventTimer::~EventTimer(){

}


/** @brief Notify when timer times out.
 *
 *  Calls the update method from the object pointer
 *  The timer is always reset. To stop the timer
 *  call TurnOff().
 */
Bool_t EventTimer::Notify(){
    Timeout();       // emit Timeout() signal
    if (m_object) m_object->UpdateConsole( m_rate );

    Reset();
    return kTRUE;
 }
