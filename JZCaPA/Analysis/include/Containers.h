/** @file Containers.h
 *  @brief Class to define containers: child classes with only public data memebers
 *
 *
 *  @author Riccardo Longo
 *  @bug No known bugs.
 */

#ifndef CONTAINERS_H
#define CONTAINERS_H

#include "Containers.h"
#include "TH1.h"

#include <iostream>
#include <vector>
#include <utility>

class Channel {

 public :
    /** Detector associated to this channel. RPD or ZDC
    *   Mapping (sitting on the impinging beam) [R,C] for ZDC:
    *   ZDC1 -> [0,1]    |
    *   ZDC2 -> [0,2]   \ /
    *   Mapping (sitting on the impinging beam) [R,C] for RPD:
    *   Jura [1,4] [1,3] [1,2] [1,1] Saleve
    *        [2,4] [2,3] [2,2] [2,1]
    *        [3,4] [3,3] [3,2] [3,1]
    *        [4,4] [4,3] [4,2] [4,1]
    */
	/** Type of detector - ZDC or RPD **/
    std::string detector;
    /** Channel name - CXX with XX in [1,20]  */
    std::string name;
    /** High voltage set for this channel - in V */
    double HV;
    /** Delay set for this channel - in ns */
    double delay;
    /** Offset set for this channel */
    double offset;
    /** Mapping in the horizontal direction [R] */
    int mapping_row;
    /** Mapping in the vertical direction [C] */
    int mapping_column;
    /** Was the channel functioning */
    bool is_on;
    /** Operating voltage of the channel*/
    int Vop;
    /** Raw waveform for a particular event */
    std::vector < float > WF;
    /** Pointer to the WF vector */
    std::vector < float > *pWF = &WF;
    /** Pointer to the DRS4 time vector */
    std::vector < float > *pTimeVec = 0;
    /** Histrogram for visualization and analysis of the waveform */
    TH1D* WF_histo;
    /** Histogram of the processed waveform */
    TH1D* PWF_histo;
    /** Histogram with first derivative **/
    TH1D* FirstDerivative;
    /** RMS value of the first derrivative of the waveform **/
    double FirstDerivativeRMS;
    /** Bin number of the peak center*/
    int Peak_center;
    /** Bin number of the derivative peak */
    int Diff_Peak_center;
    /** Calibrated time of the peak center in ns*/
    double Peak_time;
    /** Calibrated time of max slope in ns*/
    double Diff_Peak_time;
    /** Height of the peak */
    double Peak_max;
    /** Max value of the derivative */
    double Diff_max;
    /** Bin number of 1/3 peak value on the rising edge*/
    int DiscriminatorTime;
    /** Start and end of the hit window*/
    std::pair< int, int > hit_window;
    /** Flag if the channel was hit */
    bool was_hit;
    /** Flag if the waveform saturated */
    bool saturated;
    /** Pedestal of raw waveform */
    double PedMean;
    /** Pedestal RMS of raw waveform */
    double PedRMS;
    /** Integral of PWF_histo in the hit window in pC */
    double Charge;
    /** Crossing zero points - dummy - to be checked by Sheng **/
    std::vector < int > CrossZeroPoints;

};

class Alignment {

 public:

    /** Run number being analyzed **/
    int runNumber;
    /** X position of the Desy Table **/
    double x_table;
    /** Y position of the Desy Table **/
    double y_table;
    /** First detector met by the beam **/
    std::string upstream_Det;
    /** Second detector met by the beam **/
    std::string mid_Det;
    /** Third detector met by the beam **/
    std::string downstream_Det;
    /** GOLIATH magnet status **/
    bool magnet_On;
    /** Target in **/
    bool target_In;
    /** Lead absorber in **/
    bool lead_In;

};


#endif
