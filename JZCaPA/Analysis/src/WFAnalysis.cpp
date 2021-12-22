/** @ingroup ana
 *  @file WFAnalysis.cpp
 *  @brief Implementation of WFAnalysis.
 *
 *  Function definitions for WFAnalysis are provided.
 *  This class is the main  class for the waveform analysis.
 *  It initializes the histograms used for output of processed events.
 *  Also includes methods that accept all waveforms in an event.
 *  This is where the analysis should be done.
 *
 *  @author Sheng Yang, Chad Lantz, Yakov Kulinich
 *  @bug No known bugs.
 */

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TRandom.h>
#include <TTree.h>
#include <TCanvas.h>//Sheng
#include <iostream>

#include "WFAnalysis.h"

/** @brief Default Constructor for WFAnalysis.
 */
WFAnalysis::WFAnalysis( ){

}

/** @brief Destructor for WFAnalysis.
 */
WFAnalysis::~WFAnalysis( ){

    delete f;
    delete nlFit;

}

/** @brief Initialization method for WFAnalysis
 *
 *  Can add other things here that you would
 *  perhaps not put into the constructor.
 *  I.e. a TTree, some tools. Etc.
 *
 */
void WFAnalysis::Initialize( ){

     nlFit = new TF1("nl1", "(x<660)*(-0.188891+1.03623*x) + \(x>660 && x<850)*( -51210.2 + 260.781*x - 0.4916*pow(x,2) + 0.00041227*pow(x,3) - (1.29681e-7)*pow(x,4) ) + \(x>850)*( -526857 + 1844.37*x - 2.1493*pow(x,2) + 0.000834971*pow(x,3) )", 0, 1000);

}


/** @brief Historgam Setup method for WFAnalysis
 *
 *  Should instantiate any histograms you wish to output here.
 *
 *  @return none
 */
void WFAnalysis::SetupHistograms( ){

  f = new TF1("f","gaus",-50,50);
  hRMS = new TH1D("RMS","RMS",75,-75,75);
  hRMSrpd = new TH1D("RMSrpd","RMSrpd",50,-30,30);
  hPed = new TH1D("ped","ped", 1024, -512, 512);

}


/** @brief Analyze Events method for WFAnalysis
 *  @param vWF 1D vector of all waveforms for N ch
 *
 *  Here is the event-based analysis code.
 *  A const 1D vector of N TH1* is received.
 *  N is the number of channels
 *
 */
void WFAnalysis::AnalyzeEvent( const std::vector< TH1* >& vWFH ){

  // example... you can loop through all histos as follows
  for( unsigned int ch = 0; ch < vWFH.size(); ch++ ){
    TH1* h = vWFH[ch];
    for( unsigned int samp = 0; samp < h->GetNbinsX(); samp++ ){
      // will print what each sample in each channel equals
      // std::cout << ch << " " << samp << " = " << h->GetBinContent( samp + 1 ) << std::endl;
    }
  }
}


/** @brief Analyze Events method for WFAnalysis
 *  @param vWF const 2D vector for all waveforms in event
 *
 *  Here is the event-based analysis code.
 *  A const 2D vector of NxM is received.
 *  N is the number of channels
 *  M is the number of samples per channel
 *
 *  @return none
 */
void WFAnalysis::AnalyzeEvent( const std::vector< std::vector< float > >& vWF ){

  // example... you can loop through each sample in each channel
  // Can do the same in previous AnalyzeEvent( .. ) function, just
  // looking at vWFH[ch]->GetBinContent( samp + 1 );
  for( unsigned int ch = 0; ch < vWF.size(); ch++ ){
    for( unsigned int samp = 0; samp < vWF[ch].size(); samp++ ){
      // will print what each sample in each channel equals
      // std::cout << ch << " " << samp << " = " << vWF[ch][samp] << std::endl;
    }
  }
}


/**
 * @brief Analyze Event method for WF analysis
 * @param vCh A vector of pointers to Channel objects
 *
 *  Receives a const vector of Channels. Loops over all Channels within
 *  the vector for analysis. Saves all outputs to members of each Channel
 *
 */
void WFAnalysis::AnalyzeEvent( const std::vector< Channel* > vCh ){

    for( unsigned int ch = 0; ch < vCh.size(); ch++ ){
        if( !vCh.at(ch)->is_on){continue;}
      //retrieving information for each channel as a histogram
        TH1D* h = vCh.at(ch)->WF_histo;
        TH1D* hProcessed = vCh.at(ch)->PWF_histo;
        TH1D* hDiff = vCh.at(ch)->FirstDerivative;

      //If Channel belongs to an RPD, set RPD processing values
        if( vCh.at(ch)->detector.find_first_of("R") != std::string::npos ){
            m_diffSens  = m_RPDdiffSens;
            m_Tmultiple = m_RPDTmultiple;
            fCutoff     = m_RPDfCutoff;
            m_isRPD     = true;
        }else{ // If detector = ZDC, set ZDC processing values and invert the signal
            m_diffSens  = m_ZDCdiffSens;
            m_Tmultiple = m_ZDCTmultiple;
            fCutoff     = m_ZDCfCutoff;
            m_isRPD     = false;
            //The channel has an offset we set at the time of the Testbeam. Invert that,
            //then set the WF_histo with the inverted values of the raw WF vector. Finally,
            //invert the raw WF vector.
            vCh.at(ch)->offset = -1*vCh.at(ch)->offset;
            for(int bin = 0; bin < h->GetNbinsX(); bin++){
                h->SetBinContent( bin, -1*vCh.at(ch)->pWF->at(bin) );
                vCh.at(ch)->pWF->at(bin) *= -1;
            }
        }

        //Get the first derivative and determine the hit window from it
        //if a valid hit window is not found was_hit is set to false
        GetDifferential( vCh.at(ch) );

        //We set the offset of the DRS4 channels to -250 at first. This lead to clipping
        //of the baseline. If that's the case, estimate the RMS to be 5 based on good runs.
        vCh.at(ch)->FirstDerivativeRMS = ( vCh.at(ch)->offset != -250 ) ? GetRMS( vCh.at(ch) ) : 5.0;
        FindHitWindow( vCh.at(ch) );

        //If the channel was hit, proceed with processing
        if( vCh.at(ch)->was_hit ){
            //Get and subtract the pedestal from the data outside the hit window
            GetPedestal( vCh.at(ch) );
            //The DRS4 saturates around 800-900 mV. If we have those values, flag it as saturated
            if( vCh.at(ch)->PWF_histo->GetMaximum() > 890.0 ){ vCh.at(ch)->saturated = true; }

            //Calibrate out the DRS4 non-linearity and apply FFT low pass filter if desired
            if( m_DRS4){ DRS4Cal( vCh.at(ch) ); }
            if( m_FFT ){ LowPassFilter( vCh.at(ch), hProcessed ); }

            //Zero Suppress the processed waveform vector and retrieve the energy related values from it
            ZeroSuppress( vCh.at(ch) );
            GetCharge( vCh.at(ch) );
            vCh.at(ch)->Diff_Peak_time = vCh.at(ch)->pTimeVec->at( vCh.at(ch)->Diff_Peak_center );

            //If the algorithm didn't work, go with a simpler method, just to get values that make sense.
            //This should be changed to work better when I can think straight.
            if( vCh.at(ch)->Peak_center == -1 ){
                vCh.at(ch)->Peak_max  = vCh.at(ch)->PWF_histo->GetMaximum();
                vCh.at(ch)->Peak_time = vCh.at(ch)->pTimeVec->at( vCh.at(ch)->PWF_histo->GetMaximumBin() );

            }else{
                vCh.at(ch)->Peak_max  = hProcessed->GetBinContent( vCh.at(ch)->Peak_center );
                vCh.at(ch)->Peak_time =  vCh.at(ch)->pTimeVec->at( vCh.at(ch)->Peak_center );

            }
        }
    }
}

/** @brief GetDifferential method for WFAnalysis
 *  @param Ch channel who's waveform will be differentiated
 *
 * The derivative is calculated as the difference of the summation N points before and after the ith point
 *  @f[
 *       \delta_i(N)= \sum_{k=1}^{N} (s_i + k)  -  \sum_{k=1}^{N} (s_i - k)
 *  @f]
 *
 * N is given by m_diffSens and is set by the constructor or SetDiffSense()
 */
void WFAnalysis::GetDifferential( Channel* Ch ){

    // Loop over histogram
    for( unsigned int bin = m_diffSens; bin < Ch->WF_histo->GetNbinsX() - m_diffSens ; bin++ ){
      //decalare temp variable to store the sum before and after the i data point
      double sum_before = 0;
      double sum_after = 0;
      //calculate the sum before and after the i data point
      for (int i = 0; i < m_diffSens; i++  ){
        sum_before += Ch->WF.at( bin - i );
        sum_after  += Ch->WF.at( bin + i );
      }
        //set the bin to the calculated derivative value
        Ch->FirstDerivative->SetBinContent(bin,(sum_after - sum_before));

    }//end derivative loop
}

/** @brief GetRMS method for WFAnalysis
 *  @param Input Channel
 *  @return Width of gaussian fit
 *
 *  Histgrams the second derivative for a given Channel and gets the baseline RMS
 *  from a gaussian fit which ignores tails created by the peaks.
 *
 */
double WFAnalysis::GetRMS( Channel* Ch ){

    //Make a histogram with x-range to match the y-range of the input
    TH1D* hUsed = (m_isRPD) ? hRMSrpd : hRMS;
    hUsed->Reset();
    Double_t xmin,xmax;
    Ch->FirstDerivative->GetMinimumAndMaximum(xmin,xmax);
    if(xmax == 0) return 0;

    //Loop over the histogram excluding the window used for differentiating to fill hRMS
    Int_t nbins = Ch->FirstDerivative->GetNbinsX();
    for(int bin = m_diffSens; bin < nbins - m_diffSens; bin++){
        hUsed->Fill( Ch->FirstDerivative->GetBinContent( bin ) );
    }

    //Set the fit range narrower for the RPD than the ZDC. RPD data is actually smoother
    //Set all fit parameters with random numbers before fitting.
    m_isRPD ? f->SetRange(-15,15) : f->SetRange(-50,50);
    f->SetParameters( gRandom->Uniform(0,1000) , gRandom->Uniform(-300,300), gRandom->Uniform(0,100));
    hUsed->Fit("f","qR");

    //Return parameter 2. "gaus" is [0]*exp(-0.5*((x-[1])/[2])**2)
    //And reset the histogram
    return f->GetParameter(2);

}


/** @brief Defines the hit window for a given channel
 *  @param ch Channel to be processed
 *
 *  Determines if the channel was hit first, then determines the hit window using
 *  the first derivative and first derivative RMS. Also determines the peak height
 *  and peak center using the raw waveform value at the first derivative zero crossing.
 *  Saves the results to Channel members.
 */
void WFAnalysis::FindHitWindow( Channel* ch ){

    double threshold = m_Tmultiple*ch->FirstDerivativeRMS;
    ch->Diff_max = ch->FirstDerivative->GetMaximum();
    int risingEdge = ch->Diff_Peak_center = ch->FirstDerivative->GetMaximumBin();
    int fallingEdge = ch->FirstDerivative->GetMinimumBin();

    //  If the derivative maximum or minimum are below threshold, no hit
    //  Also, if the rising edge is after the falling edge (i.e. a negative pulse), no hit
    if( ch->Diff_max <= threshold || ch->FirstDerivative->GetMinimum() >= -1*threshold || risingEdge > fallingEdge){
        if( ch->is_on && m_verbose > 1){ std::cerr << std::endl << "No hit found on " << ch->name << std::endl; }
        ch->was_hit = false;

        //If the channel didn't register a hit, set all values to -1 if int and 0.0 if double
        ch->Peak_center = ch->Diff_Peak_center = ch->hit_window.first = ch->hit_window.second = -1;
        ch->Peak_max = ch->Charge = ch->Diff_max = ch->Peak_time = ch->Diff_Peak_time = 0.0;
        return;
    }

    //If it made it here, we probably have a hit and we can start processing
    ch->was_hit = true;
    int nBins = ch->FirstDerivative->GetNbinsX();

    // Find the beginning of the hit window
    for(int bin = risingEdge; bin > 0; bin--){
        if(ch->FirstDerivative->GetBinContent(bin) < threshold){
            ch->hit_window.first = bin;
            break;
        }
    }

    // Find the peak center using the derivative
    for(int bin = risingEdge; bin < fallingEdge; bin++){
        if(ch->FirstDerivative->GetBinContent(bin) < 0){
            ch->Peak_center = bin;
            break;
        }
    }

    // Find the end of the hit window
    for(int bin = fallingEdge; bin < nBins; bin++){
        if(ch->FirstDerivative->GetBinContent(bin) > -1*threshold){
            ch->hit_window.second = bin;
            break;
        }
    }
}


/** @brief Finds the pedestal and pedestal rms of the raw waveform
 *  @param ch Channel to be processed
 *
 *  Histograms data from the raw waveform excluding the hit window.
 *  Saves the mean and rms of that histogram as PedMean and PedRMS.
 *  Also saves a pedestal subtracted, zero supressed version of WF_histo to PWF_histo.
 *  Allows for a reference histogram to be supplied i.e. PWF_histo. Otherwise, defaults to WF_histo.
 */
void WFAnalysis::GetPedestal( Channel* ch ){
    int nBins = ch->WF_histo->GetNbinsX();
    //TH1D h("ped","ped", nBins, -nBins/2, nBins/2);

    if( ch->was_hit && ch->hit_window.first > ch->hit_window.second ){
        if(ch->is_on && m_verbose > 1){
            std::cerr << std::endl << "Bad hit window on " << ch->name << ": Cannot get pedestal" << std::endl;
        }
        return;
    }
    // Find PedMean and PedRMS
    for(int bin = 0; bin < nBins; bin++){
        // Skip the hit window
        if( bin == ch->hit_window.first ){ bin = ch->hit_window.second; }

        hPed->Fill( ch->pWF->at(bin) );
    }

    ch->PedMean = hPed->GetMean();
    ch->PedRMS  = hPed->GetRMS();

    // Subtract PedMean from PWF_histo
    for(int bin = 0; bin < nBins; bin++){
        double content = ch->pWF->at(bin) - ch->PedMean;
        ch->PWF_histo->SetBinContent( bin, content );
    }
}

/** @brief Supress values below PedRMS
 *
 *  Changes all values smaller than PedRMS after
 *  pedestal subtraction to zero.
 */
void WFAnalysis::ZeroSuppress( Channel* ch ){
    int nBins = ch->PWF_histo->GetNbinsX();

    for(int bin = 0; bin < nBins; bin++){
        double content = ch->PWF_histo->GetBinContent(bin);
        if( bin > ch->hit_window.first && bin < ch->hit_window.second  ){ continue; }
        if( content <= ch->PedRMS ){ ch->PWF_histo->SetBinContent( bin, 0 ); }
    }
}

/** @brief Applies a low pass filter to the processed waveform
 *  @param ch Channel to be processed
 *
 *  Applies FFT to the raw waveform, removes high frequencies
 *  and reconstructs the signal in PWF_histo.
 */
void WFAnalysis::LowPassFilter( Channel* ch, TH1D* hIn ){
   if(!hIn){ hIn = ch->WF_histo; }
   Int_t n = hIn->GetNbinsX();
   double re[n], im[n], reOut[n], imOut[n];

   TH1 *hm =0;
   hm = hIn->FFT(hm, "MAG");

   // Apply the filter
   TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
   fft->GetPointsComplex(re,im);
   for(int i = 0; i< fCutoff; i++){
       reOut[i] = re[i];
       imOut[i] = im[i];
   }
   for(int i = fCutoff; i < n; i++){
       reOut[i] = 0.0;
       imOut[i] = 0.0;
   }

   //Make a new TVirtualFFT with an inverse transform. Set the filtered points and apply the transform
   TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
   fft_back->SetPointsComplex(reOut,imOut);
   fft_back->Transform();
   ch->PWF_histo = (TH1D*)TH1::TransformHisto(fft_back,ch->PWF_histo,"MAG");

   //The transform scales the signal by sqrt(n) each time, so rescale by 1/n
   ch->PWF_histo->Scale(1.0/n);

   delete hm;
   delete fft;
   delete fft_back;

}


/** @brief Provides calibration for DRS4 voltage response
 *  @param ch Channel to be processed
 *
 *  Uses a TF1 with the DRS4 response to convert from the recorded value in mV
 *  to the actual value of the waveform in mV and reconstructs the signal in PWF_histo.
 *
 */
void WFAnalysis::DRS4Cal( Channel* ch ){

    for(int bin = 0; bin < ch->PWF_histo->GetNbinsX(); bin++){
        ch->PWF_histo->SetBinContent(bin, nlFit->Eval( ch->PWF_histo->GetBinContent(bin) ) );
    }
}

/** @brief Uses time vector and input resistance to find the actual charge detected by the DRS4
 *  @param ch Channel to be processed
 *
 *
 */
void WFAnalysis::GetCharge( Channel* ch ){
    ch->Charge = 0.0;

    for(int bin = ch->hit_window.first; bin < ch->hit_window.second; bin++){
        //Charge for a given time bin = dt*I = (t(bin+1)-t(bin))*V/R
        //Units are Charge(pC), time(ns), Voltage(mV), resistance(ohm), current(Amp = C/s)
        ch->Charge += (ch->pTimeVec->at(bin+1) - ch->pTimeVec->at(bin)) * ch->PWF_histo->GetBinContent(bin)/Rin;
    }
}

/** @brief Finalize method for WFAnalysis
 *  @return none
 *
 *  Write histograms, TTree if it exists.
 *
 */
void WFAnalysis::Finalize( ){

  // If these exist...
  // m_tree->Write();
  // m_hist->Write();
}
