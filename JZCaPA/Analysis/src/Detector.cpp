/** @ingroup ana
 *  @file Detector.cpp
 *  @brief Implementation of Detector.
 *
 *  Function definitions for Detector are provided.
 *  This is the mother class for detectors.
 *  Methods common to all detectors are implemented here.
 *
 *  @author Chad Lantz, Riccardo Longo
 *  @bug No known bugs.
 */

#include "Detector.h"
#include "Containers.h"
#include <vector>


/** @brief Default Constructor for Detector.
 */
Detector::Detector( ){

}

/** @brief Destructor for Detector.
 */
Detector::~Detector( ){

}

/** @brief Get the properties of a detector element
 *  @param row Row of element to be accessed
 *  @param column Column of element to be accessed
 *  @return Pointer to the Channel of the requested element
 *
 * Returns a pointer to the Channel stored in m_Elements with the
 * requested row and column.
 * If the requested element is not found, return a warning message and a NULL pointer.
 *
 */
Channel* Detector::GetElement(int row, int column){

    for(int i=0; i < (int)m_Element.size(); i++){
            if(row == m_Element[i]->mapping_row && column == m_Element[i]->mapping_column){
            return m_Element[i];
    }
  }
  std::cerr << " WARNING: Element (" << row << "," << column << ") not found! " << std::endl;
  return nullptr;

}

/** @brief Get the properties of a detector element
 *
 * Returns a pointer to the Channel stored in m_Elements with the
 * requested channel name.
 * If the requested element is not found, return a warning message and a NULL pointer.
 *
 */
Channel* Detector::GetElement(std::string _name){

    for(int i=0; i < (int)m_Element.size(); i++){
            if(!_name.compare(m_Element.at(i)->name)) return m_Element[i];
  }
  std::cerr << " WARNING: Element (" << _name << ") not found! " << std::endl;
  return nullptr;

}

/**
 * @brief Set the branches of the tree to the channels of the detectors (according to their name, read from the mapping)
 * @param _dataTree : processed data tree
 */
void Detector::SetBranches( TTree *_dataTree ){

    for( uint ch = 0; ch < m_Element.size(); ch++ ){
      _dataTree->SetBranchAddress( ("Raw" + m_Element.at(ch)->name).c_str(), &m_Element.at(ch)->pWF );
    }

}

/**
 * @brief Declare histograms to be filled with the raw waveform
 */
void Detector::DeclareHistograms(){

    for( uint ch = 0; ch < m_Element.size(); ch++ ){
        m_Element.at(ch)->WF_histo        = new TH1D( m_Element.at(ch)->name.c_str(), (m_Element.at(ch)->name + ", " + m_Element.at(ch)->detector).c_str(), m_nSamp, 0, m_nSamp);
        m_Element.at(ch)->PWF_histo       = new TH1D((m_Element.at(ch)->name + " Processed").c_str(),  (m_Element.at(ch)->name + ", " + m_Element.at(ch)->detector + " Processed").c_str() , m_nSamp, 0, m_nSamp);
        m_Element.at(ch)->FirstDerivative = new TH1D((m_Element.at(ch)->name + " Derivative").c_str(), (m_Element.at(ch)->name + ", " + m_Element.at(ch)->detector + " Derivative").c_str(), m_nSamp, 0, m_nSamp);
    }

}

/**
 * @brief Fill histograms with the current raw waveform
 */
void Detector::FillHistograms(){
    for( uint ch = 0; ch < m_Element.size(); ch++ ){
        m_Element.at(ch)->WF_histo->Reset();
        m_Element.at(ch)->PWF_histo->Reset();
        m_Element.at(ch)->FirstDerivative->Reset();
        m_Element.at(ch)->FirstDerivativeRMS = 0;
        m_Element.at(ch)->CrossZeroPoints.clear();
        // Loop over samples in each channel
        for( uint samp = 0; samp < m_nSamp; samp++ ){
          m_Element.at(ch)->WF_histo->SetBinContent( samp + 1, m_Element.at(ch)->WF[ samp ] );
          } // End loop over samples in each channel
        }
}
