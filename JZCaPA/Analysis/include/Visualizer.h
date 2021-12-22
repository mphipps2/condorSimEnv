/** @file Visualizer.h
 *  @brief Class to have a simple visualization of plots through JZCaPA. If you want to plot something, just define here a child method and use it recursively, w/o filling the code with canvases et al.
 *
 *  @author Riccardo Longo
 *  @bug No known bugs.
 */

#ifndef VISUALIZER_H
#define VISUALIZER_H

#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TVector.h"
#include "TGraph.h"
#include "TVirtualPad.h"
#include "TMarker.h"

#include "Containers.h"

#include "TString.h"

#include <iostream>
#include <vector>

class Visualizer {

    public :
    Visualizer();
    Visualizer( std::string _style );

    //Styles (so far only ATLAS)
    TStyle* AtlasStyle();
    //Methods to set a given style
    void SetAtlasStyle();

    /** @brief Sets the label to be printed on each plot **/
    void SetTestBeamLabel(std::string _category, std::string _year, std::string _beam) { category = _category; year = _year; beam = _beam; }
    void SetTestBeamLabel( int runNo, Alignment* alignment );

    /** @brief allow the user to define the extension of the plots once they're printed. ".pdf" by default */
    void SetPlotExtension( std::string _extension ) { m_extension = _extension; }

    /** @brief Set the output directory of Visualizer */
    void SetOutputDirectory( std::string _oFolder ){ m_oFolder = _oFolder; }

    /** @brief make the argument true to activate the debug mode. False to deactivate it*/
    void   SetDebugMode   ( bool _isDebugon ) { m_debug = _isDebugon; }

    //Special plot treatments
    void   OverlayHistos  ( TH1 *h1, TH1 *h2 , TVirtualPad* pad, double _line = 0);
    void   ScatterPlot    ( std::vector< double > _v1, std::vector< double > _v2, TVirtualPad* pad);

    //Main visualization methods
    void   ManyPadsPlot   ( std::vector< TH1* > _first_form, std::vector< TH1* > _second_form, int _ncol, int _nrow, std::string _out_name, TString _treatment );
    void   SinglePlot     ( std::vector< double > _v1, std::vector< double > _v2, std::string _out_name, std::string _treatment, double _line = 0);

    void   DrawPlot (TH1* h2, std::string _xTitle, std::string _yTitle, std::string _saveName, std::string _opt, TMarker* _marker = 0, std::string _oFolder = "" );

    private :
    /** String identifying the style to be applied */
    std::string m_style;
    /** Strings encoding the year of the test and the beam **/
    std::string category = "NO";
    std::string year = "NO";
    std::string beam = "NO";
    /** Output folder for plots. Defaults to $JZCaPA/results/ */
    std::string m_oFolder = "${JZCaPA}/results/";
    /** Debug flag */
    bool m_debug = false;

    /** String defining the extension used to print the plots */
    std::string m_extension = ".pdf";
};

#endif
