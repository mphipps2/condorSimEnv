/** @ingroup ana
 *  @file Visualizer.cpp
 *  @brief Implementation of Visualizer class.
 *
 *
 *  @author Riccardo Longo
 *  @bug No known bugs.
 */

#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TClass.h>

#include <iostream>

#include "Visualizer.h"
#include "Containers.h"

/** @brief Visualizer Constructor for Visualizer.
 */
Visualizer::Visualizer( ){

}

Visualizer::Visualizer( std::string _style ){
    m_style = _style;
    if(m_style == "ATLAS" || m_style == "atlas") Visualizer::SetAtlasStyle();
}

/**
 * @brief Function inherited from rcdaq. Creates the ATLAS TStyle object
 * @return
 */
TStyle* Visualizer::AtlasStyle(){
    TStyle *atlasStyle = new TStyle("ATLAS","Atlas style");

      // use plain black on white colors
      Int_t icol=0; // WHITE
      atlasStyle->SetFrameBorderMode(icol);
      atlasStyle->SetFrameFillColor(icol);
      atlasStyle->SetCanvasBorderMode(icol);
      atlasStyle->SetCanvasColor(icol);
      atlasStyle->SetPadBorderMode(icol);
      atlasStyle->SetPadColor(icol);
      atlasStyle->SetStatColor(icol);
      //atlasStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

      // set the paper & margin sizes
      atlasStyle->SetPaperSize(20,26);

      // set margin sizes
      atlasStyle->SetPadTopMargin(0.05);
      atlasStyle->SetPadRightMargin(0.05);
      atlasStyle->SetPadBottomMargin(0.16);
      atlasStyle->SetPadLeftMargin(0.16);

      // set title offsets (for axis label)
      atlasStyle->SetTitleXOffset(1.4);
      atlasStyle->SetTitleYOffset(1.4);

      // use large fonts
      //Int_t font=72; // Helvetica italics
      Int_t font=42; // Helvetica
      Double_t tsize=0.05;
      atlasStyle->SetTextFont(font);

      atlasStyle->SetTextSize(tsize);
      atlasStyle->SetLabelFont(font,"x");
      atlasStyle->SetTitleFont(font,"x");
      atlasStyle->SetLabelFont(font,"y");
      atlasStyle->SetTitleFont(font,"y");
      atlasStyle->SetLabelFont(font,"z");
      atlasStyle->SetTitleFont(font,"z");

      atlasStyle->SetLabelSize(tsize,"x");
      atlasStyle->SetTitleSize(tsize,"x");
      atlasStyle->SetLabelSize(tsize,"y");
      atlasStyle->SetTitleSize(tsize,"y");
      atlasStyle->SetLabelSize(tsize,"z");
      atlasStyle->SetTitleSize(tsize,"z");

      // use bold lines and markers
      atlasStyle->SetMarkerStyle(20);
      atlasStyle->SetMarkerSize(1.2);
      atlasStyle->SetHistLineWidth(2.);
      atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

      // get rid of X error bars
      //atlasStyle->SetErrorX(0.001);
      // get rid of error bar caps
      atlasStyle->SetEndErrorSize(0.);

      // do not display any of the standard histogram decorations
      atlasStyle->SetOptTitle(0);
      //atlasStyle->SetOptStat(1111);
      atlasStyle->SetOptStat(0);
      //atlasStyle->SetOptFit(1111);
      atlasStyle->SetOptFit(0);

      // put tick marks on top and RHS of plots
      atlasStyle->SetPadTickX(1);
      atlasStyle->SetPadTickY(1);

      return atlasStyle;

}

/**
 * @brief Function inherited from rcdaq. Set the ATLAS style.
 */
void Visualizer::SetAtlasStyle(){

    static TStyle* atlasStyle = 0;
    std::cout << "\nApplying ATLAS style settings...\n" << std::endl ;
    if ( atlasStyle==0 ) atlasStyle = AtlasStyle();
    gROOT->SetStyle("ATLAS");
    gROOT->ForceStyle();

}


/** @brief OverlayHistos method for Visualizer
 *
 *  Plots two input histograms on the same, specified pad with
 *  separate axis. Draws two horizontal lines at +_line and -_line
 *  on h2's axis if f _line is non-zero.
 *
 *  @param h1 - Base histogram (left y-axis)
 *  @param h2 - Overlayed histogram (right y-axis)
 *  @param pad - Address of a pad to be drawn on
 *  @param _line - y value for which a horizontal line will be drawn
 */
void Visualizer::OverlayHistos( TH1 *h1, TH1 *h2 , TVirtualPad* pad, double _line){

    // If there is no pad or no data in the histograms, return
    if( pad == nullptr ) {std::cerr<< "WARNING: No pad to overlay histos onto" << std::endl; return;}
    if( !h1->GetMinimum() && !h1->GetMaximum()) {std::cerr << "WARNING: "<< h1->GetTitle() << " is empty. Can't overlay" << std::endl; return;}

    //Remove Stat box and double the y-axis range to include negative values
    gStyle->SetOptStat( kFALSE );
    h1->DrawCopy();
    //h1->SetAxisRange( -h1->GetMaximum()*1.1, h1->GetMaximum()*1.1, "Y");
    h1->SetAxisRange( -500., 500., "Y");
    pad->Update();

   //scale h2 to the pad coordinates
   //float rightmax = 1.1*h2->GetMaximum();
   float rightmax = 1.1*2500;
   float scale = gPad->GetUymax()/rightmax;
   h2->SetLineColor( kRed );
   h2->Scale( scale );
   h2->DrawCopy( "same" );

   //draw an axis on the right side
   TGaxis axis( gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), -rightmax, rightmax, 510, "+L");
   axis.SetLineColor( kRed );
   axis.SetLabelColor( kRed );
   axis.DrawClone();

   //Draw two horizontal lines
   if(_line != 0){
     TLine lineLow (0, -_line*scale, h1->GetNbinsX(), -_line*scale );
     TLine lineHigh(0,  _line*scale, h1->GetNbinsX(),  _line*scale );
     lineLow.SetLineColor ( kGreen );
     lineHigh.SetLineColor( kGreen );
     lineLow.DrawClone( );
     lineHigh.DrawClone( );
   }

   //Rescale h2 back to the original size
   h2->Scale( 1/scale );
   if( m_debug ) pad->Print( Form( "%s_Overlay.pdf", h1->GetTitle() ) ) ;
}

/**
 * @brief Implementation of Visualizer::ManyPadsPlot. This method takes two vector of histograms that needs to be treated in a peculiar way and
 * takes care of plot those together in a canvas with a given number of columns and rows.
 * @param _first_form - vector of the first type of histograms to be plotted
 * @param _second_form - vector of the second type of histograms to be plotted s
 * @param _ncol - number of columns of the canvas
 * @param _nrow - number of rows of the canvas
 * @param _out_name - name of the plot (w/o extension, that's defined by a data member [ .pdf by default ]
 * @param _treatment - treatment of the plots (at the moment only one available, overlay)
 */
void Visualizer::ManyPadsPlot( std::vector< TH1* > _first_form, std::vector< TH1* > _second_form, int _ncol, int _nrow, std::string _out_name, TString _treatment){
    _treatment.ToLower();
    if(_treatment != "overlay"){
        std::cerr << "WARNING!!! You're looking for a treatment that has not been implemented yet! Please check it carefully" << std::endl;
        std::cerr << "Exiting w/o doing anything .." << std::endl;
        return;
    }

    if(_first_form.size() != _second_form.size())
        std::cerr << "WARNING!!! The two vectors of histograms "
                     "have different size. May result in a crash..." << std::endl;

    if( _first_form.size() < _ncol*_nrow ||  _second_form.size() < _ncol*_nrow )
        std::cerr << "WARNING!!! You have selected a vector of histrograms that will not fill all your pads "
                     "This may result in a crash..." << std::endl;

    if( _first_form.size() > _ncol*_nrow ||  _second_form.size() > _ncol*_nrow )
        std::cerr << "WARNING!!! You have selected a vector of histrograms that is bigger than the number of requested pads "
                     "This may result in histograms lost w/o plotting..." << std::endl;

    //This for the moment is hardcoded. Maybe can be moved to data member to have more general usage also for single plots.
    int squared_pad_size = 300;
    TCanvas* canv = new TCanvas( _out_name.c_str(),_out_name.c_str(),
                                 squared_pad_size*_ncol, squared_pad_size*_nrow);
    canv->Divide(_ncol,_nrow);
    for( int idraw = 0; idraw < _first_form.size(); idraw++){
        canv->cd(idraw+1);
        if( _treatment == "overlay"){
            OverlayHistos( _first_form.at(idraw), _second_form.at(idraw), gPad);
        }
    }
    canv->Print(( _out_name + m_extension ).c_str());

}

/**
 * @brief Implementation of Visualizer::SinglePadPlot. This method takes two vectors of floats and plots them on a single pad with the requested "treatment".
 * Currently supports scatter plot and overlay
 * @param _v1 - vector of floats. Can represent x values for a scatter, or y values of a histogram
 * @param _v2 - vector of the second type of histograms to be plotted s
 * @param _out_name - name of the plot (w/o extension, that's defined by a data member [ .pdf by default ]
 * @param _treatment - treatment of the plots (scatter, overlay, overlay with lines)
 * @param _line - y value for which a horizontal line will be drawn
 */
void Visualizer::SinglePlot( std::vector< double > _v1, std::vector< double > _v2, std::string _out_name, std::string _treatment, double _line){

    if(_v1.size() != _v2.size())
    std::cerr << "WARNING!!! The two vectors have different size. "
                 "May result in a crash..." << std::endl;

    int sw = -1;
    if( _treatment == "overlay" || _treatment == "OVERLAY" || _treatment == "Overlay" ){ sw = 0; }
    if( _treatment == "overlay with lines" || _treatment == "OVERLAY WITH LINES" || _treatment == "Overlay with lines" ){ sw = 1; }
    if( _treatment == "scatter" || _treatment == "SCATTER" || _treatment == "Scatter" ){ sw = 2; }

    if(sw == -1){
        std::cerr << "WARNING!!! You're looking for a treatment that has not been implemented yet! Please check it carefully" << std::endl;
        std::cerr << "Exiting w/o doing anything .." << std::endl;
        return;
    }

    int squared_pad_size = 300;
    TCanvas* canv = new TCanvas( _out_name.c_str(),_out_name.c_str(),
                                 squared_pad_size, squared_pad_size*2);
    //TVirtualPad *gpad = canv->cd();

    if(sw == 0 || sw == 1){
        TH1D* h1 = new TH1D( _out_name.c_str(), _out_name.c_str(), _v1.size(), _v1.at(0), _v1.at(_v1.size()-1));
        TH1D* h2 = new TH1D( (_out_name + "(1)").c_str(), (_out_name + "(1)").c_str(), _v1.size(), _v2.at(0), _v2.at(_v1.size()-1));

        //Fill each histogram separately in case they are different sizes
        for(int bin = 0; bin < _v1.size(); bin++){
            h1->SetBinContent(bin,_v1.at(bin));
        }
        for(int bin = 0; bin < _v2.size(); bin++){
            h2->SetBinContent(bin,_v2.at(bin));
        }

        //Draw the overlayed histos either with or without lines, depending on the selection
        if(sw == 0) OverlayHistos( h1, h2, canv->cd());
        if(sw == 1) OverlayHistos( h1, h2, canv->cd(), _line);
    }
    if(sw == 2){
        ScatterPlot(_v1, _v2, canv->cd());
    }

    canv->Print(( _out_name + m_extension ).c_str());
}

/**
 * @brief Draws a scatter plot from two vectors
 *
 * @param _vx - Vector of x values
 * @param _vy - Vector of y values
 * @param pad - Address of a pad to be drawn on
 */
void Visualizer::ScatterPlot( std::vector< double > _vx, std::vector< double > _vy, TVirtualPad* pad){
    //Declare TVectors using the input std::vectors
    TVectorD TVx;
    TVx.Use(_vx.size(),&_vx[0]);
    TVectorD TVy;
    TVy.Use(_vy.size(),&_vy[0]);

    TGraph g(TVx, TVy);
    g.DrawClone("ap");
}



/**
 * @brief Sets plot labels based on run number and alignment data
 *
 */
void Visualizer::SetTestBeamLabel( int runNo, Alignment* alignment ){
    if(runNo == 1){
      year = "";
      category = "Simulated";
    }
    if(runNo <= 405 && runNo > 10){
      year = "2018";
      category = "Test Beam";
    }
    if(alignment->magnet_On && alignment->target_In) beam = "Fragments - Magnet On";
    if(!alignment->magnet_On && alignment->target_In) beam = "Fragments - Magnet Off";
    if(!alignment->magnet_On && alignment->target_In && alignment->lead_In) beam = "Fragments - Magnet Off - Pb absorber";
    if(!alignment->magnet_On && !alignment->target_In && !alignment->lead_In) beam = "Pb ions, 150 GeV/A";
}


void Visualizer::DrawPlot(TH1 *h2, std::string _xTitle, std::string _yTitle, std::string _saveName, std::string _opt, TMarker* _marker, std::string _oFolder)
{
    TCanvas* c1 = new TCanvas(_saveName.c_str(),"c1",600,500);
    c1->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.12);
    gPad->SetTopMargin(0.10);
    h2->Draw(_opt.c_str());
    h2->SetFillColorAlpha(kAzure+1,0.3);
    h2->SetLineColor(kAzure+1);
    h2->GetXaxis()->SetTitle(_xTitle.c_str());
    h2->GetXaxis()->SetTitleOffset(1.4);
    h2->GetYaxis()->SetTitle(_yTitle.c_str());
    h2->GetYaxis()->SetTitleOffset(1.4);
    TLatex* lx = new TLatex();
    lx->SetTextFont( 62 );
    lx->SetTextSize( 0.048 );
    if (!h2->InheritsFrom("TH2")){
        lx->DrawLatexNDC(0.49,0.83,("ZDC " + category + year).c_str());
        lx->SetTextFont( 42 );
        lx->SetTextSize( 0.038 );
        lx->DrawLatexNDC(0.49,0.79,beam.c_str());
        lx->SetTextFont( 52 );
        lx->DrawLatexNDC(0.49,0.75,"Ongoing Analysis");
    }
    if (h2->InheritsFrom("TH2")){
        lx->DrawLatexNDC(0.15,0.93,("ZDC " + category + year).c_str());
        lx->SetTextFont( 42 );
        lx->SetTextSize( 0.038 );
        lx->DrawLatexNDC(0.55,0.95,beam.c_str());
        lx->SetTextFont( 52 );
        lx->DrawLatexNDC(0.55,0.91,"Ongoing Analysis");
    }
    if(_marker){
        _marker->Draw();
    }
    if( _oFolder == "" ){ _oFolder = m_oFolder; }
    c1->Print( ( _oFolder + _saveName).c_str() );
    delete c1;

}
