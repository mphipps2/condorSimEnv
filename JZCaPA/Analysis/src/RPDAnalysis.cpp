/** @ingroup ana
 *  @file RPDAnalysis.cpp
 *  @brief Implementation of RPDAnalysis.
 *
 *
 *  Function definitions for RPDAnalysis are provided.
 *  This class is the main class for the RPD analysis.
 *  The analysis is done on both RPDs simultaneously
 *
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
*/


#include "RPDAnalysis.h"
#include "RPD.h"
#include "Containers.h"
#include "TMarker.h"
#include "TH2D.h"
#include "TCanvas.h"


/** @brief Default Constructor for RPDAnalysis.
 */
RPDAnalysis::RPDAnalysis( ){

}

/** @brief Destructor for RPDAnalysis.
 */
RPDAnalysis::~RPDAnalysis( ){

}

/** @brief Initialization method for ZDCAnalysis
 *
 *  Takes a vector of detectors, picks out the RPD
 *  and channels, then assigns them to member pointers
 *
 */
void RPDAnalysis::Initialize( std::vector < Detector* > _vDet ){
  m_counter = 0;

    for( auto& det : _vDet ){
        if(det->GetChannelsVector().size() == 0) continue;
        if(det->GetChannelsVector().at(0)->detector == "RPD" ){
            m_RPD = (RPD*) det;
            m_alignment = m_RPD->GetAlignment();
            for(int row = 0; row < 4; row++){
                for(int col = 0; col < 4; col++){
                    rpd[row][col] = m_RPD->GetElement(row,col);
                }
            }
        }
    }

    int tableOffsetX, tableOffsetY;
    if(m_alignment->runNumber >=79 && m_alignment->runNumber <= 413){
      tableOffsetX = 2250;
      tableOffsetY = -500;
    }else{
      tableOffsetX = 0;
      tableOffsetY = 0;
    }

    beamPosX = -(m_alignment->x_table + tableOffsetX)/100;
    beamPosY = -(m_alignment->y_table + tableOffsetY)/100;
}

/** @brief Historgam Setup method for RPDAnalysis
 *
 *  Should instantiate any histograms you wish to output here.
 *
 *  @return none
 */
void RPDAnalysis::SetupHistograms( ){
  if(m_viz == NULL) m_viz = new Visualizer( "ATLAS" );

   hCharge     = new TH2D("RPD Charge", "average charge per tile", 4, 0, 4, 4, 0, 4);
   hPeak       = new TH2D("RPD Peak", "average peak height per tile", 4, 1, 5, 4, 1, 5);
   hCenter     = new TH2D("Shower Center", "Calculated center of mass", 200, -50, 50, 200, -50, 50);
   hChgVsPk    = new TH2D("Charge vs Peak height","Q vs Peak; Q; Peak", 300, 0, 300000, 300, 0, 1100);
   hPkVsDiffPk = new TH2D("Peak height vs Differential peak height","Peak vs Diff Peak", 200, 0, 5000, 200, 0, 300);
   hChargeSum  = new TH1D("RPD Integrated Signal","RPD Integrated Signal", 200, 0, 160000);
   hPeakSum    = new TH1D("RPD Peak Height Sum","RPD Peak Height Sum", 200, 0, 3000);
   hDiffPeakSum= new TH1D("RPD Differential Peak Height Sum","RPD Diff Peak Sum", 200, 0, 55000);

   hChargeArr.resize(4);
   hPeakArr.resize(4);
   hDPeakArr.resize(4);
   for(int row = 0; row < 4; row++){
     hChargeArr[row].resize(4);
     hPeakArr[row].resize(4);
     hDPeakArr[row].resize(4);
       for(int col = 0; col < 4; col++){
            hChargeArr[row][col] = new TH1D(Form("rpd%d_%d_Charge",row,col),Form("rpd%d_%d_Charge",row,col), 200 , 0, 15000);
            hChargeArr[row][col]->SetAxisRange(0,1500,"Y");

            hPeakArr[row][col]   = new TH1D(Form("rpd%d_%d_Peak",row,col),Form("rpd%d_%d_Peak",row,col), 200, 0, 300);
            hPeakArr[row][col]->SetAxisRange(0,1500,"Y");

            hDPeakArr[row][col]  = new TH1D(Form("rpd%d_%d_Diff_peak",row,col),Form("rpd%d_%d_Diff_peak",row,col), 200, 0, 5000);
            hDPeakArr[row][col]->SetAxisRange(0,1500,"Y");
       }
   }
}

/** @brief Branch setup method for ZDCAnalysis
 *
 *  Adds branches with data created by the analysis
 *
 */
void RPDAnalysis::SetBranches( TTree* _tree ){
    m_AnalysisTree = _tree;

    for(int row = 0; row < 4; row++){
        for(int col = 0; col < 4; col++){
            m_AnalysisTree->Branch( Form("rpd%d_%d_Charge", row, col),           &rpd[row][col]->Charge,           Form("rpd%d_%d_Charge/D", row, col) );
            m_AnalysisTree->Branch( Form("rpd%d_%d_Peak_max", row, col),         &rpd[row][col]->Peak_max,         Form("rpd%d_%d_Peak_max/D", row, col) );
            m_AnalysisTree->Branch( Form("rpd%d_%d_Diff_max", row, col),         &rpd[row][col]->Diff_max,         Form("rpd%d_%d_Diff_max/D", row, col) );
            m_AnalysisTree->Branch( Form("rpd%d_%d_Peak_center", row, col),      &rpd[row][col]->Peak_center,      Form("rpd%d_%d_Peak_center/I", row, col) );
            m_AnalysisTree->Branch( Form("rpd%d_%d_Diff_Peak_center", row, col), &rpd[row][col]->Diff_Peak_center, Form("rpd%d_%d_Diff_Peak_center/I", row, col) );
        }
    }

    m_AnalysisTree->Branch("rpd_xCoM", &xCoM, "xCoM/D" );
    m_AnalysisTree->Branch("rpd_yCoM", &yCoM, "yCoM/D" );
    m_AnalysisTree->Branch("rpd_Charge_sum", &ChargeSum, "ChargeSum/D" );
    m_AnalysisTree->Branch("rpd_Peak_sum",   &PeakSum,   "PeakSum/D" );
    m_AnalysisTree->Branch("rpd_Diff_Peak_sum",   &DiffPeakSum,   "DiffPeakSum/D" );

}

/** @brief Analyze Events method for RPDAnalysis
 *
 *
 */
void RPDAnalysis::AnalyzeEvent( ){
    ChargeSum = 0;
    PeakSum = 0;
    DiffPeakSum = 0;

    //This loop takes the running average of charge and peak height per tile
    for(int row = 0; row < 4; row++){
        for(int col = 0; col < 4; col++){
            if( rpd[row][col]->is_on ){
                ChargeSum   += rpd[row][col]->Charge;
                PeakSum     += rpd[row][col]->Peak_max;
                DiffPeakSum += rpd[row][col]->Diff_max;

                hChargeArr[row][col]->Fill(rpd[row][col]->Charge);
                hPeakArr[row][col]->Fill(rpd[row][col]->Peak_max);
                hDPeakArr[row][col]->Fill(rpd[row][col]->Diff_max);
            }
        }
    }

    //This loop finds the center of mass per event (should probably be it's own function)

    double totalCharge = 0;
    double weightedRow = 0;
    double weightedCol = 0;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if( rpd[i][j]->is_on && rpd[i][j]->was_hit ){
                totalCharge += rpd[i][j]->Charge;
                weightedRow += rpd[i][j]->Charge * yPos[i];
                weightedCol += rpd[i][j]->Charge * xPos[j];
            }
        }
    }
    xCoM = weightedCol/totalCharge;
    yCoM = weightedRow/totalCharge;
    hCenter->Fill(xCoM,yCoM);

    hChargeSum->Fill(ChargeSum);
    hPeakSum->Fill(PeakSum);
    hDiffPeakSum->Fill(DiffPeakSum);

}

/** @brief Finalize method for RPDAnalysis
 *
 *
 */
void RPDAnalysis::Finalize( ){

    if(m_viz == NULL) m_viz = new Visualizer( "ATLAS" );

    TCanvas cCharge("charge","charge",800,600);
    TCanvas cPeak("peak","peak",800,600);
    TCanvas cDpeak("dpeak","dpeak",800,600);
    cCharge.Divide(4,4);
    cPeak.Divide(4,4);
    cDpeak.Divide(4,4);
    int pad;

    for(int row = 0; row < 4; row++){
        for(int col = 0; col < 4; col++){
            pad = row*4 + col + 1;

            cCharge.cd(pad);
            hChargeArr[row][col]->Draw();

            cPeak.cd(pad);
            hPeakArr[row][col]->Draw();

            cDpeak.cd(pad);
            hDPeakArr[row][col]->Draw();

            hCharge->SetBinContent(col+1, 4-row, hChargeArr[row][col]->GetMean() );
            hPeak->SetBinContent(  col+1, 4-row, hPeakArr[row][col]->GetMean()  );
        }
    }


    std::string output =  std::getenv("JZCaPA");
    output += "/results/";

    TMarker *marker = new TMarker(beamPosX,beamPosY,5);
    marker->SetMarkerColor(kRed);
    marker->SetMarkerSize(4);

    //cCharge.Print( (output + "RPD_Charge_per_tile.png").c_str() );
    //cPeak.Print( (output + "RPD_Peak_per_tile.png").c_str() );
    //cDpeak.Print( (output + "RPD_Diff_peak_per_tile.png").c_str() );
    /*
    m_viz->DrawPlot(hChgVsPk,"RPD Q_{total}","RPD Peak_{sum}","RPD_TotalCharge.png","");
    m_viz->DrawPlot(hPkVsDiffPk,"RPD Peak_{sum}","#frac{#partial V}{#partial t}_{max}","RPD_TotalCharge.png","");
    m_viz->DrawPlot(hChargeSum,"RPD Q_{total}","Counts","RPD_TotalCharge.png","");
    m_viz->DrawPlot(hPeakSum,"RPD Peak_{sum}","Counts","RPD_PeakSum.png","");
    m_viz->DrawPlot(hDiffPeakSum,"#frac{#partial V}{#partial t}_{max}","Counts","RPD_DiffSum.png","");
    m_viz->DrawPlot(hCharge,"Col","Row","RPD_Charge.png","COLZ text",marker);
    m_viz->DrawPlot(hPeak,"Col","Row","RPD_Peak.png","COLZ text",marker);
    m_viz->DrawPlot(hCenter,"x (mm)","y (mm)","RPD_CoM.png","COLZ",marker);
    */
}
