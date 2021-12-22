/** @ingroup ana
 *  @file ZDCAnalysis.cpp
 *  @brief Implementation of ZDCAnalysis.
 *
 *
 *  Function definitions for ZDCAnalysis are provided.
 *  This class is the main class for the ZDC analysis.
 *  The analysis is done on both ZDCs simultaneously
 *
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
*/


#include "ZDCAnalysis.h"
#include "ZDC.h"
#include "Containers.h"
#include "TH2D.h"
#include "TH1.h"
#include "TCanvas.h"
#include "Visualizer.h"

/** @brief Default Constructor for ZDCAnalysis.
 */
ZDCAnalysis::ZDCAnalysis( ){

}

/** @brief Destructor for ZDCAnalysis.
 */
ZDCAnalysis::~ZDCAnalysis( ){

}

/** @brief Initialization method for ZDCAnalysis
 *
 *  Can add other things here that you would
 *  perhaps not put into the constructor.
 *  I.e. a TTree, some tools. Etc.
 *
 */
void ZDCAnalysis::Initialize( ){

}

/** @brief Initialization method for ZDCAnalysis
 *
 *  Takes a vector of detectors, picks out the ZDCs
 *  and assigns them to member pointers
 *
 */
void ZDCAnalysis::Initialize( std::vector < Detector* > _vDet ){

    for( auto& det : _vDet ){
        if( det->GetChannelsVector().at(0)->detector == "ZDC" && det->GetChannelsVector()[0]->mapping_column == 1){
            m_zdc1 = (ZDC*)det;
            zdc1 = det->GetElement(0,1);
        }
        if( det->GetChannelsVector().at(0)->detector == "ZDC" && det->GetChannelsVector()[0]->mapping_column == 2){
            m_zdc2 = (ZDC*)det;
            zdc2 = det->GetElement(0,2);
        }
    }
    m_alignment = m_zdc1->GetAlignment();
}


/** @brief Historgam Setup method for ZDCAnalysis
 *
 *  Should instantiate any histograms you wish to output here.
 *
 *  @return none
 */
void ZDCAnalysis::SetupHistograms( ){

    //TH1D
    hChargeRatio = new TH1D("ChargeZDC2_over_ChargeZDC1","Q_{ZDC2}/Q_{ZDC1}",100,0,5);
    //hChargeRatio->SetCanExtend(TH1::kXaxis);

    hChargeSum = new TH1D("Charge_sum","Q_{total} (pC)",300,0,450);

    hPeakRatio = new TH1D("PeakZDC2_over_PeakZDC1","Peak_{ZDC2}/Peak_{ZDC1}",100,0,5);
    //hPeakRatio->SetCanExtend(TH1::kXaxis);

    hCharge1 = new TH1D("ZDC1_Charge","Q_{ZDC1} (pC)",200,0,200);

    hCharge2 = new TH1D("ZDC2_Charge","Q_{ZDC2} (pC)",200,0,200);

    hPeak1   = new TH1D("ZDC1_Peak","Peak_{ZDC1}",200,0,1000);

    hPeak2   = new TH1D("ZDC2_Peak","Peak_{ZDC2}",200,0,1000);

    hDpeak1  = new TH1D("ZDC1_Diff_Peak","#frac{#partial V}{#partial t}_{max ZDC1}",200,0,4500);

    hDpeak2  = new TH1D("ZDC2_Diff_Peak","#frac{#partial V}{#partial t}_{max ZDC2}",200,0,4500);

    hArrival1  = new TH1D("ZDC1_Arrival_time","Peak Center_{ZDC1} (ns)",150,240,390);

    hArrival2  = new TH1D("ZDC2_Arrival_time","Peak Center_{ZDC2} (ns)",150,240,390);

    hToF  = new TH1D("Time_of_Flight","Peak_Center_{ZDC1-ZDC2} (ns)",30,-15,15);

    //TH2D
    hCharge = new TH2D("ZDC_Charge_Correlation", "ZDC Charge Correlation", 100, 0, 450, 100, 0, 450);
    //hCharge->SetCanExtend(TH1::kAllAxes);

    hPeak   = new TH2D("ZDC_Peak_Correlation", "ZDC Peak Correlation", 100, 0, 1000, 100, 0, 1000);
    //hPeak->SetCanExtend(TH1::kAllAxes);

    hDpeak = new TH2D("ZDC_Diff_Peak_Correlation", "ZDC Diff Peak Correlation", 200, 0, 5000, 200, 0, 5000);

    hChargePeakZDC1 = new TH2D("ZDC1_ChargePeakCorrelation","Q_{ZDC1} (pC) vs Peak_{ZDC1}",50,0,200,50,0,1000);
    //hChargePeakZDC1->SetCanExtend(TH1::kAllAxes);

    hChargePeakZDC2 = new TH2D("ZDC2_ChargePeakCorrelation","Q_{ZDC2} (pC) vs Peak_{ZDC2}",50,0,200,50,0,1000);
    //hChargePeakZDC2->SetCanExtend(TH1::kAllAxes);





}

/** @brief Branch setup method for ZDCAnalysis
 *
 *  Adds branches with data created by the analysis
 *
 */
void ZDCAnalysis::SetBranches( TTree* _tree ){
    m_AnalysisTree = _tree;

    m_AnalysisTree->Branch("zdc1_Charge",           &zdc1->Charge,           "zdc1->Charge/D" );
    m_AnalysisTree->Branch("zdc1_Peak_max",         &zdc1->Peak_max,         "zdc1->Peak_max/D" );
    m_AnalysisTree->Branch("zdc1_Diff_max",         &zdc1->Diff_max,         "zdc1->Diff_max/D" );
    m_AnalysisTree->Branch("zdc1_Peak_center",      &zdc1->Peak_center,      "zdc1->Peak_center/I" );
    m_AnalysisTree->Branch("zdc1_Peak_time",        &zdc1->Peak_time,        "zdc1->Peak_Peak_time/D" );
    m_AnalysisTree->Branch("zdc1_Diff_Peak_center", &zdc1->Diff_Peak_center, "zdc1->Diff_Peak_center/I" );
    m_AnalysisTree->Branch("zdc1_Diff_Peak_time",   &zdc1->Diff_Peak_time,   "zdc1->Diff_Peak_Peak_time/D" );

    m_AnalysisTree->Branch("zdc2_Charge",           &zdc2->Charge,           "zdc2->Charge/D" );
    m_AnalysisTree->Branch("zdc2_Peak_max",         &zdc2->Peak_max,         "zdc2->Peak_max/D" );
    m_AnalysisTree->Branch("zdc2_Diff_max",         &zdc2->Diff_max,         "zdc2->Diff_max/D" );
    m_AnalysisTree->Branch("zdc2_Peak_center",      &zdc2->Peak_center,      "zdc2->Peak_center/I" );
    m_AnalysisTree->Branch("zdc2_Peak_time",        &zdc2->Peak_time,        "zdc2->Peak_Peak_time/D" );
    m_AnalysisTree->Branch("zdc2_Diff_Peak_center", &zdc2->Diff_Peak_center, "zdc2->Diff_Peak_center/I" );
    m_AnalysisTree->Branch("zdc2_Diff_Peak_time",   &zdc2->Diff_Peak_time,   "zdc2->Diff_Peak_Peak_time/D" );

}

/** @brief Analyze Events method for ZDCAnalysis
 *
 *
 */
void ZDCAnalysis::AnalyzeEvent( ){

    //If neither ZDC was saturated
    if(zdc1->Peak_max<750.0 && zdc2->Peak_max<750.0 &&  zdc1->was_hit && zdc2->was_hit){

        hChargeRatio->Fill(zdc2->Charge/zdc1->Charge);
        hPeakRatio->Fill(zdc2->Peak_max/zdc1->Peak_max);

        hCharge->Fill( zdc1->Charge,   zdc2->Charge );
        hPeak->Fill(   zdc1->Peak_max, zdc2->Peak_max );

        hChargePeakZDC1->Fill(zdc1->Charge,zdc1->Peak_max);
        hChargePeakZDC2->Fill(zdc2->Charge,zdc2->Peak_max);

        hCharge2->Fill(zdc2->Charge);
        hChargeSum->Fill(zdc1->Charge + zdc2->Charge);

        hToF->Fill(zdc2->Diff_Peak_time - zdc1->Diff_Peak_time);
    }

    //If ZDC1 wasn't saturated
    if( zdc1->Peak_max<900.0 && zdc1->was_hit){
    //if( zdc1->was_hit){
        hCharge1->Fill(zdc1->Charge);
        hPeak1->Fill(zdc1->Peak_max);
        hDpeak1->Fill(zdc1->Diff_max);
        hArrival1->Fill(zdc1->Diff_Peak_time);
    }

    //If ZDC1 wasn't saturated
    if( zdc2->Peak_max<900.0 && zdc2->was_hit ){
    //if( zdc2->was_hit ){
        hCharge2->Fill(zdc2->Charge);
        hPeak2->Fill(zdc2->Peak_max);
        hDpeak2->Fill(zdc2->Diff_max);
        hArrival2->Fill(zdc2->Diff_Peak_time);
    }

}

/** @brief Finalize method for ZDCAnalysis
 *
 *
 */
void ZDCAnalysis::Finalize( ){

    std::string output =  std::getenv("JZCaPA");
    output += "/results/";

    TCanvas *c = new TCanvas("ZDCAnalysis","ZDCAnalysis",800,600);
    c->cd();

    if(m_viz == NULL) m_viz = new Visualizer( "ATLAS" );


    //Raw data plots
    m_viz->DrawPlot(hCharge1,"Q_{ZDC1} (pC)","Counts","ZDC1_Charge.png","");
    m_viz->DrawPlot(hCharge2,"Q_{ZDC2} (pC)","Counts","ZDC2_Charge.png","");
    m_viz->DrawPlot(hPeak1,"Peak_{max ZDC1} (mV)","Counts","ZDC1_Peak.png","");
    m_viz->DrawPlot(hPeak2,"Peak_{max ZDC1} (mV)","Counts","ZDC2_Peak.png","");
    m_viz->DrawPlot(hDpeak1,"#frac{#partial V}{#partial t}_{max ZDC1}","Counts","ZDC1_DiffPeak.png","");
    m_viz->DrawPlot(hDpeak2,"#frac{#partial V}{#partial t}_{max ZDC2}","Counts","ZDC2_DiffPeak.png","");
    m_viz->DrawPlot(hArrival1,"Arrival time_{ZDC1} (ns)","Counts","ZDC1_Arrival.png","");
    m_viz->DrawPlot(hArrival2,"Arrival time_{ZDC2} (ns)","Counts","ZDC2_Arrival.png","");

    //Correlation plots
    m_viz->DrawPlot(hCharge,"Q_{ZDC1} (pC)","Q_{ZDC2} (pC)","ZDC_charge.png","COLZ");
    m_viz->DrawPlot(hPeak,"Peak_{ZDC1} (mV)","Peak_{ZDC2} (mV)","ZDC_peak.png","COLZ");
    m_viz->DrawPlot(hChargePeakZDC1,"Q_{ZDC1} (pC)","Peak_{ZDC1} (mV)","ZDC1_ChargePeak.png","COLZ");
    m_viz->DrawPlot(hChargePeakZDC2,"Q_{ZDC2} (pC)","Peak_{ZDC2} (mV)","ZDC2_ChargePeak.png","COLZ");
    m_viz->DrawPlot(hChargeRatio,"Q_{ZDC2}/Q_{ZDC1}","Counts","ZDC_chargeRatio.png","");
    m_viz->DrawPlot(hPeakRatio,"Peak_{ZDC2}/Peak_{ZDC1}","Counts","ZDC_peakRatio.png","");
    m_viz->DrawPlot(hToF,"Time of Flight (ns)","Counts","ZDC_ToF.png","");
    m_viz->DrawPlot(hChargeSum,"Q_{total} (pC)","Counts","ZDC_Qtot.png","");
    m_viz->DrawPlot(hDpeak,"#frac{#partial V}{#partial t}_{max ZDC1}","#frac{#partial V}{#partial t}_{max ZDC2}","ZDC_Dpeak_corr.png","COLZ");

    delete c;

}
