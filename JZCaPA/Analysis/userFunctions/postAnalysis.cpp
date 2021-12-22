/** @file postAnalysis.cpp
 *  @brief Example of an analysis of the output of the primary JZCaPA analysis
 *
 *  Takes run number as an argument and searches $JZCaPA/results directory for the data.
 *  This is not for the analysis of raw data, but instead is for faster analyis of pre-processed data.
 *  Be sure to add/remove branches based on your output tree.
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
*/
#include <stdlib.h>
#include "Visualizer.h"
#include "Containers.h"
#include "DataReader.h"
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TImage.h>
#include <TText.h>
#include <TASImage.h>
#include <TH2.h>

using namespace std;

Visualizer* viz = new Visualizer( "ATLAS" );
TFile *f;
TTree *t;

//Make variables to hold ZDC data
double zdc1_Charge, zdc1_Peak_max, zdc1_Peak_time, zdc1_Diff_Peak_time, zdc2_Charge, zdc2_Peak_max, zdc2_Peak_time, zdc2_Diff_Peak_time;
int runNumber, evNo, zdc1_Peak_center, zdc1_Diff_Peak_center, zdc2_Peak_center, zdc2_Diff_Peak_center;
//Make variables to hold RPD data
double rpd_xCoM, rpd_yCoM, rpd_Charge_sum, rpd_Peak_sum, rpd_Charge[5][5], rpd_Peak_max[5][5], rpd_Diff_max[5][5];
int rpd_Peak_center[5][5], rpd_Diff_Peak_center[5][5];

// Assign branches with global variables
void SetBranches();
// Feed the visualizer run data
void SetupVisualizer( int runNo, string outputPath );


int main(int argc, char *argv[]){

    //Argument is run number
    int runNo = atoi(argv[1]);

    //Be sure to set the path to your data
    string Path = Form("/data/phenix/data/TestBeam2018/Post_processing/run%d/", runNo);

    //Open the file and exit if it fails, then get the tree and set branch addresses
    f = new TFile( Form("%soutput%d.root", Path.c_str(), runNo) );
    if( f->IsZombie() ){
        cerr << "File didn't open... closing" << endl;
        exit(-1);
    }
    t = (TTree*)f->Get("AnalysisTree");
    SetBranches();

    //Set the visualizer up with run info from config files
    SetupVisualizer( runNo, Path );

    int nEntries = t->GetEntries();
    for(int ev = 0; ev < nEntries; ev++){
      t->GetEntry(ev);


    }

    delete f;
    return 0;
}


/* Assign branches to global variables here
 *
 */
void SetBranches(){
  //Set branch addresses
  t->SetBranchAddress("evNo", &evNo);

  //ZDC1 branches
  t->SetBranchAddress("zdc1_Charge",           &zdc1_Charge);
  t->SetBranchAddress("zdc1_Peak_max",         &zdc1_Peak_max);
  t->SetBranchAddress("zdc1_Peak_center",      &zdc1_Peak_center);
  t->SetBranchAddress("zdc1_Peak_time",        &zdc1_Peak_time);
  t->SetBranchAddress("zdc1_Diff_Peak_center", &zdc1_Diff_Peak_center);
  t->SetBranchAddress("zdc1_Diff_Peak_time",   &zdc1_Diff_Peak_time);

  //ZDC2 branches
  t->SetBranchAddress("zdc2_Charge",           &zdc2_Charge);
  t->SetBranchAddress("zdc2_Peak_max",         &zdc2_Peak_max);
  t->SetBranchAddress("zdc2_Peak_center",      &zdc2_Peak_center);
  t->SetBranchAddress("zdc2_Peak_time",        &zdc2_Peak_time);
  t->SetBranchAddress("zdc2_Diff_Peak_center", &zdc2_Diff_Peak_center);
  t->SetBranchAddress("zdc2_Diff_Peak_time",   &zdc2_Diff_Peak_time);

  //RPD branches
  t->SetBranchAddress("rpd_xCoM",       &rpd_xCoM);
  t->SetBranchAddress("rpd_yCoM",       &rpd_yCoM);
  t->SetBranchAddress("rpd_Charge_sum", &rpd_Charge_sum);
  t->SetBranchAddress("rpd_Peak_sum",   &rpd_Peak_sum);
  for(int row = 1; row <= 4; row++){
    for(int col = 1; col <=4; col++){
      t->SetBranchAddress(Form("rpd%d_%d_Charge",row,col), &rpd_Charge[row][col]);
      t->SetBranchAddress(Form("rpd%d_%d_Peak_max",row,col), &rpd_Peak_max[row][col]);
      t->SetBranchAddress(Form("rpd%d_%d_Diff_max",row,col), &rpd_Diff_max[row][col]);
      t->SetBranchAddress(Form("rpd%d_%d_Peak_center",row,col), &rpd_Peak_center[row][col]);
      t->SetBranchAddress(Form("rpd%d_%d_Diff_Peak_center",row,col), &rpd_Diff_Peak_center[row][col]);
    }
  }
}

/* Feed the visualizer run data
 *
 */
void SetupVisualizer( int runNo, string outputPath ){

  //Create a DataReader to get run info and label plots with our Visualizer
  DataReader* r = new DataReader( 20, 1024, "", runNo );
  r->LoadAlignmentFile();

  //Set plot labels and output location in the visualizer
  viz->SetTestBeamLabel( runNo, r->GetAlignment() );
  viz->SetOutputDirectory( outputPath );

}
