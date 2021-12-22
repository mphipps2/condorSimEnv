/** @file RPD_Matrix_Analysis.cpp
 *  @brief Matrix elements analysis for the RPD
 *  @author Riccardo Longo, Rachel Ratvasky
 *  @bug No known bugs and overlaps.
 */

#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <regex>


#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TKey.h>
#include <TIterator.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TFile.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLatex.h>
#include <Visualizer.h>
#include <TPaveStats.h>
#include <TPad.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>

using namespace std;

int choice;
vector < vector < TH1* > > H1;
vector < vector < TGraphErrors* > > G1; //Outer: Tile number //Inner:
vector < vector < string > > N1; //Outer: Tile number //Inner:
string filename;
string fold_dir;
string beam;
TH1D* hsom;
string var_name, var_title;
int marker_col;
bool enable_gaus_fit(false);

vector < TH1* >  LoadHistogramsFromFile(string _filename){

    vector < TH1* > BFH1;
    TFile* f1 = TFile::Open(_filename.c_str());
    TIter next(f1->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
       TClass *cl = gROOT->GetClass(key->GetClassName());
       if (!cl->InheritsFrom("TH1")) continue;
       TH1 *h = (TH1*)key->ReadObj();
       BFH1.push_back(h);
    }
    cout << "Load Histograms detected " << BFH1.size() << " TH1 histograms in file " << _filename << endl;
    return BFH1;
}

void EditFitandPrint(TH1* h1){

    TCanvas* c1 = new TCanvas(h1->GetName(),"c1",600,500);
    c1->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.12);
    gPad->SetTopMargin(0.15);
    h1->Draw();
    h1->SetFillColorAlpha(kAzure+1,0.3);
    h1->SetLineColor(kAzure+1);
    h1->GetXaxis()->SetTitle(h1->GetTitle());
    h1->GetXaxis()->SetTitleOffset(1.4);
    h1->GetYaxis()->SetTitle("Counts");
    h1->GetYaxis()->SetTitleOffset(1.4);
    TLatex* lx = new TLatex();
    lx->SetTextFont( 62 );
    lx->SetTextSize( 0.048 );
    if (h1->InheritsFrom("TH1")){
        h1->GetYaxis()->SetRangeUser(0,h1->GetMaximum()*1.4);
        lx->DrawLatexNDC(0.4,0.79,"RPD - 2018 configuration");
        lx->SetTextFont( 42 );
        lx->SetTextSize( 0.038 );
        lx->DrawLatexNDC(0.4,0.75,beam.c_str());
        lx->SetTextFont( 52 );
        lx->DrawLatexNDC(0.4,0.71,"Ongoing MC Analysis");
    }
    h1->Fit("gaus");

    gPad->Update();
    TPaveStats* pv = (TPaveStats*)h1->FindObject("stats");
    pv->SetX1NDC(0.15);
    pv->SetX2NDC(0.875);
    pv->SetY1NDC(0.86);
    pv->SetY2NDC(0.995);
    pv->SetTextSize(0.0325);
    pv->SetTextColor(kBlack);
    pv->SetTextFont(42);
    gPad->Update();
    c1->Print((fold_dir + "/" + h1->GetName() + ".pdf").c_str());

    TF1 *myfunc = h1->GetFunction("gaus");
    hsom->Fill(myfunc->GetParameter(2)/myfunc->GetParameter(1));

}

void ProduceDependecyHistograms(vector < vector < TH1* > > vH1, double *_x_val, int nLines){

    for(int i = 0; i < vH1.at(0).size(); i++){
        double y_val[nLines];
        double x_val[nLines];
        double ey_val[nLines];
        double ex_val[nLines];

        for(int j = 0; j < vH1.size(); j++){
            if(!enable_gaus_fit)
            {
                y_val[j] = vH1.at(j).at(i)->GetMean();
                ey_val[j] = vH1.at(j).at(i)->GetRMS();
                ex_val[j] = 0;
            }
            else if(enable_gaus_fit && vH1.at(j).at(i)->GetMean() != 0){
                vH1.at(j).at(i)->Fit("gaus");
                TF1 *myfunc = vH1.at(j).at(i)->GetFunction("gaus");
                y_val[j] = myfunc->GetParameter(1);
                ey_val[j] = myfunc->GetParameter(2);
                ex_val[j] = 0;
            }
            else if(enable_gaus_fit && vH1.at(j).at(i)->GetMean() == 0){
                y_val[j] = 0;
                ey_val[j] = 0;
                ex_val[j] = 0;
            }
        }//End of Outer Loop

        int tile_giving;
        int tile_receiving;
        string hname = vH1.at(0).at(i)->GetName();
        sscanf( hname.c_str(), "%*[^_]_%d_%d", &tile_giving, &tile_receiving );
        //sanity check
        //cout << "a_" << tile_giving << "_" << tile_receiving << endl;
        ostringstream tg, tr;
        tg << tile_giving;
        tr << tile_receiving;
        hname = "a_{" + tg.str() + "," + tr.str() + "}";
        if(y_val[0] != 0 || y_val[1] != 0 || y_val[2] != 0) N1.at(tile_receiving-1).push_back(hname);

        double xOff = 0;
        if(N1.at(tile_receiving-1).size() > 1) xOff = (N1.at(tile_receiving-1).size()-1)*(_x_val[nLines-1]-_x_val[0])/130;
        for(int k = 0; k < nLines; k++){
           x_val[k] = _x_val[k]+xOff;
        }

        TGraphErrors* gr = new TGraphErrors(nLines, x_val, y_val, ex_val, ey_val );
        TCanvas* cxx = new TCanvas();
        cxx->cd();
        gPad->SetTopMargin(0.15);
        gr->SetMarkerStyle(21);
        gr->SetMarkerColor(marker_col);
        gr->GetXaxis()->SetTitle(var_title.c_str());
        gr->GetYaxis()->SetTitle(("<"+ (string)vH1.at(0).at(i)->GetTitle()+ "> #pm #sigma_{" + (string)vH1.at(0).at(i)->GetTitle() + "}").c_str());
        gr->Draw("AP");
        gr->GetYaxis()->SetRangeUser(0,1.6);
        if(y_val[0] != 0 || y_val[1] != 0 || y_val[2] != 0) G1.at(tile_receiving-1).push_back(gr);

        TLatex* lx = new TLatex();
        lx->SetTextFont( 62 );
        lx->SetTextSize( 0.048 );
        lx->DrawLatexNDC(0.4,0.95,"RPD - 2018 configuration");
        lx->SetTextFont( 42 );
        lx->SetTextSize( 0.038 );
        lx->DrawLatexNDC(0.4,0.91,beam.c_str());
        lx->SetTextFont( 52 );
        lx->DrawLatexNDC(0.4,0.87,"Ongoing Monte Carlo Analysis");
        cxx->Print((var_name + "/" + vH1.at(0).at(i)->GetName() + ".pdf").c_str());

    }//End of Inner Loop
}

void ProduceOverlayPlots(){
    TCanvas* cg = new TCanvas();
     gStyle->SetPalette(kRainBow);
    for(int i = 0; i < G1.size(); i++){
        cg->cd();
        TLegend* leg1 = new TLegend(0.16,0.86,0.5,0.99);
        leg1->SetNColumns(2);
        leg1->SetTextSize(0.045);
        leg1->SetTextFont(42);
        gPad->SetTopMargin(0.15);
        G1.at(i).at(0)->GetYaxis()->SetTitle("<a_{i,j}> #pm #sigma_{a_{i,j}}");
        for(int j = 0; j < G1.at(i).size(); j++){
                if( j == 0 ) G1.at(i).at(j)->Draw("PA PMC");
                else G1.at(i).at(j)->Draw("P PMC");
                leg1->AddEntry(G1.at(i).at(j),N1.at(i).at(j).c_str(),"p");
        }
        ostringstream tile;
        tile << i+1;
        TLatex* lx = new TLatex();
        lx->SetTextFont( 62 );
        lx->SetTextSize( 0.048 );
        lx->DrawLatexNDC(0.55,0.95,"RPD - 2018 configuration");
        lx->SetTextFont( 42 );
        lx->SetTextSize( 0.038 );
        lx->DrawLatexNDC(0.55,0.91,("Tile " + tile.str() + ": cross talking (input)").c_str());
        lx->SetTextFont( 52 );
        lx->DrawLatexNDC(0.55,0.87,"Ongoing Monte Carlo Analysis");
        leg1->SetBorderSize(1);
        leg1->Draw();
        cg->Print((var_name + "/same_tile/" + tile.str() + ".pdf" ).c_str());
        cg->Clear();
    }

}

int main(int argc, char *argv[]){

  //For option 1 filename must be a root file. For option 2, filename must be the name of a list of file
  //Structure of the list --> var_value root_file_name
  filename = argv[1];
  beam = argv[2];
  cout << "Working on " << filename << endl;
  cout << "File related to beam: " << beam << endl;
  cout << "Select what you want to do: " << endl;
  cout << "1) Plot coeffients histograms stored in " << filename << " file" << endl;
  cout << "2) Inspect dependency in different runs " << endl;

  cin >> choice;

  cout << "You selected " << choice  << endl;

  //Visualizer* v1 = new Visualizer("ATLAS", true);
  Visualizer* v1 = new Visualizer("ATLAS");

  if( choice == 1 ){
      if(filename.substr(filename.find(".")+1,filename.length()) != "root"){
          cout << "ERROR IN INPUT FILE" << endl;
          return 0;
      }
      //Creating folder for output
      fold_dir = filename.substr(0, filename.find("."));
      cout << "Creating folder " << fold_dir << endl;
      gSystem->Exec(("rm -r " + fold_dir).c_str());
      gSystem->Exec(("mkdir " + fold_dir).c_str());
      //Loading histograms from file filename
      H1.push_back(LoadHistogramsFromFile(filename));
      //Create a histo to handle sigma over mean ratio
      hsom = new TH1D("sigma_over_mean","#sigma_{a_{i,j}}/#mu_{a_{i,j}}",50,0,2);
      //Editing and printing those guys
      for(int i = 0; i < H1.size(); i++){
          for(int j = 0; j < H1.at(i).size(); j++){
                EditFitandPrint(H1.at(i).at(j));
          }
      }

      EditFitandPrint(hsom);
  }

  if(choice == 2){
    if(filename.substr(filename.find(".")+1,filename.length()) == "root"){
        cout << "ERROR IN INPUT FILE! NO ROOT FOR OPT 2" << endl;
        return 0;
    }
    ifstream FileList(filename);
    int nLines = count(std::istreambuf_iterator<char>(FileList),
                 std::istreambuf_iterator<char>(), '\n');
    cout << "Your list has " << nLines << " files that will go into data points " << endl;

    int sub_choice;

    cout << "Select what is your dependency now: " << endl;
    cout << "1) z vertex [mm]" << endl;
    cout << "2) Beam Energy [TeV]" << endl;
    cout << "3) Beam x-position [cm]" << endl;
    cout << "4) Beam y-position [cm]" << endl;
    cin >> sub_choice;
    if(sub_choice == 1){
            var_name = "z_vtx"; var_title = "Z^{Shower}_{Init} [mm]"; marker_col = kBlue+1;
    }
    if(sub_choice == 2){
            var_name = "E_beam"; var_title = "E_{Beam} [TeV]"; marker_col = kRed+1;
    }
    if(sub_choice == 3){
            var_name = "x_beam"; var_title = "x_{Beam} [cm]"; marker_col = kOrange+1;
    }
    if(sub_choice == 4){
            var_name = "y_beam"; var_title = "y_{Beam} [cm]"; marker_col = kGreen+3;
    }

    double x_val[nLines];
    string buffer_str;
    vector < string > analysis_files;
    int count = 0;
    //Going back to the beginning of the file
    FileList.seekg(0);
    while(FileList >> x_val[count] >> buffer_str){
        analysis_files.push_back(buffer_str);
        count++;
    }

    gSystem->Exec(("rm -r " + var_name).c_str());
    cout << "Creating folder " << var_name << endl;
    gSystem->Exec(("mkdir " + var_name).c_str());
    gSystem->Exec(("mkdir " + var_name + "/same_tile").c_str());

    for(int i = 0; i < analysis_files.size(); i++)    H1.push_back(LoadHistogramsFromFile(analysis_files.at(i)));

    cout << "Vector of histograms size: " << H1.size() << endl;
    //Creating vectors of histograms for each tile
    for(int i = 0; i < 16; i++)
    {
        vector < TGraphErrors* > GE;
        vector < string > name_string;
        G1.push_back(GE);
        N1.push_back(name_string);
    }
    cout << "Here we go " << endl;
    ProduceDependecyHistograms(H1,x_val,nLines);
    for(int i = 0; i < G1.size(); i++)
    {
        cout << "Tile " << i << " has " << G1.at(i).size() << " coefficient(s) " <<  endl;
        for(int j = 0; j < N1.at(i).size(); j++) { cout << N1.at(i).at(j) << " ::: "; }
        cout << " " << endl;
    }

    ProduceOverlayPlots();
  }
  return 0;
}
