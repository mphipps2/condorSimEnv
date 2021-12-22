#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TRandom.h"
#include <TTree.h>
#include <vector>
#include <fstream>

using namespace std;

TGraph* getPMTdata        ( string fileName );
double  evaluateSum       ( vector< TF1* >&  vec, double x );
double  GetAverageCharge ( TGraph* PMTcathodeRadSens, double gain);
void    addTimingJitter   ( vector< float >* vec, float timeInterval);
void    outputDigitization  ( string fileName, bool DRAW = false);
void    getFile           (string filename,  /*out*/ ifstream& file);
double  LTQuadStepFilt(double *t, double *amp);
double  LTQuadStepFiltKernel(double x, double a, double b, double c, double d);

/*
 * This PMT waveform function came from the ZDC_PileUpTool in Athena
 * The function was made piecewise to avoid extreme values before the
 * delay time
 * Parameter 0: Amplitude
 * Parameter 1: Delay
 * Time units are in ns
 */
Double_t PMTpulseFunction(Double_t *t, Double_t *par){
   Float_t tt = t[0];

   if(tt < par[1]) return 0.0; // cut out the time where the function behaves poorly
   //      std::cout << " amp " << par[0] << " delay " << par[1] << " tt " << tt << std::endl;
   Double_t f = par[0]*pow((tt-par[1])/10.0, 3.4)*exp(-(tt-par[1])/10.0);
   return f;
}


/*
 * Calculate the number of Cherenkov photons likely to make it out
 * of the ZDC starting from the output of the JZCaPA MC. Accounting
 * for losses due to absorption in fused silica after heavy irradiation.
 */
void cutZDCoutput(string fileName, string rod, double xOffset, double zOffset){
	int nBins = 3649;
	double energy,X,Z,Px,Py,Pz;
	TRandom rnd;

  //Create input vectors
  vector<double> *Xf=0, *Yf=0, *Zf=0, *Pxf=0, *Pyf=0, *Pzf=0, *Energyf=0, *NCherenkovs=0;
  vector<int> *PIDf=0, *IDf=0, *EventNof=0;

  //Create output vectors
  vector<double> *x=0, *z=0,*x_0=0, *z_0=0, *px=0, *py=0, *pz=0, *Energy=0;
  vector<int> *EventNo=0;


  // Set the input and output files up
  //
  TFile* inFile = new TFile(fileName.c_str());
  TTree* inTree = (TTree*)inFile->Get("ZDCtree");

  inTree->SetBranchAddress("X",&Xf);
  inTree->SetBranchAddress("Z",&Zf);
  inTree->SetBranchAddress("Px",&Pxf);
  inTree->SetBranchAddress("Py",&Pyf);
  inTree->SetBranchAddress("Pz",&Pzf);
  //inTree->SetBranchAddress("EventNo",&EventNf);
  // intree->SetBranchAddress("NCherenkovs",&NCherenkovs);

  //Output root file
  TFile* outFile = new TFile( Form("output/%s_%s.root",fileName.substr(0,fileName.find_last_of(".")).c_str(), rod.c_str() ), "RECREATE");
  TTree* outTree = new TTree( "tree", "tree" );

  // outTree->Branch("X0",&x_0);
  // outTree->Branch("Z0",&z_0);
  outTree->Branch("X",&x);
  outTree->Branch("Z",&z);
  outTree->Branch("Px",&px);
  outTree->Branch("Py",&py);
  outTree->Branch("Pz",&pz);
  outTree->Branch("Energy",&Energy);
  outTree->Branch("EventNo",&EventNo);
  // outTree->SetBranchAddress("NCherenkovs",&NCherenkovs);

  // output txt file
  ofstream outputFile( Form("output/%s_%s.txt",fileName.substr(0,fileName.find_last_of("." ) ).c_str(), rod.c_str() ) );

  //// Get the transmission data
  ////
  vector< double > transData;
  string str;
  ifstream file( Form("data/rod%s.txt", rod.c_str() ) );
  if(!file.is_open()){
    cout << Form("data/rod%s.txt didn't open", rod.c_str() ) <<  endl;
    return;
  }
  int i = 0;
  while(file){
    getline(file,str);
    transData.push_back( atof( str.c_str() ) );
    i++;
  }//End line loop
  file.close();

  // Load that data into a histogram so we can use the GetBin function
  TH1F* hTrans = new TH1F( "TransData", "TransData", nBins, 197.2888,  1027.23596);
  TH1F* hTrans2 = new TH1F( "TransData2", "TransData2", nBins, 197.2888,  1027.23596);
  float trans, rand;

  for(int bin = 0; bin < nBins; bin++){
    hTrans->SetBinContent( bin, transData[bin] );
  }

  // Apply a modified baseline to the fused quartz spectrum
  if(rod == "fusedQuartz"){
    TF1* bLine = new TF1("baseline","0.0883*log(0.013087*(x-156.8) ) + 0.802",197,1028);
    hTrans->Multiply(bLine);
  }

  //Get the number of events in the file
  int nEvents = inTree->GetEntries();
  //Find out how many events contain photons
  float wl;
  int realNevents = 0;
  int totalCut = 0, totalPhotons = 0;
  for (int ev = 0; ev < nEvents ; ev++){
    inTree->GetEntry(ev);
    if( Xf->size() > 0 ) realNevents++;
  }

  outputFile << Form("Total events = %d",realNevents) << endl;
  for (int ev = 0; ev < nEvents ; ev++){
    inTree->GetEntry(ev);

    int nCut = 0;
    int nPhotons = Xf->size();
    if(nPhotons) outputFile << Form("Event %d",ev) << endl;
    for (int k=0; k < nPhotons ; k++){	//begin loop over hits
      Px = Pxf->at(k)*1e6;
      Py = Pyf->at(k)*1e6;
      Pz = Pzf->at(k)*1e6;
      energy = sqrt( pow(Px,2) + pow(Py,2) + pow(Pz,2) );

      // If the wavelength of the photon is < 197nm, set the transmission
      // % equal to the value for the shortest wavelength we have. Otherwise,
      // get the data from the appropriate bin.
      wl = 1240./energy;
      if(wl < 197.2888){
     	trans = hTrans->GetBinContent(10);
      }else{
	trans = hTrans->GetBinContent( hTrans->FindBin( wl ) );
      }
      // Give the photon a random chance to make it weighted
      // by transmission %
      // if( rnd.Rndm() < trans ){
      if( rnd.Rndm() < 1 ){

	// If the photon is transmitted, add it to the vector
	X = Xf->at(k) + xOffset;
	Z = Zf->at(k) + zOffset;
	x->push_back( X );
	z->push_back( Z );

	//We want momentum direction, not magnitude
	Px/=energy;
	Py/=energy;
	Pz/=energy;

	px->push_back( Px );
	py->push_back( Py );
	pz->push_back( Pz );
	Energy->push_back( energy );
	outputFile << Form("V,%17.11f,%17.11f,%17.11f",X,0.0,Z) << endl;
	outputFile << Form("P,%17.14f,%17.14f,%17.14f,%17.13f",Px,Py,Pz,energy) << endl;
      }else{
	nCut++;
      }//End cuts
    }//end loop over fiber hits
    if(nPhotons){
      outputFile << Form("End event %d", ev) << endl;
      cout << Form("%7d of %7d photons cut in event %3d",nCut,nPhotons,ev) << endl;
    }
    outTree->Fill();
    totalCut += nCut;
    totalPhotons += nPhotons;

    //Clear vectors
    x->clear();
    z->clear();
    px->clear();
    py->clear();
    pz->clear();
    Energy->clear();
  }//end loop over events
  cout << Form("%5d of %5d photons cut %3.0f%%",totalCut,totalPhotons,100.*(float)totalCut/(float)totalPhotons ) << endl;
  inFile->Close();
  outFile->Close();
}





/*
 * Apply quantum efficiency cuts to the output of the lightguide MC
 */
void PMTcuts(string fileName, int PMTmodel = 6091){
  //Create input vectors
  vector<double> *Xf=0, *Yf=0, *Zf=0, *Pxf=0, *Pyf=0, *Pzf=0, *Energyf=0, *NCherenkovs=0;
  vector<int> *PIDf=0, *IDf=0, *EventNof=0;

  // Set the input and output files up
  //
  TFile* inFile = new TFile(fileName.c_str());
  if(!inFile->IsOpen()){
    cout << "Input file didn't open" << endl;
    return;
  }
  TTree* inTree = (TTree*)inFile->Get("lightGuide");

  inTree->SetBranchAddress("hitX",&Xf);
  inTree->SetBranchAddress("hitZ",&Zf);
  inTree->SetBranchAddress("energy",&Energyf);


  //// Get the quantum efficiency data
  ////

  ifstream file( Form("data/model%dQE.txt",PMTmodel) );
  if( !file.is_open() ){
    cout << Form("data/model%dQE.csv didn't open... exiting",PMTmodel) << endl;
    return;
  }

  vector<double> wavelength, efficiency;

  string str;
  size_t comma;
  getline(file,str); //Burn the header
  while( !file.eof() ){
    getline(file,str);
    if( str.length() ){
      comma = str.find_first_of(",");

      wavelength.push_back( atof( str.substr(0,comma).c_str() ) );
      efficiency.push_back( atof( str.substr( comma+1, str.length()-comma ).c_str() ) );
    }
  }

  TVectorD TvWL(wavelength.size(),&wavelength[0]);
  TVectorD TvQE(efficiency.size(),&efficiency[0]);
  TGraph* gQE = new TGraph(TvWL,TvQE);

  //// begin loop over events
  ////
  TRandom rnd;
  int nCut, nPhotons, totalCut=0, totalPhotons=0;
  double X,Z,wl,trans,rand;
  int nEvents = inTree->GetEntries();
  for (int ev = 0; ev < nEvents ; ev++){
    inTree->GetEntry(ev);

    nCut = 0;
    nPhotons = Xf->size();
    for (int k = 0; k < nPhotons; k++){	//begin loop over hits



      // If the wavelength isn't defined by the QE data QE=0,
      // so cut it and move on to the next photon
      // Else give it a chance to make it under the curve
      wl = 1240./(1e6*Energyf->at(k) );
      if(wl < wavelength.front() || wl > wavelength.back()){
        nCut++;
      } else{
        // Get transmission % from data
        trans = gQE->Eval( wl );
        // Give the photon a random chance to make it past the
        // Quantum Efficiency check
	rand = 100*rnd.Rndm();
        if( 100*rnd.Rndm() > trans ){
          nCut++;
        }
      }//end cuts






    }//end nPhotons loop
    totalCut += nCut;
    totalPhotons += nPhotons;
    cout << Form("%7d of %7d photons cut in event %3d",nCut,nPhotons,ev) << endl;
  }//end event looop
  cout << Form("%7d of %7d photons cut %3.0f%%",totalCut,totalPhotons,100.*(float)totalCut/(float)totalPhotons ) << endl;
  inFile->Close();
  delete inFile;
}//end PMTcuts


/* (experimental)
 * Generate a waveform as the sum of single photon responses
 * Pulse integral to amplitude conversion factor is a result of the following
 * integral_21^∞ A*((t - t0)/10)^3.4 exp(-1/10 (t - t0)) dt = A*101.3581
 * DRAW saves event display to .png
 */
void outputDigitization(string fileName, bool DRAW=false){
  std::cout << " running outputDigitization; input filename: " << fileName << std::endl;
  // Choose a path
  int method = 2;
  // method 1: (superposition) Create a TF1 for every photon using each photon's time and energy/wavelength
  // method 2: (per time bin) Create a TF1 for every time bin using the sum of that bin's photons with an averaged pulse amplitude
  bool isLTQuadStepFilt = 1;
  bool useCalibConstant = 0;
  bool useHighGain = 1;
  ////////////////////// Readout constants //////////////////////
  //  float sampleFrequency = 5.; // GHz
  float sampleFrequency = 0.32; // GHz
  int nSamples = 24;
  float timeBinWidth = 1./sampleFrequency; // ns
  float timeWindow = nSamples*timeBinWidth; // ns
  int nChannels = 16;
  double calibConstant;
  if (isLTQuadStepFilt) calibConstant = 0.00812013;
  else calibConstant = 0.0016127; // calibration constant
  ////////////////////// Noise constants //////////////////////
  //  float sigma = 1.0; //mV
  float sigma = 0.0; //mV
  float mean = 0.0; //mV
  ////////////////////// PMT constants //////////////////////
  // Will probably make a PMT object when appropriate
  // For R1635
  // Dark Current = 1nA typical, 50 nA max
  // gain 3e6 - 8e6
    float highGain = 400e6; // unitless
  //  float highGain = 2e6; // unitless
  //  float lowGain = 3e6; // unitless
  float lowGain = 1e6; // unitless

  float gain;
  if (useHighGain) gain = highGain;
  else gain = lowGain;
  //  float darkCurrent = 50e-6; // mA aka mC/s. Mfr claims this is 1nA typical, but up to 50nA
  float darkCurrent = 1e-6; // mA aka mC/s. Mfr claims this is 1nA typical, but up to 50nA
  
  //Get the integral of a pulse with unit amplitude
  //  float triggerDelay = 21.0; // original setting
  float triggerDelay = 0.;
  // how many ns to include before daq window -- meant to ensure falling edge of dark current shows up in daq window
  float bufferWindow = 80; 
  TF1 *testPulse;
  if (method == 1) testPulse = new TF1("testPulse", LTQuadStepFilt, 0.0, timeWindow, 2);
  else testPulse = new TF1("testPulse", PMTpulseFunction, 0.0, timeWindow, 2);
  testPulse->SetParameters(1.0,triggerDelay);
  float unitPulseIntegral = testPulse->Integral( 0.0, timeWindow, 1e-4 );
  TCanvas *c1 = new TCanvas("c1","c1",550,500);
  c1->cd(1);
  testPulse->Draw();
  //  c1->SaveAs(Form("testPulse_method%d.png",method));
  //  TGraph* g = getPMTdata( Form("model%dresponse.txt",2056) );
  TGraph* g = getPMTdata( Form("model%dresponse.txt",1635) );
  if(g == NULL) return;
  float avePhotonEnergy = (1240./g->GetMean(1))*1.6022e-19; // Joule
  float aveRadSens = g->GetMean(2); // mA/W aka mC/J
  //  float gain = 1e6; // unitless
  float outputImpedance = 50.0; // ohms
  float aveChargePerPulse = GetAverageCharge ( g, gain); // mC
  float aveNpulsePerNanoSec = 1e-9*darkCurrent/aveChargePerPulse; // 1/ns
  float aveNpulses = (timeWindow+bufferWindow)*aveNpulsePerNanoSec; // unitless
  float avePulseAmp = aveChargePerPulse*outputImpedance/(1e-9*timeBinWidth*unitPulseIntegral); // mV

  //  cout << "Average charge (ie DC amp) is  " << aveChargePerPulse << " avePulseAmp " << avePulseAmp << " impedance " << outputImpedance << " timebinWidth " << timeBinWidth << " unitPulseIntegral " << unitPulseIntegral << " aveNpulsePerNanoSec " << aveNpulsePerNanoSec << " aveNpulses " << aveNpulses << endl;
  // return;

  TH1D *h[16];
  TCanvas *c = new TCanvas("simWF","Simulated Waveforms",1200,700);
  c->Divide(4,4);

  // Input root file
  TFile* inFile = new TFile(fileName.c_str());
  if(inFile->IsZombie()){
    cout << fileName << " didn't open" << endl;
    return;
  }

  TTree* eventDatatree = (TTree*)inFile->Get("EventData");
  TTree* tree = (TTree*)inFile->Get("RPD1tree");

  vector<double> *Px=0, *Py=0, *Pz=0, *energy=0, *time=0, *LastStepZf = 0;
  vector<int> *rodNo=0, *timeHist=0, *nCherenkovs=0, *channel_rpd=0;
  vector<double> *Xf=0, *Yf=0, *Zf=0;
  TBranch *b_timeHist = 0;

  TH1D *h_time[nChannels];
  
  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("channel",&channel_rpd);
  tree->SetBranchAddress("time",&time);
  tree->SetBranchAddress("nCherenkovs",&nCherenkovs);
    
  eventDatatree->SetBranchAddress("gunPosX",&Xf);
  eventDatatree->SetBranchAddress("gunPosY",&Yf);
  eventDatatree->SetBranchAddress("gunPosZ",&Zf);
  //  eventDatatree->SetBranchAddress("lastStepZ",&LastStepZf);

  // Output root file
  string outputName = fileName;
  if(outputName.find("/") != string::npos) outputName.erase( 0, outputName.find_last_of("/") + 1 );
  if(outputName.find(".") != string::npos) outputName.erase( outputName.find_last_of(".") );
  std::cout << " creating output File: " << Form("WF_%s.root", outputName.c_str()) << std::endl;

  TFile* outFile = new TFile( Form("WF_%s.root", outputName.c_str() ), "RECREATE");
  TTree* outTree = new TTree( "tree", "tree" );
  int eventNo;
  //  outTree->Branch("EventNo",&eventNo);
  outTree->Branch("gunPosX",&Xf);
  outTree->Branch("gunPosY",&Yf);
  outTree->Branch("gunPosZ",&Zf);
  //  outTree->Branch("lastStepZ",&LastStepZf);

  vector< vector< float >* > waveforms(16,0);
  vector< vector< float >* > heights(16,0);
  vector< vector< float >* > charges(16,0);
  for(int tile = 0; tile < nChannels; tile++){
    h_time[tile] = new TH1D(Form("h_time_ch%d",tile),Form("h_time_ch%d",tile), nSamples, 0, timeWindow);
    waveforms[tile] = new vector< float >;
    outTree->Branch( Form("RawSignal%d",tile), &waveforms[tile] );
    outTree->Branch( Form("pulseHeight%d",tile), &heights[tile] );
    outTree->Branch( Form("pulseCharge%d",tile), &charges[tile] );
    if(DRAW) h[tile] = new TH1D( Form("tile%dwaveform",tile), Form("Tile %d Simulated Waveform;time (ns);Amplitude (mV)",tile), nSamples, 0, timeWindow);
  }

    int nEntries = tree->GetEntries();
    for(int eventNo = 0; eventNo < nEntries; eventNo++){
      eventDatatree->GetEntry(eventNo);
      tree->GetEntry(eventNo);

      int nHits = energy->size();
      for(int hit = 0; hit < nHits; hit++){

      }
    }
  
  ////////////////////// The hard way //////////////////////
  //
  //
  //
  // //////////////////// The hard way //////////////////////
  if(method == 1){

    double wavelength, pulseIntegral, amplitude, photonEnergy;
    vector< vector< TF1* > > pulses;
    pulses.resize(16);

    TRandom2 rnd;
    rnd.SetSeed(gSystem->Now());

    int nEntries = tree->GetEntries();
    for(int eventNo = 0; eventNo < nEntries; eventNo++){
      eventDatatree->GetEntry(eventNo);
      tree->GetEntry(eventNo);

      int nHits = energy->size();
      for(int hit = 0; hit < nHits; hit++){

	if(hit%500 == 0) cout << "\r" << std::left << Form("Event %d, Hit %d",eventNo,hit) << flush;

        // Get the photon energy and determine the pulse amplitude based on that energy
	//        photonEnergy = 1e6*energy->at(hit); //eV
        photonEnergy = energy->at(hit); //eV
	wavelength = 1240./photonEnergy; // nm
	// what unit should the pulse integral be? Charge?
	pulseIntegral = g->Eval( wavelength )*gain*photonEnergy*1.6022e-19; // mC
	// pulse integral / unit pulse integral gives relative charge -- dividing Q by t gives you current. V = IR
       	amplitude = outputImpedance*pulseIntegral/(1e-9*timeBinWidth*unitPulseIntegral); // mV
	//	amplitude = outputImpedance*pulseIntegral/(1e-9*timeBinWidth); // mV

        charges[channel_rpd->at(hit)]->push_back(pulseIntegral);
        heights[channel_rpd->at(hit)]->push_back(amplitude);


	if (isLTQuadStepFilt) {
	  pulses[channel_rpd->at(hit)].push_back(new TF1( Form("ev%d",eventNo), LTQuadStepFilt, -1*bufferWindow, timeWindow, 2) );
	  pulses[channel_rpd->at(hit)].back()->SetParameters(amplitude, time->at(hit) + triggerDelay);
	}
	else {
	  pulses[channel_rpd->at(hit)].push_back( new TF1( Form("ev%d",eventNo), PMTpulseFunction, -1*bufferWindow, timeWindow, 2) );
	  pulses[channel_rpd->at(hit)].back()->SetParameters( amplitude, time->at(hit) + triggerDelay);
	}
       

      }//end hit loop
    
      for(int tile = 0; tile < 16; tile++){

        // Add dark current pulses
        int nPulses = rnd.Poisson( aveNpulses );
	std::cout << "tile " << tile << " n dark current pulses " << nPulses << " aveNpulses " << aveNpulses << " avePulseAmp " << avePulseAmp << std::endl;
      	for(int i = 0; i < nPulses; i++){
	  // how should we get this width?
	  //	  double amp_darkcurrent = rnd.Gaus( avePulseAmp, avePulseAmp/2);
	  // if dark current is just random thermal excitations in photocathode, each exitation should see the full gain and should be one full unit in amplitude, ie) no variance in amplitude
	  double amp_darkcurrent = avePulseAmp;
	  //	  std::cout << " avePulseAmp "<< avePulseAmp  << " avePulseAmp/2 " << avePulseAmp/2 << " amp_darkcurrent " << amp_darkcurrent << std::endl;
	  double time_darkcurrent = rnd.Uniform(0,timeWindow);
	  pulses[tile].push_back(new TF1( Form("darkCurrent%d",i), PMTpulseFunction, -1*bufferWindow, timeWindow, 2) );
	  pulses[tile].back()->SetParameters( amp_darkcurrent, time_darkcurrent );
	  if (tile == 6) std::cout << " pulse " << i << " dc amp " << amp_darkcurrent << " time " << time_darkcurrent << " timeWindow " << timeWindow << " avePulseAmp " << avePulseAmp << std::endl;
      	}

        // Evaluate the TF1s for each time bin + gaussian noise
	//	for(int bin = 0; bin < 1024; bin ++){
	for(int bin = 0; bin < nSamples; bin ++){
	  waveforms[tile]->push_back(evaluateSum( pulses[tile], timeBinWidth*bin ) + rnd.Gaus(mean,sigma) );
	  //	  waveforms[tile]->push_back(evaluateSum( pulses[tile], timeBinWidth*bin ) );
	}//end waveform loop

	if(DRAW){
	  //	  for(int bin = 0; bin < 1024; bin ++){
	  for(int bin = 0; bin < nSamples; bin ++){
	    h[tile]->SetBinContent(bin, waveforms[tile]->at(bin));
	  }//end waveform loop
	  c->cd(tile+1);
	  h[tile]->Draw();
	  //					h[tile]->SetAxisRange(-10,750,"Y");
	  // h[tile]->SetAxisRange(-10,100000,"Y");
	}
      }//end tile loop
      if(DRAW) c->Print( Form("event%d_method%d_QuadStepFilt%d_noiseRMS%dmV_calib%d.png",eventNo,method,isLTQuadStepFilt,(int) sigma, useCalibConstant) );

      outTree->Fill();

      //Delete all TF1 objects, clear vectors, reset histograms
      for(int tile = 0; tile < 16; tile++){
	waveforms[tile]->clear();
	int nHits = pulses[tile].size();
	if(nHits == 0) continue;
	for(int hit = 0; hit < nHits; hit++){
	  delete pulses[tile][hit];
	}//end nHits loop
	pulses[tile].clear();
        if(DRAW)h[tile]->Reset();
      }//end tile loop

    }//end event loop
    cout << endl;



    ////////////////////// The middle way //////////////////////
    //
    //
    //
    // //////////////////// The middle way //////////////////////

  }else if( method == 2 ){
    vector< vector< TF1* > > pulses;
    pulses.resize(16);

    TRandom2 rnd;
    rnd.SetSeed(gSystem->Now());

    int nEntries = tree->GetEntries();
    for(int eventNo = 0; eventNo < nEntries; eventNo++){
      if(eventNo%50 == 0) cout << "\r" << std::left << Form("Event %d",eventNo) << flush;
      tree->GetEntry(eventNo);
      eventDatatree->GetEntry(eventNo);

      int nHits = energy->size();
      for(int hit = 0; hit < nHits; hit++){
	h_time[channel_rpd->at(hit)]->Fill(time->at(hit));
      }
      for(int channel = 0; channel < nChannels; channel++){       
	
        for(int bin = 0; bin < nSamples; bin++){
	  int nPhotons = h_time[channel]->GetBinContent(bin+1);
	  double timeBinVal = h_time[channel]->GetBinCenter(bin+1);
	  double amp;
	  if (useCalibConstant) amp = calibConstant;
	  else amp = avePulseAmp;
	  if (isLTQuadStepFilt) {
	    pulses[channel].push_back(new TF1( Form("ev%d",eventNo), LTQuadStepFilt, -1*bufferWindow, timeWindow, 2) );
	    pulses[channel].back()->SetParameters(nPhotons*amp, timeBinVal + triggerDelay);
	  }
	  else {
	    pulses[channel].push_back( new TF1( Form("ev%d",eventNo), PMTpulseFunction, -1*bufferWindow, timeWindow, 2) );
	    //	    std::cout <<" channel " << channel << " timeBinVal " << timeBinVal << " triggerDelay " << triggerDelay <<" nPhotons " << nPhotons <<  std::endl;
	    pulses[channel].back()->SetParameters( nPhotons*amp, timeBinVal + triggerDelay);
	  }
			       
        }//end bin loop
	//	std::cout << " calibConstant " << calibConstant << " avePulseAmp " << avePulseAmp << std::endl; 
        // Add dark current pulses
        int nPulses = rnd.Poisson( aveNpulses );
      	for(int i = 0; i < nPulses; i++){
	  pulses[channel].push_back(new TF1( Form("darkCurrent%d",i), PMTpulseFunction, -1*bufferWindow, timeWindow, 2) );
	  pulses[channel].back()->SetParameters( rnd.Gaus( avePulseAmp, 2*avePulseAmp), rnd.Uniform(-1*bufferWindow,timeWindow) );
      	}

        // Evaluate the TF1s for each time bin + gaussian noise
	//	for(int bin = 0; bin < 1024; bin ++){
	for(int bin = 0; bin < nSamples; bin ++){
	  waveforms[channel]->push_back(evaluateSum( pulses[channel], timeBinWidth*bin ) + rnd.Gaus(mean,sigma) );
	}//end waveform loop
	
	if(DRAW){
	  //	  for(int bin = 0; bin < 1024; bin ++){
	  for(int bin = 0; bin < nSamples; bin ++){
	    h[channel]->SetBinContent(bin, waveforms[channel]->at(bin));
	  }//end waveform loop
	  c->cd(channel+1);
	  h[channel]->Draw();
	  //	  h[channel]->SetAxisRange(-10,500,"Y");
	}
      }//end channel loop
      if(DRAW) c->Print( Form("event%d_method%d_QuadStepFilt%d_noiseRMS%dmV_calib%d_highGain%d.png",eventNo,method,isLTQuadStepFilt,(int) sigma, useCalibConstant, useHighGain) );

      outTree->Fill();
      outFile->Write();
      
      //Delete all TF1 objects, clear vectors, reset histograms
      for(int channel = 0; channel < 16; channel++){
	waveforms[channel]->clear();
	int nHits = pulses[channel].size();
	if(nHits == 0) continue;
	for(int hit = 0; hit < nHits; hit++){
	  delete pulses[channel][hit];
	}//end nHits loop
	pulses[channel].clear();
        if(DRAW)h[channel]->Reset();
      }//end channel loop
    }//end event loop
    cout << endl;
  

    // Not implemented yet

  }//end the easy way
  //  outFile->Write();
  outFile->Close();
  delete outFile; //Deletes histograms
  inFile->Close();
  delete inFile;
  for(int tile = 0; tile < 16; tile++){
    delete waveforms[tile];
  }
  if(!DRAW) delete c;
}


double GetAverageCharge( TGraph* PMTcathodeRadSens, double gain ){
  TFile* f = new TFile("Energy_dist.root","read");
  TH1D* energyHist = (TH1D*)f->Get("Energy dist");

  //Normalize the histogram so it represents a single photon
  energyHist->Scale( 1.0/energyHist->GetEntries() );
  double avgCharge = 0.0;
  double photonEnergy;
  double evToJ = 1.6022e-19;
  for(int bin = 0; bin < energyHist->GetNbinsX(); bin++){

    photonEnergy = 1e6*energyHist->GetBinCenter(bin); //eV

    // current is fraction of photons in this bin * the single photon current at this energy/wavelength
    // (mA/W) (J) -> mC
    double charge = energyHist->GetBinContent(bin) * PMTcathodeRadSens->Eval( 1240./photonEnergy )*gain*photonEnergy*evToJ; // mC 
    avgCharge += charge;
    //    std::cout << "bin " << bin << " totalBins " << energyHist->GetNbinsX() << "  bincontent " << energyHist->GetBinContent(bin) << " cathode Rad sensitivity: " << PMTcathodeRadSens->Eval( 1240./photonEnergy ) << " gain " << gain << " photon energy " << photonEnergy << " avgCharge " << avgCharge << " charge " << charge << " unnormalized " << PMTcathodeRadSens->Eval( 1240./photonEnergy )*gain*photonEnergy*evToJ << std::endl;
  }
  f->Close();
  delete f;
  return avgCharge;
}




TGraph* getPMTdata(string fileName){
  vector<double> x, y;
  ifstream file( Form("%s", fileName.c_str() ) );

  if( !file.is_open() ){
    file.open( fileName.c_str() );
    if( !file.is_open() ){
      cout << Form( "%s didn't open... exiting", fileName.c_str() ) << endl;
      return NULL;
    }
  }
  string str;
  size_t comma;
  getline(file,str); //Burn the header
  while( !file.eof() ){
    getline(file,str);
    if( str.length() ){

      comma = str.find_first_of(",");

      x.push_back( atof( str.substr(0,comma).c_str() ) );
      y.push_back( atof( str.substr( comma+1, str.length()-comma ).c_str() ) );
    }
  }
  TVectorD TvX(x.size(),&x[0]);
  TVectorD TvY(y.size(),&y[0]);
  TGraph* g = new TGraph(TvX,TvY);
  return g;
}

/*
 * Generate timing data to be used by JZCaPA
 * Argument is digitization speed in GHz
 */
void generateTimingData( float digiSpeed = 0.32 ){
  ofstream outputFile( Form("output/MC_%3.2fGHz.txt", digiSpeed ) );

  outputFile << Form("%3.2fGHz simulated timing data", digiSpeed) << endl;
  outputFile << " Module1    Module2    Module3    Module4    Module5" << endl;

  float binWidth = 1./digiSpeed;
  for(int i = 0; i < 1024; i++){
  //  for(int i = 0; i < nSamples; i++){
    outputFile << Form("%8.3f,  %8.3f,  %8.3f,  %8.3f,  %8.3f", i*binWidth, i*binWidth, i*binWidth, i*binWidth, i*binWidth) << endl;
  }
}


double evaluateSum( vector< TF1* >& funcVec, double x ){
  double value = 0.0;
  for(int i = 0; i < funcVec.size(); i++){
    value += funcVec.at(i)->Eval(x);
  }
  return value;
}

bool checkFileExistencein(const string& filename){
  ifstream f(filename.c_str());
  return f.is_open();
}


void getFile(string filename,  /*out*/ ifstream& file){

  const bool file_exists = checkFileExistencein(Form("data/%s", filename.c_str() ) );
  if (!file_exists) {
    cout << "File in " << Form("data/%s", filename.c_str() ) << " not found." << endl;
    cout << "Search " << filename << endl;
    bool file_exists = checkFileExistencein(filename);
    if (file_exists) {
      file.open(filename.c_str());
      cout <<  filename  << " is found in the same level" << endl;
    }

  }
}

double LTQuadStepFiltKernel(double x, double a, double b, double c, double d)
{
  return b*b*b*exp(-(x-a)/b)/(b-c)/(b-d) - c*c*c*exp(-(x-a)/c)/(b-c)/(c-d) - d*d*d*exp(-(x-a)/d)/(b-d)/(d-c) - (a+b+c+d) + x;
}

double LTQuadStepFilt(double *t, double *amp)
{
  double m_t0 = 0; //response time of ppm (ie nominal start time of pulse)
  double m_tau0 = 4; //response time of ppm (ie nominal start time of pulse)
  double m_tau1 = 0.5; // rise time of the nominal linear response
  double m_tau2 = 11; // fall time of the nominal linear response
  double m_tauFilt = 15; // effective time constant needed to reproduce the PPM waveform
  double deltat = t[0] - m_t0;
  double sum = 0;

  if (deltat > 0) sum += LTQuadStepFiltKernel(deltat, 0, m_tau1, m_tau2, m_tauFilt);
  if (deltat > m_tau0) sum -= 2*LTQuadStepFiltKernel(deltat, m_tau0, m_tau1, m_tau2, m_tauFilt);
  if (deltat > 2*m_tau0) sum += LTQuadStepFiltKernel(deltat, 2*m_tau0, m_tau1, m_tau2, m_tauFilt);
  return sum*amp[0];
}
