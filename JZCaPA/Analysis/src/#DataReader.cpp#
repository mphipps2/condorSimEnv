 /** @defgroup ana Analysis
 *   @ingroup ana
 *   @file DataReader.cpp
 *   @brief Implementation of DataReader.
 *
 *  Function definitions for DataReader are provided.
 *  This class reads a rootfile with raw waveforms
 *  that are processed by rcdaqAnalysis running on prdf files.
 *  Then, in the event loop, analysis classes can be called.
 *
 *  @author Yakov Kulinich, Riccardo Longo
 *  @bug No known bugs.
 */

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <fstream>

#include "DataReader.h"
#include "Analysis.h"
#include "Containers.h"
#include "RPD.h"
#include "ZDC.h"
#include "Visualizer.h"


/** @brief Default Constructor for DataReader.
 */
DataReader::DataReader() : DataReader( 0, 0, "", 0 ){

}

/** @brief Constructor for DataReader.
 *
 *  @param nCh Number of channels being read
 *  @param nSamp Number of samples per channel
 */
DataReader::DataReader( const unsigned int nCh, const unsigned int nSamp )
  : DataReader( nCh, nSamp, "", 0 ){

  // here say we will read in list of files, maybe
  // there is better way to do it, just an idea for now.
  m_readListOfFiles = true;
}

/** @brief Constructor for DataReader.
 *
 *  @param nCh Number of channels being read
 *  @param nSamp Number of samples per channel
 *  @param fNameIn Input filename.
 */
DataReader::DataReader( const uint nCh, const uint nSamp,
			const std::string& fNameIn )
  : DataReader( nCh, nSamp, fNameIn, 0 ){

}


/** @brief Constructor for DataReader.
 *
 *  @param nCh Number of channels being read
 *  @param nSamp Number of samples per channel
 *  @param4 fNameIn Output file name (custom)
 *  @param3 runNum Run number being used.

 */
DataReader::DataReader( const uint nCh, const uint nSamp,
			const std::string& fNameIn, const uint runNum )
  : m_nCh( nCh ), m_nSamp( nSamp ), m_fNameIn( fNameIn ),
    m_runNumber( runNum ), m_readListOfFiles( false ), m_fIn( NULL ){

}


/** @brief Destructor for DataReader.
 */
DataReader::~DataReader(){


  for( auto& time : m_time ){
      delete time; time = NULL;
  }
  for( auto& det : m_detectors ){
      delete det; det = NULL;
  }
  for( auto& ana : m_ana ){
    delete ana; ana = NULL;
  }
  for( auto& det_ana : m_det_ana ){
    delete det_ana; det_ana = NULL;
  }

  //This line deletes all histograms so we don't have to
  gDirectory->GetList()->Delete();

  delete m_alignment; m_alignment = NULL;
  delete m_fOut; m_fOut = NULL;
  delete m_fIn; m_fIn = NULL;

}

/** @brief Adds an analysis to vector of analysis to be executed before the detector analysis (e.g. WaveForm Analysis)
 *
 *  @param ana Pointer to an Analysis.
 *
 *  @return none
 */
void DataReader::AddPreAnalysis( Analysis* ana ){
  ana->SetVerbosity( m_verbose );
  m_ana.push_back( ana );
}

/** @brief Adds an analysis to vector of detector analysis (e.g. ZDC, RPD or combined)
 *
 *  @param ana Pointer to an Analysis.
 *
 *  @return none
 */
void DataReader::AddDetectorAnalysis( Analysis* det_ana ){
  det_ana->SetVerbosity( m_verbose );
  m_det_ana.push_back( det_ana );
}

/**
 * @brief Allows the user to select which detectors will enter the waveform analysis and the detector analysis s
 * @param _useZDC1 set it true to look to ZDC1
 * @param _useZDC2 set it true to look to ZDC2
 * @param _useRPD  set it true to look to the RPD
 */
void DataReader::SelectDetectorForAnalysis(bool _useZDC1, bool _useZDC2, bool _useRPD){

    useZDC1 = _useZDC1;
    useZDC2 = _useZDC2;
    useRPD = _useRPD;

    std::string enableDet;
    enableDet = "Detectors Selected for Analysis: ";

    if (useZDC1) enableDet += " ZDC1 " ;
    if (useZDC2) enableDet += " ZDC2 " ;
    if (useRPD) enableDet += " RPD " ;

    enableDet += ";";
    std::cout << enableDet << std::endl;
}

/** @brief Enables the read from list of files option for DataReader
 *
 *  @param listname name of the list of files to be processed (with the full path if it's not in the execution folder)
 *
 *  @return none
 */
void DataReader::ReadListOfFiles( std::string listname ){

    m_readListOfFiles = true;
    m_fListOfFiles = listname;
}


/**
 * @brief Reads the .xml configuration file and load characteristics for all the channels, immediately sorted into detectors objects
 * @param _inFile
 */
void DataReader::LoadAlignmentFile(std::string _inFile ){

    m_XMLparser = new XMLSettingsReader();

    if (!m_XMLparser->parseFile(_inFile)) {
            std::cerr << " Data Reader could not parse file : " << _inFile << std::endl;
            return;
    }

    m_alignment = new Alignment();
    m_alignment->runNumber = m_runNumber;

    std::cout << "Loading .xml Alignment File..." << std::endl;
    std::cout << "Found " << m_XMLparser->getBaseNodeCount("Alignment") << " alignment entries " << std::endl;
    std::cout << "Retrieving the information for run " << m_runNumber << std::endl;

    int run;
    for (unsigned int i = 0; i < m_XMLparser->getBaseNodeCount("Alignment"); i++) {
        m_XMLparser->getChildValue("Alignment",i,"run",run);
        if(run != m_runNumber) continue;
        std::cout << "Found Run Entry in Alignment file for run " << m_runNumber << std::endl;
        m_XMLparser->getChildValue("Alignment",i,"x_table",m_alignment->x_table);
        m_XMLparser->getChildValue("Alignment",i,"y_table",m_alignment->y_table);
        m_XMLparser->getChildValue("Alignment",i,"upstream_Det",m_alignment->upstream_Det);
        m_XMLparser->getChildValue("Alignment",i,"mid_Det",m_alignment->mid_Det);
        m_XMLparser->getChildValue("Alignment",i,"downstream_Det",m_alignment->downstream_Det);
        m_XMLparser->getChildValue("Alignment",i,"target_In",m_alignment->target_In);
        m_XMLparser->getChildValue("Alignment",i,"lead_In",m_alignment->lead_In);
        m_XMLparser->getChildValue("Alignment",i,"magnet_On",m_alignment->magnet_On);
    }

    if(m_alignment == NULL) std::cout << "WARNING: ALIGNMENT NOT FOUND!!!" << std::endl;
    delete m_XMLparser; m_XMLparser = NULL;
    return;

}

/**
 * @brief Reads the .xml configuration file and load characteristics for all the channels, immediately sorted into detectors objects
 * @param _inFile
 */
void DataReader::LoadConfigurationFile(std::string _inFile ){

    m_XMLparser = new XMLSettingsReader();

    if (!m_XMLparser->parseFile(_inFile)) {
            std::cerr << " Data Reader could not parse file : " << _inFile << std::endl;
            return;
    }

    std::cout << "Loading .xml Configuration File..." << std::endl;
    std::cout << "Found " << m_XMLparser->getBaseNodeCount("channel") << " channel entries " << std::endl;

    std::vector < Channel* > channelEntries;
    int first_run, last_run;

    for (unsigned int i = 0; i < m_XMLparser->getBaseNodeCount("channel"); i++) {
        Channel *buffer_ch = new Channel();
        m_XMLparser->getChildValue("channel",i,"start_run",first_run);
        m_XMLparser->getChildValue("channel",i,"end_run",last_run);

        //Discard entries for any channel that does not apply to our run
        if(m_runNumber < first_run || m_runNumber > last_run) continue;

        //If the entry applies, we store it in the vector
        m_XMLparser->getChildValue("channel",i,"detector",buffer_ch->detector);
        m_XMLparser->getChildValue("channel",i,"name",buffer_ch->name);
        m_XMLparser->getChildValue("channel",i,"mapping_row",buffer_ch->mapping_row);
        m_XMLparser->getChildValue("channel",i,"mapping_column",buffer_ch->mapping_column);
        m_XMLparser->getChildValue("channel",i,"delay",buffer_ch->delay);
        m_XMLparser->getChildValue("channel",i,"offset",buffer_ch->offset);
        m_XMLparser->getChildValue("channel",i,"HV",buffer_ch->HV);
        m_XMLparser->getChildValue("channel",i,"is_on",buffer_ch->is_on);
        m_XMLparser->getChildValue("channel",i,"Vop",buffer_ch->Vop);

        bool isNew(true);
        for( int k = 0; k < channelEntries.size(); k++){
            if(buffer_ch->name == channelEntries.at(k)->name){
                std::cout << "WARNING!!! Redundancy in your settings file for " << buffer_ch->name << ". Check it carefully. The second entry found will be skipped..." << std::endl;
                isNew = false;
            }
        }
        if(isNew) channelEntries.push_back(buffer_ch);
    }

    std::cout << "Loaded " << channelEntries.size() << " configuration entries " << std::endl;
    if( channelEntries.size() < 18 ) std::cout << "WARNING!!!! Number of Channels < 18. Seems that some entry is missed for this run in the config.xml. BE CAREFUL!" << std::endl;
    ZDC* zdc1 = new ZDC(channelEntries,m_runNumber,1);
    ZDC* zdc2 = new ZDC(channelEntries,m_runNumber,2);
    RPD* rpd = new RPD(channelEntries,m_runNumber);

    m_detectors.push_back(zdc1); //Position 0 goes for ZDC1
    m_detectors.push_back(zdc2); //Position 1 goes for ZDC2
    m_detectors.push_back(rpd); //Position 2 goes for the RPD

    std::cout << "Detector configuration: loading complete! " << std::endl;
    delete m_XMLparser; m_XMLparser = NULL;
    return;
}


/** @brief Loads calibrated timing information for DRS4 modules based on run number
 *  @param _inFile Optional argument for loading a custom file
 *
 *  Uses run number to determine scan number and loads the appropriate DRS4 timing
 *  information from the 2018 Testbeam. After loading, we hand each Channel a
 *  pointer to its timing vector. Must be run after LoadConfigurationFile()
 *
 */
void DataReader::LoadTimingFile(std::string _inFile){

    std::string inFile,buffer;
    char data[12];
    float dat;
    std::vector< Channel* > vChannel;
    m_time.resize(5);
    for(int i = 0; i < 5; i++){
        m_time[i] = new std::vector< float >;
    }

    int chNo;
    int start[] = {79,  152, 190, 202, 215, 258 };
    int stop[]  = {112, 171, 200, 213, 231, 413 };
    int scanNum = 0;

    //If the user has selected a file, use it and skip determining the scan number
    if(_inFile != ""){
        std::cout << "Using timing file definied by user " << _inFile << std::endl;
        inFile = _inFile;
        goto open;
    }

    //Determine scan number using ranges defined by start and stop.
    //Scan 6-13 have identical timing data, so we treat them all as scan 6
    for(int i = 0; i < 6; i++){
        if(start[i]<=m_runNumber && stop[i]>=m_runNumber){
            scanNum = i+1;
        }
    }
    if(scanNum == 0){
        std::cerr << "Scan number not defined, defaulting timing data to Scan 1" << std::endl;
        scanNum = 1;
    }
    inFile = (std::string)std::getenv("JZCaPA") + Form("/Utils/Timing_data/2018scan%d.txt",scanNum);

    open:
    std::ifstream file( inFile.c_str() );
    if( !file.is_open() ){
        std::cerr << "WARNING: Timing data file didn't open " << inFile << std::endl;
        return;
    }

    //Get the headers out
    getline(file,buffer);
    getline(file,buffer);

    while(!file.eof()){
        getline(file,buffer);
        if(buffer == "") break;

        //Loop through the DRS4 modules pushing the data into the vector
        //and chopping off the front of the string
        for(int drs = 0; drs<4; drs++){
            buffer.copy( data , 11 );
            dat = atof(data);
            m_time[drs]->push_back( dat );
            buffer.erase(buffer.front(),11);
        }
        //The spacing on the file is a bit weird, so fill the last one differently
        m_time[4]->push_back( atof( buffer.c_str() ) );
    }
    file.close();

    //Loop through the detector vector and their Channel vectors assigning the
    //time vector based on hardware channel number (Channel::name)
    for( auto& det : m_detectors ){
        vChannel = det->GetChannelsVector();
        for( Channel* ch : vChannel ){
            buffer = ch->name;
            //Remove "Signal" from name and convert to int
            buffer.erase(0,6);
            chNo = atoi( buffer.c_str() );

            //If the channel already has a time vector, print an error and continue
            if(ch->pTimeVec != 0){
                std::cerr << "WARNING: Overwriting Channel time vector" << std::endl;
            }
            ch->pTimeVec = m_time[chNo/4];
        }
    }
    std::cout << "DRS4 Timing: loading complete " << std::endl;
}

/**
 * @brief DataReader::GetDetector allows the user to access the detectors objects after loading them at the beginning of the execution
 * @param _detName can be ZDC1 (upstream module), ZDC2 (downstream module) or RPD.
 * @return
 */
Detector* DataReader::GetDetector( std::string _detName ){

    if(_detName != "ZDC1" && _detName != "ZDC2" && _detName != "zdc1" && _detName != "zdc2" && _detName != "RPD" && _detName != "rpd")
    {
        std::cout << "The detector you're looking for is not ZDC1, ZDC2 or RPD. Please check and correct your request " << std::endl;
        return NULL;
    }
    if( _detName == "ZDC1" || _detName == "zdc1" ) return m_detectors.at(0);
    if( _detName == "ZDC2" || _detName == "zdc2" ) return m_detectors.at(1);
    if( _detName == "RPD" || _detName == "rpd" ) return m_detectors.at(2);

    std::cout << "WARNING: detector recognition glitch. NULL pointer being returned..." << std::endl;
    return NULL;
}


/** @brief Console output update
 *  @return none
 *
 *  Called by the TTimer in ProcessEvents
 *  Updates the console with CPU and RAM usage,
 *  number of events processed and events/second
 *  if verbosity is not zero
 */
void DataReader::UpdateConsole( Long_t _updateRate){

    if( m_verbose == 0 ){ return; }
    if(m_event!=0){


        // Get CPU information
        gSystem->GetCpuInfo(&cpuInfo, 100);
        gSystem->GetProcInfo( &procInfo );
        // Get Memory information
        gSystem->GetMemInfo(&memInfo);
        // Get events/second
        double rate = 1000*(m_event-m_event_old)/_updateRate;
        m_event_old = m_event;

        std::cout << "\r" << std::left <<  Form("Processed %5d events, ", m_event);
        std::cout << Form( "%5.1f ev/s, ", rate);
        std::cout << Form( "CPU use/time: %3d%%/%6.1fs, ", (int)cpuInfo.fTotal, (double)procInfo.fCpuSys + procInfo.fCpuUser);
        std::cout << Form( "RAM:%4.1f/%4.1fGB, ", (double)memInfo.fMemUsed/1024, (double)memInfo.fMemTotal/1024);
        std::cout << Form( "RAM used by process: %ldMB   ", (procInfo.fMemResident + procInfo.fMemVirtual)/1024 );
        std::cout << std::flush;
    }
}

/** @brief Run method for DataReader
 *  @return none
 *
 *  Call Initialize, ProcessEvents, and Finalize
 *  Made so the driver class only has to call one method.
 *
 */
void DataReader::Run(){

  Initialize();
  if( !m_fIn ){
      std::cerr << "Input file didn't open... exiting" << std::endl;
      return;
  }
  ProcessEvents();
  Finalize();
}

/** @brief Initialization method for DataReader
 *
 *  Select which file(s) to read. For now just a single
 *  file, later this can be extended to read many and make
 *  chain of files for example. Also create and initialize
 *  the analysis that will be running.
 *
 *  @return none
 */
void DataReader::Initialize(){

  if( m_readListOfFiles ){

    // Riccardo - 21/01/2019 - TChain implementation
    // The file list name is supposed to be already set
    m_fileChain = new TChain("tree");
    std::ifstream inFile;
    inFile.open(m_fListOfFiles.c_str());
    std::string reader_buff;
    while(inFile >> reader_buff){
        //Let's push all the files in the list into the TChain
        m_fileChain->Add(reader_buff.c_str());
    }
    /** TODO - Add fileChain reading below
     */
  } else {
    m_fIn = TFile::Open( m_fNameIn.c_str() );
  }

  // If no output directory is specified by the user
  // one is created for this run in $JZCaPA/results
  if( m_outputDir == "" ){
      m_outputDir = (std::string)std::getenv("JZCaPA") + Form("/results/run%d/",m_runNumber);
      gSystem->Exec( ("mkdir -p " + m_outputDir).c_str() );
  }

  // If we are reading a list of files, or have no run number
  // make default name output.root, otherwise make it
  // outputN.root, where N is a run number of a file.
  std::string fNameOut = m_readListOfFiles ?
    (m_outputDir + "output.root") : Form( (m_outputDir + "output%d.root").c_str(), m_runNumber );

  m_fOut = new TFile( fNameOut.c_str(), "RECREATE" );
  m_tOut = new TTree("AnalysisTree","AnalysisTree");
  m_tOut->Branch("runNumber", &m_runNumber, "m_runNumber/i");
  m_tOut->Branch("evNo", &m_event, "m_event/I");
  m_tOut->Branch("x_table", &m_alignment->x_table, "x_table/D");
  m_tOut->Branch("y_table", &m_alignment->y_table, "y_table/D");

  for( uint detID = 0; detID < (int) m_detectors.size(); detID++ ){
     m_detectors.at(detID)->SetAlignment(m_alignment);
  }
  for( auto& ana : m_ana ){
    ana->Initialize();
    ana->SetupHistograms();
  }
  for( auto& det_ana : m_det_ana ){
    det_ana->Initialize( m_detectors );
    det_ana->SetupHistograms();
    det_ana->SetBranches( m_tOut );
  }
}



/** @brief Process Events method for DataReader
 *
 *  Connect to files / chain of files.
 *  Read tree event by event, and fill basic waveforms.
 *  So far (12.12.18) For each event - Read Raw data for all channels,
 *  put it into 2D vector size NxM where N = nCh, M = nSamp
 *  Also, For each event - Fill N 1D histos that have size M.
 *  Then one can send these to any Analysis function they make.
 *  Here for example (12.12.18), we send to WFAnalysis::AnalyzeEvent( .. )
 *  See below - Have fun!
 *
 *  @return none
 */
void DataReader::ProcessEvents(){


  /** TODO : add reading for list of files
    * Please note that many of the implementations are now for a single-file treatment
    */
  TTree* tree = static_cast< TTree* >( m_fIn->Get( "tree" ) );

  //Specific pointers to each detector, if needed afterwards
  ZDC* zdc1 = static_cast< ZDC* >( GetDetector("ZDC1") );
  ZDC* zdc2 = static_cast< ZDC* >( GetDetector("ZDC2") );
  RPD* rpd =  static_cast< RPD* >( GetDetector("RPD") );

  //All the raw channels addresses set for read-out
  for( uint detID = 0; detID < (int) m_detectors.size(); detID++ ){
       m_detectors.at(detID)->SetBranches(tree);
       m_detectors.at(detID)->SetNSamples(m_nSamp);
       m_detectors.at(detID)->DeclareHistograms();
  }

  std::cout << "File: " << m_fIn->GetName() << " has " << tree->GetEntries() << " events." << std::endl;

  // !! EVENT LOOP
  for( int ev = 0; ev < tree->GetEntries(); ev++ ){
    m_event = ev;
    tree->GetEntry( ev );

    // Fill the raw waveforms
    for( uint detID = 0; detID < (int) m_detectors.size(); detID++ )
        m_detectors.at(detID)->FillHistograms();

    // Now call all analysis and run their AnalyzeEvent.
    // Can either send a vector of histos, or a 2D vector
    // Of all the data, depending on what you want to do.
    // Note that at the moment none of this methods is doing anything
    for( auto& ana : m_ana ){
      //raw data analysis
        if(useZDC1) ana->AnalyzeEvent( zdc1->GetChannelsVector() );
        if(useZDC2) ana->AnalyzeEvent( zdc2->GetChannelsVector() );
        if(useRPD) ana->AnalyzeEvent( rpd->GetChannelsVector() );
        }
    for( auto& det_ana : m_det_ana ){
      //Detector level analysis
      det_ana->AnalyzeEvent(  );
    }

  m_tOut->Fill();

  } // End event loop
}

/** @brief Finalize method for DataReader
 *
 *  Close any input files
 *  Call analysis finalize methods
 *
 *  @return none
 */
void DataReader::Finalize(){

  if( m_fIn ){
    m_fIn->Close();
  }

  // enter the output file since
  // we will be writing to it now.
  m_fOut->cd();

  //Add label, if wanted
  //To understand the kind of measurement, we use alignments
  Visualizer* viz = new Visualizer( "ATLAS" );

  if(useLabel){
      viz->SetOutputDirectory( m_outputDir );
      viz->SetTestBeamLabel( m_runNumber, m_alignment);
  }

  for( auto& ana : m_ana ){
    ana->Finalize();
  }
  for( auto& det_ana : m_det_ana ){
    det_ana->AssignVisualizer( viz );
    det_ana->Finalize();
  }

  m_tOut->Write();
  m_fOut->Close();
}
