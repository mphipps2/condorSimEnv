//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
/// \ingroup mc
/// \file ModTypeRPD.cc
/// \brief RPD detector construction
/// \author Aric Tate
/// \date February 2019



#include "ModTypeRPD.hh"
#include "FiberSD.hh"
#include "FastSimModelOpFiber.hh"
//#include "FastFiberModel.hh"

#include <iostream>
#include <stdio.h>

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4UserLimits.hh"


#include "G4RunManager.hh"
#include "G4SDManager.hh"

#include "G4NistManager.hh"
#include "G4CSGSolid.hh"
#include "G4Box.hh"
#include "G4Para.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeRPD::ModTypeRPD(const G4int cn, G4LogicalVolume* mother, G4ThreeVector* pos )
  : m_modNum( cn ),
    m_pos( pos ),
    m_fiberDiam(new G4ThreeVector(.71,.68,.73)),
    m_fiber_count(0),
    m_HousingThickness(2.0*mm),
    m_rpdRotation(0.*deg),
    m_n_rows(4),
    m_n_columns(4),
    m_n_cycles_per_tile(2),
    m_fiberPitchX(0.*mm),
    m_fiberPitchZ(0.*mm),
    m_tileSize(10.*mm),
    m_minWallThickness(0.*mm),
    m_distanceToReadout(0.*mm),
    m_detType(""),
    OPTICAL(false),
    CHECK_OVERLAPS(false),
    READOUT(false),
    REDUCED_TREE(false),
    ML_REDUCED_TREE(false),
    CLAD(false),
    BUFFERED(false),
    m_logicMother( mother )
{
  materials = Materials::getInstance();
  materials->UseOpticalMaterials(true);
  materials->DefineOpticalProperties();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeRPD::ModTypeRPD(const G4int cn, ModTypeRPD* right)
  : m_modNum( cn )
{
  m_pos 						  = new G4ThreeVector(*right->m_pos);
  m_fiberDiam 			  = new G4ThreeVector(*right->m_fiberDiam);
  m_fiber_count       = 0;
  m_HousingThickness  = right->m_HousingThickness;
  m_rpdRotation       = right->m_rpdRotation;
  m_n_rows            = right->m_n_rows;
  m_n_columns         = right->m_n_columns;
  m_n_cycles_per_tile = right->m_n_cycles_per_tile;
  m_fiberPitchX 		  = right->m_fiberPitchX;
  m_fiberPitchZ 		  = right->m_fiberPitchZ;
  m_tileSize 				  = right->m_tileSize;
  m_minWallThickness  = right->m_minWallThickness;
  m_distanceToReadout = right->m_distanceToReadout;
  m_detType 				  = right->m_detType;
  OPTICAL 					  = right->OPTICAL;
  CHECK_OVERLAPS 		  = right->CHECK_OVERLAPS;
  READOUT 		        = right->READOUT;
  REDUCED_TREE 		    = right->REDUCED_TREE;
  ML_REDUCED_TREE 		    = right->ML_REDUCED_TREE;
  CLAD         		    = right->CLAD;
  BUFFERED    		    = right->BUFFERED;
  materials					  = right->materials;
  m_logicMother 		  = right->m_logicMother;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeRPD::~ModTypeRPD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeRPD::Construct(){
  printf("Before materials\n");
  DefineMaterials();
  printf("after materials\n");
  if(m_detType == "cms"){
    ConstructCMSDetector();
  }else {
    ConstructPanFluteDetector();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeRPD::DefineMaterials()
{
  //----------------------------------------------
  // Define Materials
  //----------------------------------------------

  //Quartz_UMD
  m_matQuartz = materials->pQuartz;

  //SilicaCore_UI
  m_silicaCore_UI = materials->SilicaCore_UI;

  //SilicaClad_UI
  m_silicaClad_UI = materials->SilicaClad_UI;

  //Kapton_UI
  m_kapton_UI = materials->Kapton_UI;
  
  //Aluminum
  m_Al = materials->Al;

  //Air
  m_Air = materials->Air;
  //Polyethylene
  m_Poly = materials->Polyethylene;

  //Optical Grease
  m_Grease = materials->Grease;

  //Wavelength Shifter/Fiber
  m_PMMA = materials->PMMA;

  //Aluminum optical surface
  m_opAlSurface = materials->OpAlSurface;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeRPD::ConstructPanFluteDetector()
{

  //retrieve RPD parameters
  G4double fiber_diam = m_fiberDiam->x(); // Just core for now

  if(m_fiberDiam->x() < m_fiberDiam->y()){
    CLAD = true;
    fiber_diam = m_fiberDiam->y();
  }
  if(m_fiberDiam->x() < m_fiberDiam->z() && m_fiberDiam->y() < m_fiberDiam->z()){
    BUFFERED = true;
    fiber_diam = m_fiberDiam->z();
  }
  if(!CLAD || !BUFFERED){
    G4cout << "RPD" << m_modNum << " fibers will not have cladding or buffer" << G4endl;
  }

  //If the fiber is not clad but is buffered, the buffer will have to touch the core
  G4double buffer_inner = (CLAD) ? m_fiberDiam->y() : m_fiberDiam->x();

  char name[256];

  G4Colour colors[4] = { G4Colour::Cyan(),  G4Colour::Red(), G4Colour::Green(), G4Colour::Magenta() };

  G4double pitchX;
  G4double pitchZ;

  // If the user requested option is physically possible, use it. Otherwise calculate the pitch
  if( sqrt( pow(m_fiberPitchX/2.,2) + pow(m_fiberPitchZ/2.,2) ) > fiber_diam ){

    pitchX = m_fiberPitchX;
    pitchZ = m_fiberPitchZ;
    m_n_cycles_per_tile = (m_tileSize/pitchX)/m_n_rows; // This rounds down

    // If the tile size is too small, set n_cycles to 1 before recalculating m_tileSize
    if( m_n_cycles_per_tile == 0 ){
      G4cerr << "////////////////////////////" << G4endl;
      G4cerr << "Error in RPD" << m_modNum << ": Requested tile size and pitch are incompatible" << G4endl;
      G4cerr << "Recalculating tile size based on one cycle of the requested pitch" << G4endl;
      G4cerr << "////////////////////////////" << G4endl;
      m_n_cycles_per_tile = 1;
    }
    std::cout << " tileSize " << m_tileSize << " setting n cycles per tile " << m_n_cycles_per_tile << " pitchX " << pitchX << " m_n_rows " << m_n_rows << std::endl;
    // Tile size and pitch are correlated and must match precisely, so recalculate the tile size based on the pitch
    m_tileSize = pitchX*m_n_cycles_per_tile*m_n_rows;
  }else{
    //If the user specified a fiber pitch that doesn't work, inform them
    if( m_fiberPitchX < 0.0 ){
      G4cerr << "////////////////////////////" << G4endl;
      G4cerr << "Error in RPD" << m_modNum << ": Requested fiber diameter and pitch are incompatible" << G4endl;
      G4cerr << "Recalculating tile size based on one cycle of the requested fiber diameter" << G4endl;
      G4cerr << "////////////////////////////" << G4endl;
    }

    m_n_cycles_per_tile = 1;
    //Minimum possible tile size = sqrt(2) * minimum_pitch * m_n_cycles_per_tile * m_n_rows
    G4double min_tileSize = 0.707*(fiber_diam + m_minWallThickness)*m_n_cycles_per_tile*m_n_rows;
    if( m_tileSize < min_tileSize ){
      G4cerr << "////////////////////////////" << G4endl;
      G4cerr << "Error in RPD" << m_modNum << ": Requested fiber diameter and tile size are incompatible" << G4endl;
      G4cerr << "Recalculating tile size based on one cycle of the requested fiber diameter" << G4endl;
      G4cerr << "////////////////////////////" << G4endl;
      m_tileSize = min_tileSize;
    }

    pitchX = m_tileSize/(0.707*m_n_cycles_per_tile*m_n_rows);
    pitchZ = pitchX;
  }

  // Distance in X and Z from one fiber to another diagonally from it
  G4double offsetX = pitchX/2;
  G4double offsetZ = pitchZ/2;
  // Will be filled with fiber height for each row
  G4double fiber_height;
  // Positions of the current fiber for pattern1 and pattern2
  G4double posx[2], posz[2], posy;
  // Height of the active region
  G4double activeHeight = m_n_rows*m_tileSize;
  // Width (x) of the RPD housing
  G4double housingWidth = m_n_columns*m_tileSize + 2*m_HousingThickness;
  // Height (y) of the RPD housing
  G4double housingHeight = activeHeight + 2*m_HousingThickness + m_distanceToReadout;
  // Depth (z) of the RPD housing
  G4double housingDepth = (m_n_rows - 0.5)*pitchZ + fiber_diam + 2*m_HousingThickness;
  
  //create some rotation matrices
  G4RotationMatrix* stripRotation = new G4RotationMatrix();
  stripRotation->rotateX(90.*deg);
  G4RotationMatrix* nullRotation = new G4RotationMatrix();
  G4ThreeVector     nullVector(0.,0.,0.);


  m_fastOpticalRegion = new G4Region("fastOpticalRegion");
  std::cout << " RPD length (Z): " << housingDepth << std::endl;
  // Construct the housing
  m_PFrpd_housing = new G4Box( "RPDHousing", housingWidth*mm/2.0,
			       housingHeight*mm/2.0,
			       housingDepth*mm/2.0 );

  m_PFrpd_housingLogical = new G4LogicalVolume(m_PFrpd_housing,  m_Al,   "Housing_Logical");
  // add logical skin to the aluminum cavities
  //  G4LogicalSkinSurface *alSurface = new G4LogicalSkinSurface("AlSurface",m_PFrpd_housingLogical,m_opAlSurface);
  //  alSurface->DumpInfo();
  // set the optical boundary properties
  G4OpticalSurface *opAlSurface = new G4OpticalSurface("opAlSurface",unified, ground, dielectric_metal, 1);
  G4LogicalSkinSurface *alSurface = new G4LogicalSkinSurface("AlSurface",m_PFrpd_housingLogical,opAlSurface);
  if (alSurface) alSurface->DumpInfo();
  
  
  sprintf(name,"RPD%d_Case_Physical", m_modNum);
  G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
  rotationMatrix->rotateX(m_rpdRotation);
  m_PFrpd_housingPhysical =
    //    new G4PVPlacement(nullRotation,
    new G4PVPlacement(rotationMatrix,
                      G4ThreeVector(m_pos->x(), m_pos->y() + m_distanceToReadout/2.0 , m_pos->z()),
                      m_PFrpd_housingLogical,
                      name,
                      m_logicMother,
                      false,
                      m_modNum,
                      CHECK_OVERLAPS);

  m_PFrpd_housingLogical->SetVisAttributes( G4Colour(1.0,1.0,1.0,0.7));//G4Colour(1.0,0.0,0.0,0.3)

  
  // Construct the readout above RPD tiles
  if ( READOUT ){

    //Air volume inside aluminum case
    m_PFreadout_air = new G4Box( "RPDReadoutAir", ( housingWidth - m_HousingThickness) *mm/2.0,
				 m_distanceToReadout                *mm/2.0,
				 ( housingDepth - m_HousingThickness) *mm/2.0 );

    m_PFreadout_airLogical = new G4LogicalVolume( m_PFreadout_air, m_Air, "Readout_Air_Logical");

    sprintf(name,"RPD%d_ReadoutAir_Physical", m_modNum);
    new G4PVPlacement(nullRotation,
                      G4ThreeVector( 0.0 , (activeHeight )*mm/2.0 , 0.0 ),
                      m_PFreadout_airLogical,
                      name,
                      m_PFrpd_housingLogical,
                      false,
                      m_modNum,
                      CHECK_OVERLAPS);

    m_PFreadout_airLogical->SetVisAttributes(G4Colour(0.0,0.8,1.0,0.05));

    //    double airGap = 0.2;
    //Create readout fiber core
    sprintf(name,"m_PFreadoutCore_%d", m_modNum);
    m_PFreadout_fiberCore = new G4Tubs(name,
				       0.0*mm,
				       m_fiberDiam->x()*mm/2.0,
				       m_distanceToReadout*mm/2.0,
				       0.0*deg,
				       360.0*deg);

    if(CLAD){
      //Create readout fiber cladding
      sprintf(name,"m_PFreadoutClad_%d", m_modNum);
      m_PFreadout_fiberClad = new G4Tubs(name,
					 m_fiberDiam->x()*mm/2.0,
					 m_fiberDiam->y()*mm/2.0,
					 m_distanceToReadout*mm/2.0,
					 0.0*deg,
					 360.0*deg);
    }

    if(BUFFERED){
      //Create readout fiber buffer
      sprintf(name,"m_PFreadoutBuff_%d", m_modNum);
      m_PFreadout_fiberBuff = new G4Tubs(name,
					 buffer_inner*mm/2.0,
					 m_fiberDiam->z()*mm/2.0,
					 m_distanceToReadout*mm/2.0,
					 0.0*deg,
					 360.0*deg);
    }

  }

  // Loop over rows
  for(G4int row = 0; row < m_n_rows; row++){
    // Calculate the height of fibers in this row
    fiber_height = m_tileSize*(row+1) ;//Just the portion of the fiber in the active region

    // The fiber core. One length for each row
    sprintf(name,"m_PFrpd_%d", row);

   
    m_PFrpdCore.push_back( new G4Tubs(name,
				      0.0*mm,
				      m_fiberDiam->x()*mm/2.0,
				      //				      fiber_height*mm/2.0 - airGap*mm/2.0 ,
				      fiber_height*mm/2.0 ,
				      0.0*deg,
				      360.0*deg) );

    if(CLAD){
      // The fiber cladding. One length for each row
      sprintf(name,"m_PFrpd_%d", row);
      m_PFrpdClad.push_back( new G4Tubs(name,
					m_fiberDiam->x()*mm/2.0,
					m_fiberDiam->y()*mm/2.0,
					//					fiber_height*mm/2.0 - airGap*mm/2.0 ,
					fiber_height*mm/2.0 ,
					0.0*deg,
					360.0*deg) );
    }

    if(BUFFERED){
      // The fiber buffer. One length for each row
      sprintf(name,"m_PFrpd_%d", row);
      m_PFrpdBuff.push_back( new G4Tubs(name,
					buffer_inner*mm/2.0,
					m_fiberDiam->z()*mm/2.0,
					//					fiber_height*mm/2.0 - airGap*mm/2.0,
					fiber_height*mm/2.0,
					0.0*deg,
					360.0*deg) );
    }

    // Create an air channel for the fiber to be routed
    sprintf(name,"m_PFchannel_%d", row);
    m_PFrpd_channel.push_back( new G4Tubs( name,
					   fiber_diam*mm/2.0,
					   1.1*fiber_diam*mm/2.0,
					   fiber_height*mm/2.0 ,
					   0.0*deg,
					   360.0*deg) );

    for(G4int col = 0; col < m_n_columns; col++){
      //Now we're in the realm of working on a single tile
      //we have to cycle through the two patterns until the tile
      //is filled in X
      for(G4int cycle = 0; cycle < m_n_cycles_per_tile; cycle++){

	for(G4int fiber = 0; fiber < m_n_columns; fiber++){
	  // !!!!!!!Position calculations assume an even number of rows and columns!!!!!!! //////

	  //Start at RPD center + tile width* number of tiles + cycle number * cycle width + stack number in cycle * pitch
	  posx[0] = m_tileSize*((m_n_columns/2) - col ) - (cycle*m_n_columns*pitchX) - fiber*pitchX - pitchX/4;
	  posx[1] = posx[0] - offsetX;

	  //Start at Z center - half the stack depth + moving front to back and looping back to front
	  posz[0] = - pitchZ*(m_n_columns-0.5)/2 + pitchZ*((row + fiber	   )%4);
	  posz[1] = - pitchZ*(m_n_columns-0.5)/2 + pitchZ*((row + fiber + 2)%4) + offsetZ; //Pattern is offset by 2

          //Start from the top and work down.
          posy = ((m_n_rows/2 - 1) - row)*m_tileSize + fiber_height/2 - m_distanceToReadout/2;


          //Placement is comprised of two patterns
          for(G4int pattern = 0; pattern < 2; pattern ++){

            //----------------------- Place the channel in the housing -----------------------//
            sprintf(name,"m_PFchannel_log_%d",m_fiber_count);
            m_PFrpd_channelLogical.push_back( new G4LogicalVolume(m_PFrpd_channel.back(),m_Air,name) );

            m_PFrpd_channelLogical.back()->SetVisAttributes( G4Colour(0.0,0.0,1.0,0.1) );// G4Colour(0.0,0.0,1.0,0.1) or G4VisAttributes::Invisible

            sprintf(name,"m_PFchannel_physTest_%d", m_fiber_count);
            m_PFrpd_channelPhysical.push_back(
					      new G4PVPlacement(stripRotation,
								G4ThreeVector( posx[pattern]*mm , posy*mm , posz[pattern]*mm ),
								m_PFrpd_channelLogical.back(),
								name,
								m_PFrpd_housingLogical,
								false,
								m_fiber_count,
								CHECK_OVERLAPS) );

	    
            //----------------------- Place fiber core in the housing -----------------------//
            sprintf(name,"m_PFrpd_log_%d",m_fiber_count);
            m_PFrpdCoreLogical.push_back( new G4LogicalVolume( m_PFrpdCore.back(), m_silicaCore_UI, name) );

            m_PFrpdCoreLogical.back()->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));
	    m_PFrpdCoreLogical.back()->SetVisAttributes( colors[row] );

            sprintf(name,"m_PFrpd_phys_%d", m_fiber_count);
            m_PFrpdCorePhysical.push_back(
					  new G4PVPlacement(stripRotation,
							    //							    G4ThreeVector( posx[pattern]*mm , posy*mm , posz[pattern]*mm + airGap*mm),
							    G4ThreeVector( posx[pattern]*mm , posy*mm , posz[pattern]*mm),
							    m_PFrpdCoreLogical.back(),
							    name,
							    m_PFrpd_housingLogical,
							    false,
							    m_fiber_count,
							    CHECK_OVERLAPS) );
	    // build regions where fast optical physics will be used
	    m_PFrpdCoreLogical.back()->SetRegion(m_fastOpticalRegion);	       
	    m_fastOpticalRegion->AddRootLogicalVolume(m_PFrpdCoreLogical.back());

            //----------------------- Place fiber cladding in the housing -----------------------//
            if(CLAD){
	      sprintf(name,"m_PFrpdClad_log_%d",m_fiber_count);
	      m_PFrpdCladLogical.push_back( new G4LogicalVolume( m_PFrpdClad.back(), m_silicaClad_UI, name) );

	      m_PFrpdCladLogical.back()->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));
	      m_PFrpdCladLogical.back()->SetVisAttributes( colors[row] );

	      sprintf(name,"m_PFrpdClad_phys_%d", m_fiber_count);
	      m_PFrpdCladPhysical.push_back(
					    new G4PVPlacement(stripRotation,
							      //							      G4ThreeVector( posx[pattern]*mm , posy*mm , posz[pattern]*mm + airGap*mm),
							      G4ThreeVector( posx[pattern]*mm , posy*mm , posz[pattern]*mm),
							      m_PFrpdCladLogical.back(),
							      name,
							      m_PFrpd_housingLogical,
							      false,
							      m_fiber_count,
							      CHECK_OVERLAPS) );
	      
	      // build regions where fast optical physics will be used
	      // m_PFrpdCladLogical.back()->SetRegion(m_fastOpticalRegion);	       
	      //     m_fastOpticalRegion->AddRootLogicalVolume(m_PFrpdCladLogical.back());

            }
	    // Air cladding
	    else{
	      // build regions where fast optical physics will be used
	      //	      m_PFrpd_channelLogical.back()->SetRegion(m_fastOpticalRegion);	       
	      //	      m_fastOpticalRegion->AddRootLogicalVolume(m_PFrpd_channelLogical.back());
	    }
            //----------------------- Place fiber buffer in the housing -----------------------//
            if(BUFFERED){
	      sprintf(name,"m_PFrpdBuff_log_%d",m_fiber_count);
	      m_PFrpdBuffLogical.push_back( new G4LogicalVolume( m_PFrpdBuff.back(), m_kapton_UI, name) );

	      m_PFrpdBuffLogical.back()->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));
	      m_PFrpdBuffLogical.back()->SetVisAttributes( colors[row] );

	      sprintf(name,"m_PFrpdBuff_phys_%d", m_fiber_count);
	      m_PFrpdBuffPhysical.push_back(
					    new G4PVPlacement(stripRotation,
							      //							      G4ThreeVector( posx[pattern]*mm , posy*mm , posz[pattern]*mm + airGap*mm),
							      G4ThreeVector( posx[pattern]*mm , posy*mm , posz[pattern]*mm),
							      m_PFrpdBuffLogical.back(),
							      name,
							      m_PFrpd_housingLogical,
							      false,
							      m_fiber_count,
							      CHECK_OVERLAPS) );
            }

            //----------------------- Place readout fiber in the air volume if selected-----------------------//
            if ( READOUT ){
              //Core
              sprintf(name,"m_PFreadoutCore_log_%d",m_fiber_count);
              m_PFreadout_fiberCoreLogical.push_back(
						     new G4LogicalVolume(
									 m_PFreadout_fiberCore,
									 m_silicaCore_UI,
									 name) );

              m_PFreadout_fiberCoreLogical.back()->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));
	      m_PFreadout_fiberCoreLogical.back()->SetVisAttributes(colors[row]); //G4Colour(1.0,0.0,1.0,0.9)

              sprintf(name,"phys_PFreadoutCore_%d", m_fiber_count);
              m_PFreadout_fiberCorePhysical.push_back(
						      new G4PVPlacement(stripRotation,
									G4ThreeVector( posx[pattern]*mm , 0.0 , posz[pattern]*mm ),
									m_PFreadout_fiberCoreLogical.back(),
									name,
									m_PFreadout_airLogical,
									false,
									m_fiber_count,
									CHECK_OVERLAPS) );
	      // build regions where fast optical physics will be used
	      m_PFreadout_fiberCoreLogical.back()->SetRegion(m_fastOpticalRegion);	       
	      m_fastOpticalRegion->AddRootLogicalVolume(m_PFreadout_fiberCoreLogical.back());
              //Cladding
              if(CLAD){
		sprintf(name,"m_PFreadoutClad_log_%d",m_fiber_count);
		m_PFreadout_fiberCladLogical.push_back(
						       new G4LogicalVolume(
									   m_PFreadout_fiberClad,
									   m_silicaClad_UI,
									   name) );

		m_PFreadout_fiberCladLogical.back()->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));
		m_PFreadout_fiberCladLogical.back()->SetVisAttributes(colors[row]); //G4Colour(1.0,0.0,1.0,0.9)

		sprintf(name,"phys_PFreadoutClad_%d", m_fiber_count);
		m_PFreadout_fiberCladPhysical.push_back(
							new G4PVPlacement(stripRotation,
									  G4ThreeVector( posx[pattern]*mm , 0.0 , posz[pattern]*mm ),
									  m_PFreadout_fiberCladLogical.back(),
									  name,
									  m_PFreadout_airLogical,
									  false,
									  m_fiber_count,
									  CHECK_OVERLAPS) );
		// build regions where fast optical physics will be used
		//		m_PFreadout_fiberCladLogical.back()->SetRegion(m_fastOpticalRegion);	       
		//	m_fastOpticalRegion->AddRootLogicalVolume(m_PFreadout_fiberCladLogical.back());
              }
              //Buffer
              if(BUFFERED){
		sprintf(name,"m_PFreadoutBuff_log_%d",m_fiber_count);
		m_PFreadout_fiberBuffLogical.push_back(
						       new G4LogicalVolume(
									   m_PFreadout_fiberBuff,
									   m_kapton_UI,
									   name) );

		m_PFreadout_fiberBuffLogical.back()->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));
		m_PFreadout_fiberBuffLogical.back()->SetVisAttributes(colors[row]); //G4Colour(1.0,0.0,1.0,0.9)

		sprintf(name,"phys_PFreadoutBuff_%d", m_fiber_count);
		m_PFreadout_fiberBuffPhysical.push_back(
							new G4PVPlacement(stripRotation,
									  G4ThreeVector( posx[pattern]*mm , 0.0 , posz[pattern]*mm ),
									  m_PFreadout_fiberBuffLogical.back(),
									  name,
									  m_PFreadout_airLogical,
									  false,
									  m_fiber_count,
									  CHECK_OVERLAPS) );
              }
            }//end readout
	    m_fiber_count++;
          }//end pattern
	}//end fiber
      }//end cycle
    }//end column
  }//end row
  
  // Air cladding
  if (!CLAD) {
    // build regions where fast optical physics will be used
    //    m_PFreadout_airLogical->SetRegion(m_fastOpticalRegion);	       
    //    m_fastOpticalRegion->AddRootLogicalVolume(m_PFreadout_airLogical);
  }




  
  m_topOfVolume = m_pos->y() + m_distanceToReadout +  activeHeight/2.0;
  m_bottomOfVolume = -1 * activeHeight/2.0;
  std::cout << "Prototype RPD construction finished with ";
  std::cout << m_fiber_count << " Fibers, " << m_tileSize << "mm tiles." << std::endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeRPD::ConstructCMSDetector()
{
  m_test_tile_bool = false;

  if(m_test_tile_bool){
    for(G4int i =0;i<8;i++){
      rpd_comp[i]=false;
    }
  }
  else for(G4int i =0;i<8;i++){
      rpd_comp[i]=true;
    }

  // RPD components
  bool tile_flag = rpd_comp[0];				// quartz tiles
  bool fiber_flag = rpd_comp[1];			// readout fibers
  bool cladding_flag = rpd_comp[2];		// fiber cladding
  bool grease_flag = rpd_comp[3];			// optical grease
  bool foil_flag = rpd_comp[4];				// aluminum foil
  bool fr_bck_flag = rpd_comp[5]; 		// front and back aluminum casing
  bool vertical_flag = rpd_comp[6]; 	// veritcal aluminum dividers
  bool air_detec_flag = rpd_comp[7];	// photo detectors


  //retrieve RPD parameters
  G4double tileX 		= 20;
  G4double tileY 		= 20;
  G4double tileZ 		= 10;
  G4double halfX_gap 	= (1.58/2)*mm;
  G4double fiber_diam 	= 1*mm;
  G4double foil_thickness 	= 0.016*mm;


  //manually create some RPD parameters
  G4double halfY_gap 	= (2.5*foil_thickness/2)*mm;
  G4double hole_center_offset = 0.05;
  G4double case_thickness = 0.9*2.0*halfX_gap;
  G4double core_diam;
  G4double grease_offset;
  if(grease_flag) grease_offset = 0.05;
  else grease_offset = 0.0;

  if(cladding_flag) core_diam = 0.970*(fiber_diam)*mm;
  else 							core_diam = (fiber_diam)*mm;

  G4double fiberHeightY[4];
  G4double foilHeightY[4];

  char name[256];

  //create some rotation matrices
  G4RotationMatrix* stripRotation = new G4RotationMatrix();
  stripRotation->rotateX(90.*deg);
  G4RotationMatrix* nullRotation = new G4RotationMatrix();

  //create tile that has 4 holes to be used repeatedly later
  m_tile_no_fiber_hole = new G4Box("tile_no_hole",(tileX/2)*mm, (tileY/2)*mm, (tileZ/2)*mm);
  m_fiber_subtract = new G4Tubs( "fiber_sub",
				 0.0*mm,
				 ((fiber_diam/2.0)+ grease_offset )*mm, //1mm diameter
				 tileY*mm ,
				 0.0*deg,
				 360.0*deg);

  m_tile = new G4SubtractionSolid("tile",m_tile_no_fiber_hole,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 1*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset)*mm));
  m_tile = new G4SubtractionSolid("tile",m_tile,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 2*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));
  m_tile = new G4SubtractionSolid("tile",m_tile,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 3*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));
  m_tile = new G4SubtractionSolid("tile",m_tile,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 4*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));

  m_tileLogical	 = new G4LogicalVolume(m_tile, m_matQuartz, "tile_Logical");

  G4VisAttributes* quartzColor  = new G4VisAttributes(  G4Colour::Cyan());//
  //quartzColor->SetForceSolid(true);
  m_tileLogical->SetVisAttributes( quartzColor );

  G4int cn = 0, cn_fiber = 0;

  G4double RPD_centerX = m_pos->getX();
  G4double RPD_centerY = m_pos->getY();
  G4double RPD_centerZ = m_pos->getZ();

  G4double RPD_startX = RPD_centerX + 3*halfX_gap + 1.5*tileX;
  G4double RPD_startY = RPD_centerY + 3*halfY_gap + 1.5*tileY;


  //positioning variables


  G4double factor1 = 1.5;
  G4double hole_shift		= grease_offset + hole_center_offset;

  G4double tileZcenter[4];
  for(G4int i=0; i < 4; i++){
    tileZcenter[i] = RPD_centerZ + (i * (factor1*fiber_diam + foil_thickness) );
  }
  G4double fiberZcenter[4];
  for(G4int i=0; i < 4; i++){
    fiberZcenter[i] = tileZcenter[i] + tileZ/2.0 - fiber_diam/2.0 - hole_shift;
  }

  G4double foil_gap = (fiberZcenter[1]-(fiber_diam/2.0)) - (fiberZcenter[0] + (fiber_diam/2.0));

  G4double foilZcenter[4];
  for(G4int i=0; i < 4; i++){
    foilZcenter[i] = fiberZcenter[i] + fiber_diam/2.0 + foil_gap/2.0;
  }

  G4double final_foil_pos = foilZcenter[3]+foil_thickness/2.0;
  G4double assembly_midZ = RPD_centerZ + (final_foil_pos - (RPD_centerZ + (tileZ/2.0)))/2.0;
  G4double half_Z_length = final_foil_pos-assembly_midZ;


  //create fibers/cladding/grease with correct lengths
  for(G4int k=0;k<4;k++) {
    fiberHeightY[k]=((k+1)*tileY)+(k*2*halfY_gap);

    sprintf(name,"fiber_a %d", k);
    m_fiber[k] 		= new G4Tubs( name,
				      0.0*mm,
				      (core_diam/2.0)*mm, //1mm diameter
				      fiberHeightY[k]*mm/2.0 ,
				      0.0*deg,
				      360.0*deg);
    if(grease_flag){
      sprintf(name,"fibergrease_a %d", k);
      m_fibergrease[k] 		= new G4Tubs( name,
					      (fiber_diam/2.0)*mm,
					      ((fiber_diam+grease_offset)/2.0)*mm,
					      fiberHeightY[k]*mm/2.0 ,
					      0.0*deg,
					      360.0*deg);
    }

    if(cladding_flag){
      sprintf(name,"fiberclad_a %d", k);
      m_fiberclad[k] 		= new G4Tubs( name,
					      (core_diam/2.0)*mm,
					      ((fiber_diam)/2.0)*mm,
					      fiberHeightY[k]*mm/2.0 ,
					      0.0*deg,
					      360.0*deg);
    }

    if(grease_flag){
      sprintf(name,"fibergreaseLogical_a %d", k);

      m_fibergreaseLogical[k] 	= new G4LogicalVolume(m_fibergrease[k]
						      ,m_Grease,//m_Poly
						      name);
      m_fibergreaseLogical[k]->SetVisAttributes( G4Colour(1,1,1,0.3) );//G4VisAttributes::Invisible
    }

    if(cladding_flag){
      sprintf(name,"fibercladLogical_a %d", k);

      m_fibercladLogical[k] 	= new G4LogicalVolume(m_fiberclad[k]
						      ,m_Poly,//m_Poly
						      name);
      m_fibercladLogical[k]->SetVisAttributes( G4Colour(1,0,0,0.2));//G4VisAttributes::Invisible
    }
  }

  for(G4int i=0;i<4;i++) {
    for(G4int j=0;j<4;j++) {
      for(G4int k=0;k<4;k++) {
	sprintf(name,"fiberLogical_%d_%d_%d",k,j,i);

	m_fiberLogical[m_fiber_count] 	= new G4LogicalVolume(m_fiber[k]
							      ,m_PMMA,
							      name);
	m_fiberLogical[m_fiber_count]->SetVisAttributes( G4Colour(0.0,0.0,1.0,0.2) );

	m_fiberLogical[m_fiber_count]->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));

	m_fiber_count++;
      }
    }
  }

  // create alum containment

  //create foil sections
  for(G4int k=0;k<4;k++) {
    foilHeightY[k]=(2*(k+1)-2) * halfY_gap + (k+1)*tileY;
    sprintf(name,"foil_V %d", k);

    m_foilV[k] 		= new G4Box(name,(tileX/2.0)*mm, (foilHeightY[k]/2.0)*mm, (0.75*foil_thickness/2.0)*mm);//

    sprintf(name,"foilVlog %d", k);

    m_foilVLogical[k] 	= new G4LogicalVolume(m_foilV[k]
					      ,m_Al,
					      name);
    m_foilVLogical[k]->SetVisAttributes(G4Colour(0.9,0.0,0.0,0.4));//G4VisAttributes::Invisible
  }


  m_foilVfront 		= new G4Box("m_foilVfront",(tileX/2.0)*mm, (tileY/2.0)*mm, (0.75*foil_thickness/2.0)*mm);//
  m_foilVfrontLogical 	= new G4LogicalVolume(m_foilVfront
					      ,m_Al,
					      "foilVfrntlog");
  m_foilVfrontLogical->SetVisAttributes( G4Colour(0.9,0.0,0.0,0.4) );//


  m_foilH 			= new G4Box("foil_H",(tileX/2.0)*mm, (0.75*foil_thickness/2.0)*mm, ((tileZ+(foil_gap)-hole_center_offset)/2.0) *mm) ;// foil/2.0
  m_foilHLogical 	= new G4LogicalVolume(m_foilH
					      ,m_Al,
					      "foil_H_log");
  m_foilHLogical->SetVisAttributes( G4Colour(0.9,0.0,0.0,0.4));//


  m_foilHtop 			= new G4Box("foil_H_top",(tileX/2.0)*mm, (0.75*foil_thickness/2.0)*mm, ((tileZ)/2.0) *mm) ;//+(foil_gap)-hole_center_offset
  m_foilHtop_hole = new G4SubtractionSolid("foil_H_top_hole",m_foilHtop,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 1*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset)*mm));
  m_foilHtop_hole = new G4SubtractionSolid("foil_H_top_hole",m_foilHtop_hole,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 2*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));
  m_foilHtop_hole = new G4SubtractionSolid("foil_H_top_hole",m_foilHtop_hole,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 3*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));
  m_foilHtop_hole = new G4SubtractionSolid("foil_H_top_hole",m_foilHtop_hole,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 4*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));

  m_foilHtopLogical 	= new G4LogicalVolume(m_foilHtop_hole
					      ,m_Al,
					      "foil_Htop_log");
  m_foilHtopLogical->SetVisAttributes(  G4Colour(0.9,0.0,0.0,0.4) );//

  //create alum case pieces

  //vertical dividers
  m_AlcaseV 			= new G4Box("case_V",(case_thickness/2.0)*mm, ((fiberHeightY[3]+case_thickness)/2.0)*mm, ((final_foil_pos)-(RPD_centerZ-(tileZ/2.0))+case_thickness)/2.0 *mm );
  m_AlcaseVLogical 	= new G4LogicalVolume(m_AlcaseV
					      ,m_Al,
					      "case_V_log");
  m_AlcaseVLogical->SetVisAttributes( G4Colour(0.0,0.0,0.9,0.4));

  //front and back plate
  m_Alcase 			= new G4Box("case_fr_bck",((4*tileX+8*halfX_gap+case_thickness)/2.0)*mm, ((fiberHeightY[3]+case_thickness)/2.0)*mm, ((0.99*case_thickness)/2.0)*mm );
  m_AlcaseLogical 	= new G4LogicalVolume(m_Alcase
					      ,m_Al,
					      "case_fr_bck_log");
  m_AlcaseLogical->SetVisAttributes( G4Colour(0.0,0.0,0.9,0.4));//

  ///////////////////////////////////////////////////////////////////////////////////
  // if Optical is turned on we need to build air detector above each fiber

  if ( OPTICAL ){

    G4double air_detect_thickness = 0.1;
    G4double correction = 0.0;

    m_air_detect = new G4Tubs( "air_detect",
			       0.0*mm,
			       ((core_diam/2.0)-correction)*mm,
			       ((air_detect_thickness/2.0))*mm ,
			       0.0*deg,
			       360.0*deg);


    G4int air_detecT_cnt = 0;
    for(G4int i=0;i<4;i++) {
      for(G4int j=0;j<4;j++) {
	for(G4int k=0;k<4;k++) {
	  sprintf(name,"photo_detect_log_%d_%d_%d",i, j, k);
	  m_air_detect_Logical[air_detecT_cnt]	 = new G4LogicalVolume(m_air_detect, m_PMMA, name);
	  m_air_detect_Logical[air_detecT_cnt]->SetVisAttributes( G4Colour(0.0,1.0,0.0,5.0) );

	  sprintf(name,"photo_detect_phys_%d_%d_%d",i, j, k);

	  if(air_detec_flag){

	    m_air_detectPhysical[air_detecT_cnt]		  = new G4PVPlacement(
										      nullRotation,
										      G4ThreeVector(( 0.0 ) *mm ,
												    ( 0.0 ) 	 *mm ,
												    (fiberHeightY[k]/2 - (air_detect_thickness/2.0) - correction)*mm),
										      m_air_detect_Logical[air_detecT_cnt],
										      name,
										      m_fiberLogical[air_detecT_cnt],
										      false,
										      air_detecT_cnt,
										      CHECK_OVERLAPS);

	    air_detecT_cnt++;

	  }
	}
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  //place vertical dividers
  for(G4int k=0;k<5;k++) {

    sprintf(name,"case_V %d", k);

    if(vertical_flag){
      m_AlcaseVPhysical[k] = new G4PVPlacement(
					       nullRotation,
					       G4ThreeVector((( 	(RPD_startX + tileX/2.0) + halfX_gap  ) - k*(2*halfX_gap+tileX)  ) 	*mm ,
							     (RPD_centerY )   	*mm ,
							     (assembly_midZ )			*mm),
					       m_AlcaseVLogical,
					       name,
					       m_logicMother,
					       false,
					       k,
					       CHECK_OVERLAPS);
    }
  }

  //place front/back plate
  for(G4int k=0;k<2;k++) {

    sprintf(name,"case_fr_bck %d", k);

    if(fr_bck_flag){
      m_AlcasePhysical[k] = new G4PVPlacement(
					      nullRotation,
					      G4ThreeVector(  (RPD_centerX) 	*mm ,
							      (RPD_centerY)   	*mm ,
							      (assembly_midZ - half_Z_length - case_thickness + 2.0*k*(half_Z_length+case_thickness) )			*mm),
					      m_AlcaseLogical,
					      name,
					      m_logicMother,
					      false,
					      k,
					      CHECK_OVERLAPS);
    }
  }
  //place tiles, fibers, foil
  for(G4int j=0;j<4;j++) {
    for(G4int i=0;i<4;i++) {

      if(tile_flag){
	sprintf(name,"tilephys_%d_%d", i,j);
    	m_tilePhysical[i][j] = new G4PVPlacement(
						 nullRotation,
						 G4ThreeVector( ( 	RPD_startX - (i*( tileX+(2*halfX_gap) ) ) )  	*mm ,
								( 	RPD_startY - (j*( tileY+(2*halfY_gap) ) ) )	*mm ,
								tileZcenter[j] *mm),
						 m_tileLogical,
						 name,
						 m_logicMother,
						 false,
						 cn,
						 CHECK_OVERLAPS);
      }



      if(foil_flag){
	sprintf(name,"foil_V_%d_%d",  i,j);
	m_foilVPhysical[i][j] = new G4PVPlacement(
						  nullRotation,
						  G4ThreeVector( ( 	RPD_startX - (j*( tileX+(2*halfX_gap) ) ) )  	*mm ,
								 ( RPD_startY + tileY/2 - (fiberHeightY[i]/2)  ) 		*mm ,
								 ( foilZcenter[i]) *mm),
						  m_foilVLogical[i],
						  name,
						  m_logicMother,
						  false,
						  cn,
						  CHECK_OVERLAPS);

	sprintf(name,"foil_V_front_%d_%d",  i,j);

	m_foilVfrontPhysical[i][j] = new G4PVPlacement(
						       nullRotation,
						       G4ThreeVector( ( 	RPD_startX - (i*( tileX+(2*halfX_gap) ) ) )  	*mm ,
								      ( 	RPD_startY - (j*( tileY+(2*halfY_gap) ) ) ) 		*mm ,
								      ( tileZcenter[j] - tileZ/2.0 -foil_gap/2.0) *mm),//fiberZcenter[i] + fiber_diam/2.0 + foil_gap/2.0;
						       m_foilVfrontLogical,
						       name,
						       m_logicMother,
						       false,
						       cn,
						       CHECK_OVERLAPS);


	sprintf(name,"foil_H_%d_%d", i,j);

	m_foilHPhysical[i][j] = new G4PVPlacement(
						  nullRotation,
						  G4ThreeVector( ( 	RPD_startX - (j*( tileX+(2*halfX_gap) ) ) )  	*mm ,
								 ( RPD_startY + tileY/2 + 1.3*halfY_gap - (i+1) * (tileY + 2 * halfY_gap ))   		*mm ,
								 (tileZcenter[i]  - hole_center_offset/2.0)			*mm),
						  m_foilHLogical,
						  name,
						  m_logicMother,
						  false,
						  i,
						  CHECK_OVERLAPS);

	if(i==0){

	  sprintf(name,"foil_H_top_%d", j);

	  m_foilHtopPhysical[j] = new G4PVPlacement(
						    nullRotation,
						    G4ThreeVector( ( 	RPD_startX - (j*( tileX+(2*halfX_gap) ) ) )  	*mm ,
								   (RPD_startY + tileY/2 + 1.5*foil_thickness)   		*mm ,
								   (tileZcenter[0] *mm)),
						    m_foilHtopLogical,
						    name,
						    m_logicMother,
						    false,
						    j,
						    CHECK_OVERLAPS);
	}
      }

      for(G4int k=0;k<4;k++) {
	if(fiber_flag){
	  sprintf(name,"fiberphys_%d_%d_%d",i, j, k);
	  m_fiberPhysical[cn_fiber]		  = new G4PVPlacement(
								      stripRotation,
								      G4ThreeVector( ( RPD_startX - (j*( tileX+(2*halfX_gap) ) ) ) + (tileX/2) - ((i+1)*tileX/5) *mm ,
										     ( RPD_startY + tileY/2 - (fiberHeightY[k]/2)  )	 *mm ,
										     (fiberZcenter[k]) *mm),
								      m_fiberLogical[cn_fiber],
								      name,
								      m_logicMother,
								      false,
								      cn_fiber,
								      CHECK_OVERLAPS);
	}
	if(cladding_flag){
	  sprintf(name,"cladphys_%d_%d_%d",i, j, k);
	  m_fibercladPhysical[cn_fiber]		  = new G4PVPlacement(
								      stripRotation,
								      G4ThreeVector( ( RPD_startX - (j*( tileX+(2*halfX_gap) ) ) ) + (tileX/2) - ((i+1)*tileX/5) *mm ,
										     ( RPD_startY + tileY/2 - (fiberHeightY[k]/2)  )	 *mm ,
										     (fiberZcenter[k]) *mm),
								      m_fibercladLogical[k],
								      name,
								      m_logicMother,
								      false,
								      cn_fiber,
								      CHECK_OVERLAPS);
	}
	if(grease_flag){
	  sprintf(name,"greasephys_%d_%d_%d",i, j, k);
	  m_fibergreasePhysical[cn_fiber]		  = new G4PVPlacement(
									      stripRotation,
									      G4ThreeVector( ( RPD_startX - (j*( tileX+(2*halfX_gap) ) ) ) + (tileX/2) - ((i+1)*tileX/5) *mm ,
											     ( RPD_startY + tileY/2 - (fiberHeightY[k]/2)  )	 *mm ,
											     (fiberZcenter[k]) *mm),
									      m_fibergreaseLogical[k],
									      name,
									      m_logicMother,
									      false,
									      cn_fiber,
									      CHECK_OVERLAPS);
	}

	cn_fiber++;
      }
      cn++;
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // TEST SETUP
  if(m_test_tile_bool){
    G4double testtileY = 100;
    G4double testtileZ = 50;
    G4double testgreaseY = 30;
    G4double testgreaseZ = 2.5;
    G4double testPDY = 2.5;
    G4double testcoreZ = 5;
    G4double testblckZ 	= 15;
    G4double offsetZ 	= testtileZ/2+testblckZ/2;
    G4double offsetY 	= 5;

    m_test_tile 		= new G4Box("test_tile",(100/2)*mm, (testtileY/2)*mm, (testtileZ/2)*mm);
    m_test_alum 		= new G4Box("test_alum",(100/2)*mm, (testtileY/2)*mm, (testblckZ/2)*mm);
    m_test_wls 			= new G4Box("test_wls",(100/2)*mm, (testgreaseY/2)*mm, (testcoreZ/2)*mm);
    m_test_PD 			= new G4Box("test_PD",(100/2)*mm, (testPDY/2)*mm, (testcoreZ/2)*mm);
    m_test_clad 		= new G4Box("test_clad",(100/2)*mm, (testgreaseY/2)*mm, (testgreaseZ/2)*mm);
    m_test_grease 	= new G4Box("test_grease",(100/2)*mm, (testgreaseY/2)*mm, (testgreaseZ/2)*mm);
    m_test_block 		= new G4Box("test_blck1",(100/2)*mm, ((testgreaseY)/2)*mm, (testblckZ/2)*mm);

    m_test_tileLogical	   = new G4LogicalVolume(m_test_tile, m_matQuartz, "testtile_Logical");
    m_test_alumLogical	   = new G4LogicalVolume(m_test_alum, m_Al, "testalum_Logical");
    m_test_wlsLogical	     = new G4LogicalVolume(m_test_wls, m_PMMA, "testwls_Logical");
    m_test_PDLogical	     = new G4LogicalVolume(m_test_PD, m_PMMA, "testPD_Logical");
    m_test_cladLogical	   = new G4LogicalVolume(m_test_clad, m_Poly, "testclad_Logical");
    m_test_greaseLogical	 = new G4LogicalVolume(m_test_grease, m_Grease, "testgrease_Logical");
    m_test_blockLogical	   = new G4LogicalVolume(m_test_block, m_Al, "testblck_Logical");

    m_test_tileLogical->SetVisAttributes( G4Colour(0.0,1.0,0.0,0.3) );
    m_test_alumLogical->SetVisAttributes( G4Colour(0.8,0.0,0.0,0.3) );
    m_test_wlsLogical->SetVisAttributes( G4Colour(0.0,0.0,1.0,0.3) );
    m_test_PDLogical->SetVisAttributes( G4Colour(0.0,0.0,1.0,0.3) );
    m_test_cladLogical->SetVisAttributes( G4Colour(0.3,0.5,0.0,0.3) );
    m_test_greaseLogical->SetVisAttributes( G4Colour(0.0,0.6,0.8,0.3) );
    m_test_blockLogical->SetVisAttributes( G4Colour(0.8,0.0,0.0,0.3) );

    m_test_tilePhysical = new G4PVPlacement(
					    nullRotation,
					    G4ThreeVector(( 0.0  ) 	*mm ,
							  (0.0 )   	*mm ,
							  (assembly_midZ)			*mm),
					    m_test_tileLogical,
					    "testtilephys",
					    m_logicMother,
					    false,
					    1,
					    CHECK_OVERLAPS);

    m_test_greasePhysical[0] = new G4PVPlacement(
						 nullRotation,
						 G4ThreeVector(( 0.0  ) 	*mm ,
							       (testtileY/2 + testgreaseY/2 - offsetY)   	*mm ,
							       (assembly_midZ - testcoreZ/2 - 1.5*testgreaseZ +offsetZ)			*mm),
						 m_test_greaseLogical,
						 "testgrease1phys",
						 m_logicMother,
						 false,
						 0,
						 CHECK_OVERLAPS);

    m_test_greasePhysical[1] = new G4PVPlacement(
						 nullRotation,
						 G4ThreeVector(( 0.0  ) 	*mm ,
							       (testtileY/2 + testgreaseY/2- offsetY )   	*mm ,
							       (assembly_midZ + testcoreZ/2 + 1.5*testgreaseZ+offsetZ )			*mm),
						 m_test_greaseLogical,
						 "testgrease2phys",
						 m_logicMother,
						 false,
						 1,
						 CHECK_OVERLAPS);

    m_test_cladPhysical[0] = new G4PVPlacement(
					       nullRotation,
					       G4ThreeVector(( 0.0  ) 	*mm ,
							     (testtileY/2 + testgreaseY/2- offsetY )   	*mm ,
							     (assembly_midZ - testcoreZ/2 - testgreaseZ/2+offsetZ )			*mm),
					       m_test_cladLogical,
					       "testclad1phys",
					       m_logicMother,
					       false,
					       0,
					       CHECK_OVERLAPS);

    m_test_cladPhysical[1] = new G4PVPlacement(
					       nullRotation,
					       G4ThreeVector(( 0.0  ) 	*mm ,
							     (testtileY/2 + testgreaseY/2- offsetY )   	*mm ,
							     (assembly_midZ+ testcoreZ/2 + testgreaseZ/2+offsetZ )			*mm),
					       m_test_cladLogical,
					       "testclad2phys",
					       m_logicMother,
					       false,
					       1,
					       CHECK_OVERLAPS);

    m_test_wlsPhysical = new G4PVPlacement(
					   nullRotation,
					   G4ThreeVector(( 0.0  ) 	*mm ,
							 (testtileY/2 + testgreaseY/2- offsetY )   	*mm,
							 (assembly_midZ + offsetZ)			*mm),
					   m_test_wlsLogical,
					   "testtilewls",
					   m_logicMother,
					   false,
					   1,
					   CHECK_OVERLAPS);

    m_test_PDPhysical = new G4PVPlacement(
					  nullRotation,
					  G4ThreeVector(( 0.0  ) 	*mm ,
							(testgreaseY/2 - testPDY/2  )   	*mm ,
							(0.0)			*mm),
					  m_test_PDLogical,
					  "photo_detect",
					  m_test_wlsLogical,
					  false,
					  1,
					  CHECK_OVERLAPS);

    m_test_alumPhysical = new G4PVPlacement(
					    nullRotation,
					    G4ThreeVector(( 0.0  ) 	*mm ,
							  (0.0 - offsetY )   	*mm ,
							  (assembly_midZ + testblckZ/2 + testtileZ/2 )			*mm),
					    m_test_alumLogical,
					    "testalumphys",
					    m_logicMother,
					    false,
					    1,
					    CHECK_OVERLAPS);

    m_test_blockPhysical[0] = new G4PVPlacement(
						nullRotation,
						G4ThreeVector(( 0.0  ) 	*mm ,
							      (testtileY/2 + testgreaseY/2 )   	*mm ,
							      (assembly_midZ + testblckZ/2 + 2.0* testgreaseZ + testcoreZ/2 +offsetZ)			*mm),
						m_test_blockLogical,
						"testblck1phys",
						m_logicMother,
						false,
						1,
						CHECK_OVERLAPS);

    m_test_blockPhysical[1] = new G4PVPlacement(
						nullRotation,
						G4ThreeVector(( 0.0  ) 	*mm ,
							      (testtileY/2 + testgreaseY/2   )   	*mm ,
							      (assembly_midZ - testblckZ/2 - 2.0* testgreaseZ - testcoreZ/2 +offsetZ)			*mm),
						m_test_blockLogical,
						"testblck2phys",
						m_logicMother,
						false,
						1,
						CHECK_OVERLAPS);


  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeRPD::ConstructSDandField(){
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  //Create and initialize the SD
  char fiberSDname[256];
  sprintf( fiberSDname, "RPD%d_SD", m_modNum);
  FiberSD* aFiberSD = new FiberSD( fiberSDname, m_modNum, OPTICAL );
  aFiberSD->HistInitialize();

  aFiberSD->SetTopOfVolume( m_topOfVolume );
  aFiberSD->SetBottomOfVolume( m_bottomOfVolume );
  aFiberSD->SetnFibers( m_fiber_count );
  SDman->AddNewDetector( aFiberSD );
  if(REDUCED_TREE) aFiberSD->SetReducedTree( m_fiber_count, GetnChannels() );
  if(ML_REDUCED_TREE) {
    aFiberSD->SetMLReducedTree( m_fiber_count, GetnChannels() );
  }
  // Assign fibers as SD volumes
  if(m_detType == "cms"){
    if(m_test_tile_bool){
      m_test_PDLogical->SetSensitiveDetector( aFiberSD );
    }else{
      if (OPTICAL){
	for(G4int i=0;i<64;i++){
	  m_air_detect_Logical[i]->SetSensitiveDetector( aFiberSD );}
      }else{
	for(G4int k=0;k<64;k++){
	  m_fiberLogical[k]->SetSensitiveDetector( aFiberSD );}}
    }
  }else{// not CMS i.e. Pan Flute
    for(int i = 0; i < m_fiber_count; i++){
      m_PFrpdCoreLogical.at(i)->SetSensitiveDetector( aFiberSD );
      if(CLAD) m_PFrpdCladLogical.at(i)->SetSensitiveDetector( aFiberSD );
      if(BUFFERED) m_PFrpdBuffLogical.at(i)->SetSensitiveDetector( aFiberSD );
      if( READOUT ){
	m_PFreadout_fiberCoreLogical.at(i)->SetSensitiveDetector( aFiberSD );
	if(CLAD) m_PFreadout_fiberCladLogical.at(i)->SetSensitiveDetector( aFiberSD );
	if(BUFFERED) m_PFreadout_fiberBuffLogical.at(i)->SetSensitiveDetector( aFiberSD );
      }
    }// end fiber loop
    m_fastOptical = new FastSimModelOpFiber("FastSimModelOpFiber",m_fastOpticalRegion, m_distanceToReadout, m_fiber_count, GetnChannels());
    //    if (FASTOPTICAL) m_fastOptical = new FastFiberModel("FastFiberModel",m_fastOpticalRegion);
  }// end else CMS

 
  std::cout << "RPD SD construction finished: SD name " << fiberSDname << std::endl;

}
