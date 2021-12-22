/// \ingroup mc
/// \file Materials.cc
/// \author Riccardo Longo

#include "Materials.hh"
#include "G4SystemOfUnits.hh"

Materials* Materials::materials = NULL;

Materials* Materials::getInstance(void)
{
    if (materials == NULL) {
        materials = new Materials();
    }
    return materials;
}

/* Function that constructs all the materials - to be then retrieved via pointer */
Materials::Materials(void){

  H  = new G4Element("Hydrogen", "H",  1.,   1.01   *g/mole);
  N  = new G4Element("Nitrogen", "N",  7.,   14.01  *g/mole);
  Ni = new G4Element("Ni",       "Ni", 28.0, 58.6934*g/mole);
  W  = new G4Element("W",        "W",  74.0, 183.84 *g/mole);
  Fe = new G4Element("Fe",       "Fe", 26.0, 55.845 *g/mole);
  O  = new G4Element("Oxygen",   "O",  8.,   16.00  *g/mole);
  Si = new G4Element("Silicon",  "Si", 14.,  28.08  *g/mole);
  C  = new G4Element("Carbon",   "C",  6.,   12.00  *g/mole);

  //Materials definition
  Xenon = new G4Material("XenonGas", 54., 131.29 * g/mole, 5.458 * mg/cm3,
                                     kStateGas, 93.15 * kelvin, 1 * atmosphere);


  // Absorber composition:  savannah.cern.ch/task/download.php?file_id=22925
  NiW = new G4Material("Tungsten/Nickel Composite",18.155*g/cm3,3);
  NiW->AddElement(W,  0.948);
  NiW->AddElement(Ni, 0.037);
  NiW->AddElement(Fe, 0.015);

  Steel = new G4Material("Steel", 7.9*gram/cm3,2);
  Steel->AddElement(Fe, 0.98);
  Steel->AddElement(C,  0.02);
  
  SilicaCore_UI = new G4Material("SilicaCore_UI", 2.200 * g/cm3, 2);
  SilicaCore_UI->AddElement(Si,1);
  SilicaCore_UI->AddElement(O, 2);

  SilicaClad_UI = new G4Material("SilicaClad_UI", 2.200 * g/cm3, 2);
  SilicaClad_UI->AddElement(Si,1);
  SilicaClad_UI->AddElement(O, 2);

  Kapton_UI = new G4Material("Kapton_UI", 1.42 * g/cm3, 4);
  Kapton_UI->AddElement(C,22);
  Kapton_UI->AddElement(H,10);
  Kapton_UI->AddElement(N,2);
  Kapton_UI->AddElement(O,5);
  
  pQuartz = new G4Material("Quartz", 2.200 * g/cm3, 2);
  pQuartz->AddElement(Si,1);
  pQuartz->AddElement(O, 2);

  EM_Quartz = new G4Material("EMQuartz", 2.200 * g/cm3, 2);
  EM_Quartz->AddElement(Si,1);
  EM_Quartz->AddElement(O, 2);

  nist_manager = G4NistManager::Instance();
  nist_manager->SetVerbose(1);

  Air   = nist_manager->FindOrBuildMaterial("G4_AIR");
  Al    = nist_manager->FindOrBuildMaterial("G4_Al");   	//G4_Al
  Cu    = nist_manager->FindOrBuildMaterial("G4_Cu");
  Pb    = nist_manager->FindOrBuildMaterial("G4_Pb");
  pureW = nist_manager->FindOrBuildMaterial("G4_W");

  Polyethylene = new G4Material("Polyethylene", 0.975 * g/cm3, 2);
  Polyethylene->AddElement(H,4);
  Polyethylene->AddElement(C,2);
  
  PMMA =  new G4Material("Polymethylmethacrylate", 1.18 * g/cm3, 3);
  PMMA->AddElement(C,5);
  PMMA->AddElement(H,8);
  PMMA->AddElement(O,2);

  Grease =  new G4Material("Grease", 1.0 * g/cm3, 3);
  Grease->AddElement(C,1);
  Grease->AddElement(H,1);
  Grease->AddElement(O,1);

}

void Materials::DefineOpticalProperties(void){




  // UIUC RPD Materials
  
  // Silica core optical properties
  // 2 eV == 620 nm
  // 3.5 eV == 354 nm
  // in nm: {650, 600, 500, 400, 300}

  G4double photonEnergy_UI[nEntries_UI] = {1.91*eV, 2.07*eV, 2.48*eV, 3.1*eV, 4.13*eV};  
  //    G4double photonEnergy_UI[nEntries_UI] = {1.91*eV, 2.07*eV, 2.48*eV, 3.1*eV, 3.54*eV};  
  // https://www.content.molex.com/dxdam/literature/987650-8936.pdf
  G4double silica_RIND_UI[nEntries_UI] = {1.6, 1.6, 1.6, 1.6, 1.6};
  G4double silica_ABSL_UI[nEntries_UI] = {5329.*cm, 3525.*cm, 2322.*cm, 885.*cm, 271.*cm};
  //  G4double silica_ABSL_UI[nEntries_UI] = {5329.*cm, 5329.*cm, 5329.*cm, 5329.*cm, 5329.*cm,};
  //G4double silica_ABSL_UI[nEntries_UI_small] = {300.00*cm, 1720.*cm};
  //G4double silica_ABSL_UI[nEntries_UI_small] = {5329.*cm, 271.*cm};
  
  MPT_Array.push_back(new G4MaterialPropertiesTable());
  MPT_Array.back()->AddProperty("RINDEX",photonEnergy_UI,silica_RIND_UI,nEntries_UI);//index of refraction
  MPT_Array.back()->AddProperty("ABSLENGTH",photonEnergy_UI,silica_ABSL_UI,nEntries_UI);//absorption length
  SilicaCore_UI->SetMaterialPropertiesTable(MPT_Array.back());


  // Silica cladding optical properties
  // Numerical aperture is given by data sheet as 0.22 and NA = sqrt( n1^2 - n2^2 ), so n2 = sqrt( n1^2 - NA^2 )
  G4double silica_clad_RIND_UI[nEntries_UI] = {1.585, 1.585, 1.585, 1.585, 1.585}; // Only refractive index changes

  //  G4double silica_clad_RIND_UI[nEntries_UI] = {1.4575, 1.4575, 1.4575, 1.4575, 1.4575}; // Only refractive index changes
 
  //  G4double silica_clad_RIND_UI[nEntries_UI] = {1., 1., 1., 1., 1.}; // Only refractive index changes
  MPT_Array.push_back(new G4MaterialPropertiesTable());
  MPT_Array.back()->AddProperty("RINDEX",photonEnergy_UI,silica_clad_RIND_UI,nEntries_UI);//index of refraction
  MPT_Array.back()->AddProperty("ABSLENGTH",photonEnergy_UI,silica_ABSL_UI,nEntries_UI);//absorption length
  SilicaClad_UI->SetMaterialPropertiesTable(MPT_Array.back());

  
  // Kapton optical properties:
  // https://engineering.case.edu/centers/sdle/sites/engineering.case.edu.centers.sdle/files/optical_properties_of_materials.pdf
  // https://www.sciencedirect.com/science/article/pii/S0257897205010881
  G4double kapton_RIND_UI[nEntries_UI] = {1.8, 1.82, 1.86, 1.95, 2.18};
  G4double kapton_ABSL_UI[nEntries_UI] = {0.011*cm, 0.0022*cm, 0.0056*cm, 0.04*cm, 0.063*cm};
  G4double kapton_REFL_UI[nEntries_UI] = {0.35, 0.2, 0.02, 0.03, 0.04};

  MPT_Array.push_back(new G4MaterialPropertiesTable());
  MPT_Array.back()->AddProperty("RINDEX",photonEnergy_UI,kapton_RIND_UI,nEntries_UI);
  MPT_Array.back()->AddProperty("ABSLENGTH",photonEnergy_UI,kapton_ABSL_UI,nEntries_UI);
  MPT_Array.back()->AddProperty("REFLECTIVITY",photonEnergy_UI,kapton_REFL_UI,nEntries_UI);
  Kapton_UI->SetMaterialPropertiesTable(MPT_Array.back());
  

  // https://refractiveindex.info/?shelf=main&book=W&page=Werner
  //  G4double NiW_RIND_UI[nEntries_UI] = {0.92, 1.05, 2.02, 1.49, 1.43};
  //  G4double NiW_ABSL_UI[nEntries_UI_small] = {5.8e-7*cm, 5.8e-7*cm};
  // https://www.tungsten.com/materials/tungsten/
  // only defining reflectivity since we only want optical light created in the quartz cavitiies
  G4double NiW_REFL_UI[nEntries_UI] = {0.62, 0.62, 0.62, 0.62, 0.62};
  MPT_Array.push_back(new G4MaterialPropertiesTable());
  MPT_Array.back()->AddProperty("REFLECTIVITY",photonEnergy_UI,NiW_REFL_UI,nEntries_UI);
  NiW->SetMaterialPropertiesTable(MPT_Array.back());

  // https://refractiveindex.info/?shelf=main&book=W&page=Werner
  //  G4double W_RIND_UI[nEntries_UI] = {0.92, 1.05, 2.02, 1.49, 1.43};
  //  G4double W_ABSL_UI[nEntries_UI_small] = {5.8e-7*cm, 5.8e-7*cm};
  // https://www.tungsten.com/materials/tungsten/
  // only defining reflectivity since we only want optical light created in the quartz cavitiies
  G4double W_REFL_UI[nEntries_UI] = {0.62, 0.62, 0.62, 0.62, 0.62};
  MPT_Array.push_back(new G4MaterialPropertiesTable());
  MPT_Array.back()->AddProperty("REFLECTIVITY",photonEnergy_UI,W_REFL_UI,nEntries_UI);
  pureW->SetMaterialPropertiesTable(MPT_Array.back());
  
  //Air optical properties
  G4double RefractiveIndexAir[nEntries_UI] = {1.0, 1.0, 1.0, 1.0, 1.0};
  MPT_Array.push_back(new G4MaterialPropertiesTable());
  MPT_Array.back()->AddProperty("RINDEX", photonEnergy_UI, RefractiveIndexAir, nEntries_UI);
  Air->SetMaterialPropertiesTable(MPT_Array.back());


  //Aluminum optical properties
  G4double Al_refl[nEntries_UI] = {0.89, 0.89, 0.89, 0.89, 0.89};//"Hard cut"
  MPT_Array.push_back(new G4MaterialPropertiesTable());
  MPT_Array.back()->AddProperty("REFLECTIVITY",photonEnergy_UI,Al_refl,nEntries_UI);
  Al->SetMaterialPropertiesTable(MPT_Array.back());
  
  
  /*
      //quartz to air surface
  MPT_Array.push_back(new G4MaterialPropertiesTable());
  // set the optical boundary properties
  G4OpticalSurface *opQuartzToAir = new G4OpticalSurface("opQuartzToAir",unified, ground, dielectric_metal, 1);
  opAlSurface->SetMaterialPropertiesTable(MPT_Array.back());
  */









  /*


     // Maryland RPD Materials
  
  const float grease_RI = 1.46;
  const float clad_RI = 1.49;
  const float core_RI = 1.6;
  //  const float tile_RI = 1.4585;

  G4double quartz_RIND[nEntries_UMD], quartz_ABSL[nEntries_UMD], PhotonEnergy[nEntries_UMD];
  //G4double quartz_RFLT[nEntries_UMD];
  G4double quartz_EFIC[nEntries_UMD];
     for(int i = 0; i < nEntries_UMD; i++){
         PhotonEnergy[i] = 2.00*eV + i*0.03*eV;
         quartz_RIND[i] = core_RI; //Refractive Index - constants
         quartz_ABSL[i] = 300.00*cm + i*20*cm; //Attenuation length

	 //         quartz_RFLT[i] = 0.5;
	 //         quartz_EFIC[i] = 0.5;

     }

   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("RINDEX",PhotonEnergy,quartz_RIND,nEntries_UMD);//index of refraction
   MPT_Array.back()->AddProperty("ABSLENGTH",PhotonEnergy,quartz_ABSL,nEntries_UMD);//absorption length
  // MPT_Array.back()->AddProperty("REFLECTIVITY",PhotonEnergy,quartz_RFLT,nEntries_UMD);//refelectivity
  // MPT_Array.back()->AddProperty("EFFICIENCY",PhotonEnergy,quartz_EFIC,nEntries_UMD);//efficiency
   pQuartz->SetMaterialPropertiesTable(MPT_Array.back());


   //Quartz Cladding optical properties
   G4double clad_RIND[nEntries_UMD]; // Only refractive index changes
   for(int i = 0; i < nEntries_UMD; i++){
     // Numerical aperture is given by data sheet as 0.22 and NA = sqrt( n1^2 - n2^2 ), so n2 = sqrt( n1^2 - NA^2 )
     clad_RIND[i] = sqrt( pow(core_RI,2.0) - pow(0.22,2.0) ); //Refractive Index
     //clad_RIND[i] = 1.0; //Refractive Index
   }

   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("RINDEX",PhotonEnergy,clad_RIND,nEntries_UMD);//index of refraction
   MPT_Array.back()->AddProperty("ABSLENGTH",PhotonEnergy,quartz_ABSL,nEntries_UMD);//absorption length
   // MPT_Array.back()->AddProperty("REFLECTIVITY",PhotonEnergy,quartz_RFLT,nEntries_UMD);//refelectivity
   // MPT_Array.back()->AddProperty("EFFICIENCY",PhotonEnergy,quartz_EFIC,nEntries_UMD);//efficiency
   fiberClad->SetMaterialPropertiesTable(MPT_Array.back());
   
   

   //Aluminium optical properties
   G4double AllPhotonEnergies[50] = {2.00*eV, 2.04*eV, 2.07*eV, 2.11*eV, 2.15*eV, 2.18*eV, 2.22*eV, 2.26*eV, 2.29*eV,
                                              2.33*eV, 2.37*eV, 2.40*eV, 2.44*eV, 2.48*eV, 2.51*eV, 2.55*eV, 2.59*eV, 2.62*eV,
                                              2.66*eV, 2.70*eV, 2.74*eV, 2.77*eV, 2.81*eV, 2.85*eV, 2.88*eV, 2.92*eV, 2.96*eV,
                                              2.99*eV, 3.03*eV, 3.07*eV, 3.10*eV, 3.14*eV, 3.18*eV, 3.21*eV, 3.25*eV, 3.29*eV,
                                              3.32*eV, 3.36*eV, 3.40*eV, 3.43*eV, 3.47*eV, 5.0*eV,  25.0*eV,100.0*eV, 1000.00*eV,
                                              10000.0*eV, 25000.0*eV, 50000.0*eV, 250000.*eV, 1000000.*eV};
   G4double Al_refl[50] = {0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89,
                                    0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89,
                                    0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89,
                                    0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89,
                                    0.89, 0.89, 0.89, 0.89, 0.79, 0.69, 0.49, 0.29, 0.19, 0.09};	//"Hard cut"
   G4double Al_ABSL[50] = {0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm,
                              0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm,
                              0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm,
                              0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm,
                              0.01*mm, 0.10*mm, 1.00*mm,0.01*m,  0.01*m,  0.10*m, 0.10*m,  1.00*m,  10.0*m, 10.0*m};

   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("REFLECTIVITY",AllPhotonEnergies,Al_refl,50);
   MPT_Array.back()->AddProperty("ABSLENGTH",AllPhotonEnergies,Al_ABSL,50);
   Al->SetMaterialPropertiesTable(MPT_Array.back());

   //Air optical properties
   G4double RefractiveIndexAir[nEntries_UMD];
   for (int i = 0; i < nEntries_UMD; i++) RefractiveIndexAir[i] = 1.0;//

   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("RINDEX", PhotonEnergy, RefractiveIndexAir, nEntries_UMD);
   Air->SetMaterialPropertiesTable(MPT_Array.back());

   //Maryland Polyethilene optical properties

   G4double RefractiveIndexClad1[nEntries_UMD];
   G4double AbsClad[nEntries_UMD];
   for (int i = 0; i < nEntries_UMD; i++){
     RefractiveIndexClad1[i] = clad_RI;//
     AbsClad[i] = 20.0*m;
   }

   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexClad1,nEntries_UMD);
   MPT_Array.back()->AddProperty("ABSLENGTH",PhotonEnergy,AbsClad,nEntries_UMD);
   Polyethylene->SetMaterialPropertiesTable(MPT_Array.back());

   //PMMA optical properties
   G4double RefractiveIndexWLSfiber[nEntries_UMD];
   for (int i = 0; i < nEntries_UMD; i++) RefractiveIndexWLSfiber[i] = core_RI;

   G4double AbsWLSfiber[50] = { 5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
                                         5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
                                         5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,1.10*m,
                                         1.10*m,1.10*m,1.10*m,1.10*m,1.10*m,1.10*m, 1.*mm, 1.*mm, 1.*mm, 1.*mm,
                                         1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm};
   G4double EmissionWLSfiber[50] = {0.05, 0.10, 0.30, 0.50, 0.75, 1.00, 1.50, 1.85, 2.30, 2.75,
                                             3.25, 3.80, 4.50, 5.20, 6.00, 7.00, 8.50, 9.50, 11.1, 12.4,
                                             12.9, 13.0, 12.8, 12.3, 11.1, 11.0, 12.0, 11.0, 17.0, 16.9,
                                             15.0, 9.00, 2.50, 1.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00,
                                             0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

   
   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("WLSABSLENGTH", PhotonEnergy, AbsWLSfiber, nEntries_UMD);
   MPT_Array.back()->AddProperty("WLSCOMPONENT", PhotonEnergy, EmissionWLSfiber, nEntries_UMD);
   MPT_Array.back()->AddProperty("RINDEX", PhotonEnergy, RefractiveIndexWLSfiber, nEntries_UMD);
   MPT_Array.back()->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
   PMMA->SetMaterialPropertiesTable(MPT_Array.back());

    //Kapton optical properties
    // OPTICAL PROPERTIES OF MATERIALS FOR CONCENTRATOR PHOTOVOLTAIC SYSTEMS
    // https://engineering.case.edu/centers/sdle/sites/engineering.case.edu.centers.sdle/files/optical_properties_of_materials.pdf


   G4double kapton_RIND[50] = {1.71, 1.72, 1.72, 1.72, 1.72, 1.73, 1.73, 1.73, 1.73, 1.74, 1.74, 1.74,
			       1.74, 1.74, 1.75, 1.75, 1.75, 1.76, 1.76, 1.77, 1.77, 1.77, 1.78, 1.78,
			       1.78, 1.79, 1.79, 1.80, 1.80, 1.81, 1.81, 1.82, 1.82, 1.82, 1.83, 1.83,
			       1.84, 1.84, 1.85, 1.85, 1.86, 1.86, 1.87, 1.87, 1.88, 1.89, 1.90, 1.90,
			       1.92, 1.93};

   // Optical Characterization of Commonly Used Thermal Control Paints in a Simulated GEO Environment
   // https://amostech.com/TechnicalPapers/2018/Poster/Bengtson.pdf
   G4double kapton_RFLT[50] = {0.44, 0.42, 0.40, 0.37, 0.35, 0.32, 0.29, 0.26, 0.22, 0.19, 0.15, 0.11, 0.09, 0.06, 0.04,
			       0.03, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
			       0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
			       0.01, 0.01, 0.01, 0.01, 0.01};

   G4double um = micrometer;

   // Synthesis and properties of soluble aromatic polyimides from novel 4,5-diazafluorene-containing dianhydride
   // https://pubs.rsc.org/fa/content/articlehtml/2018/ra/c7ra12101f
   G4double kapton_ABSL[50] = {65.39*um, 64.28*um, 62.72*um, 61.96*um, 60.49*um, 59.43*um, 58.31*um, 57.07*um,
			       56.10*um, 54.85*um, 53.47*um, 52.21*um, 51.09*um, 49.68*um, 48.15*um, 46.51*um,
			       45.24*um, 44.06*um, 43.02*um, 42.01*um, 40.87*um, 40.07*um, 39.16*um, 38.45*um,
			       37.64*um, 37.01*um, 36.26*um, 35.71*um, 35.05*um, 34.42*um, 33.84*um, 33.14*um,
			       32.40*um, 31.73*um, 30.64*um, 29.49*um, 28.41*um, 26.86*um, 25.22*um, 23.35*um,
			       21.15*um, 19.18*um, 16.93*um, 14.84*um, 12.67*um, 10.63*um, 9.07*um, 7.72*um,
			       6.60*um, 5.64*um};
   
    MPT_Array.push_back(new G4MaterialPropertiesTable());
    MPT_Array.back()->AddProperty("RINDEX", PhotonEnergy, kapton_RIND, nEntries_UMD);
    MPT_Array.back()->AddProperty("REFLECTIVITY",AllPhotonEnergies,kapton_RFLT,50);
    MPT_Array.back()->AddProperty("ABSLENGTH",AllPhotonEnergies,kapton_ABSL,50);
    Kapton->SetMaterialPropertiesTable(MPT_Array.back());

    //Grease (silicone) optical properties
    G4double RefractiveIndexGrease[50];
    for (int i = 0; i < nEntries_UMD; i++) RefractiveIndexGrease[i] = grease_RI;
    
    MPT_Array.push_back(new G4MaterialPropertiesTable());
    MPT_Array.back()->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexGrease,nEntries_UMD);
    MPT_Array.back()->AddProperty("ABSLENGTH",PhotonEnergy,AbsClad,nEntries_UMD);
    Grease->SetMaterialPropertiesTable(MPT_Array.back());


   // set the optical boundary properties

  AlSurface = new G4OpticalSurface("AlSurface",unified, polished, dielectric_metal, 1);

   G4MaterialPropertiesTable* AlSurfaceProperty = new G4MaterialPropertiesTable();
   //AlSurfaceProperty->AddProperty("RINDEX",AllPhotonEnergies,RefractiveIndexGrease,nEntries_UMD);
   AlSurfaceProperty->AddProperty("REFLECTIVITY",AllPhotonEnergies,Al_refl,nEntries_UMD);
   AlSurfaceProperty->AddProperty("EFFICIENCY",AllPhotonEnergies,quartz_EFIC,nEntries_UMD);
   AlSurface->SetMaterialPropertiesTable(AlSurfaceProperty);

  TileSurface = new G4OpticalSurface("TileSurface",unified, polished, dielectric_dielectric,1);

   G4MaterialPropertiesTable* TileSurfaceProperty = new G4MaterialPropertiesTable();
   TileSurfaceProperty->AddProperty("RINDEX",AllPhotonEnergies,quartz_RIND,nEntries_UMD);
   //TileSurfaceProperty->AddProperty("REFLECTIVITY",AllPhotonEnergies,reflectivity,nEntries_UMD);
   TileSurfaceProperty->AddProperty("EFFICIENCY",AllPhotonEnergies,quartz_EFIC,nEntries_UMD);
   TileSurface->SetMaterialPropertiesTable(TileSurfaceProperty);

   photonDetSurface = new G4OpticalSurface("PhotonDetSurface", glisur, ground, dielectric_metal,1);

   G4MaterialPropertiesTable* photonDetSurfaceProperty = new G4MaterialPropertiesTable();

   G4double p_mppc[] = {2.00*eV, 3.47*eV};
   //   G4double p_mppc[] = {2.00*eV, 4.13*eV};
   const G4int nbins = sizeof(p_mppc)/sizeof(G4double);
   G4double refl_mppc[] = {1, 1};
   assert(sizeof(refl_mppc) == sizeof(p_mppc));
   G4double effi_mppc[] = {1, 1};
   assert(sizeof(effi_mppc) == sizeof(p_mppc));

   photonDetSurfaceProperty->AddProperty("REFLECTIVITY",p_mppc,refl_mppc,nbins);
   photonDetSurfaceProperty->AddProperty("EFFICIENCY",p_mppc,effi_mppc,nbins);

   photonDetSurface->SetMaterialPropertiesTable(photonDetSurfaceProperty);
  */

}

Materials::~Materials(void)
{
  delete H;
  delete N;
  delete O;
  delete Si;
  delete C;

  delete fiberClad;
  delete Xenon;
  delete pQuartz;
  delete SilicaCore_UI;
  delete SilicaClad_UI;
  delete Kapton_UI;
  delete EM_Quartz;
  delete Air;
  delete Al;
  delete Cu;
  delete Pb;
  delete pureW;
  delete NiW;
  delete Steel;
  delete Polyethylene;
  delete PMMA;
  delete Kapton;
  delete Grease;

  for (unsigned int i = 0; i < MPT_Array.size(); i++)
      delete MPT_Array.at(i);

  MPT_Array.clear();
}
