#ifndef MATERIALS_HH_
#define MATERIALS_HH_


#include "G4Element.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include <vector>

class Materials
{
    public:
            /*! \brief Get the instance of Materials class.*/
            static Materials* getInstance(void);
            /*! \brief Default destructor.*/
            ~Materials(void);
            /*! \enable optical materials */
            void UseOpticalMaterials(bool _useOptical) { useOptical = _useOptical; }
            /*! \brief Define indexes of interest.*/
            void DefineOpticalProperties(void);

            //Elements
            /*  \brief Hydrogen */
            G4Element* H;
            /*  \brief Nitrogen */
            G4Element* N;
            /*  \brief Oxygen */
            G4Element* O;
            /*  \brief Silicon */
            G4Element* Si;
            /*  \brief Nickel */
            G4Element* Ni;
            /*  \brief Tungsten */
            G4Element* W;
            /*  \brief Iron */
            G4Element* Fe;
            /*  \brief Carbon */
            G4Element* C;

            //Materials
            /* \brief Xenon */
            G4Material* Xenon;
            /* \brief Quartz */
            G4Material* pQuartz;
              /* \brief SilicaClad */
            G4Material* SilicaClad_UI;
                /* \brief SilicaCore */
            G4Material* SilicaCore_UI;
            /* \brief Fiber cladding */
            G4Material* fiberClad;
            /* \brief EM module Quartz (?) */
            G4Material* EM_Quartz;
            /* \brief Air */
            G4Material* Air;
            /* \brief Al no optical */
            G4Material* Al;
            /* \brief Copper */
            G4Material* Cu;
            /* \brief Lead */
            G4Material* Pb;
            /* \brief Tungsten */
            G4Material* pureW;
            /* \brief Tungsten Nickel composite */
            G4Material* NiW;
            /* \brief Tungsten Nickel composite */
            G4Material* Steel;
            /* \brief Kapton_UI */
            G4Material* Kapton_UI;
              /* \brief Polyethilene_UMD */
            G4Material* Polyethylene;
            /* \brief Polymethylmethacrylate */
            G4Material* PMMA;
            /* \brief Kapton/Polyimide */
            G4Material* Kapton;
            /* \brief Grease */
            G4Material* Grease;
            /* \brief Al Optical Surface */
            G4OpticalSurface* AlSurface;
            /* \brief Tile Optical Surface */
            G4OpticalSurface* TileSurface;
            /* \brief PhotonDetect Optical Surface */
            G4OpticalSurface* photonDetSurface;
            G4OpticalSurface* OpAlSurface;

    private:
            /*! \brief Private constructor.*/
            Materials(void);
            /*! \brief Private class pointer.*/
            static Materials* materials;
            /*! \brief True if FULL physics activated.*/
            bool useOptical = false;
            /*! \brief Entries of WLS response parametrizations */
  //              static const G4int nEntries_UMD = 50;
                static const G4int nEntries_UMD = 51;
  //              static const G4int nEntries_UI = 30;
              static const G4int nEntries_UI = 5;
              static const G4int nEntries_UI_small = 2;
  //            static const G4int nEntries_UI = 19;
            /*! \brief G4NistManager pointer.*/
            G4NistManager* nist_manager;
            /*! \brief G4MaterialPropertiesTable array */
            std::vector <G4MaterialPropertiesTable*> MPT_Array;


};

#endif /* MATERIALS_HH_ */
