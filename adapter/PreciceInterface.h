/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#ifndef PRECICEINTERFACE_H
#define PRECICEINTERFACE_H

#include <string.h>
#include "ConfigReader.h"
#include "CCXHelpers.h"

/*
 * PreciceInterface: Structure with all the information of a coupled surface
 * Includes data regarding the surface mesh(es) and the coupling data
 */
typedef struct PreciceInterface {

	char * name;
	int dim;
	int dimCCX;

	// Interface nodes
	int numNodes;
	int num2DNodes; // Nodes in a single plane in case of quasi 2D-3D coupling
	int * nodeIDs;
	int * quasiMapping; // Node IDs to filter out 2D place in quasi 2D-3D coupling
	double * nodeCoordinates;
	double * node2DCoordinates; // 2D coordinates for quasi 2D-3D coupling
	int nodeSetID;
	int * preciceNodeIDs;
	int nodesMeshID;
	char * nodesMeshName;

	// Interface face elements
	int numElements;
	int * elementIDs;
	int * faceIDs;
	double * faceCenterCoordinates;
	int faceSetID;
	int faceCentersMeshID;
	char * faceCentersMeshName;
	int * preciceFaceCenterIDs;
	int * triangles;

	// Arrays to store the coupling data
	double * nodeScalarData;
	double * nodeVectorData; //Forces, displacements, velocities, positions and displacementDeltas are vector quantities
	double * node2DVectorData; // Vector quantities in 2D in case quasi 2D-3D coupling is done
	double * faceCenterData;

	// preCICE Data IDs
	int temperatureDataID;
	int fluxDataID;
	int kDeltaWriteDataID;
	int kDeltaTemperatureWriteDataID;
	int kDeltaReadDataID;
	int kDeltaTemperatureReadDataID;
	int displacementsDataID; //New data ID for displacements
	int displacementDeltasDataID; //New data ID for displacementDeltas
	int positionsDataID; //New data ID for positions
	int velocitiesDataID; //New data ID for velocities
	int forcesDataID; //New data ID for forces

	// Indices that indicate where to apply the boundary conditions / forces
	int * xloadIndices;
	int * xbounIndices;
	int * xforcIndices;

	// Mapping type if nearest-projection mapping
	int mapNPType;

	// Indicates if pseudo 2D-3D coupling is implemented
	int quasi2D3D;

	int numReadData;
	int numWriteData;
	enum CouplingDataType *readData;
	enum CouplingDataType *writeData;

} PreciceInterface;

/*
 * SimulationData: Structure with all the CalculiX variables
 * that need to be accessed by the adapter in order to do the coupling.
 * A list of variables and their meaning is available in the documentation
 * ccx_2.10.pdf (page 518)
 */
typedef struct SimulationData {

	// CalculiX data
	ITG * ialset;
	ITG * ielmat;
	ITG * istartset;
	ITG * iendset;
	char ** lakon;
	ITG * kon;
	ITG * ipkon;
	ITG nset;
	char * set;
	double * co;
	ITG nboun;
	ITG nforc; //total number of forces
	ITG * ikboun;
	ITG * ikforc; //the DoFs are all stored here in an array in numerical order
	ITG * ilboun;
	ITG * ilforc; //number of the force is stored here
	ITG * nelemload;
	int nload;
	char * sideload;
	double nk;
	ITG mt;
	double * theta;
	double * dtheta;
	double * tper;
	ITG * nmethod;
	double * xload;
	double * xforc; //scalar value of the force in one direction
	double * xboun;
	ITG * ntmat_;
	double * vold;
	double * veold;
	double * fn;//values of forces read from calculix
	double * cocon;
	ITG * ncocon;
	ITG * mi;

	// Interfaces
	int numPreciceInterfaces;
	PreciceInterface ** preciceInterfaces;

	// Coupling data
	double * coupling_init_v;
	double coupling_init_theta;
	double coupling_init_dtheta;
	double precice_dt;
	double solver_dt;

} SimulationData;



/**
 * @brief Configures and initializes preCICE and the interfaces
 * @param configFilename: YAML config file
 * @param participantName
 * @param sim
 * @param preciceInterfaces
 * @param numPreciceInterfaces
 */
void Precice_Setup( char * configFilename, char * participantName, SimulationData * sim );

/**
 * @brief Initializes the coupling data (does an initial exchange) if necessary
 * @param sim
 * @param preciceInterfaces
 * @param numInterfaces
 */
void Precice_InitializeData( SimulationData * sim );

/**
 * @brief Adjusts the solver time step based on the coupling time step and the solver time step
 * @param sim
 */
void Precice_AdjustSolverTimestep( SimulationData * sim );

/**
 * @brief Advances the coupling
 * @param sim
 */
void Precice_Advance( SimulationData * sim );

/**
 * @brief Returns true if coupling is still ongoing
 * @return
 */
bool Precice_IsCouplingOngoing();

/**
 * @brief Returns true if checkpoint must be read
 * @return
 */
bool Precice_IsReadCheckpointRequired();

/**
 * @brief Returns true if checkpoint must be written
 * @return
 */
bool Precice_IsWriteCheckpointRequired();

/**
 * @brief Tells preCICE that the checkpoint has been read
 */
void Precice_FulfilledReadCheckpoint();

/**
 * @brief Tells preCICE that the checkpoint has been written
 */
void Precice_FulfilledWriteCheckpoint();

/**
 * @brief Reads iteration checkpoint
 * @param sim: Structure with CalculiX data
 * @param v: CalculiX array with the temperature and displacement values
 */
void Precice_ReadIterationCheckpoint( SimulationData * sim, double * v );

/**
 * @brief Writes iteration checkpoint
 * @param sim: Structure with CalculiX data
 * @param v: CalculiX array with the temperature and displacement values
 */
void Precice_WriteIterationCheckpoint( SimulationData * sim, double * v );

/**
 * @brief Reads the coupling data for all interfaces
 * @param sim
 * @param preciceInterfaces
 * @param numInterfaces
 */
void Precice_ReadCouplingData( SimulationData * sim );

/**
 * @brief Writes the coupling data of all interfaces
 * @param sim
 * @param preciceInterfaces
 * @param numInterfaces
 */
void Precice_WriteCouplingData( SimulationData * sim );

/**
 * @brief Frees the memory
 * @param sim
 * @param preciceInterfaces
 * @param numInterfaces
 */
void Precice_FreeData( SimulationData * sim );

/**
 * @brief Creates an interface that is coupled with preCICE
 * @param interface
 * @param sim
 * @param config
 */
void PreciceInterface_Create( PreciceInterface * interface, SimulationData * sim, InterfaceConfig const * config );

/**
 * @brief Configures the face centers mesh and calls setMeshVertices on preCICE
 * @param interface
 * @param sim
 */
void PreciceInterface_ConfigureFaceCentersMesh( PreciceInterface * interface, SimulationData * sim );

/**
 * @brief Configures the nodes mesh
 * @param interface
 * @param sim: Structure with CalculiX data
 */
void PreciceInterface_ConfigureNodesMesh( PreciceInterface * interface, SimulationData * sim );

/**
 * @brief Terminate execution if the nodes mesh ID is not valid
 * @param interface
 */
void PreciceInterface_EnsureValidNodesMeshID( PreciceInterface * interface );

/**
 * @brief Configures the faces mesh (for tetrahedral elements only)
 * @param interface
 * @param sim
 */
void PreciceInterface_ConfigureTetraFaces( PreciceInterface * interface, SimulationData * sim );

/**
 * @brief Configures the node connectivity for nearest-projection mapping
 * @param interface
 * @param sim
 * @param config
 */
void PreciceInterface_NodeConnectivity( PreciceInterface * interface, SimulationData * sim );

/**
 * @brief Configures the coupling data for CHT or FSI
 * @param interface
 * @param sim
 * @param config
 */
void PreciceInterface_ConfigureCouplingData( PreciceInterface * interface, SimulationData * sim, InterfaceConfig const * config );

/**
 * @brief Frees the memory
 * @param preciceInterface
 */
void PreciceInterface_FreeData( PreciceInterface * preciceInterface );

#endif // PRECICEINTERFACE_H
