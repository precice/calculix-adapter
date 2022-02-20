/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#include "PreciceInterface.h"
#include <assert.h>
#include <stdlib.h>
#include "CCXHelpers.h"
#include "ConfigReader.h"
#include "precice/SolverInterfaceC.h"

void Precice_Setup(char *configFilename, char *participantName, SimulationData *sim)
{
  assert(sim != NULL);
  assert(configFilename != NULL);
  assert(participantName != NULL);

  printf("Setting up preCICE participant %s, using config file: %s\n", participantName, configFilename);
  fflush(stdout);

  // Read the YAML config file
  AdapterConfig adapterConfig;
  ConfigReader_Read(configFilename, participantName, &adapterConfig);

  assert(adapterConfig.interfaces != NULL);
  assert(adapterConfig.preciceConfigFilename != NULL);
  assert(adapterConfig.numInterfaces > 0);

  sim->numPreciceInterfaces = adapterConfig.numInterfaces;

  // Create the solver interface and configure it - Alex: Calculix is always a serial participant (MPI size 1, rank 0)
  precicec_createSolverInterface(participantName, adapterConfig.preciceConfigFilename, 0, 1);

  // Create interfaces as specified in the config file
  sim->preciceInterfaces = (struct PreciceInterface **) calloc(adapterConfig.numInterfaces, sizeof(PreciceInterface *));

  int i;
  for (i = 0; i < adapterConfig.numInterfaces; i++) {
    InterfaceConfig *config   = adapterConfig.interfaces + i;
    sim->preciceInterfaces[i] = malloc(sizeof(PreciceInterface));
    PreciceInterface_Create(sim->preciceInterfaces[i], sim, config);
  }

  // At this point we are done with the configuration
  AdapterConfig_Free(&adapterConfig);

  // Initialize variables needed for the coupling
  NNEW(sim->coupling_init_v, double, sim->mt * sim->nk);

  // Initialize preCICE
  sim->precice_dt = precicec_initialize();

  // Initialize coupling data
  Precice_InitializeData(sim);
}

void Precice_InitializeData(SimulationData *sim)
{
  printf("Initializing coupling data\n");
  fflush(stdout);

  Precice_WriteCouplingData(sim);
  precicec_initialize_data();
  Precice_ReadCouplingData(sim);
}

void Precice_AdjustSolverTimestep(SimulationData *sim)
{
  if (isSteadyStateSimulation(sim->nmethod)) {
    printf("Adjusting time step for steady-state step\n");
    fflush(stdout);

    // For steady-state simulations, we will always compute the converged steady-state solution in one coupling step
    *sim->theta  = 0;
    *sim->tper   = 1;
    *sim->dtheta = 1;

    // Set the solver time step to be the same as the coupling time step
    sim->solver_dt = sim->precice_dt;
  } else {
    printf("Adjusting time step for transient step\n");
    printf("precice_dt dtheta = %f, dtheta = %f, solver_dt = %f\n", sim->precice_dt / *sim->tper, *sim->dtheta, fmin(sim->precice_dt, *sim->dtheta * *sim->tper));
    fflush(stdout);

    // Compute the normalized time step used by CalculiX
    *sim->dtheta = fmin(sim->precice_dt / *sim->tper, *sim->dtheta);

    // Compute the non-normalized time step used by preCICE
    sim->solver_dt = (*sim->dtheta) * (*sim->tper);
  }
}

void Precice_Advance(SimulationData *sim)
{
  printf("Adapter calling advance()...\n");
  fflush(stdout);

  sim->precice_dt = precicec_advance(sim->solver_dt);
}

bool Precice_IsCouplingOngoing()
{
  return precicec_isCouplingOngoing();
}

bool Precice_IsReadCheckpointRequired()
{
  return precicec_isActionRequired("read-iteration-checkpoint");
}

bool Precice_IsWriteCheckpointRequired()
{
  return precicec_isActionRequired("write-iteration-checkpoint");
}

void Precice_FulfilledReadCheckpoint()
{
  precicec_markActionFulfilled("read-iteration-checkpoint");
}

void Precice_FulfilledWriteCheckpoint()
{
  precicec_markActionFulfilled("write-iteration-checkpoint");
}

void Precice_ReadIterationCheckpoint(SimulationData *sim, double *v)
{

  printf("Adapter reading checkpoint...\n");
  fflush(stdout);

  // Reload time
  *(sim->theta) = sim->coupling_init_theta;

  // Reload step size
  *(sim->dtheta) = sim->coupling_init_dtheta;

  // Reload solution vector v
  memcpy(v, sim->coupling_init_v, sizeof(double) * sim->mt * sim->nk);
}

void Precice_WriteIterationCheckpoint(SimulationData *sim, double *v)
{

  printf("Adapter writing checkpoint...\n");
  fflush(stdout);

  // Save time
  sim->coupling_init_theta = *(sim->theta);

  // Save step size
  sim->coupling_init_dtheta = *(sim->dtheta);

  // Save solution vector v
  memcpy(sim->coupling_init_v, v, sizeof(double) * sim->mt * sim->nk);
}

void Precice_ReadCouplingData(SimulationData *sim)
{

  printf("Adapter reading coupling data...\n");
  fflush(stdout);

  PreciceInterface **interfaces    = sim->preciceInterfaces;
  int                numInterfaces = sim->numPreciceInterfaces;
  int                i, j;

  if (precicec_isReadDataAvailable()) {
    for (i = 0; i < numInterfaces; i++) {

      for (j = 0; j < interfaces[i]->numReadData; j++) {

        switch (interfaces[i]->readData[j]) {
        case TEMPERATURE:
          // Read and set temperature BC
          if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
            precicec_readBlockScalarData(interfaces[i]->temperatureDataID, interfaces[i]->num2DNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->node2DScalarData);
            mapData2Dto3DScalar(interfaces[i]->node2DScalarData, interfaces[i]->mapping2D3D, interfaces[i]->numNodes, interfaces[i]->nodeScalarData);
          } else {
            precicec_readBlockScalarData(interfaces[i]->temperatureDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeScalarData);
          }
          setNodeTemperatures(interfaces[i]->nodeScalarData, interfaces[i]->numNodes, interfaces[i]->xbounIndices, sim->xboun);
          printf("Reading TEMPERATURE coupling data with ID '%d'. \n", interfaces[i]->temperatureDataID);
          break;
        case HEAT_FLUX:
          // Read and set heat flux BC
          precicec_readBlockScalarData(interfaces[i]->fluxDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData);
          setFaceFluxes(interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload);
          printf("Reading HEAT_FLUX coupling data with ID '%d'. \n", interfaces[i]->fluxDataID);
          break;
        case SINK_TEMPERATURE:
          // Read and set sink temperature in convective film BC
          precicec_readBlockScalarData(interfaces[i]->kDeltaTemperatureReadDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData);
          setFaceSinkTemperatures(interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload);
          printf("Reading SINK_TEMPERATURE coupling data with ID '%d'. \n", interfaces[i]->kDeltaTemperatureReadDataID);
          break;
        case HEAT_TRANSFER_COEFF:
          // Read and set heat transfer coefficient in convective film BC
          precicec_readBlockScalarData(interfaces[i]->kDeltaReadDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData);
          setFaceHeatTransferCoefficients(interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload);
          printf("Reading HEAT_TRANSFER_COEFF coupling data with ID '%d'. \n", interfaces[i]->kDeltaReadDataID);
          break;
        case FORCES:
          // Read and set forces as concentrated loads (Neumann BC)
          if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
            precicec_readBlockVectorData(interfaces[i]->forcesDataID, interfaces[i]->num2DNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->node2DVectorData);
            mapData2Dto3DVector(interfaces[i]->node2DVectorData, interfaces[i]->mapping2D3D, interfaces[i]->numNodes, interfaces[i]->nodeVectorData);
          } else {
            precicec_readBlockVectorData(interfaces[i]->forcesDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
          }
          setNodeForces(interfaces[i]->nodeVectorData, interfaces[i]->numNodes, interfaces[i]->dimCCX, interfaces[i]->xforcIndices, sim->xforc);
          printf("Reading FORCES coupling data with ID '%d'. \n", interfaces[i]->forcesDataID);
          break;
        case DISPLACEMENTS:
          // Read and set displacements as single point constraints (Dirichlet BC)
          if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
            precicec_readBlockVectorData(interfaces[i]->displacementsDataID, interfaces[i]->num2DNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->node2DVectorData);
            mapData2Dto3DVector(interfaces[i]->node2DVectorData, interfaces[i]->mapping2D3D, interfaces[i]->numNodes, interfaces[i]->nodeVectorData);
          } else {
            precicec_readBlockVectorData(interfaces[i]->displacementsDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
          }
          setNodeDisplacements(interfaces[i]->nodeVectorData, interfaces[i]->numNodes, interfaces[i]->dimCCX, interfaces[i]->xbounIndices, sim->xboun);
          printf("Reading DISPLACEMENTS coupling data with ID '%d'. \n", interfaces[i]->displacementsDataID);
          break;
        case DISPLACEMENTDELTAS:
          printf("DisplacementDeltas cannot be used as read data\n");
          fflush(stdout);
          exit(EXIT_FAILURE);
          break;
        case VELOCITIES:
          printf("Velocities cannot be used as read data\n");
          fflush(stdout);
          exit(EXIT_FAILURE);
          break;
        case POSITIONS:
          printf("Positions cannot be used as read data.\n");
          fflush(stdout);
          exit(EXIT_FAILURE);
          break;
        }
      }
    }
  }
}

void Precice_WriteCouplingData(SimulationData *sim)
{

  printf("Adapter writing coupling data...\n");
  fflush(stdout);

  PreciceInterface **interfaces    = sim->preciceInterfaces;
  int                numInterfaces = sim->numPreciceInterfaces;
  int                i, j;
  int                iset;

  if (precicec_isWriteDataRequired(sim->solver_dt) || precicec_isActionRequired("write-initial-data")) {
    for (i = 0; i < numInterfaces; i++) {
      // Prepare data
      double *KDelta = NULL;
      double *T      = NULL;
      for (j = 0; j < interfaces[i]->numWriteData; j++) {
        enum CouplingDataType type = interfaces[i]->writeData[j];
        if (type == SINK_TEMPERATURE || type == HEAT_TRANSFER_COEFF) {
          if (KDelta == NULL) {
            int iset = interfaces[i]->faceSetID + 1; // Adjust index before calling Fortran function
            KDelta   = malloc(interfaces[i]->numElements * sizeof(double));
            T        = malloc(interfaces[i]->numElements * sizeof(double));
            FORTRAN(getkdeltatemp, (sim->co,
                                    sim->ntmat_,
                                    sim->vold,
                                    sim->cocon,
                                    sim->ncocon,
                                    &iset,
                                    sim->istartset,
                                    sim->iendset,
                                    sim->ipkon,
                                    *sim->lakon,
                                    sim->kon,
                                    sim->ialset,
                                    sim->ielmat,
                                    sim->mi,
                                    KDelta,
                                    T));
          }
        }
      }

      // Write data
      for (j = 0; j < interfaces[i]->numWriteData; j++) {

        switch (interfaces[i]->writeData[j]) {
        case TEMPERATURE:
          getNodeTemperatures(interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->vold, sim->mt, interfaces[i]->nodeScalarData);
          if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
            setDoubleArrayZero(interfaces[i]->node2DScalarData, interfaces[i]->num2DNodes, 1);
            mapData3Dto2DScalar(interfaces[i]->nodeScalarData, interfaces[i]->mapping2D3D, interfaces[i]->numNodes, interfaces[i]->node2DScalarData);
            precicec_writeBlockScalarData(interfaces[i]->temperatureDataID, interfaces[i]->num2DNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->node2DScalarData);
          } else {
            precicec_writeBlockScalarData(interfaces[i]->temperatureDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeScalarData);
          }
          printf("Writing TEMPERATURE coupling data with ID '%d'. \n", interfaces[i]->temperatureDataID);
          break;
        case HEAT_FLUX:
          iset = interfaces[i]->faceSetID + 1; // Adjust index before calling Fortran function
          FORTRAN(getflux, (sim->co,
                            sim->ntmat_,
                            sim->vold,
                            sim->cocon,
                            sim->ncocon,
                            &iset,
                            sim->istartset,
                            sim->iendset,
                            sim->ipkon,
                            *sim->lakon,
                            sim->kon,
                            sim->ialset,
                            sim->ielmat,
                            sim->mi,
                            interfaces[i]->faceCenterData));
          precicec_writeBlockScalarData(interfaces[i]->fluxDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData);
          printf("Writing HEAT_FLUX coupling data with ID '%d'. \n", interfaces[i]->fluxDataID);
          break;
        case SINK_TEMPERATURE:
          precicec_writeBlockScalarData(interfaces[i]->kDeltaTemperatureWriteDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, T);
          printf("Writing SINK_TEMPERATURE coupling data with ID '%d'. \n", interfaces[i]->kDeltaTemperatureWriteDataID);
          break;
        case HEAT_TRANSFER_COEFF:
          precicec_writeBlockScalarData(interfaces[i]->kDeltaWriteDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, KDelta);
          printf("Writing HEAT_TRANSFER_COEFF coupling data with ID '%d'. \n", interfaces[i]->kDeltaWriteDataID);
          break;
        case DISPLACEMENTS:
          getNodeDisplacements(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dimCCX, sim->vold, sim->mt, interfaces[i]->nodeVectorData);
          if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
            setDoubleArrayZero(interfaces[i]->node2DVectorData, interfaces[i]->num2DNodes, interfaces[i]->dim);
            mapData3Dto2DVector(interfaces[i]->nodeVectorData, interfaces[i]->mapping2D3D, interfaces[i]->numNodes, interfaces[i]->node2DVectorData);
            precicec_writeBlockVectorData(interfaces[i]->displacementsDataID, interfaces[i]->num2DNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->node2DVectorData);
          } else {
            precicec_writeBlockVectorData(interfaces[i]->displacementsDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
          }
          printf("Writing DISPLACEMENTS coupling data with ID '%d'. \n", interfaces[i]->displacementsDataID);
          break;
        case DISPLACEMENTDELTAS:
          getNodeDisplacementDeltas(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dimCCX, sim->vold, sim->coupling_init_v, sim->mt, interfaces[i]->nodeVectorData);
          if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
            setDoubleArrayZero(interfaces[i]->node2DVectorData, interfaces[i]->num2DNodes, interfaces[i]->dim);
            mapData3Dto2DVector(interfaces[i]->nodeVectorData, interfaces[i]->mapping2D3D, interfaces[i]->numNodes, interfaces[i]->node2DVectorData);
            precicec_writeBlockVectorData(interfaces[i]->displacementDeltasDataID, interfaces[i]->num2DNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->node2DVectorData);
          } else {
            precicec_writeBlockVectorData(interfaces[i]->displacementDeltasDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
          }
          printf("Writing DISPLACEMENTDELTAS coupling data with ID '%d'. \n", interfaces[i]->displacementDeltasDataID);
          break;
        case VELOCITIES:
          getNodeVelocities(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dim, sim->veold, sim->mt, interfaces[i]->nodeVectorData);
          precicec_writeBlockVectorData(interfaces[i]->velocitiesDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
          printf("Writing VELOCITIES coupling data with ID '%d'. \n", interfaces[i]->velocitiesDataID);
          break;
        case POSITIONS:
          getNodeCoordinates(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dim, sim->co, sim->vold, sim->mt, interfaces[i]->nodeVectorData);
          precicec_writeBlockVectorData(interfaces[i]->positionsDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
          printf("Writing POSITIONS coupling data with ID '%d'. \n", interfaces[i]->positionsDataID);
          break;
        case FORCES:
          getNodeForces(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dim, sim->fn, sim->mt, interfaces[i]->nodeVectorData);
          precicec_writeBlockVectorData(interfaces[i]->forcesDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
          printf("Writing FORCES coupling data with ID '%d'. \n", interfaces[i]->forcesDataID);
          break;
        }
      }
      // Cleanup data
      free(T);
      free(KDelta);
    }
    if (precicec_isActionRequired("write-initial-data")) {
      precicec_markActionFulfilled("write-initial-data");
    }
  }
}

void Precice_FreeData(SimulationData *sim)
{
  int i;

  if (sim->coupling_init_v != NULL) {
    free(sim->coupling_init_v);
  }

  for (i = 0; i < sim->numPreciceInterfaces; i++) {
    PreciceInterface_FreeData(sim->preciceInterfaces[i]);
    free(sim->preciceInterfaces[i]);
  }

  free(sim->preciceInterfaces);
  precicec_finalize();
}

void PreciceInterface_Create(PreciceInterface *interface, SimulationData *sim, InterfaceConfig const *config)
{
  interface->dim = precicec_getDimensions();

  // Initialize pointers as NULL
  interface->elementIDs            = NULL;
  interface->faceIDs               = NULL;
  interface->faceCenterCoordinates = NULL;
  interface->preciceFaceCenterIDs  = NULL;
  interface->nodeCoordinates       = NULL;
  interface->node2DCoordinates     = NULL;
  interface->nodeIDs               = NULL;
  interface->mapping2D3D           = NULL;
  interface->preciceNodeIDs        = NULL;
  interface->nodeScalarData        = NULL;
  interface->node2DScalarData      = NULL;
  interface->nodeVectorData        = NULL;
  interface->node2DVectorData      = NULL;
  interface->faceCenterData        = NULL;
  interface->xbounIndices          = NULL;
  interface->xloadIndices          = NULL;
  interface->xforcIndices          = NULL;

  // Initialize data ids to -1
  interface->temperatureDataID            = -1;
  interface->fluxDataID                   = -1;
  interface->kDeltaWriteDataID            = -1;
  interface->kDeltaTemperatureWriteDataID = -1;
  interface->kDeltaReadDataID             = -1;
  interface->kDeltaTemperatureReadDataID  = -1;
  interface->displacementsDataID          = -1;
  interface->displacementDeltasDataID     = -1;
  interface->positionsDataID              = -1;
  interface->velocitiesDataID             = -1;
  interface->forcesDataID                 = -1;

  // Check if quasi 2D-3D coupling needs to be implemented
  if (interface->dim == 2) {
    interface->quasi2D3D = 1;
    printf("Using quasi 2D-3D coupling\n");
    interface->dimCCX = 3; // CalculiX always solves a 3D problem
  } else {
    interface->quasi2D3D = 0;
    interface->dimCCX    = interface->dim;
  }

  //Mapping Type
  // The patch identifies the set used as interface in Calculix
  interface->name = strdup(config->patchName);
  // Calculix needs to know if nearest-projection mapping is implemented. config->map = 1 is for nearest-projection, config->map = 0 is for everything else
  interface->mapNPType = config->map;

  // Nodes mesh
  interface->nodesMeshID   = -1;
  interface->nodesMeshName = NULL;
  if (config->nodesMeshName) {
    interface->nodesMeshName = strdup(config->nodesMeshName);
    PreciceInterface_ConfigureNodesMesh(interface, sim);
  }

  // Face centers mesh
  interface->faceCentersMeshID   = -1;
  interface->faceCentersMeshName = NULL;
  if (config->facesMeshName) {
    interface->faceCentersMeshName = strdup(config->facesMeshName);
    // Only configure a face center mesh if necesary; i.e. do not configure it for FSI simulations, also do not configure tetra faces if no face center mesh is used (as in FSI simulations)
    PreciceInterface_ConfigureFaceCentersMesh(interface, sim);
    // Triangles of the nodes mesh (needs to be called after the face centers mesh is configured!)
    PreciceInterface_ConfigureTetraFaces(interface, sim);
  }

  PreciceInterface_ConfigureCouplingData(interface, sim, config);
}

void PreciceInterface_ConfigureFaceCentersMesh(PreciceInterface *interface, SimulationData *sim)
{
  //printf("Entering ConfigureFaceCentersMesh \n");
  char *faceSetName      = toFaceSetName(interface->name);
  interface->faceSetID   = getSetID(faceSetName, sim->set, sim->nset);
  interface->numElements = getNumSetElements(interface->faceSetID, sim->istartset, sim->iendset);

  interface->elementIDs = malloc(interface->numElements * sizeof(ITG));
  interface->faceIDs    = malloc(interface->numElements * sizeof(ITG));
  getSurfaceElementsAndFaces(interface->faceSetID, sim->ialset, sim->istartset, sim->iendset, interface->elementIDs, interface->faceIDs);

  interface->faceCenterCoordinates = malloc(interface->numElements * 3 * sizeof(double));
  interface->preciceFaceCenterIDs  = malloc(interface->numElements * 3 * sizeof(int));
  getTetraFaceCenters(interface->elementIDs, interface->faceIDs, interface->numElements, sim->kon, sim->ipkon, sim->co, interface->faceCenterCoordinates);

  interface->faceCentersMeshID    = precicec_getMeshID(interface->faceCentersMeshName);
  interface->preciceFaceCenterIDs = malloc(interface->numElements * sizeof(int));

  precicec_setMeshVertices(interface->faceCentersMeshID, interface->numElements, interface->faceCenterCoordinates, interface->preciceFaceCenterIDs);
}

void PreciceInterface_ConfigureNodesMesh(PreciceInterface *interface, SimulationData *sim)
{
  //printf("Entering configureNodesMesh \n");
  char *nodeSetName    = toNodeSetName(interface->name);
  interface->nodeSetID = getSetID(nodeSetName, sim->set, sim->nset);
  interface->numNodes  = getNumSetElements(interface->nodeSetID, sim->istartset, sim->iendset);
  interface->nodeIDs   = &sim->ialset[sim->istartset[interface->nodeSetID] - 1]; //Lucia: make a copy

  interface->nodeCoordinates = malloc(interface->numNodes * interface->dimCCX * sizeof(double));
  getNodeCoordinates(interface->nodeIDs, interface->numNodes, interface->dimCCX, sim->co, sim->vold, sim->mt, interface->nodeCoordinates);

  // Extract 2D coordinates from 3D coordinates for quasi 2D-3D coupling
  if (isQuasi2D3D(interface->quasi2D3D)) {
    int count                    = 0;
    int dim                      = interface->dim;
    int dimCCX                   = interface->dimCCX;
    interface->num2DNodes        = interface->numNodes / 2;
    interface->node2DCoordinates = malloc(interface->num2DNodes * dim * sizeof(double));
    interface->mapping2D3D       = malloc(interface->numNodes * sizeof(int));
    for (int i = 0; i < interface->numNodes; i++) {
      // Filter out nodes which are in the XY plane (Z = 0) for getting 2D mesh
      if (isDoubleEqual(interface->nodeCoordinates[i * dimCCX + 2], 0.0)) {
        interface->node2DCoordinates[count * dim]     = interface->nodeCoordinates[i * dimCCX];
        interface->node2DCoordinates[count * dim + 1] = interface->nodeCoordinates[i * dimCCX + 1];
        count += 1;
      }
    }
    assert(count == interface->num2DNodes && "Filtering of 2D nodes not done properly. Please make sure the out-of-plane axis is Z-axis");

    count = 0;
    for (int i = 0; i < interface->numNodes; i++) {
      for (int ii = 0; ii < interface->numNodes; ii++) {
        // Compare each node with every other node to find nodes with matching X and Y coordinates
        if (isDoubleEqual(interface->nodeCoordinates[ii * dimCCX], interface->nodeCoordinates[i * dimCCX]) &&
            isDoubleEqual(interface->nodeCoordinates[ii * dimCCX + 1], interface->nodeCoordinates[i * dimCCX + 1]) &&
            !isDoubleEqual(interface->nodeCoordinates[ii * dimCCX + 2], interface->nodeCoordinates[i * dimCCX + 2])) {
          if (!isDoubleEqual(interface->nodeCoordinates[i * dimCCX + 2], 0.0)) {
            interface->mapping2D3D[i]  = count;
            interface->mapping2D3D[ii] = count;
            count += 1;
          }
        }
      }
    }
  }

  if (interface->nodesMeshName != NULL) {
    //printf("nodesMeshName is not null \n");
    interface->nodesMeshID = precicec_getMeshID(interface->nodesMeshName);
    if (isQuasi2D3D(interface->quasi2D3D)) {
      interface->preciceNodeIDs = malloc(interface->num2DNodes * sizeof(int));
      precicec_setMeshVertices(interface->nodesMeshID, interface->num2DNodes, interface->node2DCoordinates, interface->preciceNodeIDs);
    } else {
      interface->preciceNodeIDs = malloc(interface->numNodes * sizeof(int));
      precicec_setMeshVertices(interface->nodesMeshID, interface->numNodes, interface->nodeCoordinates, interface->preciceNodeIDs);
    }
  }

  if (interface->mapNPType == 1) {
    assert(!interface->quasi2D3D && "Quasi 2D - 3D configuration does not work for nearest-projection mapping");
    PreciceInterface_NodeConnectivity(interface, sim);
  }
}

void PreciceInterface_NodeConnectivity(PreciceInterface *interface, SimulationData *sim)
{
  int   numElements;
  char *faceSetName                = toFaceSetName(interface->name);
  interface->faceSetID             = getSetID(faceSetName, sim->set, sim->nset);
  numElements                      = getNumSetElements(interface->faceSetID, sim->istartset, sim->iendset);
  interface->elementIDs            = malloc(numElements * sizeof(ITG));
  interface->faceIDs               = malloc(numElements * sizeof(ITG));
  interface->faceCenterCoordinates = malloc(numElements * 3 * sizeof(double));
  getSurfaceElementsAndFaces(interface->faceSetID, sim->ialset, sim->istartset, sim->iendset, interface->elementIDs, interface->faceIDs);
  interface->numElements = numElements;
  PreciceInterface_ConfigureTetraFaces(interface, sim);
}

void PreciceInterface_EnsureValidNodesMeshID(PreciceInterface *interface, const char* type)
{
  if (interface->nodesMeshID < 0) {
    printf("Nodes mesh not provided in YAML config file. They are required for writing/reading the data %s.\n", type);
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
}

void PreciceInterface_ConfigureTetraFaces(PreciceInterface *interface, SimulationData *sim)
{
  int i;
  printf("Setting node connectivity for nearest projection mapping: \n");
  if (interface->nodesMeshName != NULL) {
    int *triangles = malloc(interface->numElements * 3 * sizeof(ITG));
    getTetraFaceNodes(interface->elementIDs, interface->faceIDs, interface->nodeIDs, interface->numElements, interface->numNodes, sim->kon, sim->ipkon, triangles);

    for (i = 0; i < interface->numElements; i++) {
      precicec_setMeshTriangleWithEdges(interface->nodesMeshID, triangles[3 * i], triangles[3 * i + 1], triangles[3 * i + 2]);
    }
    free(triangles);
  }
}

void PreciceInterface_EnsureValidFacesMeshID(PreciceInterface* interface, const char *type) {
  if (interface->faceCentersMeshID < 0) {
    printf("Face centers mesh not provided in YAML config file. They are required for writing/reading the data %s.\n", type);
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
}

void PreciceInterface_ConfigureCouplingData(PreciceInterface *interface, SimulationData *sim, InterfaceConfig const *config)
{
  interface->nodeScalarData = malloc(interface->numNodes * sizeof(double));
  interface->nodeVectorData = malloc(interface->numNodes * 3 * sizeof(double));

  if (isQuasi2D3D(interface->quasi2D3D)) {
    interface->node2DScalarData = malloc(interface->num2DNodes * sizeof(double));
    interface->node2DVectorData = malloc(interface->num2DNodes * 2 * sizeof(double));

    int dim = interface->dim;
    for (int i = 0; i < interface->num2DNodes; i++) {
      interface->node2DScalarData[i]           = 0.0;
      interface->node2DVectorData[i * dim]     = 0.0;
      interface->node2DVectorData[i * dim + 1] = 0.0;
    }
  }

  interface->faceCenterData = malloc(interface->numElements * sizeof(double));

  int i;
  interface->numReadData = config->numReadData;
  if (config->numReadData > 0)
    interface->readData = malloc(config->numReadData * sizeof(int));
  for (i = 0; i < config->numReadData; i++) {
    if (isEqual(config->readDataNames[i], "Temperature")) {
      PreciceInterface_EnsureValidNodesMeshID(interface, "Temperature");
      interface->readData[i]       = TEMPERATURE;
      interface->xbounIndices      = malloc(interface->numNodes * sizeof(int));
      interface->temperatureDataID = precicec_getDataID("Temperature", interface->nodesMeshID);
      getXbounIndices(interface->nodeIDs, interface->numNodes, sim->nboun, sim->ikboun, sim->ilboun, interface->xbounIndices, TEMPERATURE);
      printf("Read data '%s' found with ID # '%d'.\n", config->readDataNames[i], interface->temperatureDataID);
    } else if (isEqual(config->readDataNames[i], "Heat-Flux")) {
      interface->readData[i]  = HEAT_FLUX;
      interface->xloadIndices = malloc(interface->numElements * sizeof(int));
      getXloadIndices("DFLUX", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices);
      PreciceInterface_EnsureValidFacesMeshID(interface, "Heat-Flux");
      interface->fluxDataID = precicec_getDataID("Heat-Flux", interface->faceCentersMeshID);
      printf("Read data '%s' found with ID # '%d'.\n", config->readDataNames[i], interface->fluxDataID);
    } else if (startsWith(config->readDataNames[i], "Sink-Temperature-")) {
      interface->readData[i]  = SINK_TEMPERATURE;
      interface->xloadIndices = malloc(interface->numElements * sizeof(int));
      getXloadIndices("FILM", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices);
      PreciceInterface_EnsureValidFacesMeshID(interface, "Sink Temperature");
      interface->kDeltaTemperatureReadDataID = precicec_getDataID(config->readDataNames[i], interface->faceCentersMeshID);
      printf("Read data '%s' found with ID # '%d'.\n", config->readDataNames[i], interface->kDeltaTemperatureReadDataID);
    } else if (startsWith(config->readDataNames[i], "Heat-Transfer-Coefficient-")) {
      interface->readData[i]      = HEAT_TRANSFER_COEFF;
      PreciceInterface_EnsureValidFacesMeshID(interface, "Heat Transfer Coefficient");
      interface->kDeltaReadDataID = precicec_getDataID(config->readDataNames[i], interface->faceCentersMeshID);
      printf("Read data '%s' found with ID # '%d'.\n", config->readDataNames[i], interface->kDeltaReadDataID);
    } else if (startsWith(config->readDataNames[i], "Force")) {
      PreciceInterface_EnsureValidNodesMeshID(interface, "Force");
      interface->readData[i]  = FORCES;
      interface->xforcIndices = malloc(interface->numNodes * 3 * sizeof(int));
      interface->forcesDataID = precicec_getDataID(config->readDataNames[i], interface->nodesMeshID);
      getXforcIndices(interface->nodeIDs, interface->numNodes, sim->nforc, sim->ikforc, sim->ilforc, interface->xforcIndices);
      printf("Read data '%s' found with ID # '%d'.\n", config->readDataNames[i], interface->forcesDataID);
    } else if (startsWith(config->readDataNames[i], "Displacement")) {
      PreciceInterface_EnsureValidNodesMeshID(interface, "Displacement");
      interface->readData[i]         = DISPLACEMENTS;
      interface->xbounIndices        = malloc(interface->numNodes * 3 * sizeof(int));
      interface->displacementsDataID = precicec_getDataID(config->readDataNames[i], interface->nodesMeshID);
      getXbounIndices(interface->nodeIDs, interface->numNodes, sim->nboun, sim->ikboun, sim->ilboun, interface->xbounIndices, DISPLACEMENTS);
      printf("Read data '%s' found with ID # '%d'.\n", config->readDataNames[i], interface->displacementsDataID);
    } else {
      printf("ERROR: Read data '%s' does not exist!\n", config->readDataNames[i]);
      exit(EXIT_FAILURE);
    }
  }

  interface->numWriteData = config->numWriteData;
  if (config->numWriteData > 0)
    interface->writeData = malloc(config->numWriteData * sizeof(int));

  for (i = 0; i < config->numWriteData; i++) {
    if (isEqual(config->writeDataNames[i], "Temperature")) {
      PreciceInterface_EnsureValidNodesMeshID(interface, "Temperature");
      interface->writeData[i]      = TEMPERATURE;
      interface->temperatureDataID = precicec_getDataID("Temperature", interface->nodesMeshID);
      printf("Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i], interface->temperatureDataID);
    } else if (isEqual(config->writeDataNames[i], "Heat-Flux")) {
      interface->writeData[i] = HEAT_FLUX;
      PreciceInterface_EnsureValidFacesMeshID(interface, "Heat Flux");
      interface->fluxDataID   = precicec_getDataID("Heat-Flux", interface->faceCentersMeshID);
      printf("Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i], interface->fluxDataID);
    } else if (startsWith(config->writeDataNames[i], "Sink-Temperature-")) {
      interface->writeData[i]                 = SINK_TEMPERATURE;
      PreciceInterface_EnsureValidFacesMeshID(interface, "Sink temperature");
      interface->kDeltaTemperatureWriteDataID = precicec_getDataID(config->writeDataNames[i], interface->faceCentersMeshID);
      printf("Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i], interface->kDeltaTemperatureWriteDataID);
    } else if (startsWith(config->writeDataNames[i], "Heat-Transfer-Coefficient-")) {
      interface->writeData[i]      = HEAT_TRANSFER_COEFF;
      PreciceInterface_EnsureValidFacesMeshID(interface, "Heat Transfer Coefficient");
      interface->kDeltaWriteDataID = precicec_getDataID(config->writeDataNames[i], interface->faceCentersMeshID);
      printf("Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i], interface->kDeltaWriteDataID);
    } else if (startsWith(config->writeDataNames[i], "DisplacementDelta")) {
      PreciceInterface_EnsureValidNodesMeshID(interface, "DisplacementDeltas");
      interface->writeData[i]             = DISPLACEMENTDELTAS;
      interface->displacementDeltasDataID = precicec_getDataID(config->writeDataNames[i], interface->nodesMeshID);
      printf("Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i], interface->displacementDeltasDataID);
    } else if (startsWith(config->writeDataNames[i], "Displacement")) {
      PreciceInterface_EnsureValidNodesMeshID(interface, "Displacement");
      interface->writeData[i]        = DISPLACEMENTS;
      interface->displacementsDataID = precicec_getDataID(config->writeDataNames[i], interface->nodesMeshID);
      printf("Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i], interface->displacementsDataID);
    } else if (startsWith(config->writeDataNames[i], "Position")) {
      PreciceInterface_EnsureValidNodesMeshID(interface, "Position");
      interface->writeData[i]    = POSITIONS;
      interface->positionsDataID = precicec_getDataID(config->writeDataNames[i], interface->nodesMeshID);
      printf("Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i], interface->positionsDataID);
    } else if (startsWith(config->writeDataNames[i], "Velocity")) {
      PreciceInterface_EnsureValidNodesMeshID(interface, "Velocity");
      interface->writeData[i]     = VELOCITIES;
      interface->velocitiesDataID = precicec_getDataID(config->writeDataNames[i], interface->nodesMeshID);
      printf("Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i], interface->velocitiesDataID);
    } else if (startsWith(config->writeDataNames[i], "Force")) {
      PreciceInterface_EnsureValidNodesMeshID(interface, "Force");
      interface->writeData[i] = FORCES;
      interface->forcesDataID = precicec_getDataID(config->writeDataNames[i], interface->nodesMeshID);
      printf("Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i], interface->forcesDataID);
    } else {
      printf("ERROR: Write data '%s' does not exist!\n", config->writeDataNames[i]);
      exit(EXIT_FAILURE);
    }
  }
}

void PreciceInterface_FreeData(PreciceInterface *preciceInterface)
{
  if (preciceInterface->readData != NULL) {
    free(preciceInterface->readData);
  }

  if (preciceInterface->writeData != NULL) {
    free(preciceInterface->writeData);
  }

  if (preciceInterface->elementIDs != NULL) {
    free(preciceInterface->elementIDs);
  }

  if (preciceInterface->faceIDs != NULL) {
    free(preciceInterface->faceIDs);
  }

  if (preciceInterface->faceCenterCoordinates != NULL) {
    free(preciceInterface->faceCenterCoordinates);
  }

  if (preciceInterface->preciceFaceCenterIDs != NULL) {
    free(preciceInterface->preciceFaceCenterIDs);
  }

  if (preciceInterface->nodeCoordinates != NULL) {
    free(preciceInterface->nodeCoordinates);
  }

  if (preciceInterface->node2DCoordinates != NULL) {
    free(preciceInterface->node2DCoordinates);
  }

  if (preciceInterface->preciceNodeIDs != NULL) {
    free(preciceInterface->preciceNodeIDs);
  }

  if (preciceInterface->nodeScalarData != NULL) {
    free(preciceInterface->nodeScalarData);
  }

  if (preciceInterface->node2DScalarData != NULL) {
    free(preciceInterface->node2DScalarData);
  }

  if (preciceInterface->nodeVectorData != NULL) {
    free(preciceInterface->nodeVectorData);
  }

  if (preciceInterface->node2DVectorData != NULL) {
    free(preciceInterface->node2DVectorData);
  }

  if (preciceInterface->faceCenterData != NULL) {
    free(preciceInterface->faceCenterData);
  }

  if (preciceInterface->xbounIndices != NULL) {
    free(preciceInterface->xbounIndices);
  }

  if (preciceInterface->xloadIndices != NULL) {
    free(preciceInterface->xloadIndices);
  }

  if (preciceInterface->xforcIndices != NULL) {
    free(preciceInterface->xforcIndices);
  }

  if (preciceInterface->mapping2D3D != NULL) {
    free(preciceInterface->mapping2D3D);
  }
}
