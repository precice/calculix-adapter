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
#include "precice/preciceC.h"

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
  precicec_createParticipant(participantName, adapterConfig.preciceConfigFilename, 0, 1);

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
  if (precicec_requiresInitialData()) {
    Precice_WriteCouplingData(sim);
  }

  // Initialize preCICE
  precicec_initialize();

  // Initialize coupling data
  printf("Initializing coupling data\n");
  fflush(stdout);
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
    sim->solver_dt = precicec_getMaxTimeStepSize();
  } else {
    double precice_dt = precicec_getMaxTimeStepSize();

    // Compute the time step size of CalculiX
    double solver_dt = (*sim->dtheta) * (*sim->tper);

    // Synchronize CalculiX time step with preCICE time window end
    double dt = fmin(precice_dt, solver_dt);

    // Normalize the agreed-on time step size
    double new_dtheta = dt / (*sim->tper);

    printf("Adjusting time step for transient step\n");
    printf("precice_dt = %f, ccx_dt = %f (dtheta = %f, tper = %f) -> dt = %f (dtheta = %f)\n", precice_dt, solver_dt, *sim->dtheta, *sim->tper, dt, new_dtheta);
    fflush(stdout);

    // Update dt and dtheta
    *sim->dtheta   = new_dtheta;
    sim->solver_dt = dt;
  }
}

void Precice_Advance(SimulationData *sim)
{
  printf("Adapter calling advance()...\n");
  fflush(stdout);

  precicec_advance(sim->solver_dt);
}

bool Precice_IsCouplingOngoing()
{
  return precicec_isCouplingOngoing();
}

bool Precice_requiresReadingCheckpoint()
{
  return precicec_requiresReadingCheckpoint();
}

bool Precice_requiresWritingCheckpoint()
{
  return precicec_requiresWritingCheckpoint();
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

void Precice_ReadIterationCheckpointModal(SimulationData *sim, double *dofs, double *derivatives, int nev)
{

  printf("Adapter reading checkpoint...\n");
  fflush(stdout);

  // Reload time
  *(sim->theta) = sim->coupling_init_theta;

  // Reload step size
  *(sim->dtheta) = sim->coupling_init_dtheta;

  // Reload DOFs in eigenmodes space & counters
  memcpy(dofs, sim->eigenDOFs, sizeof(double) * nev);
  memcpy(derivatives, sim->eigenDOFsDerivatives, sizeof(double) * nev);
}

void Precice_WriteIterationCheckpointModal(SimulationData *sim, const double *dofs, const double *derivatives, int nev)
{

  printf("Adapter writing checkpoint...\n");
  fflush(stdout);

  // Save time
  sim->coupling_init_theta = *(sim->theta);

  // Save step size
  sim->coupling_init_dtheta = *(sim->dtheta);

  // Save DOFs in eigenmodes space & counters
  memcpy(sim->eigenDOFs, dofs, sizeof(double) * nev);
  memcpy(sim->eigenDOFsDerivatives, derivatives, sizeof(double) * nev);
}

void Precice_ReadCouplingData(SimulationData *sim)
{

  printf("Adapter reading coupling data...\n");
  fflush(stdout);

  PreciceInterface **interfaces    = sim->preciceInterfaces;
  int                numInterfaces = sim->numPreciceInterfaces;
  int                i, j;

  for (i = 0; i < numInterfaces; i++) {

    for (j = 0; j < interfaces[i]->numReadData; j++) {

      switch (interfaces[i]->readData[j]) {
      case TEMPERATURE:
        // Read and set temperature BC
        if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
          consistentScalarRead(interfaces[i]->mappingQuasi2D3D, interfaces[i]->couplingMeshName, interfaces[i]->temperature, sim->solver_dt);
          setNodeTemperatures(interfaces[i]->mappingQuasi2D3D->bufferScalar3D, interfaces[i]->numNodes, interfaces[i]->xbounIndices, sim->xboun);
        } else {
          precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->temperature, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, sim->solver_dt, interfaces[i]->nodeScalarData);
          setNodeTemperatures(interfaces[i]->nodeScalarData, interfaces[i]->numNodes, interfaces[i]->xbounIndices, sim->xboun);
        }
        printf("Reading TEMPERATURE coupling data.\n");
        break;
      case HEAT_FLUX:
        // Read and set heat flux BC
        // Not working in 2D-3D now
        precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->flux, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, sim->solver_dt, interfaces[i]->faceCenterData);
        setFaceFluxes(interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload);
        printf("Reading HEAT_FLUX coupling data.\n");
        break;
      case SINK_TEMPERATURE:
        // Read and set sink temperature in convective film BC
        // Not working in 2D-3D now

        precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->kDeltaTemperatureRead, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, sim->solver_dt, interfaces[i]->faceCenterData);
        setFaceSinkTemperatures(interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload);
        printf("Reading SINK_TEMPERATURE coupling data.\n");
        break;
      case HEAT_TRANSFER_COEFF:
        // Read and set heat transfer coefficient in convective film BC
        // Not working in 2D-3D now

        precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->kDeltaRead, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, sim->solver_dt, interfaces[i]->faceCenterData);
        setFaceHeatTransferCoefficients(interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload);
        printf("Reading HEAT_TRANSFER_COEFF coupling data.\n");
        break;
      case FORCES:
        // Read and set forces as concentrated loads (Neumann BC)
        if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
          conservativeVectorRead(interfaces[i]->mappingQuasi2D3D, interfaces[i]->couplingMeshName, interfaces[i]->forces, sim->solver_dt);
          setNodeForces(interfaces[i]->mappingQuasi2D3D->bufferVector3D, interfaces[i]->numNodes, interfaces[i]->dimCCX, interfaces[i]->xforcIndices, sim->xforc);
        } else {
          precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->forces, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, sim->solver_dt, interfaces[i]->nodeVectorData);
          setNodeForces(interfaces[i]->nodeVectorData, interfaces[i]->numNodes, interfaces[i]->dimCCX, interfaces[i]->xforcIndices, sim->xforc);
        }
        printf("Reading FORCES coupling data.\n");
        break;
      case PRESSURE:
        precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->pressure, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, sim->solver_dt, interfaces[i]->faceCenterData);
        setFacePressure(interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload);
        printf("Reading PRESSURE coupling data.\n");
        break;
      case DISPLACEMENTS:
        // Read and set displacements as single point constraints (Dirichlet BC)
        if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
          conservativeVectorRead(interfaces[i]->mappingQuasi2D3D, interfaces[i]->couplingMeshName, interfaces[i]->displacements, sim->solver_dt);
          setNodeDisplacements(interfaces[i]->mappingQuasi2D3D->bufferVector3D, interfaces[i]->numNodes, interfaces[i]->dimCCX, interfaces[i]->xbounIndices, sim->xboun);
        } else {
          precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->displacements, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, sim->solver_dt, interfaces[i]->nodeVectorData);
          setNodeDisplacements(interfaces[i]->nodeVectorData, interfaces[i]->numNodes, interfaces[i]->dimCCX, interfaces[i]->xbounIndices, sim->xboun);
        }
        printf("Reading DISPLACEMENTS coupling data.\n");
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

void Precice_WriteCouplingData(SimulationData *sim)
{

  printf("Adapter writing coupling data...\n");
  fflush(stdout);

  PreciceInterface **interfaces    = sim->preciceInterfaces;
  int                numInterfaces = sim->numPreciceInterfaces;
  int                i, j;
  int                iset;

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
                                  sim->lakon,
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
        if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
          getNodeTemperatures(interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->vold, sim->mt, interfaces[i]->mappingQuasi2D3D->bufferScalar3D);
          consistentScalarWrite(interfaces[i]->mappingQuasi2D3D, interfaces[i]->couplingMeshName, interfaces[i]->temperature);
        } else {
          getNodeTemperatures(interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->vold, sim->mt, interfaces[i]->nodeScalarData);
          precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->temperature, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeScalarData);
        }
        printf("Writing TEMPERATURE coupling data.\n");
        break;
      case HEAT_FLUX:
        // Not implemented: 2D-3D
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
                          sim->lakon,
                          sim->kon,
                          sim->ialset,
                          sim->ielmat,
                          sim->mi,
                          interfaces[i]->faceCenterData));
        precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->flux, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData);
        printf("Writing HEAT_FLUX coupling data.\n");
        break;
      case SINK_TEMPERATURE:
        // Not implemented: 2D-3D

        precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->kDeltaTemperatureWrite, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, T);
        printf("Writing SINK_TEMPERATURE coupling data.\n");
        break;
      case HEAT_TRANSFER_COEFF:
        // Not implemented: 2D-3D

        precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->kDeltaWrite, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, KDelta);
        printf("Writing HEAT_TRANSFER_COEFF coupling data.\n");
        break;
      case DISPLACEMENTS:
        if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
          getNodeDisplacements(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dimCCX, sim->vold, sim->mt, interfaces[i]->mappingQuasi2D3D->bufferVector3D);
          consistentVectorWrite(interfaces[i]->mappingQuasi2D3D, interfaces[i]->couplingMeshName, interfaces[i]->displacements);
        } else {
          getNodeDisplacements(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dimCCX, sim->vold, sim->mt, interfaces[i]->nodeVectorData);
          precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->displacements, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
        }
        printf("Writing DISPLACEMENTS coupling data.\n");
        break;
      case DISPLACEMENTDELTAS:
        if (isQuasi2D3D(interfaces[i]->quasi2D3D)) {
          getNodeDisplacementDeltas(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dimCCX, sim->vold, sim->coupling_init_v, sim->mt, interfaces[i]->mappingQuasi2D3D->bufferVector3D);
          consistentVectorWrite(interfaces[i]->mappingQuasi2D3D, interfaces[i]->couplingMeshName, interfaces[i]->displacementDeltas);
        } else {
          getNodeDisplacementDeltas(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dimCCX, sim->vold, sim->coupling_init_v, sim->mt, interfaces[i]->nodeVectorData);
          precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->displacementDeltas, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
        }
        printf("Writing DISPLACEMENTDELTAS coupling data.\n");
        break;
      case VELOCITIES:
        getNodeVelocities(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dim, sim->veold, sim->mt, interfaces[i]->nodeVectorData);
        precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->velocities, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
        printf("Writing VELOCITIES coupling data.\n");
        break;
      case POSITIONS:
        getNodeCoordinates(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dim, sim->co, sim->vold, sim->mt, interfaces[i]->nodeVectorData);
        precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->positions, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
        printf("Writing POSITIONS coupling data.\n");
        break;
      case FORCES:
        getNodeForces(interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dim, sim->fn, sim->mt, interfaces[i]->nodeVectorData);
        precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->forces, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData);
        printf("Writing FORCES coupling data.\n");
        break;
      }
    }
    // Cleanup data
    free(T);
    free(KDelta);
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

  // Clean up checkpointing buffers
  if (sim->eigenDOFs != NULL) {
    free(sim->eigenDOFs);
  }
  if (sim->eigenDOFsDerivatives != NULL) {
    free(sim->eigenDOFsDerivatives);
  }

  precicec_finalize();
}

void PreciceInterface_Create(PreciceInterface *interface, SimulationData *sim, InterfaceConfig const *config)
{
  if (config->nodesMeshName) {
    interface->dim = precicec_getMeshDimensions(config->nodesMeshName);
  } else {
    if (config->facesMeshName) {
      interface->dim = precicec_getMeshDimensions(config->facesMeshName);
    } else {
      printf("ERROR: You need to define either a face or a nodes mesh. Check the adapter configuration file.\n");
      exit(EXIT_FAILURE);
    }
  }

  // Initialize pointers as NULL
  interface->elementIDs            = NULL;
  interface->faceIDs               = NULL;
  interface->faceCenterCoordinates = NULL;
  interface->preciceFaceCenterIDs  = NULL;
  interface->nodeCoordinates       = NULL;
  interface->nodeIDs               = NULL;
  interface->mappingQuasi2D3D      = NULL;
  interface->preciceNodeIDs        = NULL;
  interface->nodeScalarData        = NULL;
  interface->node2DScalarData      = NULL;
  interface->nodeVectorData        = NULL;
  interface->node2DVectorData      = NULL;
  interface->faceCenterData        = NULL;
  interface->xbounIndices          = NULL;
  interface->xloadIndices          = NULL;
  interface->xforcIndices          = NULL;
  interface->writeData             = NULL;
  interface->readData              = NULL;

  // Initialize preCICE mesh name as NULL
  interface->couplingMeshName = NULL;

  // Initialize data names as NULL
  interface->temperature            = NULL;
  interface->flux                   = NULL;
  interface->kDeltaWrite            = NULL;
  interface->kDeltaTemperatureWrite = NULL;
  interface->kDeltaRead             = NULL;
  interface->kDeltaTemperatureRead  = NULL;
  interface->displacements          = NULL;
  interface->displacementDeltas     = NULL;
  interface->positions              = NULL;
  interface->velocities             = NULL;
  interface->forces                 = NULL;
  interface->pressure               = NULL;

  // Check if quasi 2D-3D coupling needs to be implemented
  if (interface->dim == 2) {
    interface->quasi2D3D = 1;
    printf("Using quasi 2D-3D coupling\n");
    interface->dimCCX = 3; // CalculiX always solves a 3D problem
  } else {
    interface->quasi2D3D = 0;
    interface->dimCCX    = interface->dim;
  }

  // Mapping Type
  //  The patch identifies the set used as interface in Calculix
  interface->name = strdup(config->patchName);
  // Calculix needs to know if nearest-projection mapping is implemented. config->map = 1 is for nearest-projection, config->map = 0 is for everything else
  interface->mapNPType = config->map;

  // Nodes mesh
  interface->nodesMeshName = NULL;
  if (config->nodesMeshName) {
    interface->nodesMeshName = strdup(config->nodesMeshName);
    PreciceInterface_ConfigureNodesMesh(interface, sim);

    interface->couplingMeshName = interface->nodesMeshName;
  }

  // Face centers mesh
  interface->faceCentersMeshName = NULL;
  if (config->facesMeshName) {
    interface->faceCentersMeshName = strdup(config->facesMeshName);
    // Only configure a face center mesh if necesary; i.e. do not configure it for FSI simulations, also do not configure tetra faces if no face center mesh is used (as in FSI simulations)
    PreciceInterface_ConfigureFaceCentersMesh(interface, sim);
    // Triangles of the nodes mesh (needs to be called after the face centers mesh is configured!)
    PreciceInterface_ConfigureTetraFaces(interface, sim);

    interface->couplingMeshName = interface->faceCentersMeshName;
  }

  PreciceInterface_ConfigureCouplingData(interface, sim, config);
}

static enum ElemType findSimulationMeshType(SimulationData *sim)
{
  // Assuming only tetrahedra are used, or only hexaedral, for faces meshes.
  // Return first non-zero mesh type.
  // lakon tab takes 8 chars per element

  const char *lakon_ptr = sim->lakon;
  for (int i = 0; i < sim->ne; ++i) {
    if (startsWith(lakon_ptr, "C3D4") || startsWith(lakon_ptr, "C3D10")) {
      return TETRAHEDRA;
    } else if (startsWith(lakon_ptr, "C3D8") || startsWith(lakon_ptr, "C3D20")) {
      return HEXAHEDRA;
    }
    lakon_ptr += 8;
  }

  return INVALID_ELEMENT;
}

void PreciceInterface_ConfigureFaceCentersMesh(PreciceInterface *interface, SimulationData *sim)
{
  // printf("Entering ConfigureFaceCentersMesh \n");
  char *faceSetName      = toFaceSetName(interface->name);
  interface->faceSetID   = getSetID(faceSetName, sim->set, sim->nset);
  interface->numElements = getNumSetElements(interface->faceSetID, sim->istartset, sim->iendset);

  interface->elementIDs = malloc(interface->numElements * sizeof(ITG));
  interface->faceIDs    = malloc(interface->numElements * sizeof(ITG));
  getSurfaceElementsAndFaces(interface->faceSetID, sim->ialset, sim->istartset, sim->iendset, interface->elementIDs, interface->faceIDs);

  interface->faceCenterCoordinates = malloc(interface->numElements * 3 * sizeof(double));
  interface->preciceFaceCenterIDs  = malloc(interface->numElements * 3 * sizeof(int));

  enum ElemType elemType = findSimulationMeshType(sim);

  if (elemType == TETRAHEDRA) {
    printf("Configuring faces mesh with tetrahedra.\n");
    getTetraFaceCenters(interface->elementIDs, interface->faceIDs, interface->numElements, sim->kon, sim->ipkon, sim->co, interface->faceCenterCoordinates);
  } else if (elemType == HEXAHEDRA) {
    printf("Configuring faces mesh with hexahedra.\n");
    getHexaFaceCenters(interface->elementIDs, interface->faceIDs, interface->numElements, sim->kon, sim->ipkon, sim->co, interface->faceCenterCoordinates);
  } else {
    supportedElementError();
  }

  interface->preciceFaceCenterIDs = malloc(interface->numElements * sizeof(int));

  sendFaceCentersVertices(interface);
}

void sendFaceCentersVertices(PreciceInterface *interface)
{
  // Send the data of face centers to preCICE. If 2D coupling is used, skip the z component!

  if (!isQuasi2D3D(interface->quasi2D3D)) {
    precicec_setMeshVertices(interface->faceCentersMeshName, interface->numElements, interface->faceCenterCoordinates, interface->preciceFaceCenterIDs);
  } else {
    double *coordinates2D = malloc(interface->numElements * 2 * sizeof(double));
    for (int i = 0; i < interface->numElements; ++i) {
      coordinates2D[2 * i]     = interface->faceCenterCoordinates[3 * i];
      coordinates2D[2 * i + 1] = interface->faceCenterCoordinates[3 * i + 1];
    }
    precicec_setMeshVertices(interface->faceCentersMeshName, interface->numElements, coordinates2D, interface->preciceFaceCenterIDs);
    free(coordinates2D);
  }
}

void PreciceInterface_ConfigureNodesMesh(PreciceInterface *interface, SimulationData *sim)
{
  // printf("Entering configureNodesMesh \n");
  char *nodeSetName    = toNodeSetName(interface->name);
  interface->nodeSetID = getSetID(nodeSetName, sim->set, sim->nset);
  interface->numNodes  = getNumSetElements(interface->nodeSetID, sim->istartset, sim->iendset);
  interface->nodeIDs   = &sim->ialset[sim->istartset[interface->nodeSetID] - 1]; // Lucia: make a copy

  interface->nodeCoordinates = malloc(interface->numNodes * interface->dimCCX * sizeof(double));
  getNodeCoordinates(interface->nodeIDs, interface->numNodes, interface->dimCCX, sim->co, sim->vold, sim->mt, interface->nodeCoordinates);

  // If 2D-3Q coupling is used (for a node mesh) delegate this to the specialized data structure.
  if (interface->nodesMeshName != NULL) {

    if (isQuasi2D3D(interface->quasi2D3D)) {
      interface->mappingQuasi2D3D = createMapping(interface->nodeCoordinates, interface->numNodes, interface->nodesMeshName);
    } else {
      interface->preciceNodeIDs = malloc(interface->numNodes * sizeof(int));
      precicec_setMeshVertices(interface->nodesMeshName, interface->numNodes, interface->nodeCoordinates, interface->preciceNodeIDs);
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

void PreciceInterface_EnsureValidRead(SimulationData *sim, enum CouplingDataType type)
{
  // Forbidden read data in modal dynamic simulations
  if (sim->isModalDynamic) {
    if (type == TEMPERATURE || type == DISPLACEMENTS || type == DISPLACEMENTDELTAS || type == VELOCITIES || type == POSITIONS) {
      printf("Error: in modal dynamic simulations, only loads (forces, pressures, heat fluxes) can be read.\n"
             "Degrees of freedom (positions, velocities, temperatures) cannot be read from preCICE.");
      fflush(stdout);
      exit(EXIT_FAILURE);
    }
  }
}

void PreciceInterface_ConfigureTetraFaces(PreciceInterface *interface, SimulationData *sim)
{
  int i;
  printf("Setting node connectivity for nearest projection mapping: \n");
  if (interface->nodesMeshName != NULL) {
    int *triangles = malloc(interface->numElements * 3 * sizeof(ITG));
    getTetraFaceNodes(interface->elementIDs, interface->faceIDs, interface->nodeIDs, interface->numElements, interface->numNodes, sim->kon, sim->ipkon, triangles);

    precicec_setMeshTriangles(interface->nodesMeshName, interface->numElements, triangles);
    free(triangles);
  }
}

void PreciceInterface_ConfigureCouplingData(PreciceInterface *interface, SimulationData *sim, InterfaceConfig const *config)
{
  interface->nodeScalarData = malloc(interface->numNodes * sizeof(double));
  interface->nodeVectorData = malloc(interface->numNodes * 3 * sizeof(double));
  interface->faceCenterData = malloc(interface->numElements * sizeof(double));

  // Configure all the read data, then the write data
  int i;
  interface->numReadData = config->numReadData;
  if (config->numReadData > 0)
    interface->readData = malloc(config->numReadData * sizeof(int));
  for (i = 0; i < config->numReadData; i++) {
    if (startsWith(config->readDataNames[i], "Temperature")) {
      PreciceInterface_EnsureValidRead(sim, TEMPERATURE);
      interface->readData[i]  = TEMPERATURE;
      interface->xbounIndices = malloc(interface->numNodes * sizeof(int));
      interface->temperature  = strdup(config->readDataNames[i]);
      getXbounIndices(interface->nodeIDs, interface->numNodes, sim->nboun, sim->ikboun, sim->ilboun, interface->xbounIndices, TEMPERATURE);
      printf("Read data '%s' found.\n", interface->temperature);
    } else if (startsWith(config->readDataNames[i], "Heat-Flux")) {
      interface->readData[i]  = HEAT_FLUX;
      interface->xloadIndices = malloc(interface->numElements * sizeof(int));
      getXloadIndices("DFLUX", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices);
      PreciceInterface_EnsureValidRead(sim, HEAT_FLUX);
      interface->flux = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", interface->flux);
    } else if (startsWith(config->readDataNames[i], "Sink-Temperature")) {
      interface->readData[i]  = SINK_TEMPERATURE;
      interface->xloadIndices = malloc(interface->numElements * sizeof(int));
      getXloadIndices("FILM", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices);
      PreciceInterface_EnsureValidRead(sim, SINK_TEMPERATURE);
      interface->kDeltaTemperatureRead = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", interface->kDeltaTemperatureRead);
    } else if (startsWith(config->readDataNames[i], "Heat-Transfer-Coefficient")) {
      interface->readData[i] = HEAT_TRANSFER_COEFF;
      PreciceInterface_EnsureValidRead(sim, HEAT_TRANSFER_COEFF);
      interface->kDeltaRead = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", interface->kDeltaRead);
    } else if (startsWith(config->readDataNames[i], "Pressure")) {
      interface->readData[i]  = PRESSURE;
      interface->xloadIndices = malloc(interface->numElements * sizeof(int));
      PreciceInterface_EnsureValidRead(sim, PRESSURE);
      getXloadIndices("PRESSUREDLOAD", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices);
      interface->pressure = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", interface->pressure);
    } else if (startsWith(config->readDataNames[i], "Force")) {
      PreciceInterface_EnsureValidRead(sim, FORCES);
      interface->readData[i]  = FORCES;
      interface->xforcIndices = malloc(interface->numNodes * 3 * sizeof(int));
      interface->forces       = strdup(config->readDataNames[i]);
      getXforcIndices(interface->nodeIDs, interface->numNodes, sim->nforc, sim->ikforc, sim->ilforc, interface->xforcIndices);
      printf("Read data '%s' found.\n", interface->forces);
    } else if (startsWith(config->readDataNames[i], "Displacement")) {
      PreciceInterface_EnsureValidRead(sim, DISPLACEMENTS);
      interface->readData[i]   = DISPLACEMENTS;
      interface->xbounIndices  = malloc(interface->numNodes * 3 * sizeof(int));
      interface->displacements = strdup(config->readDataNames[i]);
      getXbounIndices(interface->nodeIDs, interface->numNodes, sim->nboun, sim->ikboun, sim->ilboun, interface->xbounIndices, DISPLACEMENTS);
      printf("Read data '%s' found.\n", interface->displacements);
    } else {
      printf("ERROR: Read data '%s' is not of a known type for the CalculiX-preCICE adapter. Check the adapter configuration file.\n", config->readDataNames[i]);
      exit(EXIT_FAILURE);
    }
  }

  interface->numWriteData = config->numWriteData;
  if (config->numWriteData > 0)
    interface->writeData = malloc(config->numWriteData * sizeof(int));

  for (i = 0; i < config->numWriteData; i++) {
    if (startsWith(config->writeDataNames[i], "Temperature")) {
      interface->writeData[i] = TEMPERATURE;
      interface->temperature  = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found.\n", interface->temperature);
    } else if (startsWith(config->writeDataNames[i], "Heat-Flux")) {
      interface->writeData[i] = HEAT_FLUX;
      interface->flux         = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found'.\n", interface->flux);
    } else if (startsWith(config->writeDataNames[i], "Sink-Temperature")) {
      interface->writeData[i]           = SINK_TEMPERATURE;
      interface->kDeltaTemperatureWrite = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found.\n", interface->kDeltaTemperatureWrite);
    } else if (startsWith(config->writeDataNames[i], "Heat-Transfer-Coefficient")) {
      interface->writeData[i] = HEAT_TRANSFER_COEFF;
      interface->kDeltaWrite  = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found.\n", interface->kDeltaWrite);
    } else if (startsWith(config->writeDataNames[i], "DisplacementDelta")) {
      interface->writeData[i]       = DISPLACEMENTDELTAS;
      interface->displacementDeltas = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found.\n", interface->displacementDeltas);
    } else if (startsWith(config->writeDataNames[i], "Displacement")) {
      interface->writeData[i]  = DISPLACEMENTS;
      interface->displacements = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found.\n", interface->displacements);
    } else if (startsWith(config->writeDataNames[i], "Position")) {
      interface->writeData[i] = POSITIONS;
      interface->positions    = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found.\n", interface->positions);
    }
    /* Both "Velocities" and "Velocity" are valid, so we accept Velocit[...]*/
    else if (startsWith(config->writeDataNames[i], "Velocit")) {
      interface->writeData[i] = VELOCITIES;
      interface->velocities   = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found.\n", interface->velocities);
    } else if (startsWith(config->writeDataNames[i], "Force")) {
      interface->writeData[i] = FORCES;
      interface->forces       = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found.\n", interface->forces);
    } else {
      printf("ERROR: Write data '%s' is not of a known type for the CalculiX-preCICE adapter. Check the adapter configuration file.\n", config->writeDataNames[i]);
      exit(EXIT_FAILURE);
    }
  }
}

void PreciceInterface_FreeData(PreciceInterface *preciceInterface)
{
  free(preciceInterface->readData);
  free(preciceInterface->writeData);
  free(preciceInterface->elementIDs);
  free(preciceInterface->faceIDs);
  free(preciceInterface->preciceFaceCenterIDs);
  free(preciceInterface->faceCenterCoordinates);
  free(preciceInterface->nodeCoordinates);
  free(preciceInterface->preciceNodeIDs);
  free(preciceInterface->nodeScalarData);
  free(preciceInterface->node2DScalarData);
  free(preciceInterface->nodeVectorData);
  free(preciceInterface->node2DVectorData);
  free(preciceInterface->faceCenterData);
  free(preciceInterface->xbounIndices);
  free(preciceInterface->xloadIndices);
  free(preciceInterface->xforcIndices);

  freeMapping(preciceInterface->mappingQuasi2D3D);

  // Mesh names
  free(preciceInterface->faceCentersMeshName);
  free(preciceInterface->nodesMeshName);

  // Data names
  free(preciceInterface->displacementDeltas);
  free(preciceInterface->displacements);
  free(preciceInterface->flux);
  free(preciceInterface->forces);
  free(preciceInterface->kDeltaRead);
  free(preciceInterface->kDeltaTemperatureRead);
  free(preciceInterface->kDeltaTemperatureWrite);
  free(preciceInterface->kDeltaWrite);
  free(preciceInterface->positions);
  free(preciceInterface->pressure);
  free(preciceInterface->temperature);
  free(preciceInterface->velocities);
}
