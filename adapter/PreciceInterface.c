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
  // Precice_ReadCouplingData(sim);
}

void Precice_AdjustSolverTimestep(SimulationData *sim)
{
  double precice_dt = precicec_getMaxTimeStepSize();

  if (isSteadyStateSimulation(sim->nmethod)) {
    printf("Adjusting time step for steady-state step\n");
    fflush(stdout);

    // For steady-state simulations, we will always compute the converged steady-state solution in one coupling step
    sim->theta  = 0;
    sim->tper   = 1;
    sim->dtheta = 1;

    // Set the solver time step to be the same as the coupling time step
    sim->solver_dt = precice_dt;
  } else {
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

  printf("*(sim->theta): %f\n", *(sim->theta));
  printf("*(sim->dtheta): %f\n", *(sim->dtheta));
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

  printf("Adapter reading coupling data....\n");
  // printf("precicec_isReadDataAvailable()  %d \n,",precicec_isReadDataAvailable());
  fflush(stdout);

  PreciceInterface **interfaces    = sim->preciceInterfaces;
  int                numInterfaces = sim->numPreciceInterfaces;
  int                i, j, idx;

  for (i = 0; i < numInterfaces; i++) {

    for (j = 0; j < interfaces[i]->numReadData; j++) {

      printf("Read data: %d\n", interfaces[i]->readData[j]);

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
        // precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->pressure, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, sim->solver_dt, interfaces[i]->faceCenterData);
        // setFacePressure(interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload);
        // printf("Reading PRESSURE coupling data.\n");
        // break;
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

      // VOLUMETRIC COUPLING - MULTISCALE
      case CONV_FLAG:
      // READ CONVERGENCE FLAG
        // precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->macroInputData, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, sim->solver_dt, interfaces[i]->elementIPScalarData);
        printf("Reading CONVERGENCE FLAG coupling data.\n");
        break;

      case CMAT1:
        // READ MATERIAL MATRIX COMPONENTS - C11, C12, C13        
        idx=1;
        precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->materialTangent1Data, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, sim->solver_dt, interfaces[i]->elementIPVectorData);
        // for (int k = 0; k < interfaces[i]->numIPTotal; k++) {
        //     printf( " %i, %e, %e, %e \n", k, interfaces[i]->elementIPVectorData[k*3], interfaces[i]->elementIPVectorData[k*3+1], interfaces[i]->elementIPVectorData[k*3+2]);
        // }
        FORTRAN(precice_multiscale_set_xstiff, (sim->mi,
                                  idx, 
                                  interfaces[i]->numElements,
                                  interfaces[i]->elementIPVectorData,
                                  sim->xstiff));
        printf("Reading MATERIAL TANGENT 1 coupling data.\n");
        break;

      case CMAT2:
        // READ MATERIAL MATRIX COMPONENTS -  C14, C15, C16
        // idx=4;
        // precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->materialTangent2Data, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, sim->solver_dt, interfaces[i]->elementIPVectorData);
        // FORTRAN(precice_multiscale_set_xstiff, (sim->mi,
        //                           &idx, 
        //                           &interfaces[i]->numElements,
        //                           interfaces[i]->elementIPVectorData,
        //                           sim->xstiff));
        printf("Reading MATERIAL TANGENT 2 coupling data.\n");
        break;

      case CMAT3:
        // READ MATERIAL MATRIX COMPONENTS -  C22, C23, C24
        // idx=7;
        // precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->materialTangent3Data, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, sim->solver_dt, interfaces[i]->elementIPVectorData);
        // FORTRAN(precice_multiscale_set_xstiff, (sim->mi,
        //                           &idx, 
        //                           &interfaces[i]->numElements,
        //                           interfaces[i]->elementIPVectorData,
        //                           sim->xstiff));
        printf("Reading MATERIAL TANGENT 3 coupling data.\n");
        break;

      case CMAT4:
        // READ MATERIAL MATRIX COMPONENTS -  C25, C26, C33
        // idx=10;
        // precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->materialTangent4Data, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, sim->solver_dt, interfaces[i]->elementIPVectorData);
        // FORTRAN(precice_multiscale_set_xstiff, (sim->mi,
        //                           &idx, 
        //                           &interfaces[i]->numElements,
        //                           interfaces[i]->elementIPVectorData,
        //                           sim->xstiff));
        printf("Reading MATERIAL TANGENT 4 coupling data.\n");
        break;

      case CMAT5:
        // READ MATERIAL MATRIX COMPONENTS -  C34, C35, C36
        // idx=13;
        // precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->materialTangent5Data, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, sim->solver_dt, interfaces[i]->elementIPVectorData);
        // FORTRAN(precice_multiscale_set_xstiff, (sim->mi,
        //                           &idx, 
        //                           &interfaces[i]->numElements,
        //                           interfaces[i]->elementIPVectorData,
        //                           sim->xstiff));
        printf("Reading MATERIAL TANGENT 5 coupling data.\n");
        break;

      case CMAT6:
        // READ MATERIAL MATRIX COMPONENTS -  C44, C45, C46
        // idx=16;
        // precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->materialTangent6Data, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, sim->solver_dt, interfaces[i]->elementIPVectorData);
        // FORTRAN(precice_multiscale_set_xstiff, (sim->mi,
        //                           &idx, 
        //                           &interfaces[i]->numElements,
        //                           interfaces[i]->elementIPVectorData,
        //                           sim->xstiff));
        printf("Reading MATERIAL TANGENT 6 coupling data.\n");
        break;

      case CMAT7:
        // READ MATERIAL MATRIX COMPONENTS -  C55, C56, C66
        // idx=19;
        // precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->materialTangent7Data, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, sim->solver_dt, interfaces[i]->elementIPVectorData);
        // FORTRAN(precice_multiscale_set_xstiff, (sim->mi,
        //                           &idx, 
        //                           &interfaces[i]->numElements,
        //                           interfaces[i]->elementIPVectorData,
        //                           sim->xstiff));
        printf("Reading MATERIAL TANGENT 7 coupling data.\n");
        break;

      case STRESS1TO3:
      // READ STRESS COMPONENTS - S11, S22, S33
        idx=1;
        // precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->stress1to3Data, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, sim->solver_dt, interfaces[i]->elementIPVectorData);
        // FORTRAN(precice_multiscale_set_stx, (sim->mi,
        //                           &idx,
        //                           &interfaces[i]->numElements,
        //                           interfaces[i]->elementIPVectorData,
        //                           sim->stx));
        printf("Reading STRESS1TO3 coupling data.\n");
        break;

      case STRESS4TO6:
        // READ STRESS COMPONENTS - S23, S13, S12
        // idx=4;
        // precicec_readData(interfaces[i]->couplingMeshName, interfaces[i]->stress4to6Data, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, sim->solver_dt, interfaces[i]->elementIPVectorData);
        // FORTRAN(precice_multiscale_set_stx, (sim->mi,
        //                           &idx,
        //                           &interfaces[i]->numElements,
        //                           interfaces[i]->elementIPVectorData,
        //                           sim->stx));
        printf("Reading STRESS4TO6 coupling data.\n");
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
      case MACRO_IP_ID:
        printf("MACRO IP ID  cannot be used as read data.\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
        break;
      case INPUT_ID:
        printf("Input ID cannot be used as read data.\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
        break;
      case STRAIN1TO3:
        printf("Strain 1to3 cannot be used as read data.\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
        break;
      case STRAIN4TO6:
        printf("Strain 4to6 cannot be used as read data.\n");
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
  int                i, j, idx;
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
      printf("Write data: %d\n", interfaces[i]->writeData[j]);
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
      /* VOLUMETRIC COUPLING - MULTISCALE */
      case INPUT_ID:
        // WRITE INPUT ID
        for (int k = 0; k < interfaces[i]->numIPTotal; k++) {
          interfaces[i]->elementIPScalarData[k] = 1;
        };
        precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->macroInputData, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, interfaces[i]->elementIPScalarData);
        printf("Writing INPUT ID coupling data.\n");
        break;
      case STRAIN1TO3:
        idx=0;
        getElementStrain(idx, sim->mi, interfaces[i]->numElements, sim->eei, interfaces[i]->elementIPVectorData);
        precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->strain1to3Data, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, interfaces[i]->elementIPVectorData);
        printf("Writing STRAIN1TO3 coupling data.\n");
        break;
      case STRAIN4TO6:
        idx=3;
        getElementStrain(idx, sim->mi, interfaces[i]->numElements, sim->eei, interfaces[i]->elementIPVectorData);
        precicec_writeData(interfaces[i]->couplingMeshName, interfaces[i]->strain4to6Data, interfaces[i]->numIPTotal, interfaces[i]->elemIPID, interfaces[i]->elementIPVectorData);
        printf("Writing STRAIN4TO6 coupling data.\n");
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
  precicec_finalize();
}

void PreciceInterface_Create(PreciceInterface *interface, SimulationData *sim, InterfaceConfig const *config)
{
  // Deduce configured dimensions
  if (config->nodesMeshName == NULL && config->facesMeshName == NULL && config->elementsMeshName == NULL) {
    printf("ERROR: You need to define a face mesh, nodes mesh or element mesh. Check the adapter configuration file.\n");
    exit(EXIT_FAILURE);
  }
  if (config->nodesMeshName && config->facesMeshName) {
    int facesDims = precicec_getMeshDimensions(config->facesMeshName);
    int nodesDims = precicec_getMeshDimensions(config->nodesMeshName);
    if (facesDims != nodesDims) {
      printf("ERROR: You configured a %dD face mesh '%s' and a %dD nodes mesh '%s'. These meshes needs to configured with the same dimensionality .\n", facesDims, config->facesMeshName, nodesDims, config->nodesMeshName);
      exit(EXIT_FAILURE);
    }
    interface->dim = nodesDims;
  } else {
    if (config->nodesMeshName) {
      interface->dim = precicec_getMeshDimensions(config->nodesMeshName);
    }
    if (config->facesMeshName) {
      interface->dim = precicec_getMeshDimensions(config->facesMeshName);
    }
  }

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
  interface->macroInputData         = NULL;
  interface->strain1to3Data         = NULL;
  interface->strain4to6Data         = NULL;
  interface->stress1to3Data         = NULL;
  interface->stress4to6Data         = NULL;
  interface->materialTangent1Data   = NULL;
  interface->materialTangent2Data   = NULL;
  interface->materialTangent3Data   = NULL;
  interface->materialTangent4Data   = NULL;
  interface->materialTangent5Data   = NULL;
  interface->materialTangent6Data   = NULL;
  interface->materialTangent7Data   = NULL;

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

  // Element mesh
  interface->elementMeshName = NULL;
  if (config->elementsMeshName) {
    interface->elementMeshName = strdup(config->elementsMeshName);
    // Configuring element mesh
    PreciceInterface_ConfigureElementsMesh(interface, sim);
    interface->couplingMeshName = interface->elementMeshName;
  }

  PreciceInterface_ConfigureCouplingData(interface, sim, config);
}

void PreciceInterface_ConfigureElementsMesh(PreciceInterface *interface, SimulationData *sim)
{
  printf("Entering ConfigureElementsMesh \n");
  char *elementSetName    = interface->name; //toFaceSetName(interface->name);
  interface->elementSetID = getSetID(elementSetName, sim->set, sim->nset);
  interface->numElements  = getNumSetElements(interface->elementSetID, sim->istartset, sim->iendset);

  printf("element set id : %d\n", interface->elementSetID);
  printf("num elements : %d\n", interface->numElements);

  interface->elementIDs = malloc(interface->numElements * sizeof(ITG));
  getElementsIDs(interface->elementSetID, sim->ialset, sim->istartset, sim->iendset, interface->elementIDs);

  for (int j = 0; j < interface->numElements; j++) {
    printf(" %d, element id: %d \n", j, interface->elementIDs[j]);
  }

  // Find guass point coordinates of the element -> Serves as mesh for data transfer
  int numElt                   = interface->numElements;
  interface->numIPTotal        = 8 * interface->numElements; // Gauss point mesh coordinate -Each element 8 gauss points
  interface->elemIPCoordinates = malloc(interface->numIPTotal * 3 * sizeof(double));
  interface->elemIPID          = malloc(interface->numIPTotal * sizeof(int));

  for (int j = 0; j < interface->numIPTotal; j++) {
    interface->elemIPID[j] = j;
    interface->elemIPCoordinates[j * 3]     = j;
    interface->elemIPCoordinates[j * 3 + 1] = 0.0;
    interface->elemIPCoordinates[j * 3 + 2] = 0.0;
  }

  // getElementGaussPointCoordinates(interface->numElements, interface->numIPTotal, interface->elementIDs, sim->co,
  //                                 sim->kon, sim->lakon, sim->ipkon, interface->elemIPID, interface->elemIPCoordinates);

  precicec_setMeshVertices(interface->elementMeshName, interface->numIPTotal, interface->elemIPCoordinates, interface->elemIPID);
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
  //printf("Entering configureNodesMesh \n");
  char *nodeSetName    = toNodeSetName(interface->name);
  interface->nodeSetID = getSetID(nodeSetName, sim->set, sim->nset);
  interface->numNodes  = getNumSetElements(interface->nodeSetID, sim->istartset, sim->iendset);
  interface->nodeIDs   = &sim->ialset[sim->istartset[interface->nodeSetID] - 1]; //Lucia: make a copy

  interface->nodeCoordinates = malloc(interface->numNodes * interface->dimCCX * sizeof(double));
  getNodeCoordinates(interface->nodeIDs, interface->numNodes, interface->dimCCX, sim->co, sim->vold, sim->mt, interface->nodeCoordinates);

  // If 2D-3Q coupling is used (for a node mesh) delegate this to the specialized data structure.
  if (interface->nodesMeshName != NULL) {

    int count = 0;
    int dimCCX = interface->dimCCX;
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
    // interface->nodesMeshID = precicec_getMeshID(interface->nodesMeshName);
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

void PreciceInterface_EnsureValidRead(PreciceInterface *interface, enum CouplingDataType type)
{

  if (interface->elementMeshName == NULL) {
    printf("Element mesh not provided in YAML config file\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
}

void PreciceInterface_ConfigureTetraFaces(PreciceInterface *interface, SimulationData *sim)
{
  // int i;
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

  /* Allocating and initilizing memory for multiscale coupling */
  interface->elementIPScalarData = malloc(interface->numIPTotal * sizeof(double));
  interface->elementIPVectorData = malloc(interface->numIPTotal * 3 * sizeof(double));
  for (int i = 0; i < interface->numIPTotal; i++) {
    interface->elementIPScalarData[i] = 0.0;
    interface->elementIPVectorData[i * 3]     = 0.0;
    interface->elementIPVectorData[i * 3 + 1] = 0.0;
    interface->elementIPVectorData[i * 3 + 2] = 0.0;
  }

  int i;
  interface->numReadData = config->numReadData;
  if (config->numReadData > 0)
    interface->readData = malloc(config->numReadData * sizeof(int));
  for (i = 0; i < config->numReadData; i++) {
    if (startsWith(config->readDataNames[i], "Temperature")) {
      PreciceInterface_EnsureValidRead(interface, TEMPERATURE);
      interface->readData[i]  = TEMPERATURE;
      interface->xbounIndices = malloc(interface->numNodes * sizeof(int));
      interface->temperature  = strdup(config->readDataNames[i]);
      getXbounIndices(interface->nodeIDs, interface->numNodes, sim->nboun, sim->ikboun, sim->ilboun, interface->xbounIndices, TEMPERATURE);
      printf("Read data '%s' found.\n", interface->temperature);
    } else if (startsWith(config->readDataNames[i], "Heat-Flux")) {
      interface->readData[i]  = HEAT_FLUX;
      interface->xloadIndices = malloc(interface->numElements * sizeof(int));
      getXloadIndices("DFLUX", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices);
      PreciceInterface_EnsureValidRead(interface, HEAT_FLUX);
      interface->flux = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", interface->flux);
    } else if (startsWith(config->readDataNames[i], "Sink-Temperature")) {
      interface->readData[i]  = SINK_TEMPERATURE;
      interface->xloadIndices = malloc(interface->numElements * sizeof(int));
      getXloadIndices("FILM", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices);
      PreciceInterface_EnsureValidRead(interface, SINK_TEMPERATURE);
      interface->kDeltaTemperatureRead = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", interface->kDeltaTemperatureRead);
    } else if (startsWith(config->readDataNames[i], "Heat-Transfer-Coefficient")) {
      interface->readData[i] = HEAT_TRANSFER_COEFF;
      PreciceInterface_EnsureValidRead(interface, HEAT_TRANSFER_COEFF);
      interface->kDeltaRead = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", interface->kDeltaRead);
    } else if (startsWith(config->readDataNames[i], "Pressure")) {
      interface->readData[i]  = PRESSURE;
      interface->xloadIndices = malloc(interface->numElements * sizeof(int));
      PreciceInterface_EnsureValidRead(interface, PRESSURE);
      getXloadIndices("PRESSUREDLOAD", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices);
      interface->pressure = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", interface->pressure);
    } else if (startsWith(config->readDataNames[i], "Force")) {
      PreciceInterface_EnsureValidRead(interface, FORCES);
      interface->readData[i]  = FORCES;
      interface->xforcIndices = malloc(interface->numNodes * 3 * sizeof(int));
      interface->forces       = strdup(config->readDataNames[i]);
      getXforcIndices(interface->nodeIDs, interface->numNodes, sim->nforc, sim->ikforc, sim->ilforc, interface->xforcIndices);
      printf("Read data '%s' found.\n", interface->forces);
    } else if (startsWith(config->readDataNames[i], "Displacement")) {
      PreciceInterface_EnsureValidRead(interface, DISPLACEMENTS);
      interface->readData[i]       = DISPLACEMENTS;
      interface->xbounIndices      = malloc(interface->numNodes * 3 * sizeof(int));
      interface->displacementsData = strdup(config->readDataNames[i]);
      getXbounIndices(interface->nodeIDs, interface->numNodes, sim->nboun, sim->ikboun, sim->ilboun, interface->xbounIndices, DISPLACEMENTS);
      printf("Read data '%s' found.\n", config->readDataNames[i]);
      /* MICROMANAGER COUPLING */
    } else if (startsWith(config->readDataNames[i], "cmat1")) {
      PreciceInterface_EnsureValidRead(interface, CMAT1);
      interface->readData[i]          = CMAT1;
      interface->materialTangent1Data = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", config->readDataNames[i]);
    } else if (startsWith(config->readDataNames[i], "cmat2")) {
      PreciceInterface_EnsureValidRead(interface, CMAT2);
      interface->readData[i] = CMAT2;
      interface->materialTangent2Data = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", config->readDataNames[i]);
    } else if (startsWith(config->readDataNames[i], "cmat3")) {
      PreciceInterface_EnsureValidRead(interface, CMAT3);
      interface->readData[i] = CMAT3;
      interface->materialTangent3Data = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", config->readDataNames[i]);
    } else if (startsWith(config->readDataNames[i], "cmat4")) {
      PreciceInterface_EnsureValidRead(interface, CMAT4);
      interface->readData[i] = CMAT4;
      interface->materialTangent4Data = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", config->readDataNames[i]);
    } else if (startsWith(config->readDataNames[i], "cmat5")) {
      PreciceInterface_EnsureValidRead(interface, CMAT5);
      interface->readData[i] = CMAT5;
      interface->materialTangent5Data = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", config->readDataNames[i]);
    } else if (startsWith(config->readDataNames[i], "cmat6")) {
      PreciceInterface_EnsureValidRead(interface, CMAT6);
      interface->readData[i] = CMAT6;
      interface->materialTangent6Data = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", config->readDataNames[i]);
    } else if (startsWith(config->readDataNames[i], "cmat7")) {
      PreciceInterface_EnsureValidRead(interface, CMAT7);
      interface->readData[i] = CMAT7;
      interface->materialTangent7Data = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", config->readDataNames[i]);
    } else if (startsWith(config->readDataNames[i], "conv_flag")) {
      PreciceInterface_EnsureValidRead(interface, CONV_FLAG);
      interface->readData[i] = CONV_FLAG;
      printf("Read data '%s' found.\n", config->readDataNames[i]);
    } else if (startsWith(config->readDataNames[i], "stress1to3")) {
      PreciceInterface_EnsureValidRead(interface, STRESS1TO3);
      interface->readData[i] = STRESS1TO3;
      interface->stress1to3Data = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", config->readDataNames[i]);
    } else if (startsWith(config->readDataNames[i], "stress4to6")) {
      PreciceInterface_EnsureValidRead(interface, STRESS4TO6);
      interface->readData[i] = STRESS4TO6;
      interface->stress4to6Data = strdup(config->readDataNames[i]);
      printf("Read data '%s' found.\n", config->readDataNames[i]);
    } else {
      printf("ERROR: Read data '%s' does not exist!\n", config->readDataNames[i]);
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
    } else if (isEqual(config->writeDataNames[i], "input_id")) {
      interface->writeData[i] = INPUT_ID;
      interface->macroInputData    = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found.\n", config->writeDataNames[i]);
    } else if (isEqual(config->writeDataNames[i], "strain1to3")) {
      interface->writeData[i]   = STRAIN1TO3;
      interface->strain1to3Data = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found.\n", config->writeDataNames[i]);
    } else if (isEqual(config->writeDataNames[i], "strain4to6")) {
      interface->writeData[i]      = STRAIN4TO6;
      interface->strain4to6Data = strdup(config->writeDataNames[i]);
      printf("Write data '%s' found.\n", config->writeDataNames[i]);
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
  free(preciceInterface->strain1to3Data);
  free(preciceInterface->strain4to6Data);
  free(preciceInterface->stress1to3Data);
  free(preciceInterface->stress4to6Data);
  free(preciceInterface->materialTangent1Data);
  free(preciceInterface->materialTangent2Data);
  free(preciceInterface->materialTangent3Data);
  free(preciceInterface->materialTangent4Data);
  free(preciceInterface->materialTangent5Data);
  free(preciceInterface->materialTangent6Data);
  free(preciceInterface->materialTangent7Data);
}

void PreciceInterface_MultiscaleCheckpoint(SimulationData *sim){

  if( Precice_IsCouplingOngoing()){

  // Write strain data
  Precice_WriteCouplingData(sim);

  // Advance time
  Precice_Advance(sim);

  // Read stress, material tangent data
  Precice_ReadCouplingData(sim);

  }
}
