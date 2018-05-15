/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#include <stdlib.h>
#include "PreciceInterface.h"
#include "ConfigReader.h"
#include "precice/adapters/c/SolverInterfaceC.h"

void Precice_Setup( char * configFilename, char * participantName, SimulationData * sim )
{

	printf( "Setting up preCICE participant %s, using config file: %s\n", participantName, configFilename );
	fflush( stdout );

	int i;
	char * preciceConfigFilename;
	InterfaceConfig * interfaces;

	// Read the YAML config file
	ConfigReader_Read( configFilename, participantName, &preciceConfigFilename, &interfaces, &sim->numPreciceInterfaces );

	// Create the solver interface and configure it - Alex: Calculix is always a serial participant (MPI size 1, rank 0)
	precicec_createSolverInterface( participantName, preciceConfigFilename, 0, 1 );

	// Create interfaces as specified in the config file
	sim->preciceInterfaces = (struct PreciceInterface**) malloc( sim->numPreciceInterfaces * sizeof( PreciceInterface* ) );

	for( i = 0 ; i < sim->numPreciceInterfaces ; i++ )
	{
		sim->preciceInterfaces[i] = malloc( sizeof( PreciceInterface ) );
		PreciceInterface_Create( sim->preciceInterfaces[i], sim, &interfaces[i] );
	}
	// Initialize variables needed for the coupling
	NNEW( sim->coupling_init_v, double, sim->mt * sim->nk );

	// Initialize preCICE
	sim->precice_dt = precicec_initialize();

	// Initialize coupling data
	Precice_InitializeData( sim );

}

void Precice_InitializeData( SimulationData * sim )
{
	printf( "Initializing coupling data\n" );
	fflush( stdout );

	Precice_WriteCouplingData( sim );
	precicec_initialize_data();
	Precice_ReadCouplingData( sim );
}

void Precice_AdjustSolverTimestep( SimulationData * sim )
{
	if( isSteadyStateSimulation( sim->nmethod ) )
	{
		printf( "Adjusting time step for steady-state step\n" );
		fflush( stdout );

		// For steady-state simulations, we will always compute the converged steady-state solution in one coupling step
		*sim->theta = 0;
		*sim->tper = 1;
		*sim->dtheta = 1;

		// Set the solver time step to be the same as the coupling time step
		sim->solver_dt = sim->precice_dt;
	}
	else
	{
		printf( "Adjusting time step for transient step\n" );
		printf( "precice_dt dtheta = %f, dtheta = %f, solver_dt = %f\n", sim->precice_dt / *sim->tper, *sim->dtheta, fmin( sim->precice_dt, *sim->dtheta * *sim->tper ) );
		fflush( stdout );

		// Compute the normalized time step used by CalculiX
		*sim->dtheta = fmin( sim->precice_dt / *sim->tper, *sim->dtheta );

		// Compute the non-normalized time step used by preCICE
		sim->solver_dt = ( *sim->dtheta ) * ( *sim->tper );
	}
}

void Precice_Advance( SimulationData * sim )
{
	printf( "Adapter calling advance()...\n" );
	fflush( stdout );

	sim->precice_dt = precicec_advance( sim->solver_dt );
}

bool Precice_IsCouplingOngoing()
{
	return precicec_isCouplingOngoing();
}

bool Precice_IsReadCheckpointRequired()
{
	return precicec_isActionRequired( "read-iteration-checkpoint" );
}

bool Precice_IsWriteCheckpointRequired()
{
	return precicec_isActionRequired( "write-iteration-checkpoint" );
}

void Precice_FulfilledReadCheckpoint()
{
	precicec_fulfilledAction( "read-iteration-checkpoint" );
}

void Precice_FulfilledWriteCheckpoint()
{
	precicec_fulfilledAction( "write-iteration-checkpoint" );
}

void Precice_ReadIterationCheckpoint( SimulationData * sim, double * v )
{

	printf( "Adapter reading checkpoint...\n" );
	fflush( stdout );

	// Reload time
	*( sim->theta ) = sim->coupling_init_theta;

	// Reload step size
	*( sim->dtheta ) = sim->coupling_init_dtheta;

	// Reload solution vector v
	memcpy( v, sim->coupling_init_v, sizeof( double ) * sim->mt * sim->nk );

}

void Precice_WriteIterationCheckpoint( SimulationData * sim, double * v )
{

	printf( "Adapter writing checkpoint...\n" );
	fflush( stdout );

	// Save time
	sim->coupling_init_theta = *( sim->theta );

	// Save step size
	sim->coupling_init_dtheta = *( sim->dtheta );

	// Save solution vector v
	memcpy( sim->coupling_init_v, v, sizeof( double ) * sim->mt * sim->nk );

}

void Precice_ReadCouplingData( SimulationData * sim )
{

	printf( "Adapter reading coupling data...\n" );
	fflush( stdout );

	PreciceInterface ** interfaces = sim->preciceInterfaces;
	int numInterfaces = sim->numPreciceInterfaces;
	int i;
	
	if( precicec_isReadDataAvailable() )
	{
		for( i = 0 ; i < numInterfaces ; i++ )
		{
			switch( interfaces[i]->readData )
			{
			case TEMPERATURE:
				// Read and set temperature BC
				precicec_readBlockScalarData( interfaces[i]->temperatureDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeScalarData );
				setNodeTemperatures( interfaces[i]->nodeScalarData, interfaces[i]->numNodes, interfaces[i]->xbounIndices, sim->xboun );
				break;
			case HEAT_FLUX:
				// Read and set heat flux BC
				precicec_readBlockScalarData( interfaces[i]->fluxDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData );
				setFaceFluxes( interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload );
				break;
			case CONVECTION:
				// Read and set sink temperature in convective film BC
				precicec_readBlockScalarData( interfaces[i]->kDeltaTemperatureReadDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData );
				setFaceSinkTemperatures( interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload );
				// Read and set heat transfer coefficient in convective film BC
				precicec_readBlockScalarData( interfaces[i]->kDeltaReadDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData );
				setFaceHeatTransferCoefficients( interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload );
				break;
            case FORCES:
                //Read and set forces as point loads
                precicec_readBlockVectorData( interfaces[i]->forcesDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
                setNodeForces( interfaces[i]->nodeVectorData, interfaces[i]->numNodes, interfaces[i]->xforcIndices, sim->xforc );
                break;
            case DISPLACEMENTS:
                printf( "Displacements cannot be used as read data\n" );
                fflush( stdout );
                exit( EXIT_FAILURE );
                break;
	    case DISPLACEMENTDELTAS:
                printf( "DisplacementDeltas cannot be used as read data\n" );
                fflush( stdout );
                exit( EXIT_FAILURE );
                break;
			}
		}
	}
}

void Precice_WriteCouplingData( SimulationData * sim )
{

	printf( "Adapter writing coupling data...\n" );
	fflush( stdout );

	PreciceInterface ** interfaces = sim->preciceInterfaces;
	int numInterfaces = sim->numPreciceInterfaces;
	int i;
	int iset;

	if( precicec_isWriteDataRequired( sim->solver_dt ) || precicec_isActionRequired( "write-initial-data" ) )
	{
		for( i = 0 ; i < numInterfaces ; i++ )
		{
			switch( interfaces[i]->writeData )
			{
			case TEMPERATURE:
				getNodeTemperatures( interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->vold, sim->mt, interfaces[i]->nodeScalarData );
				precicec_writeBlockScalarData( interfaces[i]->temperatureDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeScalarData );
				break;
			case HEAT_FLUX:
				iset = interfaces[i]->faceSetID + 1; // Adjust index before calling Fortran function
				FORTRAN( getflux, ( sim->co,
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
									interfaces[i]->faceCenterData
									)
						 );
				precicec_writeBlockScalarData( interfaces[i]->fluxDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData );
				break;
			case CONVECTION:
				iset = interfaces[i]->faceSetID + 1; // Adjust index before calling Fortran function
				double * myKDelta = malloc( interfaces[i]->numElements * sizeof( double ) );
				double * T = malloc( interfaces[i]->numElements * sizeof( double ) );
				FORTRAN( getkdeltatemp, ( sim->co,
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
										  myKDelta,
										  T
										  )
						 );
				precicec_writeBlockScalarData( interfaces[i]->kDeltaWriteDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, myKDelta );
				precicec_writeBlockScalarData( interfaces[i]->kDeltaTemperatureWriteDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, T );
				free( myKDelta );
				free( T );
				break;
            case DISPLACEMENTS:
                getNodeDisplacements( interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->vold, sim->mt, interfaces[i]->nodeVectorData );
                precicec_writeBlockVectorData( interfaces[i]->displacementsDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
                break;
	    case DISPLACEMENTDELTAS:
                getNodeDisplacementDeltas( interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->vold, sim->coupling_init_v, sim->mt, interfaces[i]->nodeVectorData );
                precicec_writeBlockVectorData( interfaces[i]->displacementDeltasDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
                break;
            case FORCES:
                printf( "Forces cannot be used as write data\n" );
                fflush( stdout );
                exit( EXIT_FAILURE );
                break;
			}
		}

		if( precicec_isActionRequired( "write-initial-data" ) )
		{
			precicec_fulfilledAction( "write-initial-data" );
		}
	}
}

void Precice_FreeData( SimulationData * sim )
{
	int i;

	if( sim->coupling_init_v != NULL ){
		free( sim->coupling_init_v );
	}

	for( i = 0 ; i < sim->numPreciceInterfaces ; i++ )
	{
		PreciceInterface_FreeData( sim->preciceInterfaces[i] );
		if( sim->preciceInterfaces[i] != NULL ){
			free( sim->preciceInterfaces[i] );
		}
	}
	
	precicec_finalize();
}

void PreciceInterface_Create( PreciceInterface * interface, SimulationData * sim, InterfaceConfig * config )
{

	// Initialize pointers as NULL
	interface->elementIDs = NULL;
	interface->faceIDs = NULL;
	interface->faceCenterCoordinates = NULL;
	interface->preciceFaceCenterIDs = NULL;
	interface->nodeCoordinates = NULL;
	interface->preciceNodeIDs = NULL;
	interface->triangles = NULL;
	interface->nodeScalarData = NULL;
        interface->nodeVectorData = NULL;
	interface->faceCenterData = NULL;
	interface->xbounIndices = NULL;
	interface->xloadIndices = NULL;
    	interface->xforcIndices = NULL;

    	//The patch identifies the set used as interface in Calculix
	interface->name = config->patchName;

	// Nodes mesh
	interface->nodesMeshID = -1;
	interface->nodesMeshName = config->nodesMeshName;
	PreciceInterface_ConfigureNodesMesh( interface, sim );

	// Face centers mesh
	interface->faceCentersMeshID = -1;
	interface->faceCentersMeshName = config->facesMeshName;
    	//Only configure a face center mesh if necesary; i.e. do not configure it for FSI simulations, also do not configure tetra faces if no face center mesh is used (as in FSI simulations)
    	if ( interface->faceCentersMeshName != NULL) {
        	PreciceInterface_ConfigureFaceCentersMesh( interface, sim );
        	// Triangles of the nodes mesh (needs to be called after the face centers mesh is configured!)
        	PreciceInterface_ConfigureTetraFaces( interface, sim );
    	}

	PreciceInterface_ConfigureCouplingData( interface, sim, config );

}

void PreciceInterface_ConfigureFaceCentersMesh( PreciceInterface * interface, SimulationData * sim )
{

	char * faceSetName = toFaceSetName( interface->name );
	interface->faceSetID = getSetID( faceSetName, sim->set, sim->nset );
	interface->numElements = getNumSetElements( interface->faceSetID, sim->istartset, sim->iendset );

	interface->elementIDs = malloc( interface->numElements * sizeof( ITG ) );
	interface->faceIDs = malloc( interface->numElements * sizeof( ITG ) );
	getSurfaceElementsAndFaces( interface->faceSetID, sim->ialset, sim->istartset, sim->iendset, interface->elementIDs, interface->faceIDs );

	interface->faceCenterCoordinates = malloc( interface->numElements * 3 * sizeof( double ) );
	getTetraFaceCenters( interface->elementIDs, interface->faceIDs, interface->numElements, sim->kon, sim->ipkon, sim->co, interface->faceCenterCoordinates );

	interface->faceCentersMeshID = precicec_getMeshID( interface->faceCentersMeshName );
	interface->preciceFaceCenterIDs = malloc( interface->numElements * sizeof( int ) );
	precicec_setMeshVertices( interface->faceCentersMeshID, interface->numElements, interface->faceCenterCoordinates, interface->preciceFaceCenterIDs );

}

void PreciceInterface_ConfigureNodesMesh( PreciceInterface * interface, SimulationData * sim )
{

	char * nodeSetName = toNodeSetName( interface->name );
	interface->nodeSetID = getSetID( nodeSetName, sim->set, sim->nset );
	interface->numNodes = getNumSetElements( interface->nodeSetID, sim->istartset, sim->iendset );
	interface->nodeIDs = &sim->ialset[sim->istartset[interface->nodeSetID] - 1]; //Lucia: make a copy

	interface->nodeCoordinates = malloc( interface->numNodes * 3 * sizeof( double ) );
	getNodeCoordinates( interface->nodeIDs, interface->numNodes, sim->co, sim->vold, sim->mt, interface->nodeCoordinates );

	if( interface->nodesMeshName != NULL )
	{
		interface->nodesMeshID = precicec_getMeshID( interface->nodesMeshName );
		interface->preciceNodeIDs = malloc( interface->numNodes * sizeof( int ) );
		precicec_setMeshVertices( interface->nodesMeshID, interface->numNodes, interface->nodeCoordinates, interface->preciceNodeIDs );
	}

}

void PreciceInterface_EnsureValidNodesMeshID( PreciceInterface * interface )
{
	if( interface->nodesMeshID < 0 )
	{
		printf( "Nodes mesh not provided in YAML config file\n" );
		fflush( stdout );
		exit( EXIT_FAILURE );
	}
}

void PreciceInterface_ConfigureTetraFaces( PreciceInterface * interface, SimulationData * sim )
{
	int i;

	if( interface->nodesMeshName != NULL )
	{
		interface->triangles = malloc( interface->numElements * 3 * sizeof( ITG ) );
		getTetraFaceNodes( interface->elementIDs, interface->faceIDs,  interface->nodeIDs, interface->numElements, interface->numNodes, sim->kon, sim->ipkon, interface->triangles );

		for( i = 0 ; i < interface->numElements ; i++ )
		{
			precicec_setMeshTriangleWithEdges( interface->nodesMeshID, interface->triangles[3*i], interface->triangles[3*i+1], interface->triangles[3*i+2] );
		}
	}
}

void PreciceInterface_ConfigureCouplingData( PreciceInterface * interface, SimulationData * sim, InterfaceConfig * config )
{

	interface->nodeScalarData = malloc( interface->numNodes * sizeof( double ) );
    	interface->nodeVectorData = malloc( interface->numNodes * 3 * sizeof( double ) );
	interface->faceCenterData = malloc( interface->numElements * sizeof( double ) );

	int i;

	for( i = 0 ; i < config->numReadData ; i++ )
	{
        
		if( strcmp( config->readDataNames[i], "Temperature" ) == 0 )
		{

			PreciceInterface_EnsureValidNodesMeshID( interface );
			interface->readData = TEMPERATURE;
			interface->xbounIndices = malloc( interface->numNodes * sizeof( int ) );
			interface->temperatureDataID = precicec_getDataID( "Temperature", interface->nodesMeshID );
			getXbounIndices( interface->nodeIDs, interface->numNodes, sim->nboun, sim->ikboun, sim->ilboun, interface->xbounIndices );
			printf( "Read data '%s' found.\n", config->readDataNames[i] );
			break;
		}
		else if ( strcmp( config->readDataNames[i], "Heat-Flux" ) == 0 )
		{
			interface->readData = HEAT_FLUX;
			interface->xloadIndices = malloc( interface->numElements * sizeof( int ) );
			getXloadIndices( "DFLUX", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices );
			interface->fluxDataID = precicec_getDataID( "Heat-Flux", interface->faceCentersMeshID );
			printf( "Read data '%s' found.\n", config->readDataNames[i] );
			break;
		}
		else if ( strcmp1( config->readDataNames[i], "Sink-Temperature-" ) == 0 )
		{
			interface->readData = CONVECTION;
			interface->xloadIndices = malloc( interface->numElements * sizeof( int ) );
			getXloadIndices( "FILM", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices );
			interface->kDeltaTemperatureReadDataID = precicec_getDataID( config->readDataNames[i], interface->faceCentersMeshID );
			printf( "Read data '%s' found.\n", config->readDataNames[i] );
		}
		else if ( strcmp1( config->readDataNames[i], "Heat-Transfer-Coefficient-" ) == 0 )
		{
			interface->kDeltaReadDataID = precicec_getDataID( config->readDataNames[i], interface->faceCentersMeshID );
			printf( "Read data '%s' found.\n", config->readDataNames[i] );
		}
        	else if ( strcmp1( config->readDataNames[i], "Forces" + i ) == 0 )
        	{
            		PreciceInterface_EnsureValidNodesMeshID( interface );
            		interface->readData = FORCES;
            		interface->xforcIndices = malloc( interface->numNodes * 3 * sizeof( int ) );
            		interface->forcesDataID = precicec_getDataID( config->readDataNames[i], interface->nodesMeshID );
            		getXforcIndices( interface->nodeIDs, interface->numNodes, sim->nforc, sim->ikforc, sim->ilforc, interface->xforcIndices );
            		printf( "Read data '%s' found.\n", config->readDataNames[i] );
            		break;
            
        	}
		else
		{
			printf( "ERROR: Read data '%s' does not exist!\n", config->readDataNames[i] );
			exit( EXIT_FAILURE );
		}
	}
    
	for( i = 0 ; i < config->numWriteData ; i++ )
	{
		if( strcmp( config->writeDataNames[i], "Temperature" ) == 0 )
		{
			PreciceInterface_EnsureValidNodesMeshID( interface );
			interface->writeData = TEMPERATURE;
			interface->temperatureDataID = precicec_getDataID( "Temperature", interface->nodesMeshID );
			printf( "Write data '%s' found.\n", config->writeDataNames[i] );
			break;
		}
		else if ( strcmp( config->writeDataNames[i], "Heat-Flux" ) == 0 )
		{
			interface->writeData = HEAT_FLUX;
			interface->fluxDataID = precicec_getDataID( "Heat-Flux", interface->faceCentersMeshID );
			printf( "Write data '%s' found.\n", config->writeDataNames[i] );
			break;
		}
		else if ( strcmp1( config->writeDataNames[i], "Sink-Temperature-" ) == 0 )
		{
			interface->writeData = CONVECTION;
			interface->kDeltaTemperatureWriteDataID = precicec_getDataID( config->writeDataNames[i], interface->faceCentersMeshID );
			printf( "Write data '%s' found.\n", config->writeDataNames[i] );
		}
		else if ( strcmp1( config->writeDataNames[i], "Heat-Transfer-Coefficient-" ) == 0 )
		{
			interface->kDeltaWriteDataID = precicec_getDataID( config->writeDataNames[i], interface->faceCentersMeshID );
			printf( "Write data '%s' found.\n", config->writeDataNames[i] );
		}
                else if ( strcmp1( config->writeDataNames[i], "Displacements" + i ) == 0 )
                {
                        PreciceInterface_EnsureValidNodesMeshID( interface );
            		interface->writeData = DISPLACEMENTS;
                        interface->displacementsDataID = precicec_getDataID( config->writeDataNames[i], interface->nodesMeshID );
                        printf( "Write data '%s' found.\n", config->writeDataNames[i] );
                }
		else if ( strcmp1( config->writeDataNames[i], "DisplacementDeltas" + i ) == 0 )
                {
                        PreciceInterface_EnsureValidNodesMeshID( interface );
                        interface->writeData = DISPLACEMENTDELTAS;
                        interface->displacementDeltasDataID = precicec_getDataID( config->writeDataNames[i], interface->nodesMeshID );
                        printf( "Write data '%s' found.\n", config->writeDataNames[i] );
                }
		else
		{
			printf( "ERROR: Write data '%s' does not exist!\n", config->writeDataNames[i] );
			exit( EXIT_FAILURE );
		}
	}
}

void PreciceInterface_FreeData( PreciceInterface * preciceInterface )
{
	if( preciceInterface->elementIDs != NULL ){
		free( preciceInterface->elementIDs );
	}
	
	if( preciceInterface->faceIDs != NULL ){
		free( preciceInterface->faceIDs );
	}
	
	if( preciceInterface->faceCenterCoordinates != NULL ){
		free( preciceInterface->faceCenterCoordinates );
	}
	
	if( preciceInterface->preciceFaceCenterIDs != NULL ){
		free( preciceInterface->preciceFaceCenterIDs );
	}
	
	if( preciceInterface->nodeCoordinates != NULL ){
		free( preciceInterface->nodeCoordinates );
	}

	if( preciceInterface->preciceNodeIDs != NULL ){
		free( preciceInterface->preciceNodeIDs );
	}

	if( preciceInterface->triangles != NULL ){
		free( preciceInterface->triangles );
	}

	if( preciceInterface->nodeScalarData != NULL ){
		free( preciceInterface->nodeScalarData );
	}
	
	if( preciceInterface->nodeVectorData != NULL ){
    		free( preciceInterface->nodeVectorData );
	}
	
	if( preciceInterface->faceCenterData != NULL ){
		free( preciceInterface->faceCenterData );
	}

	if( preciceInterface->xbounIndices != NULL ){
		free( preciceInterface->xbounIndices );
	}

	if( preciceInterface->xloadIndices != NULL ){
		free( preciceInterface->xloadIndices );
	}
    
    	if ( preciceInterface->xforcIndices != NULL ) {
        	free( preciceInterface->xforcIndices );
    	}

}

