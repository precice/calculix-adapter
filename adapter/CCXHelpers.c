/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#include "CCXHelpers.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

char* toNodeSetName( char const * name )
{
	char * prefix = "N";
	char * suffix = "N";
	return concat( prefix, name, suffix );
}

char* toFaceSetName( char const * name )
{
	char * prefix = "S";
	char * suffix = "T";
	return concat( prefix, name, suffix );
}

ITG getSetID( char const * setName, char const * set, ITG nset )
{

	ITG i;
	ITG nameLength = 81;

	for( i = 0 ; i < nset ; i++ )
	{
		if( strcmp1( &set[i * nameLength], setName ) == 0 )
		{
			printf("Set ID Found \n");
			return i;
		}
	}
	// Set not found:

	if( setName[0] == (char) 'N' )
	{
		printf("Set ID NOT Found \n");
		nodeSetNotFoundError( setName );
	}
	else if ( setName[0] == (char) 'S' )
	{
		printf("Set ID NOT Found \n");
		faceSetNotFoundError( setName );
	}
  unreachableError();
  return -1;
}

ITG getNumSetElements( ITG setID, ITG * istartset, ITG * iendset )
{
	return iendset[setID] - istartset[setID] + 1;
}

void getSurfaceElementsAndFaces( ITG setID, ITG * ialset, ITG * istartset, ITG * iendset, ITG * elements, ITG * faces )
{

	ITG i, k = 0;

	for( i = istartset[setID]-1 ; i < iendset[setID] ; i++ )
	{
		elements[k] = ialset[i] / 10;
		faces[k] = ialset[i] % 10;
		k++;
	}
}

void getNodeCoordinates( ITG * nodes, ITG numNodes, int dim, double * co, double * v, int mt, double * coordinates )
{

	ITG i, j;

	for( i = 0 ; i < numNodes ; i++ )
	{
		int nodeIdx = nodes[i] - 1;
		//The displacements are added to the coordinates such that in case of a simulation restart the displaced coordinates are used for initializing the coupling interface instead of the initial coordinates
		for( j = 0 ; j < dim ; j++ )
    {
      coordinates[i * dim + j] = co[nodeIdx * 3 + j] + v[nodeIdx * mt + j + 1];
    }
	}
}

void getNodeTemperatures( ITG * nodes, ITG numNodes, double * v, int mt, double * temperatures )
{

	// CalculiX variable mt = 4 : temperature + 3 displacements (depends on the type of analysis)
	ITG i;

	for( i = 0 ; i < numNodes ; i++ )
	{
		int nodeIdx = nodes[i] - 1;
		temperatures[i] = v[nodeIdx * mt];
	}
}

void getNodeForces( ITG * nodes, ITG numNodes, int dim, double * fn, ITG mt, double * forces )
{
	ITG i, j;

	for ( i = 0 ; i < numNodes ; i++ ) 
	{
		int nodeIdx = nodes[i] - 1;
		for( j = 0 ; j < dim ; j++ )
    {
      forces[dim * i + j] = fn[nodeIdx * mt + j + 1];
    }
	}
}

void getNodeDisplacements( ITG * nodes, ITG numNodes, int dim, double * v, ITG mt, double * displacements )
{

	// CalculiX variable mt = 4 : temperature + 3 displacements (depends on the type of analysis)
	// where 0 index corresponds to temp; 1, 2, 3 indices correspond to the displacements, respectively
	ITG i, j;

	for( i = 0 ; i < numNodes ; i++ )
	{
		int nodeIdx = nodes[i] - 1; //The node Id starts with 1, not with 0, therefore, decrement is necessary
		for( j = 0 ; j < dim ; j++ )
  	{
      displacements[dim * i + j] = v[nodeIdx * mt + j + 1];
    }
	}
}

void getNodeDisplacementDeltas( ITG * nodes, ITG numNodes, int dim, double * v, double * v_init, ITG mt, double * displacementDeltas )
{

	// CalculiX variable mt = 4 : temperature + 3 displacements (depends on the type of analysis)
	// where 0 index corresponds to temp; 1, 2, 3 indices correspond to the displacements, respectively
	ITG i, j;

	for( i = 0 ; i < numNodes ; i++ )
	{
		int nodeIdx = nodes[i] - 1; //The node Id starts with 1, not with 0, therefore, decrement is necessary
		for( j = 0 ; j < dim ; j++ )
    {
      displacementDeltas[dim * i + j] = v[nodeIdx * mt + j + 1] - v_init[nodeIdx * mt + j + 1];
    }
	}
}

void getNodeVelocities( ITG * nodes, ITG numNodes, int dim, double * ve, ITG mt, double * velocities )
{

	// CalculiX variable mt = 4 : temperature rate + 3 velocities (depends on the type of analysis)
	// where 0 index corresponds to temp rate; 1, 2, 3 indices correspond to the velocities, respectively
	ITG i, j;

	for( i = 0 ; i < numNodes ; i++ )
	{
		int nodeIdx = nodes[i] - 1; //The node Id starts with 1, not with 0, therefore, decrement is necessary
		for( j = 0 ; j < dim ; j++ )
    {
      velocities[dim * i + j] = ve[nodeIdx * mt + j + 1];
    }
	}
}

/*
   int getNodesPerFace(char * lakon, int elementIdx) {

		int nodesPerFace;
		if(strcmp1(&lakon[elementIdx * 8], "C3D4") == 0) {
				nodesPerFace = 3;
		} else if(strcmp1(&lakon[elementIdx * 8], "C3D10") == 0) {
				nodesPerFace = 6;
		}
		return nodesPerFace;

   }
 */

void getTetraFaceCenters( ITG * elements, ITG * faces, ITG numElements, ITG * kon, ITG * ipkon, double * co, double * faceCenters )
{

	// Assume all tetra elements -- maybe implement checking later...

	// Node numbering for faces of tetrahedral elements (in the documentation the number is + 1)
	// Numbering is the same for first and second order elements
	int faceNodes[4][3] = { { 0,1,2 }, { 0,3,1 }, { 1,3,2 }, { 2,3,0 } };

	ITG i, j;

	for( i = 0 ; i < numElements ; i++ )
	{

		ITG faceIdx = faces[i] - 1;
		ITG elementIdx = elements[i] - 1;
		double x = 0, y = 0, z = 0;

		for( j = 0 ; j < 3 ; j++ )
		{

			ITG nodeNum = faceNodes[faceIdx][j];
			ITG nodeID = kon[ipkon[elementIdx] + nodeNum];
			ITG nodeIdx = ( nodeID - 1 ) * 3;
			// The nodeIdx is already multiplied by 3, therefore it must be divided by 3 ONLY when checking if coordinates match getNodeCoordinates

			x += co[nodeIdx + 0];
			y += co[nodeIdx + 1];
			z += co[nodeIdx + 2];

		}
		faceCenters[i * 3 + 0] = x / 3;
		faceCenters[i * 3 + 1] = y / 3;
		faceCenters[i * 3 + 2] = z / 3;

	}

	
}

/*
   void getSurfaceGaussPoints(int setID, ITG * co, ITG istartset, ITG iendset, ITG * ipkon, ITG * lakon, ITG * kon, ITG * ialset, double * coords) {

        int iset = setID + 1; // plus one because of fortran indices
        FORTRAN(getgausspointscoords,(co,&iset,istartset,iendset,ipkon,lakon,kon,ialset, coords));

   }
 */


void getTetraFaceNodes( ITG * elements, ITG * faces, ITG * nodes, ITG numElements, ITG numNodes, ITG * kon, ITG * ipkon, int * tetraFaceNodes )
{
	// Assume all tetra elements -- maybe implement checking later...

	// Node numbering for faces of tetrahedral elements (in the documentation the number is + 1)
	int faceNodes[4][3] = { { 0,1,2 }, { 0,3,1 }, { 1,3,2 }, { 2,3,0 } };

	ITG i, j, k;

	

	for( i = 0 ; i < numElements ; i++ )
	{

		ITG faceIdx = faces[i] - 1;
		ITG elementIdx = elements[i] - 1;

		for( j = 0 ; j < 3 ; j++ )
		{

			ITG nodeNum = faceNodes[faceIdx][j];
			ITG nodeID = kon[ipkon[elementIdx] + nodeNum];

			for( k = 0 ; k < numNodes ; k++ )
			{

				if( nodes[k] == nodeID )
				{
					tetraFaceNodes[i*3 + j] = k;
				}
			}
			
		}
	}
}

void getXloadIndices( char const * loadType, ITG * elementIDs, ITG * faceIDs, ITG numElements, ITG nload, ITG * nelemload, char const * sideload, ITG * xloadIndices )
{

	ITG i, k;
	int nameLength = 20;
	char faceLabel[] = { 'x', 'x', '\0' };

	/* Face number is prefixed with 'S' if it is DFLUX boundary condition
	 * and with 'F' if it is a FILM boundary condition */
	if( strcmp( loadType, "DFLUX" ) == 0 )
	{
		faceLabel[0] = (char) 'S';
	}
	else if( strcmp( loadType, "FILM" ) == 0 )
	{
		faceLabel[0] = (char) 'F';
	}

	for( k = 0 ; k < numElements ; k++ )
	{

		ITG faceID = faceIDs[k];
		ITG elementID = elementIDs[k];
		faceLabel[1] = faceID + '0';
		int found = 0;

		for( i = 0 ; i < nload ; i++ )
		{
			if( elementID == nelemload[i * 2] && strcmp1( &sideload[i * nameLength], faceLabel ) == 0 )
			{
				xloadIndices[k] = 2 * i;
				found = 1;
				break;
			}
		}

		// xload index not found:
		if( !found && strcmp( loadType, "DFLUX" ) == 0 )
		{
			missingDfluxBCError();
		}
		else if ( !found && strcmp( loadType, "FILM" ) == 0 )
		{
			missingFilmBCError();
		}


	}
}

// Get the indices for the xboun array, corresponding to the temperature or displacement DOF of the nodes passed to the function
void getXbounIndices( ITG * nodes, ITG numNodes, int nboun, int * ikboun, int * ilboun, int * xbounIndices, enum CouplingDataType couplDataType )
{
	ITG i;

	switch( couplDataType )
	{
	case TEMPERATURE:
		for( i = 0 ; i < numNodes ; i++ )
		{
			int idof = 8 * ( nodes[i] - 1 ) + 0; // 0 for temperature DOF
			int k;
			FORTRAN( nident, ( ikboun, &idof, &nboun, &k ) );
			k -= 1; // Adjust because of FORTRAN indices
			int m = ilboun[k] - 1; // Adjust because of FORTRAN indices
			xbounIndices[i] = m;
		}
		// See documentation ccx_2.13.pdf for the definition of ikboun and ilboun

		for( i = 0 ; i < numNodes ; i++ )
		{
			if( xbounIndices[i] < 0 )
			{
				missingTemperatureBCError();
			}
		}
    break;
	case DISPLACEMENTS:
		for( i = 0 ; i < numNodes ; i++ )
		{
			int idof_x = 8 * ( nodes[i] - 1 ) + 1; // x-component of displacement
			int idof_y = 8 * ( nodes[i] - 1 ) + 2; // y-component of displacement
			int idof_z = 8 * ( nodes[i] - 1 ) + 3; // z-component of displacement
			int kx, ky, kz;
			FORTRAN( nident, ( ikboun, &idof_x, &nboun, &kx ) );
			FORTRAN( nident, ( ikboun, &idof_y, &nboun, &ky ) );
			FORTRAN( nident, ( ikboun, &idof_z, &nboun, &kz ) );
			xbounIndices[3 * i] = ilboun[kx - 1] - 1;
			xbounIndices[3 * i + 1] = ilboun[ky - 1] - 1;
			xbounIndices[3 * i + 2] = ilboun[kz - 1] - 1;
		}

		// TODO This warning is not triggered if there are other BC (e.g. fixed nodes)
		for( i = 0 ; i < 3 * numNodes; i++ )
		{
			if( xbounIndices[i] < 0 )
			{
				missingDisplacementBCError();
			}
		}
    break;
  default:
    unreachableError();
	}
}

// Get the indices for the xforc array, corresponding to the force DOFs of the nodes passed to the function
void getXforcIndices( ITG * nodes, ITG numNodes, int nforc, int * ikforc, int * ilforc, int * xforcIndices )
{
	ITG i;

	for( i = 0 ; i < numNodes ; i++ )
	{
		//x-direction
		int idof = 8 * ( nodes[i] - 1 ) + 1; // 1 for x force DOF
		int k;
		FORTRAN( nident, ( ikforc, &idof, &nforc, &k ) );
		k -= 1; // Adjust because of FORTRAN indices
		int m = ilforc[k] - 1; // Adjust because of FORTRAN indices
		xforcIndices[3 * i] = m;

		//y-direction
		idof = 8 * ( nodes[i] - 1 ) + 2; // 2 for y force DOF
		FORTRAN( nident, ( ikforc, &idof, &nforc, &k ) );
		k -= 1; // Adjust because of FORTRAN indices
		m = ilforc[k] - 1; // Adjust because of FORTRAN indices
		xforcIndices[3 * i + 1] = m;

		//z-direction
		idof = 8 * ( nodes[i] - 1 ) + 3; // 3 for z force DOF
		FORTRAN( nident, ( ikforc, &idof, &nforc, &k ) );
		k -= 1; // Adjust because of FORTRAN indices
		m = ilforc[k] - 1; // Adjust because of FORTRAN indices
		xforcIndices[3 * i + 2] = m;
	}
	// See documentation ccx_2.10.pdf for the definition of ikforc and ilforc

	for( i = 0 ; i < numNodes ; i++ )
	{
		if( xforcIndices[i] < 0 )
		{
			missingForceError();
		}
	}
}


int getXloadIndexOffset( enum xloadVariable xloadVar )
{
	/*
	 * xload is the CalculiX array where the DFLUX and FILM boundary conditions are stored
	 * the array has two components:
	 * - the first component corresponds to the flux value and the heat transfer coefficient
	 * - the second component corresponds to the sink temperature
	 * */
	switch( xloadVar )
	{
	case DFLUX:
		return 0;
	case FILM_H:
		return 0;
	case FILM_T:
		return 1;
  default:
    unreachableError();
    return -1;
	}
}

void setXload( double * xload, int * xloadIndices, double * values, int numValues, enum xloadVariable xloadVar )
{
	ITG i;
	int indexOffset = getXloadIndexOffset( xloadVar );

	for( i = 0 ; i < numValues ; i++ )
	{
		xload[xloadIndices[i] + indexOffset] = values[i];
	}
}

void setFaceFluxes( double * fluxes, ITG numFaces, int * xloadIndices, double * xload )
{
	setXload( xload, xloadIndices, fluxes, numFaces, DFLUX );
}

void setFaceHeatTransferCoefficients( double * coefficients, ITG numFaces, int * xloadIndices, double * xload )
{
	setXload( xload, xloadIndices, coefficients, numFaces, FILM_H );
}

void setFaceSinkTemperatures( double * sinkTemperatures, ITG numFaces, int * xloadIndices, double * xload )
{
	setXload( xload, xloadIndices, sinkTemperatures, numFaces, FILM_T );
}

void setNodeTemperatures( double * temperatures, ITG numNodes, int * xbounIndices, double * xboun )
{
	ITG i;

	for( i = 0 ; i < numNodes ; i++ )
	{
		xboun[xbounIndices[i]] = temperatures[i];
	}
}

void setNodeForces( double * forces, ITG numNodes, int dim, int * xforcIndices, double * xforc )
{
	ITG i, j;

	for ( i = 0 ; i < numNodes ; i++ ) 
	{
		for( j = 0 ; j < dim ; j++ )
    {
      xforc[xforcIndices[3 * i + j]] = forces[dim * i + j];
    }
	}
}

void setNodeDisplacements( double * displacements, ITG numNodes, int dim, int * xbounIndices, double * xboun )
{
	ITG i, j;

	for( i = 0 ; i < numNodes ; i++ )
	{
		for( j = 0 ; j < dim ; j++ ) 
    {
      xboun[xbounIndices[3 * i + j]] = displacements[dim * i + j];
    }
	}
}

bool isSteadyStateSimulation( ITG * nmethod )
{
	return *nmethod == 1;
}

char* concat( char const * prefix, char const * string, char const * suffix )
{
	int nameLength = strlen( string ) + strlen( prefix ) + strlen( suffix ) + 1;
	char * result = malloc( nameLength );
	strcpy( result, prefix );
	strcat( result, string );
	strcat( result, suffix );
	return result;
}

bool startsWith(const char * string, const char * prefix)
{
    char sc, pc;
    do {
        sc =*string++;
        pc =*prefix++;
        if (pc == '\0') return true;
    } while (sc == pc);
    return false;
}

bool isEqual(const char * lhs, const char * rhs)
{
  return strcmp(lhs, rhs) == 0;
}

/* Errors messages */

void nodeSetNotFoundError( char const * setName )
{
	printf( "ERROR: Set %s does not exist! Please check that the interface names are correct and that .nam file is provided.\n", setName );
	fflush( stdout );
	exit( EXIT_FAILURE );
}

void faceSetNotFoundError( char const * setName )
{
	printf( "ERROR: Set %s does not exist! Please check the following: \n 1) If nearest projection mapping is required, check that the interface names are correct and that .sur file is provided.\n 2) If nearest-projection mapping is not required, remove 'nodes-mesh-mesh-connectivity' and replace with 'nodes-mesh' in the config.yml file", setName );
	fflush( stdout );
	exit( EXIT_FAILURE );
}

void missingTemperatureBCError()
{
	printf( "ERROR: Cannot apply temperature BC to one or more interface nodes.  Please make sure that a temperature boundary condition is set for the interface, when using a Dirichlet coupling BC.\n" );
	exit( EXIT_FAILURE );
}

void missingForceError()
{
	printf( "ERROR: Cannot apply forces to one or more interface nodes.\n" );
	exit( EXIT_FAILURE );
}

void missingDisplacementBCError()
{
	printf( "ERROR: Cannot apply displacement to one or more interface nodes. Please make sure that a single point constraint in each direction is set for all interface nodes.\n" );
	exit( EXIT_FAILURE );
}

void missingDfluxBCError()
{
	printf( "ERROR: Cannot apply DFLUX BC to one or more interface elements.  Please make sure that a .dfl file is provided for the interface, when using a Neumann coupling BC.\n" );
	exit( EXIT_FAILURE );
}

void missingFilmBCError()
{
	printf( "ERROR: Cannot apply FILM BC to one or more interface elements.  Please make sure that a .flm file is provided for the interface, when using a Robin coupling BC.\n" );
	exit( EXIT_FAILURE );
}

void unreachableError()
{
	printf( "ERROR: The preCICE adapter just entered an unreachable state. Something is very wrong!\n" );
	exit( EXIT_FAILURE );
}

