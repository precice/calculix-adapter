/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#ifndef CCXHELPERS_H
#define CCXHELPERS_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../CalculiX.h"

/*
 * These are some helper functions for handling CalculiX data structures.
 * The original names of the CalculiX variables were kept.
 * Please refer to the CalculiX documentation ccx_2.10.pdf on page 518
 * for a "List of variables and their meaning"
 * */

/**
 * @brief Possible variables in xload (CalculiX variable for the load)
 */
enum xloadVariable { DFLUX, FILM_H, FILM_T };

/**
 * @brief Type of coupling data
 *  Temperature - Dirichlet
 *  Heat Flux - Neumann
 *  Sink Temperature & Heat Transfer Coefficient - Robin
 *  Forces - dynamics data to be read/written (by the Calculix adapter)
 *  Displacements - dynamics data to be read/written (by the Calculix adapter)
 *  DisplacementDeltas - FSI data to be written (by the Calculix adapter)
 *  Velocities - FSI data to be written (by the Calculix adapter)
 *  Positions - FSI data to be written (by the Calculix adapter)
 */
enum CouplingDataType { TEMPERATURE, HEAT_FLUX, SINK_TEMPERATURE, HEAT_TRANSFER_COEFF, FORCES, DISPLACEMENTS, DISPLACEMENTDELTAS, VELOCITIES, POSITIONS };

/**
 * @brief Returns node set name with internal CalculiX format
 * Prepends and appends an N: e.g. If the input name is "interface",
 * it returns NinterfaceN
 */
char* toNodeSetName( char const * name );

/**
 * @brief Returns face set name with internal CalculiX format
 * Prepends an S and appends a T: e.g. If the input name is "interface",
 * it returns SinterfaceT
 */
char* toFaceSetName( char const * name );

/**
 * @brief Returns id of a set given its name
 * @param setName: set name as returned by toNodeSetName or toFaceSetName
 * @param set: CalculiX array for all the set names
 * @param nset: CalculiX variable for the number of sets
 */
ITG getSetID( char const * setName, char const * set, ITG nset );

/**
 * @brief Returns number of elements in a set
 */
ITG getNumSetElements( ITG setID, ITG * istartset, ITG * iendset );

/**
 * @brief Gets the element and face IDs given a set ID
 * @param setID: input set id
 * @param ialset: CalculiX variable
 * @param istartset: CalculiX variable
 * @param iendset: CalculiX variable
 * @param elements: output element IDs
 * @param faces: output face IDs (local IDs: e.g. 1, 2, 3, 4 for tetrahedral elements)
 */
void getSurfaceElementsAndFaces( ITG setID, ITG * ialset, ITG * istartset, ITG * iendset, ITG * elements, ITG * faces );

/**
 * @brief Gets the coordinates of a list of input node IDs
 * @param nodes: input node IDs
 * @param numNodes: number of input nodes
 * @param co: CalculiX array with all the node coordinates
 * @param v: CalculiX array containing the solution variables (displacements needed here)
 * @param mt: CalculiX variable describing the number of solution variables in the solution array v
 * @param coordinates: output array with the coordinates of the input nodes
 */
void getNodeCoordinates( ITG * nodes, ITG numNodes, int dim, double * co, double * v, int mt, double * coordinates );

/**
 * @brief getNodeTemperatures
 * @param nodes: input node IDs
 * @param numNodes: number of input nodes
 * @param dim: coupling dimension
 * @param v: CalculiX array containing the temperatures
 * @param mt: CalculiX variable describing the number of solution variables in the solution array v
 * @param temperatures: output array with the temperatures of the input nodes
 */
void getNodeTemperatures( ITG * nodes, ITG numNodes, double * v, int mt, double * temperatures );

/**
 * @brief getNodeForces
 * @param nodes: input node IDs
 * @param numNodes: number of nodes
 * @param fn: Calculix array containing the actual forces
 * @param mt: CalculiX variable describing the number of force variables in fn
 * @param forces: output array with the forces of the input nodes
 */
void getNodeForces( ITG * nodes, ITG numNodes, int dim, double * fn, ITG mt, double * forces );

/**
 * @brief getNodeDisplacements
 * @param nodes: input node IDs
 * @param numNodes: number of input nodes
 * @param dim: coupling dimension
 * @param v: CalculiX solution array containing the displacements
 * @param mt: CalculiX variable describing the number of solution variables in the solution array v
 * @param displacements: output array with the displacements in preCICE-conform order of the input nodes
 */
void getNodeDisplacements( ITG * nodes, ITG numNodes, int dim, double * v, ITG mt, double * displacements );

/**
 * @brief getNodeVelocities
 * @param nodes: input node IDs
 * @param numNodes: number of input nodes
 * @param dim: coupling dimension
 * @param ve: CalculiX solution array containing the velocities
 * @param mt: CalculiX variable describing the number of solution variables in the solution array v
 * @param velocities: output array with the velocities in preCICE-conform order of the input nodes
 */
void getNodeVelocities( ITG * nodes, ITG numNodes, int dim, double * ve, ITG mt, double * velocities );

/**
 * @brief getNodeDisplacementDeltas
 * @param nodes: input node IDs
 * @param numNodes: number of input nodes
 * @param dim: coupling dimension
 * @param v: CalculiX solution array containing the displacements
 * @param v_init: Calculix solution array containing the displacements of the last converged step
 * @param mt: CalculiX variable describing the number of solution variables in the solution array v
 * @param displacementDeltas: output array with the displacementDeltas in preCICE-conform order of the input nodes
 */
void getNodeDisplacementDeltas( ITG * nodes, ITG numNodes, int dim, double * v, double * v_init, ITG mt, double * displacementDeltas );

/**
 * @brief Computes the center of one of the faces of a tetrahedral element
 * @param elements: input list of tetrahedral elements
 * @param faces: input list of local face IDs
 * @param numElements: number of input elements
 * @param kon: CalculiX variable
 * @param ipkon: CalculiX variable
 * @param co: CalculiX array with the coordinates of all the nodes
 * @param faceCenters: output array with the face centers of the input element faces
 */
void getTetraFaceCenters( ITG * elements, ITG * faces, ITG numElements, ITG * kon, ITG * ipkon, double * co, double * faceCenters );

/**
 * @brief Gets a list of node IDs from a list of input element faces
 * @param elements: list of element IDs
 * @param faces: list of local face IDs
 * @param nodes: list of node IDs
 * @param numElements: number of input elements
 * @param numNodes: number of input nodes
 * @param kon: CalculiX array with the connectivity information
 * @param ipkon: CalculiX array (see description in ccx_2.10.pdf)
 * @param tetraFaceNodes: output list of node IDs that belong to the input element faces
 */
void getTetraFaceNodes( ITG * elements, ITG * faces, ITG * nodes, ITG numElements, ITG numNodes, ITG * kon, ITG * ipkon, int * tetraFaceNodes );

/**
 * @brief Gets the indices of the xload where the DFLUX and FILM boundary conditions must be applied
 * @param loadType: DFLUX or FILM
 * @param elementIDs: list of IDs of elements on which the boundary conditions must be applied
 * @param faceIDs: list of local face IDs
 * @param numElements: number of elements
 * @param nload: CalculiX variable that indicates the size of the xload array
 * @param nelemload: CalculiX array containing the element IDs with DFLUX and FILM boundary conditions
 * @param sideload: CalculiX array containing the faces to which the DFLUX or FILM boundary conditions are applied
 * @param xloadIndices: output list of indices of the xload array
 */
void getXloadIndices( char const * loadType, ITG * elementIDs, ITG * faceIDs, ITG numElements, ITG nload, ITG * nelemload, char const * sideload, ITG * xloadIndices );

/**
 * @brief Gets the indices of the xboun array where the boundary conditions must be applied
 * @param nodes: list of node IDs
 * @param numNodes: number of nodes
 * @param nboun: CalculiX variable for the number of SPCs (single point constraints)
 * @param ikboun: CalculiX ordered array of the DOFs corresponding to the SPCs
 * @param ilboun: CalculiX auxiliary array for dealing with the indices of the SPCs
 * @param xbounIndices: output list of indices of the xboun array
 * @param couplDataType: indicates whether indices must be for temperature or displacement BC
 */
void getXbounIndices( ITG * nodes, ITG numNodes, int nboun, int * ikboun, int * ilboun, int * xbounIndices, enum CouplingDataType couplDataType );

/**
 * @brief Gets the indices of the xforc array where the forces must be applied
 * @param nodes: list of node IDs
 * @param numNodes: number of nodes
 * @param nforc: CalculiX variable for the number of point loads (point forces)
 * @param ikforc: CalculiX ordered array of the DOFs corresponding to the point loads
 * @param ilforc: CalculiX auxiliary array for dealing with the indices of the point loads
 * @param xforcIndices: output list of indices of the xforc array
 */
void getXforcIndices( ITG * nodes, ITG numNodes, int nforc, int * ikforc, int * ilforc, int * xforcIndices );

/**
 * @brief Modifies the values of a DFLUX or FILM boundary condition
 * @param xload: CalculiX array for the loads
 * @param xloadIndices: list of indices where the values must be set in the xload array
 * @param values: boundary values to apply
 * @param numValues: number of boundary values provided
 * @param xloadVar: variable that is actually modified: DFLUX for heat flux, FILM_H for heat transfer coeff, FILM_T for sink temperature
 */
void setXload( double * xload, int * xloadIndices, double * values, int numValues, enum xloadVariable xloadVar );

/**
 * @brief Calls setXload to update the flux values at the specified indices
 * @param fluxes: flux values
 * @param numFaces: number of faces
 * @param xloadIndices: indices of the xload array where the values must be updated
 * @param xload: CalculiX array for the loads
 */
void setFaceFluxes( double * fluxes, ITG numFaces, int * xloadIndices, double * xload );

/**
 * @brief Calls setXload to update the heat transfer coefficients at the specified indices
 * @param coefficients: values for the heat transfer coefficient
 * @param numFaces: number of faces
 * @param xloadIndices: indices of the xload array where the values must be updated
 * @param xload: CalculiX array for the loads
 */
void setFaceHeatTransferCoefficients( double * coefficients, ITG numFaces, int * xloadIndices, double * xload );

/**
 * @brief Calls setXload to update the sink temperature at the specified indices
 * @param sinkTemperatures: values for the sink temperature
 * @param numFaces: number of faces
 * @param xloadIndices: indices of the xload array where the values must be updated
 * @param xload: CalculiX array for the loads
 */
void setFaceSinkTemperatures( double * sinkTemperatures, ITG numFaces, int * xloadIndices, double * xload );

/**
 * @brief Modifies the values of the temperature boundary condition
 * @param temperatures: temperature values to apply
 * @param numNodes: number of nodes
 * @param xbounIndices: indices of the xboun array to modify
 * @param xboun: CalculiX array containing temperature and displacement boundary values
 */
void setNodeTemperatures( double * temperatures, ITG numNodes, int * xbounIndices, double * xboun );

/**
 * @brief Modifies the values of the concentrated loads (point forces) applied to the structure
 * @param nodes: list of node IDs
 * @param forces: force values to apply
 * @param numNodes: number of nodes
 * @param dim: coupling dimension
 * @param xforcIndices: indices of the xforc array to modify
 * @param xforc: CalculiX array containing the (componentwise) assigned force values
 */
void setNodeForces( double * forces, ITG numNodes, int dim, int * xforcIndices, double * xforc );

/**
 * @brief Modifies the values of the displacements at the interface, as a Dirichlet boundary condition
 * @param displacements: displacement values to apply
 * @param numNodes: number of nodes
 * @param dim: coupling dimension
 * @param xbounIndices: indices of the xboun array to modify
 * @param xboun: CalculiX array containing temperature and displacement boundary values
 */
void setNodeDisplacements( double * displacements, ITG numNodes, int dim, int * xbounIndices, double * xboun );

/**
 * @brief Returns whether it is a steady-state simulation based on the value of nmethod
 * @param nmethod: CalculiX variable with information regarding the type of analysis
 */
bool isSteadyStateSimulation( ITG * nmethod );

/**
 * @brief Concatenates prefix + string + suffix
 * @param prefix
 * @param string
 * @param suffix
 */
char * concat(char const * prefix, char const * string, char const * suffix);

/**
 * @brief Checks wheather one zero-terminated string is prefixed by another
 * @param string the string to inspect
 * @param prefix the prefix to look for
 */
bool startsWith(const char * string, const char * prefix);

/**
 * @brief Checks wheather two zero-terminated strings are identical
 * @param lhs the left-hand string
 * @param rhs the right-hand string
 */
bool isEqual(const char * lhs, const char * rhs);

/* Error messages */

/**
 * @brief Terminate program if a node set is not defined for the interface (e.g. missing interface.nam file)
 * @param setName
 */
void nodeSetNotFoundError( char const * setName );

/**
 * @brief Terminate program if a face set is not defined for the interface (e.g. missing interface.sur file)
 * @param setName
 */
void faceSetNotFoundError( char const * setName );

/**
 * @brief Terminate program if a temperature BC is not defined when using Dirichlet BC for coupling (e.g. missing line under *BOUNDARY)
 */
void missingTemperatureBCError();

/**
 * @brief Terminate program if forces are not defined
 */
void missingForceError();

/**
 * @brief Terminate program if displacement Dirichlet BC are not defined for the interface in .inp file. (e.g. missing line under *BOUNDARY)
 */
void missingDisplacementBCError();

/**
 * @brief Terminate program if a DFLUX BC is not defined when using Neumann BC for coupling (e.g. missing interface.dfl file)
 */
void missingDfluxBCError();

/**
 * @brief Terminate program if a FILM BC is not defined when using Robin BC for coupling (e.g. missing interface.flm file)
 */
void missingFilmBCError(); 

/**
 * @brief Terminate program if the adapter reaches an unreachable state. This should never occur.
 */
void unreachableError(); 

#endif // CCXHELPERS_H
