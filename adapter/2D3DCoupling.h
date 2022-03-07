/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#ifndef CCX_2D3D_H
#define CCX_2D3D_H

/**
 * @brief Structure for holding data of a 2D-3D (nodal) mapping. It also handles the 2D preCICE mesh.
 * Usage: initialize with geometry data and mesh ID. This will compute the mapping and create the mesh.
 * Then:
 * - to write data: fill the scalar or vector 3D buffer and call the relevant function with the data preCICE ID.
 * - to read data: Call the relevant function with data ID. Read data will be stored in the scalar of vector 3D buffer and can be recovered.
 *  Don't forget the cleanup call.
 */

typedef struct {
  double *pos2D;
  int *   preciceNodesIDs;
  int *   numParentNodes;
  int *   mapping3D2D;
  int     num2DNodes;
  int     num3DNodes;
  double *bufferVector2D;
  double *bufferScalar2D;
  double *bufferVector3D;
  double *bufferScalar3D;

} Mapping2D3D;

Mapping2D3D *createMapping(const double *nodeCoordinates, int num3Dnodes, int nodesMeshID);
void         freeMapping(Mapping2D3D *);
void         consistentScalarRead(Mapping2D3D *map, int dataID);
void         consistentVectorRead(Mapping2D3D *map, int dataID);
void         conservativeScalarRead(Mapping2D3D *map, int dataID);
void         conservativeVectorRead(Mapping2D3D *map, int dataID);
void         consistentScalarWrite(Mapping2D3D *map, int dataID);
void         consistentVectorWrite(Mapping2D3D *map, int dataID);
void         conservativeScalarWrite(Mapping2D3D *map, int dataID);
void         conservativeVectorWrite(Mapping2D3D *map, int dataID);

#endif
