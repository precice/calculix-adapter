/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#include "2D3DCoupling.hpp"
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <vector>
#include "precice/SolverInterfaceC.h"

namespace MappingHelper {

struct Point2D {
  double x;
  double y;

  friend bool operator<(const Point2D &lhs, const Point2D &rhs)
  {
    if (((lhs.x - rhs.x) * (lhs.x - rhs.x)) + (lhs.y - rhs.y) * (lhs.y - rhs.y) < 1e-8) {
      // Assume equality up to tolerance
      return false;
    }

    // Lexicographical sort
    if (lhs.x < rhs.x) {
      return true;
    } else if (lhs.x > rhs.x) {
      return false;
    } else {
      //Exactly same x => compare y
      return lhs.y < rhs.y;
    }
  }
};

using Bijection = std::map<Point2D, std::vector<int>>;

} // namespace MappingHelper

static void setDoubleArrayZero(double *values, const int length, const int dim)
{
  int i, j;

  for (i = 0; i < length; i++) {
    for (j = 0; j < dim; j++) {
      values[i * dim + j] = 0.0;
    }
  }
}

Mapping2D3D *createMapping(const double *nodeCoordinates, int num3Dnodes, int nodesMeshID)
{

  Mapping2D3D *map = (Mapping2D3D *) malloc(sizeof(Mapping2D3D));

  map->num3DNodes = num3Dnodes;

  // For each 2D point, keep track of relevant 3D points
  MappingHelper::Bijection helper;

  for (int i = 0; i < num3Dnodes; ++i) {
    MappingHelper::Point2D pos{.x = nodeCoordinates[3 * i], .y = nodeCoordinates[3 * i + 1]};
    helper[pos].emplace_back(i);
  }

  map->num2DNodes = helper.size();
  printf("2D-3D mapping results in a 2D mesh of size %d from the %d points in 3D space.\n", map->num2DNodes, num3Dnodes);

  // Setup 3D to 2D map, 2D ancestors count and 2D mesh
  map->mapping3D2D    = (int *) malloc(num3Dnodes * sizeof(int));
  map->numParentNodes = (int *) malloc(map->num2DNodes * sizeof(int));
  map->pos2D          = (double *) malloc(map->num2DNodes * 2 * sizeof(double));

  int counter = 0;

  for (const auto &kv : helper) {
    MappingHelper::Point2D pos       = kv.first;
    const auto &           indices3D = kv.second;

    map->pos2D[counter * 2]     = pos.x;
    map->pos2D[counter * 2 + 1] = pos.y;

    for (int i : indices3D) {
      map->mapping3D2D[i] = counter;
    }

    map->numParentNodes[counter] = indices3D.size();

    ++counter;
  }

  assert(counter == map->num2DNodes);

  // Setup preCICE mesh
  map->preciceNodesIDs = (int *) malloc(map->num2DNodes * sizeof(int));
  precicec_setMeshVertices(nodesMeshID, map->num2DNodes, map->pos2D, map->preciceNodesIDs);

  // Initialize buffers
  map->bufferScalar2D = (double *) malloc(map->num2DNodes * sizeof(double));
  map->bufferVector2D = (double *) malloc(2 * map->num2DNodes * sizeof(double));
  map->bufferScalar3D = (double *) malloc(map->num3DNodes * sizeof(double));
  map->bufferVector3D = (double *) malloc(3 * map->num3DNodes * sizeof(double));

  return map;
}

void freeMapping(Mapping2D3D *map)
{
  if (map == nullptr)
    return;

  free(map->pos2D);
  free(map->preciceNodesIDs);
  free(map->mapping3D2D);
  free(map->numParentNodes);
  free(map->bufferVector2D);
  free(map->bufferScalar2D);
  free(map->bufferVector3D);
  free(map->bufferScalar3D);
  free(map);
}

void consistentScalarRead(Mapping2D3D *map, int dataID)
{
  precicec_readBlockScalarData(dataID,
                               map->num2DNodes,
                               map->preciceNodesIDs,
                               map->bufferScalar2D);

  // For each 3D point, copy the 2D counterpart
  for (int i = 0; i < map->num3DNodes; ++i) {
    map->bufferScalar3D[i] = map->bufferScalar2D[map->mapping3D2D[i]];
  }
}

void consistentVectorRead(Mapping2D3D *map, int dataID)
{
  precicec_readBlockVectorData(dataID,
                               map->num2DNodes,
                               map->preciceNodesIDs,
                               map->bufferVector2D);

  // For each 3D point, copy the 2D counterpart, component-wise
  for (int i = 0; i < map->num3DNodes; ++i) {
    map->bufferVector3D[3 * i]     = map->bufferVector2D[map->mapping3D2D[2 * i]];
    map->bufferVector3D[3 * i + 1] = map->bufferVector2D[map->mapping3D2D[2 * i + 1]];
  }
}

void conservativeScalarRead(Mapping2D3D *map, int dataID)
{
  precicec_readBlockScalarData(dataID,
                               map->num2DNodes,
                               map->preciceNodesIDs,
                               map->bufferScalar2D);

  // For each 3D point, take a fraction of the value in 2D
  for (int i = 0; i < map->num3DNodes; ++i) {
    map->bufferScalar3D[i] = map->bufferScalar2D[map->mapping3D2D[i]] / map->numParentNodes[map->mapping3D2D[i]];
  }
}

void conservativeVectorRead(Mapping2D3D *map, int dataID)
{
  precicec_readBlockVectorData(dataID,
                               map->num2DNodes,
                               map->preciceNodesIDs,
                               map->bufferVector2D);

  // For each 3D point, take a fraction of the value in 2D
  for (int i = 0; i < map->num3DNodes; ++i) {
    map->bufferVector3D[3 * i]     = map->bufferVector2D[map->mapping3D2D[i] * 2] / map->numParentNodes[map->mapping3D2D[i]];
    map->bufferVector3D[3 * i + 1] = map->bufferVector2D[map->mapping3D2D[i] * 2 + 1] / map->numParentNodes[map->mapping3D2D[i]];
  }
}

void consistentScalarWrite(Mapping2D3D *map, int dataID)
{
  // For each 2D point, write the average of the relevant 3D points,
  // in two steps: sum then divide
  setDoubleArrayZero(map->bufferScalar2D, map->num2DNodes, 1);

  // Assumption: CalculiX data written in the buffer 3D

  for (int i = 0; i < map->num3DNodes; ++i) {
    map->bufferScalar2D[map->mapping3D2D[i]] += map->bufferScalar3D[i];
  }

  for (int i2d = 0; i2d < map->num2DNodes; ++i2d) {
    map->bufferScalar2D[i2d] /= map->numParentNodes[i2d];
  }

  precicec_writeBlockScalarData(dataID,
                                map->num2DNodes,
                                map->preciceNodesIDs,
                                map->bufferScalar2D);
}

void consistentVectorWrite(Mapping2D3D *map, int dataID)
{
  // For each 2D point, write the average of the relevant 3D points,
  // in two steps: sum then divide
  setDoubleArrayZero(map->bufferVector2D, map->num2DNodes, 2);

  // Assumption: CalculiX data written in the buffer 3D
  for (int i = 0; i < map->num3DNodes; ++i) {
    map->bufferVector2D[2 * map->mapping3D2D[i]] += map->bufferVector3D[3 * i];
    map->bufferVector2D[2 * map->mapping3D2D[i] + 1] += map->bufferVector3D[3 * i + 1];
  }

  for (int i2d = 0; i2d < map->num2DNodes; ++i2d) {
    map->bufferVector2D[2 * i2d] /= map->numParentNodes[i2d];
    map->bufferVector2D[2 * i2d + 1] /= map->numParentNodes[i2d];
  }

  precicec_writeBlockVectorData(dataID,
                                map->num2DNodes,
                                map->preciceNodesIDs,
                                map->bufferVector2D);
}

void conservativeScalarWrite(Mapping2D3D *map, int dataID)
{
  // For each 2D point, write the sum of the relevant 3D points,
  setDoubleArrayZero(map->bufferScalar2D, map->num2DNodes, 1);

  // Assumption: CalculiX data written in the buffer 3D

  for (int i = 0; i < map->num3DNodes; ++i) {
    map->bufferScalar2D[map->mapping3D2D[i]] += map->bufferScalar3D[i];
  }

  precicec_writeBlockScalarData(dataID,
                                map->num2DNodes,
                                map->preciceNodesIDs,
                                map->bufferScalar2D);
}

void consservativeVectorWrite(Mapping2D3D *map, int dataID)
{
  // For each 2D point, write the sum of the relevant 3D points,
  setDoubleArrayZero(map->bufferVector2D, map->num2DNodes, 2);

  // Assumption: CalculiX data written in the buffer 3D
  for (int i = 0; i < map->num3DNodes; ++i) {
    map->bufferVector2D[2 * map->mapping3D2D[i]] += map->bufferVector3D[3 * i];
    map->bufferVector2D[2 * map->mapping3D2D[i] + 1] += map->bufferVector3D[3 * i + 1];
  }

  precicec_writeBlockVectorData(dataID,
                                map->num2DNodes,
                                map->preciceNodesIDs,
                                map->bufferVector2D);
}
