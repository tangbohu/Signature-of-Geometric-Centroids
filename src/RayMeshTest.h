///////////////////////////////////////////////////////////////
//
// RayMeshTest.h
//
//   Intersection Test between Axis-aligned Ray and Triangles
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 02/Apr/2015
//
///////////////////////////////////////////////////////////////

#ifndef _RAY_MESH_TEST_H
#define _RAY_MESH_TEST_H

#include "trimesh\include\Vec.h"
#include <Vector>

using namespace trimesh;
using namespace std;

struct Triangle;

struct HitPoint 
{
	double dist;     // Distance to ray origin
	Vector3f pos;    // Hit point position
	Vector3f nor;    // Normal of hit triangle
};


///////////////////////////////////////////////////////////////
// Function Declaration
///////////////////////////////////////////////////////////////

// Intersection Test for General Rays
bool LineMeshIntersect(Vector3f rayOrg, Vector3f rayDir, vector<Triangle*> _triList, HitPoint &hitPoint);

// Intersection Test for Axis-aligned Rays
vector<HitPoint> RayMeshIntersect(Vector3f rayOrg, Vector3f rayDir, vector<Triangle*> _triList);
void RemoveDuplicatedPoints(vector<HitPoint> &hitPtList);
bool IsDepthInList(double testDepth, vector<double> depthList);
void GetTriangleBBox(Triangle tri, Vector3f &boxMinPt, Vector3f &boxMaxPt);
bool RayBoxIntersect(Vector3f rayOrg, Vector3f rayDir, Vector3f boxMinPt, Vector3f boxMaxPt);


#endif