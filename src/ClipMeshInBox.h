///////////////////////////////////////////////////////////////
//
// ClipMeshInBox.h
//
//   Clip Mesh that is Inside a Box
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 07/May/2015
//
//
///////////////////////////////////////////////////////////////

#ifndef _CLIP_MESH_IN_BOX_H
#define _CLIP_MESH_IN_BOX_H

#include "Vec.h"
#include <Vector>
#include "HelpStructs.h"

using namespace std;
using namespace trimesh;

struct TempPoint
{
	Vector3f point;    // Point position
	float dist;        // Distance between the point and a plane
};

struct Plane
{
	Vector3f normal;
	Vector3f point;
};


///////////////////////////////////////////////////////////////
// Function Declaration
///////////////////////////////////////////////////////////////

vector<Triangle*> ClipMeshInBox(vector<Triangle*> origTris, Vector3f boxCen, Vector3f boxSize);
vector<Triangle*> CopyTriangles(vector<Triangle*> triList);
vector<Plane> CreateClipPlanes(Vector3f center, Vector3f size);

vector<Triangle*> ClipTrisWithPlane(vector<Triangle*> triList, Vector3f planeNor, Vector3f planePt);
vector<Triangle*> PlaneTriTest(Vector3f planeNor, Vector3f planePt, Triangle *triangle);
Vector3f LerpPoint(Vector3f pointA, float distA, Vector3f pointB, float distB);

#endif